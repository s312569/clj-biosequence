(ns clj-biosequence.core
  (:require [clojure.java.io :refer [writer reader]]
            [fs.core :refer [file? absolute-path]]
            [clj-http.client :refer [post]]
            [clojure.string :refer [trim split]]
            [taoensso.nippy :refer [freeze thaw]]
            [clj-biosequence.alphabet :as ala]
            [clojure.edn :as edn]
            [clojure.data.xml :as xml]
            [miner.tagged :as tag]
            [taoensso.nippy :as nip]
            [iota :as iot]
            [clojure.core.reducers :as r])
  (:import [java.io RandomAccessFile]))

(declare init-fasta-store init-fasta-sequence translate init-indexed-fasta)

(defprotocol biosequenceIO
  (bs-reader [this]
    "Returns a reader for a file containing biosequences. Use with
    `with-open'"))

(defprotocol biosequenceReader
  (biosequence-seq [this]
    "Returns a lazy sequence of biosequence objects.")
  (biosequence-seq-lazy [this]
    "Returns a lazy sequence of biosequence objects with a lazy stream
    of sequences.")
  (parameters [this]
    "Returns parameters, if any.")
  (get-biosequence [this accession]
    "Returns the biosequence object with the corresponding
    accession."))

(defprotocol Biosequence
  (accession [this]
    "Returns the accession of a biosequence object.")
  (accessions [this]
    "Returns a list of accessions for a biosequence object.")
  (def-line [this]
    "Returns the description of a biosequence object.")
  (bs-seq [this]
    "Returns the sequence of a biosequence as a vector.")
  (fasta-string [this]
    "Returns the biosequence as a string in fasta format.")
  (protein? [this]
    "Returns true if a protein and false otherwise.")
  (alphabet [this]
    "Returns the alphabet of a biosequence.")
  (feature-seq [this]
    "Returns a lazy list of features in a sequence."))

(defprotocol biosequenceFile
  (bs-path [this]
    "Returns the path of the file as string.")
  (index-file [this] [this ofile]))

(defprotocol biosequenceCitation
  (ref-type [this]
    "Returns the citation type from a citation object.")
  (title [this]
    "Returns the title of a citation object.")
  (journal [this]
    "Returns the journal of a citation object.")
  (year [this]
    "Returns the year of a citation object.")
  (volume [this]
    "Returns the volume of a citation object.")
  (pstart [this]
    "Returns the start page of a citation object.")
  (pend [this]
    "Returns the end page of a citation object.")
  (authors [this]
    "Returns the authors from a citation object."))

(defprotocol biosequenceFeature
  (interval-seq [this]
    "Returns a lazy list of intervals in a sequence.")
  (feature-type [this]
    "Returns the feature key."))

(defprotocol biosequenceInterval
  (start [this]
    "Returns the start position of an interval as an integer.")
  (end [this]
    "Returns the end position of an interval as an integer.")
  (comp? [this]
    "Is the interval complementary to the biosequence sequence. Boolean"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn bioseq->string
  "Returns the sequence of a biosequence as a string."
  [bs]
  (apply str (interpose "\n" (map #(apply str %) (partition-all 80 (bs-seq bs))))))


(defn sub-bioseq
  "Returns a new fasta sequence object with the sequence corresponding
   to 'beg' (inclusive) and 'end' (exclusive) of 'bs'. If no 'end'
   argument returns from 'start' to the end of the sequence. Zero
   based index."
  ([bs beg] (sub-bioseq bs beg nil))
  ([bs beg end]
     (init-fasta-sequence (accession bs)
                          (str (def-line bs)
                               " [" beg " - "
                               (if end end "End") "]")
                          (alphabet bs)
                          (if end
                            (subvec (bs-seq bs) beg end)
                            (subvec (bs-seq bs) beg)))))

(defn reverse-comp [this]
  "Returns a new fastaSequence with the reverse complement sequence."
  (if (protein? this)
    (throw (IllegalArgumentException. "Can't reverse/complement a protein sequence."))
    (init-fasta-sequence (accession this)
                         (str (def-line this) " - Reverse-comp")
                         (alphabet this)
                         (ala/revcom (bs-seq this)))))

(defn reverse-seq [this]
  "Returns a new fastaSequence with the reverse sequence."
  (init-fasta-sequence (accession this)
                       (str (def-line this) " - Reversed")
                       (alphabet this) 
                       (vec (reverse (bs-seq this)))))

(defn translate
  "Returns a fastaSequence sequence with a sequence translated in the
   specified frame."
  [bs frame & {:keys [table id-alter] :or {table (ala/codon-tables 1) id-alter true}}]
  (let [f ({1 1 2 2 3 3 4 4 5 5 6 6 -1 4 -2 5 -3 6} frame)]
    (cond (protein? bs)
          (throw (IllegalArgumentException. "Can't translate a protein sequence!"))
          (not (#{1 2 3 4 5 6} f))
          (throw (IllegalArgumentException. (str "Invalid frame: " frame)))
          :else
          (init-fasta-sequence (if id-alter (str (accession bs) "-" f) (accession bs))
                               (str (def-line bs) " - Translated frame: " f)
                               :iupacAminoAcids
                               (let [v (cond (#{1 2 3} f)
                                             (sub-bioseq bs (- f 1))
                                             (#{4 5 6} f)
                                             (-> (reverse-comp bs)
                                                 (sub-bioseq ( - f 4))))]
                                 (vec (map #(ala/codon->aa % table)
                                           (partition-all 3 (bs-seq v)))))))))

(defn six-frame-translation
  "Returns a lazy list of fastaSequence objects representing translations of
   a nucleotide biosequence object in six frames."
  ([nucleotide] (six-frame-translation nucleotide (ala/codon-tables 1)))
  ([nucleotide table]
     (map #(translate nucleotide % :table table)
          '(1 2 3 -1 -2 -3))))

(defn get-feature-sequence
  "Returns a fastaSequence object containing the sequence specified in a 
   genbankFeature object from a genbankSequence object. Designed for applying
   intervals to the sequence entry they originate from."
  [feature bs]
  (let [intervals (interval-seq feature)]
    (init-fasta-sequence
     (accession bs)
     (str (def-line bs) " - Feature: " (feature-type feature)
          " - [" (start (first intervals)) "-" (end (last intervals)) "]")
     (alphabet bs)
     (vec (mapcat #(if (comp? %)
                     (subvec (ala/revcom (bs-seq bs))
                             (- (end %) 1)
                             (start %))
                     (subvec (bs-seq bs)
                             (- (start %) 1)
                             (end %))) intervals)))))

(defn get-interval-sequence
  "Returns a fasta sequence of an interval."
  [interval bs]
  (let [dna (bs-seq bs)
        start (start interval)
        end (end interval)]
    (init-fasta-sequence
     (accession bs)
     (str (def-line bs) " [" start "-" end "]")
     (if (protein? bs) :iupacAminoAcids :iupacNucleicAcids)
     (if (false? (comp? interval))
       (subvec dna (- start 1) end)
       (subvec (ala/revcom dna) (- end 1) start)))))

;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;

(defn biosequence->file
  "Takes a collection of biosequences and prints them to file. To
  append to an existing file use `:append true` and the `:func`
  argument can be used to pass a function that will be used to prepare
  the printed output, the default is `fasta-string` which will print
  the biosequences to the file in fasta format. Returns the file."
  [bs file & {:keys [append func] :or {append true func fasta-string}}]
  (with-open [w (writer file :append append)]
    (dorun (map #(let [n (func %)]
                   (if n
                     (.write w (str n "\n"))))
                bs)))
  file)

;; sequence

(defn clean-sequence
  "Removes spaces and newlines and checks that all characters are
   legal characters for the supplied alphabet. Replaces non-valid
   characters with \\X. If `a' is not a defined alphabet throws an
   exception."
  [s a]
  (let [k (ala/get-alphabet a)
        w #{\space \newline}]
    (loop [l s a []]
      (if-not (seq l)
        a
        (let [c (first l)]
          (if (not (w c))
            (recur (rest l) (conj a (if (k c) c \X)))
            (recur (rest l) a)))))))

;; utilities

(defn protein-charge
  "Calculates the theoretical protein charge at the specified
  pH (default 7). Uses pKa values set out in the protein alphabets
  from cl-biosequence.alphabet. Considers Lys, His, Arg, Glu, Asp, Tyr
  and Cys residues only and ignores all other amino acids. The number
  of disulfides can be specified and 2 times this figure will be
  deducted from the number of Cys residues used in the calculation.
  Values used for the pKa of the N-term and C-term are 9.69 and 2.34
  respectively."
  [p & {:keys [ph disulfides] :or {ph 7 disulfides 0}}]
  (if-not (protein? p)
    (throw (Throwable. "Protein charge calculations can only be performed on protein objects.")))
  (let [a (ala/get-alphabet :signalpAminoAcids)
        ncalc (fn [x n]
                (* n
                   (/ (Math/pow 10 x) (+ (Math/pow 10 ph) (Math/pow 10 x)))))
        ccalc (fn [x n]
                (* n
                   (/ (Math/pow 10 ph) (+ (Math/pow 10 ph) (Math/pow 10 x)))))
        freq (frequencies (bs-seq p))
        cys (let [c (freq \C)
                  f (if c (- c (* 2 disulfides)) 0)]
              (if (>= f 0)
                f
                (throw (Throwable. "More disulfides specified than Cys residues."))))]
    (- (reduce + (cons (ncalc 9.69 1) (map (fn [[x n]] (ncalc (:pka (a x)) n))
                                           (select-keys freq [\K \H \R]))))
       (+ (ccalc (:pka (a \C)) cys)
          (reduce + (cons (ccalc 2.34 1) (map (fn [[x n]] (ccalc (:pka (a x)) n))
                                              (select-keys freq [\E \D \Y]))))))))

;; serialising

(defn- bs-freeze
  [this]
  (freeze this))

(defn- bs-thaw
  [this]
  (thaw this))

(defrecord indexWriter [strm]

  java.io.Closeable
  
  (close [this]
    (.close ^java.io.RandomAccessFile (:strm this))))

(defn init-index-writer
  [file]
  (->indexWriter (RandomAccessFile. (str (bs-path file) ".bin") "rw")))

(defn print-tagged
  [o w]
  (tag/pr-tagged-record-on o w))

(defprotocol indexFileIO
  (bs-writer [this]))

(defn write-and-position
  [obj writer acc]
  (let [o (nip/freeze obj)
        off (.getFilePointer (:strm writer))]
    (.write (:strm writer) o)
    (vector acc (list off (count o)))))

(defn index-biosequence-file
  [file]
  (let [ifile (index-file file)]
    (with-open [o (bs-writer ifile)]
      (let [i (assoc ifile :index
                (with-open [r (bs-reader file)]
                  (->> (biosequence-seq r)
                       (map #(write-and-position % o (accession %)))
                       (into {}))))]
        (spit (str (bs-path ifile) ".idx") (pr-str i))
        i))))

(defn merge-files
  [files outfile]
  (let [ifile (index-file (first files) outfile)]
    (with-open [o (bs-writer ifile)]
      (let [i (assoc ifile :index
                     (merge (mapcat (fn [x] (with-open [r (bs-reader x)]
                                             (->> (biosequence-seq r)
                                                  (map #(write-and-position % o (accession %)))
                                                  (into {}))))
                                    files)))]
           (spit (str (bs-path ifile) ".idx") (pr-str i))
           i))))

(defn read-one
  [off len file]
  (let [bb (byte-array len)]
    (with-open [r (RandomAccessFile. file "r")]
      (.seek r off)
      (.read r bb)
      (bs-thaw bb))))

(defn- my-record-tag-reader
  [tag val]
  (when-let [factory (and (map? val)
                          (tag/tag->factory tag))]
    (factory val)))

(def my-tagged-default-reader 
  (tag/some-tag-reader-fn my-record-tag-reader))

(defn load-indexed-file
  [path]
  (edn/read-string {:default my-tagged-default-reader
                    :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
                   (slurp (str path ".idx"))))

;; helper files

(load "mapping")
(load "fasta")

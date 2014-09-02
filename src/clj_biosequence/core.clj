(ns clj-biosequence.core
  (:require [clojure.java.io :refer [writer reader]]
            [fs.core :refer [file? absolute-path]]
            [clj-http.client :refer [post]]
            [clojure.string :refer [trim split upper-case]]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.indexing :as ind]
            [clojure.edn :as edn]
            [clojure.data.xml :as xml]
            [miner.tagged :as tag]
            [iota :as iot]
            [clojure.core.reducers :as r]))

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
    "Returns para
meters, if any.")
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
    "Returns a lazy list of features in a sequence.")
  (references [this]
    "Returns a collection of references in a sequence record."))

(defprotocol biosequenceFile
  (bs-path [this]
    "Returns the path of the file as string.")
  (index-file [this] [this ofile])
  (empty-instance [this path]))

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
                            (apply str (subvec (bs-seq bs) beg end))
                            (apply str (subvec (bs-seq bs) beg))))))

(defn reverse-comp
  "Returns a new fastaSequence with the reverse complement sequence."
  [this]
  (if (protein? this)
    (throw (IllegalArgumentException. "Can't reverse/complement a protein sequence."))
    (init-fasta-sequence (accession this)
                         (str (def-line this) " - Reverse-comp")
                         (alphabet this)
                         (apply str (ala/revcom (bs-seq this))))))

(defn reverse-seq
  "Returns a new fastaSequence with the reverse sequence."
  [this]
  (init-fasta-sequence (accession this)
                       (str (def-line this) " - Reversed")
                       (alphabet this) 
                       (apply str (reverse (bs-seq this)))))

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
                                 (apply str (map #(ala/codon->aa % table)
                                                 (partition-all 3 (bs-seq v)))))))))

(defn six-frame-translation
  "Returns a lazy list of fastaSequence objects representing translations of
   a nucleotide biosequence object in six frames."
  ([nucleotide] (six-frame-translation nucleotide (ala/codon-tables 1)))
  ([nucleotide table]
     (map #(translate nucleotide % :table table)
          '(1 2 3 -1 -2 -3))))

(defn n50
  "Takes anything that can have `biosequence-seq' called on it and
  returns the N50 of the sequences therein."
  [reader]
  (let [sa (sort > (pmap #(count (bs-seq %)) (biosequence-seq reader)))
        t (/ (reduce + sa) 2)
        n50 (atom 0)]
    (loop [l sa]
      (if (seq l)
        (do (swap! n50 + (first l))
            (if (>= @n50 t)
              (first l)
              (recur (rest l))))
        (throw (Throwable. "N50 calculation failed!"))))))

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
                     (apply str (subvec (ala/revcom (bs-seq bs))
                                        (- (end %) 1)
                                        (start %)))
                     (apply str (subvec (bs-seq bs)
                                        (- (start %) 1)
                                        (end %)))) intervals)))))

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
       (apply str (subvec dna (- start 1) end))
       (apply str (subvec (ala/revcom dna) (- end 1) start))))))

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
        w #{\space \newline}
        t (upper-case s)]
    (loop [l t a []]
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

(defn index-biosequence-file
  [file]
  (try
    (let [index (index-file file)
          i (with-open [w (ind/index-writer (bs-path index))]
              (assoc index :index
                     (with-open [r (bs-reader file)]
                       (ind/index-objects w (biosequence-seq r) accession))))]
      (ind/save-index (bs-path file) i)
      i)
    (catch Exception e
      (ind/delete-index (bs-path file))
      (println (str "Exception: " (.getMessage e))))))

(defn load-biosequence-index
  [path]
  (assoc (ind/load-indexed-file path) :path path))

(defn merge-biosequence-indexes
  [indexes outfile]
  (try
    (let [o (empty-instance (first indexes) outfile)
          i (with-open [w (ind/index-writer (bs-path o))]
              (assoc o :index
                     (apply merge
                            (map (fn [x] (with-open [r (bs-reader x)]
                                          (ind/index-objects w (biosequence-seq r) accession)))
                                 indexes))))]
      (ind/save-index (bs-path o) i)
      i)
    (catch Exception e
      (ind/delete-index outfile)
      (println (str "Exception: " (.getMessage e))))))


(defn delete-indexed-biosequence
  [index-file]
  (ind/delete-index (bs-path index-file)))

(defn get-object
  [reader key func]
  (let [i (get (:index reader) key)]
    (if-not (nil? i)
      (first (ind/object-seq (:strm reader) (list i) func)))))

(defn indexed-seq
  [index func]
  (ind/object-seq (:strm index) (vals (:index index)) func))

(defn close-index-reader
  [reader]
  (ind/close-index-reader (:strm reader)))

(defn open-index-reader
  [path]
  (ind/index-reader path))

(defn print-tagged-index
  [index w]
  (ind/print-tagged index w))

;; helper files

(load "mapping")
(load "fasta")

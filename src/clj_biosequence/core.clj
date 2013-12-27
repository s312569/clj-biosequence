(ns clj-biosequence.core
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-http.client :as client]
            [clojure.string :as string]
            [taoensso.nippy :refer [freeze thaw]]
            [clj-biosequence.alphabet :as ala]))

(declare init-fasta-store init-fasta-sequence translate)

(defprotocol biosequenceIO
  (bs-reader [this]
    "Returns a reader for a file containing biosequences. Use with
    `with-open'"))

(defprotocol biosequenceReader
  (biosequence-seq [this]
    "Returns a lazy sequence of biosequence objects.")
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
    "Returns the alphabet of a biosequence."))

(defprotocol biosequenceFile
  (bs-path [this]
    "Returns the path of the file as string."))

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn bioseq->string
  "Returns the sequence of a biosequence as a string."
  [bs]
  (apply str (bs-seq bs)))

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
     (map #(translate nucleotide % table)
          '(1 2 3 -1 -2 -3))))

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
  (with-open [w (io/writer file :append append)]
    (dorun (map #(let [n (func %)]
                   (if n
                     (.write w n)))
                bs)))
  file)

(defn if-string-int
  "If a string an integer is parsed, if not returns e. Will throw an
   expception if no integer can be parsed if `error?' is true. Used
   only when parsing optional fields in files where value could be nil
   or a string representation of an integer."
  ([e] (if-string-int e true))
  ([e error?]
     (try
       (if (string? e) (Integer/parseInt e) e)
       (catch NumberFormatException e
         (if error?
           (throw e))))))

;; sequence

(defn clean-sequence
  "Removes spaces and newlines and checks that all characters are
   legal characters for the supplied alphabet. Replaces non-valid
   characters with \\X. If `a' is not a defined alphabet throws an
   exception."
  [s a]
  (let [k (complement (ala/alphabet-chars a))
        w #{\space \newline}]
    (vec (remove nil? (map #(cond (k %) \X
                                  (w %) nil
                                  :else %) (vec s))))))

;; stats

(defn std-dev [samples]
  (let [n (count samples)
        mean (/ (reduce + samples) n)
        intermediate (map #(Math/pow (- %1 mean) 2) samples)]
    (Math/sqrt 
     (/ (reduce + intermediate) n))))

;; serialising

(defn bs-freeze
  [this]
  (freeze this))

(defn bs-thaw
  [this]
  (thaw this))

;; helper files

(load "mapping")
(load "fasta")

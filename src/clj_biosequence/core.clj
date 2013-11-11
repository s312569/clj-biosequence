(ns clj-biosequence.core
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-http.client :as client]
            [clojure.string :as string]
            [clj-biosequence.write :refer [print-tagged]]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.store :as st]))

(declare init-fasta-store init-fasta-sequence translate)

(defprotocol biosequenceIO
  (bs-reader [this] "Returns a reader for a file containing
  biosequences. Use with `with-open'"))

(defprotocol biosequenceReader
  (biosequence-seq [this]))

(defprotocol Biosequence
  (accession [this]
    "Returns the accession of a biosequence.")
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
  (save-rep [this]))

(defprotocol biosequenceFile
  (bs-path [this]
    "Returns the path of the file as string."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn bioseq->string
  "Returns the sequence of a biosequence as a string."
  [bs]
  (apply str (bs-seq bs)))

(defn residue-frequencies
  "Returns a map of residues in a biosequence and their frequency."
  [bs]
  (frequencies (bs-seq bs)))

(defn sub-bioseq
  "Returns a new fasta sequence object with the sequence corresponding
   to 'beg' (inclusive) and 'end' (exclusive) of 'bs'. If no 'end'
   argument returns from 'start' to the end of the sequence. Indexes
   start at zero."
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

(defn partition-bioseq
  "Partitions a sequence into a lazy list of lists of 'n' size.
   Default partition size is 3."
  ([bs] (partition-bioseq bs 3))
  ([bs n]
     (partition-all n (bs-seq bs))))

(defn reverse-comp [this]
  "Returns a new fastaSequence with the reverse complement sequence."
  (if (protein? this)
    (throw (Throwable. "Can't reverse/complement a protein sequence."))
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
  ([bs frame] (translate bs frame (ala/codon-tables 1)))
  ([bs frame table]
     (let [f ({1 1 2 2 3 3 4 4 5 5 6 6 -1 4 -2 5 -3 6} frame)]
      (cond (protein? bs)
            (throw (IllegalArgumentException. "Can't translate a protein sequence!"))
            (not (#{1 2 3 4 5 6} f))
            (throw (IllegalArgumentException. (str "Invalid frame: " frame)))
            :else
            (init-fasta-sequence (str (accession bs) "-" f)
                                 (str (def-line bs) " - Translated frame: " f)
                                 :iupacAminoAcids
                                 (let [v (cond (#{1 2 3} f)
                                               (sub-bioseq bs (- f 1))
                                               (#{4 5 6} f)
                                               (-> (reverse-comp bs)
                                                   (sub-bioseq ( - f 4))))]
                                   (vec (map #(ala/codon->aa % table)
                                             (partition-bioseq v 3)))))))))

(defn six-frame-translation
  "Returns a lazy list of fastaSequence objects representing translations of
   a nucleotide biosequence object in six frames."
  ([nucleotide] (six-frame-translation nucleotide (ala/codon-tables 1)))
  ([nucleotide table]
     (map #(translate nucleotide % table)
          '(1 2 3 -1 -2 -3))))

(defn fasta->file
  "Takes a collection of biosequences and prints them to file. To
  append to an existing file use :append true and the :func argument
  can be used to pass a function that will be used to prepare the
  printed output, the default is fasta-string which will print the
  biosequences to the file in fasta format."
  [bs file & {:keys [append func] :or {append true func fasta-string}}]
  (if (not append) (fs/delete (fs/absolute-path file)))
  (with-open [w (io/writer file)]
    (dorun (map #(let [n (func %)]
                   (if n
                     (.write w n)))
                bs)))
  file)

;; helper files

(load "mapping")
(load "utilities")
(load "fasta")

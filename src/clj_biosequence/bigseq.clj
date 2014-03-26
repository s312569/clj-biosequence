(ns clj-biosequence.bigseq
  (:require [clojure.string :refer [split]]
            [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.core :as bs])
  (:import [java.io RandomAccessFile]))

(declare read-til-defline lazy-helper get-lazy-seq)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; an implementation of biosequence for big fasta sequences
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; sequence IO

(defprotocol bigseqSeq
  (sequence-seq [this]))

(defrecord bigseqSeqReader [strm alphabet line]

  java.io.Closeable

  (close [this]
    (.close (:strm this)))

  bigseqSeq

  (sequence-seq [this]
    (map #(bs/clean-sequence % (:alphabet this))
         (take-while #(not (= \> (first %)))
                     (rest (drop-while #(not (= (:line this) %))
                                       (line-seq (:strm this))))))))

(defprotocol bigseqAccess
  (sequence-reader [this]))

;; biosequence

(defrecord bigSequence [line alphabet file]

  bs/Biosequence

  (accession [this]
    (:accession this))

  (accessions [this]
    (list (bs/accession this)))

  (def-line [this]
    (:description this))

  (protein? [this]
    (ala/alphabet-is-protein (:alphabet this)))

  (fasta-string [this])

  (alphabet [this]
    (:alphabet this))

  bigseqAccess

  (sequence-reader [this]
    (let [s (io/reader (:file this))]
      (->bigseqSeqReader s (:alphabet this) (:line this)))))

;; IO

(defrecord bigseqReader [strm alphabet file]

  bs/biosequenceReader

  (biosequence-seq [this]
    (lazy-helper (:strm this) (:alphabet this) (bs/bs-path this)))

  (parameters [this])

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this)))

  java.io.Closeable

  (close [this]
    (.close (:strm this))))

(defn init-bigseq-reader
  [strm alphabet file]
  (->bigseqReader strm alphabet file))

(defrecord bigseqFile [file alphabet]

  bs/biosequenceIO

  (bs-reader [this]
    (init-bigseq-reader (io/reader (:file this))
                        (:alphabet this)
                        (:file this)))

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this))))

(defn init-bigseq-file
  "Initialises fasta protein file. Accession numbers and description are 
   processed by splitting the string on the first space, the accession 
   being the first value and description the second."
  [path alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (if (fs/file? path)
      (->bigseqFile path alphabet)
      (throw (Throwable. (str "File not found: " path))))))

;; utilities

(defn- read-til-defline [strm]
  (loop [s strm]
    (if-let [l (.readLine s)]
      (if (= \> (first l))
        l
        (recur s)))))

(defn- get-lazy-seq [strm]
  (let [l (.read strm)
        c (char l)]
    (if (and (not (= -1 l)) (not (= \> c)))
      (cons c (lazy-seq (get-lazy-seq strm))))))

(defn- lazy-helper [strm a file]
  (let [l (read-til-defline strm)]
    (if l
      (cons (->bigSequence l a file)
            (lazy-seq (lazy-helper strm a file))))))


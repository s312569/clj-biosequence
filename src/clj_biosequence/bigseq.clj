(ns clj-biosequence.bigseq
  (:require [clojure.string :refer [split]]
            [clojure.java.io :refer [reader]]
            [fs.core :as fs]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.core :as bs]
            [iota :as iot]
            [clojure.core.reducers :as r]))

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

(defrecord bigSequence [accession description s e alphabet vec]

  bs/Biosequence

  (accession [this]
    (:accession this))

  (accessions [this]
    (list (bs/accession this)))

  (bs-seq [this]
    (pmap #(second (split (first (iot/subvec (:vec this) % (+ 1 %)))
                         #"\t"))
         (range (+ 1 (:s this)) (:e this))))

  (def-line [this]
    (:description this))

  (protein? [this]
    (ala/alphabet-is-protein (:alphabet this)))

  (fasta-string [this])

  (alphabet [this]
    (:alphabet this)))

(defn bs-reducer
  [file func]
  (->> (iot/seq (:file file))
       (r/filter #(not (= \> (first %))))
       (r/map #(bs/clean-sequence % :iupacAminoAcids))
       (r/map func)
       (r/fold h2)))

(defn h2
  ([] {})
  ([a b] (merge-with + a b)))

;; IO

(defrecord bigseqReader [vec index alphabet file]

  bs/biosequenceReader

  (biosequence-seq [this]
    (let [l (iot/numbered-vec (:file this))]
      (->> (filter #(re-find #">" %) l)
           (partition-all 2 1)
           (map #(->bigSequence (second (re-find #"^>([^\s]+)" (first %)))
                                (second (re-find #">[^\s]+\s+(.+)" (first %)))
                                (Integer/parseInt (first (split (first %) #"\t")))
                                (Integer/parseInt (first (split (second %) #"\t")))
                                (:alphabet this)
                                l)))))

  (parameters [this])

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this)))

  java.io.Closeable

  (close [this]
    (.close (:vec this))))

(defrecord bigseqFile [file alphabet]

  bs/biosequenceIO

  (bs-reader [this]
    (->bigseqReader (reader (:file this)) nil (:alphabet this) (:file this)))

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

(defn- lazy-helper [strm a file]
  (let [l (read-til-defline strm)]
    (if l
      (cons (->bigSequence l a file)
            (lazy-seq (lazy-helper strm a file))))))


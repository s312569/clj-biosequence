(ns clj-biosequence.fastq
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [fs.core :as fs])
  (:import (org.apache.commons.compress.compressors.gzip GzipCompressorInputStream)
           (org.apache.commons.compress.compressors.bzip2 BZip2CompressorInputStream)))

(declare init-indexed-fastq)

(defprotocol fastqSequenceData

  (qualities [this]))

(defrecord fastqSequence [description sequence quality]

  bs/Biosequence

  (accession [this]
    (:description this))

  (accessions [this]
    (list (bs/accession this)))

  (bs-seq [this]
    (vec (:sequence this)))

  (def-line [this]
    (bs/accession this))

  (protein? [this]
    false)

  (fasta-string [this]
    (str ">" (bs/accession this) "\n" (bs/bioseq->string this) "\n"))

  (alphabet [this]
    :iupacNucleicAcids)

  fastqSequenceData

  (qualities [this]
    (:quality this)))

(defn init-fastq-sequence
  [description sequence quality]
  (->fastqSequence description sequence quality))

(defn fastq->string
  [bs]
  (str (bs/accession bs) "\n" (:sequence bs) "\n"
       "+\n" (:quality bs) "\n"))

;; IO

(defrecord fastqReader [strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (map (fn [[d s d1 q]]
           (if (and (= \@ (first d))
                    (= \+ (first d1))
                    (= (count s) (count q)))
             (init-fastq-sequence d s q)
             (throw (Throwable.
                     (str "Data corrupted at: " d)))))
         (partition 4 (line-seq (:strm this)))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord fastqFile [file]

  bs/biosequenceIO
  
  (bs-reader [this]
    (condp = (fs/extension (bs/bs-path this))
      ".gz" (->fastqReader
             (-> (bs/bs-path this)
                 io/file io/input-stream GzipCompressorInputStream. io/reader))
      ".bz2" (->fastqReader
              (-> (bs/bs-path this)
                  io/file io/input-stream BZip2CompressorInputStream. io/reader))
      (->fastqReader (io/reader (bs/bs-path this)))))
  
  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this)))
  
  (index-file [this]
    (init-indexed-fastq (bs/bs-path this)))

  (index-file [this ofile]
    (init-indexed-fastq (fs/absolute-path ofile))))

(defrecord fastqString [str]

  bs/biosequenceIO

  (bs-reader [this]
    (->fastqReader (java.io.BufferedReader. (java.io.StringReader. (:str this))))))

(defn init-fastq-file
  [path]
  (if (fs/file? path)
    (->fastqFile path)
    (throw (Throwable. (str "File not found: " path)))))

(defn init-fastq-string
  [str]
  (->fastqString str))

;; functions

(defn shuffle-fq
  "Takes forward and reverse reads from paired end sequencing data and
  interleaves them in a single file."
  [^String forward ^String reverse ^String out]
  (with-open [^java.io.BufferedReader f (io/reader (bs/bs-path forward))
              ^java.io.BufferedReader r (io/reader (bs/bs-path reverse))
              ^java.io.BufferedWriter o (io/writer out)]
    (loop [fl (.readLine f)]
      (when fl
        (.write o ^String fl)
        (.write o ^String (.readLine f))
        (.write o ^String (.readLine f))
        (.write o ^String (.readLine f))
        (.write o ^String (.readLine r))
        (.write o ^String (.readLine r))
        (.write o ^String (.readLine r))
        (.write o ^String (.readLine r))
        (recur (.readLine f))))))

;; (defn shuffle-fq
;;   "Takes forward and reverse reads from paired end sequencing data and
;;   interleaves them in a single file."
;;   [forward reverse out]
;;   (with-open [^java.io.BufferedReader f (io/reader (bs/bs-path forward))
;;               ^java.io.BufferedReader r (io/reader (bs/bs-path reverse))
;;               ^java.io.BufferedWriter o (io/writer out)]
;;     (binding [*out* o]
;;       (dorun (map (fn [x y]
;;                     (map prn x)
;;                     (map prn y))
;;                   (partition-all 4 (line-seq f)) (partition-all 4 (line-seq r)))))))

;; (defn unshuffle-fq
;;   "Takes a fastq file with interleaved paired-end sequences and
;;   separates forward and reverse reads into separate files."
;;   [forward-out reverse-out in]
;;   {:pre [(fs/exists? (bs/bs-path in))]}
;;   (with-open [f (io/writer forward-out)
;;               r (io/writer reverse-out)
;;               i (bs/bs-reader in)]
;;     (dorun (map (fn [[x y]]
;;                   (.write f (fastq->string x))
;;                   (.write r (fastq->string y)))
;;                 (partition-all 2 (bs/biosequence-seq i))))))

;; indexing

(defrecord indexedFastqFile [index path]

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:path this)))

  bs/indexFileIO

  (bs-writer [this]
    (bs/init-index-writer this))

  bs/biosequenceReader

  (biosequence-seq [this]
    (map (fn [[o l]]
           (map->fastqSequence  (bs/read-one o l (str (bs/bs-path this) ".bin"))))
         (vals (:index this))))

  (get-biosequence [this accession]
    (let [[o l] (get (:index this) accession)]
      (if o
        (map->fastqSequence (bs/read-one o l (str (bs/bs-path this) ".bin")))))))

(defn init-indexed-fastq
  [file]
  (->indexedFastqFile {} file))

(defmethod print-method clj_biosequence.fastq.indexedFastqFile
  [this w]
  (bs/print-tagged this w))


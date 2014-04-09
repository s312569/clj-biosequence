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
    (let [ifile (init-indexed-fastq this)]
      (bs/index-entries this ifile))))

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

(defn init-indexed-fastq [fastqfile]
  (->indexedFastqFile {} (bs/bs-path fastqfile)))

(defmethod print-method clj_biosequence.fastq.indexedFastqFile
  [this w]
  (bs/print-tagged this w))


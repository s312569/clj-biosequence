(ns clj-biosequence.fastq
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-biosequence.core :as bios]
            [fs.core :as fs])
  (:import (org.apache.commons.compress.compressors.gzip GzipCompressorInputStream)
           (org.apache.commons.compress.compressors.bzip2 BZip2CompressorInputStream)))

(defprotocol fastqSequenceData

  (qualities [this]))

(defrecord fastqSequence [description sequence quality]

  bios/Biosequence

  (accession [this]
    (:description this))

  (accessions [this]
    (list (bios/accession this)))

  (bs-seq [this]
    (:sequence this))

  (def-line [this]
    (bios/accession this))

  (protein? [this]
    false)

  (fasta-string [this]
    (str ">" (bios/accession this) "\n" (bios/bioseq->string this) "\n"))

  (alphabet [this]
    :iupacNucleicAcids)

  (bs-save [this]
    (let [s (pr-str (dissoc this :_id))]
      (merge {:src s}
             (dissoc this :description :sequence :quality))))

  fastqSequenceData

  (qualities [this]
    (:quality this)))

(defmethod print-method clj_biosequence.core.fastaSequence
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

(defn check-fastq
  [fq]
  (and (= \@ (first (bios/accession fq)))
       (= (count (bios/bs-seq fq)) (count (qualities fq)))
       fq))

(defn init-fastq-sequence
  [description sequence quality]
  (check-fastq (->fastqSequence description sequence quality)))

(defn fastq->string
  [bs]
  (str (bios/accession bs) "\n" (bios/bioseq->string bs) "\n"
       "+\n" (qualities bs) "\n"))

;; IO

(defrecord fastqReader [strm]

  bios/biosequenceReader

  (biosequence-seq [this]
    (map (fn [[d s d1 q]]
           (if (and (= \@ (first d))
                    (= \+ (first d1))
                    (= (count s) (count q)))
             (init-fastq-sequence d s q)
             (throw (Throwable.
                     (str "Data corrupted at: " d)))))
         (partition-all 4 (line-seq (:strm this)))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord fastqFile [file]

  bios/biosequenceIO
  
  (bs-reader [this]
    (condp = (fs/extension (:file this))
      ".gz" (->fastqReader
             (-> (:file this) io/file io/input-stream GzipCompressorInputStream. io/reader))
      ".bz2" (->fastqReader
              (-> (:file this) io/file io/input-stream BZip2CompressorInputStream. io/reader))
      (->fastqReader (io/reader (:file this)))))
  
  bios/biosequenceFile

  (bs-path [this]
    (:file this)))

(defrecord fastqString [str]

  bios/biosequenceIO

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

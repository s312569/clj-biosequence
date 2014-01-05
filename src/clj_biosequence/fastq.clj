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
    (vec (:sequence this)))

  (def-line [this]
    (bios/accession this))

  (protein? [this]
    false)

  (fasta-string [this]
    (str ">" (bios/accession this) "\n" (bios/bioseq->string this) "\n"))

  (alphabet [this]
    :iupacNucleicAcids)

  fastqSequenceData

  (qualities [this]
    (:quality this)))

(defn init-fastq-sequence
  [description sequence quality]
  (check-fastq (->fastqSequence description sequence quality)))

(defn fastq->string
  [bs]
  (str (bios/accession bs) "\n" (:sequence bs) "\n"
       "+\n" (:quality bs) "\n"))

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
         (partition 4 (line-seq (:strm this)))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord fastqFile [file]

  bios/biosequenceIO
  
  (bs-reader [this]
    (println (fs/extension (bios/bs-path this)))
    (condp = (fs/extension (bios/bs-path this))
      ".gz" (->fastqReader
             (-> (bios/bs-path this)
                 io/file io/input-stream GzipCompressorInputStream. io/reader))
      ".bz2" (->fastqReader
              (-> (bios/bs-path this)
                  io/file io/input-stream BZip2CompressorInputStream. io/reader))
      (->fastqReader (io/reader (bios/bs-path this)))))
  
  bios/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this))))

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

(ns clj-biosequence.fastq
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [fs.core :as fs])
  (:import
   (org.apache.commons.compress.compressors.gzip GzipCompressorInputStream)
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
  ;; ugly but fastest on very large files
  "Takes forward and reverse reads from paired end sequencing data and
  interleaves them in a single file."
  [^fastqFile forward ^fastqFile reverse ^String out]
  {:pre [(fs/exists? (bs/bs-path forward)) (fs/exists? (bs/bs-path reverse))]}
  (with-open [^java.io.BufferedReader f (io/reader (bs/bs-path forward))
              ^java.io.BufferedReader r (io/reader (bs/bs-path reverse))
              ^java.io.BufferedWriter o (io/writer out)]
    (loop [fl (.readLine f)]
      (when fl
        (.write o fl)
        (.write o "\n")
        (.write o (.readLine f))
        (.write o "\n")
        (.write o (.readLine f))
        (.write o "\n")
        (.write o (.readLine f))
        (.write o "\n")
        (.write o (.readLine r))
        (.write o "\n")
        (.write o (.readLine r))
        (.write o "\n")
        (.write o (.readLine r))
        (.write o "\n")
        (.write o (.readLine r))
        (.write o "\n")
        (recur (.readLine f))))))

(defn unshuffle-fq
  "Takes a fastq file with interleaved paired-end sequences and
  separates forward and reverse reads into separate files."
  [^String forward-out ^String reverse-out ^fastqFile in]
  {:pre [(fs/exists? (bs/bs-path in))]}
  (with-open [^java.io.BufferedWriter f (io/writer forward-out)
              ^java.io.BufferedWriter r (io/writer reverse-out)
              ^java.io.BufferedReader i (io/reader (bs/bs-path in))]
    (loop [fl (.readLine i)]
      (when fl
        (.write f fl)
        (.write f "\n")
        (.write f (.readLine i))
        (.write f "\n")
        (.write f (.readLine i))
        (.write f "\n")
        (.write f (.readLine i))
        (.write f "\n")
        (.write r (.readLine i))
        (.write r "\n")
        (.write r (.readLine i))
        (.write r "\n")
        (.write r (.readLine i))
        (.write r "\n")
        (.write r (.readLine i))
        (.write r "\n")
        (recur (.readLine i))))))

;; indexing

(defrecord indexedFastqReader [index strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (bs/indexed-seq this map->fastqSequence))

  (get-biosequence [this accession]
    (bs/get-object this accession map->fastqSequence))

  java.io.Closeable

  (close [this]
    (bs/close-index-reader this)))

(defrecord indexedFastqFile [index path]

  bs/biosequenceIO

  (bs-reader [this]
    (->indexedFastqReader (:index this)
                          (bs/open-index-reader (:path this))))

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:path this)))

  (empty-instance [this path]
    (init-indexed-fastq path)))

(defn init-indexed-fastq
  [file]
  (->indexedFastqFile {} file))

(defmethod print-method clj_biosequence.fastq.indexedFastqFile
  [this w]
  (bs/print-tagged-index this w))


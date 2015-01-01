(ns clj-biosequence.fastq
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [fs.core :as fs]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; sequence
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastqSequence [description sequence quality])

(extend fastqSequence
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession (fn [this] (:description this))
    :accessions (fn [this] (list (bs/accession this))))
  bs/biosequenceDescription
  (assoc bs/default-biosequence-description
    :description (fn [this] (bs/accession this)))
  bs/Biosequence
  (assoc bs/default-biosequence-biosequence
    :bs-seq (fn [this] (vec (:sequence this)))
    :protein? (fn [this] false)
    :alphabet (fn [this] :iupacNucleicAcids)))

(defn qualities
  "Returns the quality string from a fastq sequence record."
  [this]
  (:quality this))

(defn init-fastq-sequence
  "Returns a fastq sequence record."
  [description sequence quality]
  (->fastqSequence description sequence quality))

(defn fastq->string
  "Returns a fastq sequence string."
  [bs]
  (str (bs/accession bs) "\n" (:sequence bs) "\n"
       "+\n" (:quality bs) "\n"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; IO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- parse-fastq
  [r]
  (map (fn [[d s d1 q]]
         (if (and (= \@ (first d))
                  (= \+ (first d1))
                  (= (count s) (count q)))
           (init-fastq-sequence d s q)
           (throw (Throwable.
                   (str "Data corrupted at: " d)))))
       (partition 4 (line-seq (:strm r)))))

(defrecord fastqReader [strm]
  bs/biosequenceReader
  (biosequence-seq [this] (parse-fastq this))
  java.io.Closeable
  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastqFile [file opts])

(extend fastqFile
  bs/biosequenceIO
  {:bs-reader (fn [this]
                (->fastqReader
                 (apply bs/bioreader
                        (bs/bs-path this) (:opts this))))}
  bs/biosequenceFile
  bs/default-biosequence-file)

(defn init-fastq-file
  "Returns a fastqFile record."
  [^String path & opts]
  {:pre [(fs/file? path)]}
  (->fastqFile path opts))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; string
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastqString [str]
  bs/biosequenceIO
  (bs-reader [this]
    (->fastqReader (java.io.BufferedReader.
                    (java.io.StringReader. (:str this))))))

(defn init-fastq-string
  "Returns a fastqString record."
  [str]
  (->fastqString str))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn shuffle-fq
  ;; ugly but fastest on very large files
  "Takes forward and reverse reads from paired end sequencing data and
  interleaves them in a single file."
  [^fastqFile forward ^fastqFile reverse ^String out]
  {:pre [(fs/exists? (bs/bs-path forward))
         (fs/exists? (bs/bs-path reverse))]}
  (with-open [^java.io.BufferedReader f
              (io/reader (bs/bs-path forward))
              ^java.io.BufferedReader r
              (io/reader (bs/bs-path reverse))
              ^java.io.BufferedWriter o
              (io/writer out)]
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
  (with-open [^java.io.BufferedWriter f
              (io/writer forward-out)
              ^java.io.BufferedWriter r
              (io/writer reverse-out)
              ^java.io.BufferedReader i
              (io/reader (bs/bs-path in))]
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

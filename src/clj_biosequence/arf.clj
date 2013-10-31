(ns clj-biosequence.arf
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-biosequence.core :as bios]
            [fs.core :as fs]))

;; entry

(defrecord arfSequence [accession length start-read end-read seq-read genome-id length-genome start-genome end-genome seq-genome strand mismatches match-string]

  bios/Biosequence

  (accession [this]
    (:accession this))

  (accessions [this]
    (list (:accession this)))

  (bs-seq [this]
    (:seq-read this))

  (fasta-string [this]
    (str ">" (:accession this) " " (:genome-id this) " " (:start-genome this) " " (:end-genome this) \newline
         (:seq-read this) \newline))

  (protein? [this]
    false)

  (alphabet [this]
    :iupacNucleicAcids)

  (bs-save [this]
    (let [s (pr-str (dissoc this :_id))]
      (assoc {:src s}
        (:_id this)))))

(defmethod print-method clj_biosequence.arf.arfSequence
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

(defn init-arf-sequence [lst]
  (let [[accession length start-read end-read seq-read genome-id length-genome start-genome end-genome seq-genome strand mismatches match-string]
        lst]
    (->arfSequence accession
                   (bios/if-string-int length)
                   (bios/if-string-int start-read)
                   (bios/if-string-int end-read)
                   seq-read
                   genome-id
                   (bios/if-string-int length-genome)
                   (bios/if-string-int start-genome)
                   (bios/if-string-int end-genome)
                   seq-genome
                   (if (#{"+" "-"} strand) strand
                       (throw (Throwable. (str  "Disallowed strand value: " strand))))
                   (bios/if-string-int mismatches)
                   match-string)))

;; IO

(defrecord arfReader [strm]

  bios/biosequenceReader

  (biosequence-seq [this]
    (map #(init-arf-sequence (string/split % #"\t"))
         (line-seq strm)))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord arfFile [file]

  bios/biosequenceIO

  (bs-reader [this]
    (->arfReader (io/reader (:file this))))

  bios/biosequenceFile

  (bs-path [this]
    (:file this)))

(defrecord arfString [str]

  bios/biosequenceIO

  (bs-reader [this]
    (->arfReader (java.io.BufferedReader. (java.io.StringReader. (:str this))))))

(defn init-arf-file [path]
  [path]
  (if (fs/file? path)
    (->arfFile path)
    (throw (Throwable. (str "File not found: " path)))))

(defn init-arf-string [str]
  [str]
  (->arfString str))

;; functions

(defn arf->bed [arfentry]
  (str (:genome-id arfentry) "\t" (:start-genome arfentry) "\t" (:end-genome arfentry) "\n"))

;; Format details:
;; read identifier
;; length of read sequence
;; start position in read sequence that is mapped
;; end position in read sequence that is mapped
;; read sequence
;; identifier of the genome-part to which a read is mapped to. This is either a scaffold id or a chromosome name
;; length of the genome sequence a read is mapped to
;; start position in the genome where a read is mapped to
;; end position in the genome where a read is mapped to
;; genome sequence to which a read is mapped
;; genome strand information. Plus means the read is aligned to the sense-strand of the genome. Minus means it is aligned to the antisense-strand of the genome.
;; Number of mismatches in the read mapping
;; Edit string that indicates matches by lowercase 'm' and mismatches by uppercase 'M'

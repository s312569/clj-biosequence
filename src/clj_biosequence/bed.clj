(ns clj-biosequence.bed
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [clj-biosequence.store :as st]
            [fs.core :as fs]))

;; entry

(defrecord bedSequence [chrom start end name score strand thickstart thickend itemrgb blockcount blocksizes blockstarts]

  bs/Biosequence

  (accession [this]
    (:name this))

  st/mongoBSRecordIO

  (mongo-bs-save [this pname cname]
    (let [s (hash-map :acc (bs/accession this) :element "sequence"
                      :pname pname :cname cname
                      :type "biosequence/bed"
                      :src (bs/bs-freeze this))]
      (if (:_id this)
        (assoc s :_id (:_id this))
        s))))

(defn init-bed-sequence
  [chrom start end name score strand thickstart thickend itemrgb blockcount blocksizes blockstarts]
  (->bedSequence chrom
                 (bs/if-string-int start false)
                 (bs/if-string-int end false)
                 name score
                 (if (#{"+" "-"} strand) strand
                     (throw (Throwable. (str  "Disallowed strand value: " strand))))
                 (bs/if-string-int thickstart false)
                 (bs/if-string-int thickend false)
                 itemrgb
                 (bs/if-string-int blockcount false)
                 blocksizes
                 blockstarts))

;; IO

(defrecord bedReader [strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (map #(apply init-bed-sequence (first (partition 12 12 (repeat "") %)))
         (filter (fn [x] (not (#{"browser" "track"} (first x))))
                 (map #(string/split % #"\s+")
                      (line-seq (:strm this))))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord bedFile [file]

  bs/biosequenceIO

  (bs-reader [this]
    (->bedReader (io/reader (:file this))))

  bs/biosequenceFile

  (bs-path [this]
    (:file this)))

(defrecord bedString [str]

  bs/biosequenceIO

  (bs-reader [this]
    (->bedReader (java.io.BufferedReader. (java.io.StringReader. (:str this))))))

(defn init-bed-file [path]
  [path]
  (if (fs/file? path)
    (->bedFile path)
    (throw (Throwable. (str "File not found: " path)))))

(defn init-bed-string [str]
  [str]
  (->bedString str))


;; Format details:

;; chrom - The name of the chromosome on which the genome feature exists.
;; Any string can be used. For example, “chr1”, “III”, “myChrom”, “contig1112.23”.
;; This column is required.

;; start - The zero-based starting position of the feature in the chromosome.
;; The first base in a chromosome is numbered 0.
;; The start position in each BED feature is therefore interpreted to be 1 greater than the start position listed in the feature. For example, start=9, end=20 is interpreted to span bases 10 through 20,inclusive.
;; This column is required.

;; end - The one-based ending position of the feature in the chromosome.
;; The end position in each BED feature is one-based. See example above.
;; This column is required.

;; name - Defines the name of the BED feature.
;; Any string can be used. For example, “LINE”, “Exon3”, “HWIEAS_0001:3:1:0:266#0/1”, or “my_Feature”.
;; This column is optional.

;; score - The UCSC definition requires that a BED score range from 0 to 1000, inclusive. However, bedtools allows any string to be stored in this field in order to allow greater flexibility in annotation features. For example, strings allow scientific notation for p-values, mean enrichment values, etc. It should be noted that this flexibility could prevent such annotations from being correctly displayed on the UCSC browser.
;; Any string can be used. For example, 7.31E-05 (p-value), 0.33456 (mean enrichment value), “up”, “down”, etc.
;; This column is optional.

;; strand - Defines the strand - either ‘+’ or ‘-‘.
;; This column is optional.

;; thickStart - The starting position at which the feature is drawn thickly.
;; Allowed yet ignored by bedtools.

;; thickEnd - The ending position at which the feature is drawn thickly.
;; Allowed yet ignored by bedtools.

;; itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0).
;; Allowed yet ignored by bedtools.

;; blockCount - The number of blocks (exons) in the BED line.
;; Allowed yet ignored by bedtools.

;; blockSizes - A comma-separated list of the block sizes.
;; Allowed yet ignored by bedtools.

;; blockStarts - A comma-separated list of block starts.
;; Allowed yet ignored by bedtools.

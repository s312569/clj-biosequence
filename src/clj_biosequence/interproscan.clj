(ns clj-biosequence.interproscan
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :as io]
            [clj-commons-exec :as exec]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clj-biosequence.core :as bios]
            [fs.core :as fs]))

(declare run-ips)

;; ips entry

(defrecord interproscanEntry [src])

(defmethod print-method clj_biosequence.interproscan.interproscanEntry
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

(defn init-ips-entry
  [src]
  (->interproscanEntry src))

(defn ips-entry-type
  [e]
  (zf/xml1-> (zip/xml-zip (:src e))
             (zf/attr "type")))

;; ips go entry

(defrecord interproscanGO [src])

(defmethod print-method clj_biosequence.interproscan.interproscanGO
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

(defn ips-entry-go
  [e]
  (map #(->interproscanGO (zip/node %))
       (zf/xml-> (zip/xml-zip (:src e))
                 :interpro
                 :classification
                 (zf/attr= :class_type "GO"))))

;; ips protein

(defrecord interproscanProtein [src]

  bios/Biosequence

  (accession [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               (zf/attr :id))))

(defmethod print-method clj_biosequence.interproscan.interproscanProtein
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

(defn init-ips-protein
  [src]
  (->interproscanProtein src))

;; ips search

(defrecord interproscanReader [strm]

  bios/biosequenceReader

  (biosequence-seq [this]
    (->> (:content (xml/parse (:strm this)))
         (filter #(= :protein (:tag %)))
         (map #(->interproscanProtein %))))

  java.io.Closeable

  (close [this]
    (.close (:strm this))))

(defrecord interproscanSearch [file]

  bios/biosequenceIO

  (bs-reader [this]
    (->interproscanReader (io/reader (:file this)))))

(defn init-ips-search
  [file]
  (->interproscanSearch file))

;; ips run

(defn ips
  [bs & {:keys [outfile] :or {outfile (fs/temp-file "ips")}}]
  (let [i (bios/fasta->file bs (fs/temp-file "seq-") :append false)]
    (try
      (run-ips (fs/absolute-path i) (fs/absolute-path outfile))
      (finally (fs/delete i)))))

;; utilities

(defn- run-ips
  [in out]
  (let [ips @(exec/sh ["iprscan" "-cli" "-i" in "-o" out "-appl"
                       "hmmpfam" "-iprlookup" "-goterms" "-seqtype" "p"])]
    (if (= 0 (:exit ips))
      (->interproscanSearch out)
      (if (:err ips)
        (throw (Throwable. (str "Interproscan error: " (:err ips))))
        (throw (Throwable. (str "Exception: " (:exception ips))))))))

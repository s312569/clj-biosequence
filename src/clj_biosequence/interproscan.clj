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
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord interproscanSearch [file]

  bios/biosequenceIO

  (bs-reader [this]
    (->interproscanReader (io/reader (:file this)))))

(defn init-ips-search
  [file]
  (->interproscanSearch file))

;; ips run

(defn ips
  [bs seqtype & {:keys [outfile appl lookup goterms] :or {outfile (fs/temp-file "ips")
                                                          appl '("hmmpfam")
                                                          lookup true
                                                          goterms true}}]
  (if (#{"p" "n"} seqtype)
    (let [i (bios/fasta->file bs (fs/temp-file "seq-") :append false
                              :func (fn [x] (if (not (> (count (bios/bs-seq x)) 10000))
                                             (bios/fasta-string x))))]
      (try
        (run-ips (fs/absolute-path i) (fs/absolute-path outfile)
                 seqtype appl lookup goterms)
        (finally (fs/delete i))))
    (throw (Throwable. (str seqtype " Not supported. Seqtype can be 'p' or 'n' only")))))

;; utilities

(defn ips-command
  [i o s a l g]
  (vec
   (remove nil? (-> (list
                     "iprscan" "-cli" "-i" i "-o" o "-seqtype" s
                     (if (or l g) "-iprlookup")
                     (if g "-goterms"))
                    (concat (interleave (iterate identity "-appl ") a))))))

(defn- run-ips
  [i o s a l g]
  (try
    (let [ips @(exec/sh (ips-command i o s a l g))]
      (if (= 0 (:exit ips))
        (->interproscanSearch o)
        (throw (Throwable. (str "Interproscan error: " (:err ips))))))
    (catch Exception e (str "Exception: " (:exception ips)))
    (finally (fs/delete o))))

(ns clj-biosequence.interproscan
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :as io]
            [clj-commons-exec :as exec]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [fs.core :as fs]))

(declare run-ips)

;; ips go entry

(defrecord interproscanGo [src])

(defn go-component
  [go]
  (zf/xml1-> (zip/xml-zip (:src go))
             :category
             zf/text))

(defn go-accession [go]
  (zf/xml1-> (zip/xml-zip (:src go))
             (zf/attr :id)))

(defn go-description [go]
  (zf/xml1-> (zip/xml-zip (:src go))
             :description
             zf/text))

;; ips entry

(defrecord interproscanEntry [src])

(defn init-ips-entry
  [src]
  (->interproscanEntry src))

(defn ips-entry-type
  [e]
  (zf/xml1-> (zip/xml-zip (:src e))
             (zf/attr :type)))

(defn ips-entry-name
  [e]
  (zf/xml1-> (zip/xml-zip (:src e))
             (zf/attr :name)))

(defn ips-entry-accession
  [e]
  (zf/xml1-> (zip/xml-zip (:src e))
             (zf/attr :id)))

(defn ips-entry-start
  [e]
  (zf/xml1-> (zip/xml-zip (:src e))
             :match
             :location
             (zf/attr :start)))

(defn ips-entry-end
  [e]
  (zf/xml1-> (zip/xml-zip (:src e))
             :match
             :location
             (zf/attr :end)))

(defn ips-go-seq
  [e]
  (map #(->interproscanGo (zip/node %))
       (zf/xml-> (zip/xml-zip (:src e))
                 :classification
                 (zf/attr= :class_type "GO"))))

;; ips protein

(defrecord interproscanProtein [src]

  bs/Biosequence

  (accession [this]
    (apply str
           (interpose "_"
                      (drop-last 2
                                 (string/split (zf/xml1-> (zip/xml-zip (:src this))
                                                          (zf/attr :id))
                                               #"_"))))))

(defn init-ips-protein
  [src]
  (->interproscanProtein src))

(defn ips-entry-seq
  [ip]
  (map #(init-ips-entry (zip/node %)) (zf/xml-> (zip/xml-zip (:src ip)) :interpro)))

(defn prot-go-terms
  [ip]
  (flatten (pmap #(ips-go-seq %)
                 (ips-entry-seq ip))))

;; ips search

(defrecord interproscanReader [strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (->> (:content (xml/parse (:strm this)))
         (filter #(= :protein (:tag %)))
         (map #(->interproscanProtein %))))

  (parameters [this]
    ())

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord interproscanFile [file]

  bs/biosequenceIO

  (bs-reader [this]
    (->interproscanReader (io/reader (:file this))))

  bs/biosequenceFile

  (bs-path [this]
    (:file this)))

(defn init-ips-result-file
  [file]
  (->interproscanFile (fs/absolute-path file)))

;; ips run

(defn ips
  [bs seqtype & {:keys [outfile appl lookup goterms] :or {outfile (fs/temp-file "ips")
                                                          appl '("hmmpfam")
                                                          lookup true
                                                          goterms true}}]
  (if (#{"p" "n"} seqtype)
    (let [i (bs/biosequence->file bs (fs/temp-file "seq-") :append false
                              :func (fn [x] (if (not (> (count (bs/bs-seq x)) 10000))
                                             (bs/fasta-string x))))]
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
        (->interproscanFile o)
        (throw (Throwable. (str "Interproscan error: " (:err ips))))))
    (catch Exception e
      (str "Exception: " (:exception ips))
      (fs/delete o))))

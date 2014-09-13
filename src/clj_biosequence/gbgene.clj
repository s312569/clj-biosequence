(ns clj-biosequence.gbgene
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :refer [split]]
            [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.eutilities :as eu]
            [clj-biosequence.core :as bs]
            [clj-biosequence.citation :as ci]))

;; gene

(defrecord gbGene [src]

  bs/Biosequence

  (accession [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :Entrezgene_track-info
               :Gene-track
               :Gene-track_geneid
               zf/text))

  (def-line [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :Entrezgene_gene
               :Gene-ref
               :Gene-ref_desc
               zf/text))

  (protein?
    [this]
    false))

(defn locus-name
  [gbgene]
  (zf/xml1-> (zip/xml-zip (:src gbgene))
               :Entrezgene_gene
               :Gene-ref
               :Gene-ref_locus
               zf/text))

;; reader

(defrecord gbgeneReader [strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (let [xml (xml/parse (:strm this) :support-dtd false)]
      (map (fn [x]
             (vec x) ;; realising all the laziness
             (->gbGene x))
           (filter #(= (:tag %) :Entrezgene)
                   (:content xml)))))

  java.io.Closeable

  (close [this]
    (.close (:strm this))))

(defn init-gene-reader
  [strm]
  (->gbgeneReader strm))

;; file

(defrecord gbgeneFile [file]
  
  bs/biosequenceIO

  (bs-reader [this]
    (init-gene-reader (io/reader (:file this))))
  
  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this)))

  ;; (index-file [this]
  ;;   (init-indexed-gbgene (bs/bs-path this)))

  ;; (index-file [this ofile]
  ;;   (init-indexed-gbgene (fs/absolute-path ofile)))

  )

(defn init-gbgene-file
  [file]
  (->gbgeneFile file))

;; string

(defrecord gbgeneString [str]

  bs/biosequenceIO

  (bs-reader [this]
    (init-gene-reader (io/reader (:str this)))))

;; connection

(defrecord gbgeneConnection [acc-list]

  bs/biosequenceIO

  (bs-reader [this]
    (init-gene-reader (io/reader
                       (eu/e-fetch (:acc-list this) "gene" nil "xml")))))

(defn init-gbgene-connection
  [acc-list]
  (->gbgeneConnection acc-list))

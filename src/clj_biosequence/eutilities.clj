(ns clj-biosequence.eutilities
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :refer [split]]
            [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.core :as bs]))

(defn- search-helper
  ([term db retstart] (search-helper term db retstart nil))
  ([term db retstart key]
     (xml/parse-str
      (:body
       (bs/get-req
        (str "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db="
             (java.net.URLEncoder/encode (name db))
             "&term="
             (java.net.URLEncoder/encode term)
             "&retmax=" 1000
             "&retstart=" retstart
             (if key (str "&WebEnv=" key))
             "&usehistory=y"))))))

(defn e-search
  "Takes a term and a database and returns a list of accession numbers
  matching the search term."
  ([term db] (e-search term db 0 nil))
  ([term db restart key]
     (let [r (search-helper term db restart key)
           k (zf/xml1-> (zip/xml-zip r) :WebEnv zf/text)
           c (Integer/parseInt (zf/xml1-> (zip/xml-zip r) :Count zf/text))]
       (if (> restart c)
         nil
         (lazy-cat (zf/xml-> (zip/xml-zip r) :IdList :Id zf/text)
                   (e-search term db (+ restart 1000) k))))))

(defn e-fetch
  "Retrieves a stream of database entries from GenBank corresponding
  to the list of accession numbers. Sequences returned in the format
  specified by the combination of rettype and retmode as set out in
  the eutilities documentation."
  [a-list db rettype retmode]
  (if (empty? a-list)
    nil
    (let [r (bs/post-req "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                         {:query-params
                          (merge {:db (name db)
                                  :id (apply str
                                             (interpose "," a-list))}
                                 (if rettype {:rettype rettype} {})
                                 (if retmode {:retmode retmode} {}))
                          :as :stream})]
      (:body r))))

(ns clj-biosequence.eutilities
  (:require [clojure.data.xml :refer [parse-str]]
            [clojure.data.zip.xml :refer [xml1-> xml-> text]]
            [clojure.zip :refer [xml-zip]]
            [clojure.string :refer [split]]
            [clj-biosequence.core :refer [get-req post-req]]))

(defn- search-helper
  ([term db retstart] (search-helper term db retstart nil))
  ([term db retstart key]
     (parse-str
      (:body
       (get-req
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
           k (xml1-> (xml-zip r) :WebEnv text)
           c (Integer/parseInt (xml1-> (xml-zip r) :Count text))]
       (if (> restart c)
         nil
         (lazy-cat (xml-> (xml-zip r) :IdList :Id text)
                   (e-search term db (+ restart 1000) k))))))

(defn e-fetch
  "Retrieves a stream of database entries from GenBank corresponding
  to the list of accession numbers. Sequences returned in the format
  specified by the combination of rettype and retmode as set out in
  the eutilities documentation."
  [a-list db rettype retmode]
  (if (empty? a-list)
    nil
    (let [r (post-req "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                         {:query-params
                          (merge {:db (name db)
                                  :id (apply str
                                             (interpose "," a-list))}
                                 (if rettype {:rettype rettype} {})
                                 (if retmode {:retmode retmode} {}))
                          :as :stream})]
      (:body r))))

(defn- e-link-helper
  [acc idb tdb]
  (let [r (post-req "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
                    {:query-params
                     (merge {:dbfrom (name idb)
                             :db (name tdb)
                             :id acc})})]
    (:body r)))

(defn e-link
  [acc idb tdb]
  (let [x (parse-str (e-link-helper acc idb tdb))
        p (xml-> (xml-zip x) :LinkSet :LinkSetDb)]
    (into {}
          (map #(vector (xml1-> % :LinkName text)
                        (xml-> % :Link :Id text))
               p))))

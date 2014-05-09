(ns clj-biosequence.interproscan
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :as io]
            [clj-commons-exec :as exec]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [fs.core :as fs]))

(declare run-ips init-indexed-ips)

;; ips protocol

(defprotocol ipsSignature
  (library [this] "Returns a hashmap with the keys :library and :version")
  (entry [this] "Returns the IPR entry from a signature."))

(defprotocol ipsMatch
  (score [this])
  (evalue [this])
  (signature [this]))

(defprotocol ipsEntry
  (entry-type [this])
  (entry-name [this])
  (ips-go-terms [this])
  (pathways [this]))

(defprotocol ipsGoterm
  (category [this])
  (go-name [this])
  (go-string [this])
  (bp? [this])
  (cc? [this])
  (mp? [this]))

(defprotocol ipsPathway
  (pathway-database [this])
  (pathway-name [this]))

;; pathway

(defrecord interproscanPathway [src]

  ipsPathway

  (pathway-database [this]
    (:db (:attrs (:src this))))

  (pathway-name [this]
    (:name (:attrs (:src this))))

  bs/Biosequence

  (accession [this]
    (:id (:attrs (:src this)))))

;; go terms

(defrecord interproscanGoTerm [src]

  ipsGoterm

  (category [this]
    (:category (:attrs (:src this))))

  (go-name [this]
    (:name (:attrs (:src this))))

  (go-string [this]
    (let [c (category this)
          nc (condp = c
               "CELLULAR_COMPONENT" "Cellular Component"
               "BIOLOGICAL_PROCESS" "Biological Process"
               "MOLECULAR_FUNCTION" "Molecular Function"
               (str "Unexpected value: " c))]
      (str (bs/accession this) "," nc "," (go-name this))))

  (cc? [this]
    (= "CELLULAR_COMPONENT" (category this)))

  (bp? [this]
    (= "BIOLOGICAL_PROCESS" (category this)))

  (mp? [this]
    (= "MOLECULAR_FUNCTION" (category this)))

  bs/Biosequence

  (accession [this]
    (:id (:attrs (:src this)))))

;; entry

(defrecord interproscanEntry [src]

  ipsEntry

  (entry-type [this]
    (:type (:attrs (:src this))))

  (entry-name [this]
    (:name (:attrs (:src this))))

  (ips-go-terms [this]
    (doall (map #(->interproscanGoTerm (zip/node %))
                (zf/xml-> (zip/xml-zip (:src this))
                          :go-xref))))

  (pathways [this]
    (doall (map #(->interproscanPathway (zip/node %))
                (zf/xml-> (zip/xml-zip (:src this))
                          :pathway-xref))))

  bs/Biosequence

  (accession [this]
    (:ac (:attrs (:src this))))

  (def-line [this]
    (:desc (:attrs (:src this)))))

;; signature

(defrecord interproscanSignature [src]

  ipsSignature

  (library [this]
    (let [l (zf/xml1-> (zip/xml-zip (:src this))
                       :signature-library-release)]
      {:library (zf/attr l :library) :version (zf/attr l :version)}))

  (entry [this]
    (let [r (zf/xml1-> (zip/xml-zip (:src this)) :entry)]
      (if (not (nil? r))
        (->interproscanEntry (zip/node r))
        (->interproscanEntry ()))))

  bs/Biosequence

  (accession [this]
    (:ac (:attrs (:src this))))

  (def-line [this]
    (:desc (:attrs (:src this)))))

;; matches

(defrecord interproscanHmmThree [src]

  ipsMatch

  (score [this]
    (Float/parseFloat (:score (:attrs (:src this)))))

  (evalue [this]
    (Float/parseFloat (:evalue (:attrs (:src this)))))

  (signature [this]
    (->interproscanSignature
     (zip/node (zf/xml1-> (zip/xml-zip (:src this)) :signature)))))

;; ips protein

(defrecord interproscanProtein [accession src]

  bs/Biosequence

  (accession [this]
    (:accession this)))

(defn hmmer-3-seq
  ""
  [protein]
  (map #(->interproscanHmmThree (zip/node %))
       (zf/xml-> (zip/xml-zip (:src protein))
                 :matches
                 :hmmer3-match)))

(defn get-gos
  [ips]
  (->> (hmmer-3-seq ips)
       (map signature)
       (map entry)
       (mapcat ips-go-terms)))

;; ips search

(defrecord interproscanReader [strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (let [el (->> (:content (xml/parse (:strm this)))
                  (filter (fn [x] (= :protein (:tag x)))))]
      (mapcat #(->> (zf/xml-> (zip/xml-zip %)
                              :xref (zf/attr :id))
                    (map (fn [x]
                           (->interproscanProtein x %))))
              el)))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord interproscanResult [file]

  bs/biosequenceIO

  (bs-reader [this]
    (->interproscanReader (io/reader (:file this))))

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this)))

  (index-file [this]
    (init-indexed-ips (bs/bs-path this)))

  (index-file [this ofile]
    (init-indexed-ips (fs/absolute-path ofile))))

(defn init-ips-result [file]
  (->interproscanResult (fs/absolute-path file)))

;; ips run

(defn ips
  "Runs interproscan on a list of biosequences."
  [bsl outfile & {:keys [appl lookup goterms precalc pathways]
                  :or {appl '("tigrfam")
                       lookup true
                       goterms true
                       precalc false
                       pathways false}}]
  (->interproscanResult
   (let [i (bs/biosequence->file bsl (fs/temp-file "ips-") :append false)]
     (try
       (run-ips :infile (fs/absolute-path i)
                :outfile outfile
                :appl appl
                :precalc precalc
                :pathways pathways
                :lookup lookup
                :goterms goterms)
       (finally (fs/delete i))))))

;; utilities

(defn- ips-command
  [& {:keys [infile outfile appl precalc pathways lookup goterms]}]
  (vec
   (remove nil? (-> (list
                     "interproscan.sh" "-i" infile "-o" outfile "-seqtype" "p" "-f" "XML"
                     "-appl" (apply str (interpose "," appl))
                     (if (not precalc) "-dp")
                     (if (not pathways) "-pa")
                     (if (or lookup goterms) "-iprlookup")
                     (if goterms "-goterms"))))))

(defn- run-ips
  [& {:keys [infile outfile appl precalc pathways lookup goterms]}]
  (try
    (let [ips @(exec/sh (ips-command :infile infile
                                     :outfile outfile
                                     :appl appl
                                     :precalc precalc
                                     :lookup lookup
                                     :goterms goterms
                                     :pathways pathways))]
      (if (= 0 (:exit ips))
        outfile
        (do (println ips)
            (throw (Throwable. (str "Interproscan error: " (:err ips)))))))
    (catch Exception e
      (throw e)
      (fs/delete outfile))))

;; indexing

(defrecord indexedInterproscanFile [index path]

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:path this)))

  bs/indexFileIO

  (bs-writer [this]
    (bs/init-index-writer this))

  bs/biosequenceReader

  (biosequence-seq [this]
    (map (fn [[o l]]
           (map->interproscanProtein (bs/read-one o l (str (bs/bs-path this) ".bin"))))
         (vals (:index this))))

  (get-biosequence [this accession]
    (let [[o l] (get (:index this) accession)]
      (if o
        (map->interproscanProtein (bs/read-one o l (str (bs/bs-path this) ".bin")))))))

(defn init-indexed-ips [file]
  (->indexedInterproscanFile {} file))

(defmethod print-method clj_biosequence.interproscan.indexedInterproscanFile
  [this w]
  (bs/print-tagged this w))

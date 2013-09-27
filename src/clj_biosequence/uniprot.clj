(ns clj-biosequence.uniprot
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :as io]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clj-biosequence.core :as bios]
            [clj-biosequence.persistence :as ps]
            [clj-http.client :as client]
            [fs.core :as fs]))

(declare prot-name meta-data amino-acids nomenclature uniprot-process-request get-uniprot-stream)

;; protein

(defrecord uniprotProtein [src]

  bios/Biosequence
  
  (accessions [uniprot]
    (zf/xml-> (zip/xml-zip (:src uniprot)) :accession zf/text))

  (accession
    [this]
    (first (bios/accessions this)))

  (def-line [this]
    (let [nom (zip/xml-zip (nomenclature this))]
      (str (zf/xml1-> nom
                      :recommendedName
                      :fullName
                      zf/text)
           " | "
           (zf/xml1-> nom
                      :alternativeName
                      :fullName
                      zf/text)
           " ["
           (zf/xml1-> (zip/xml-zip (organism this))
                      :name
                      (zf/attr= :type "scientific")
                      zf/text)
           "]")))

  (bs-seq [this]
    (vec (apply str
                (remove {\space\newline}
                        (zf/xml1-> (zip/xml-zip (amino-acids this))
                                   zf/text)))))

  (fasta-string [this]
    (let [db (condp = (:dataset (meta-data this))
               "Swiss-Prot" "sp"
               "TrEMBL" "tr"
               "bs")]
      (str ">" db "|" (bios/accession this) "|"
           (prot-name this)
           " " (bios/def-line this)
           \newline
           (apply str (bios/bs-seq this))
           \newline)))

  (protein? [this] true)

  (reverse-seq [this]
    (bios/init-fasta-sequence (bios/accession this)
                              (bios/def-line this)
                              :iupacAminoAcids
                              (reverse (bios/bs-seq this))))

  (reverse-comp [this]
    (throw (Throwable. "Action not defined for protein sequence.")))

  (alphabet [this]
    :iupacAminoAcids))

;; IO

(defrecord uniprotReader [strm]

  bios/biosequenceReader

  (bios/biosequence-seq [this]
    (let [xml (xml/parse (:strm this))]
      (map (fn [x]
             (->uniprotProtein x))
           (filter #(= (:tag %) :entry)
                   (:content xml)))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defn init-uniprot-reader
  [strm]
  (->uniprotReader strm))

(defrecord uniprotFile [file]

  bios/biosequenceIO

  (bs-reader [this]
    (init-uniprot-reader (io/reader (:file this)))))

(defrecord uniprotString [str]

  bios/biosequenceIO

  (bs-reader [this]
    (init-uniprot-reader (io/reader (:str this)))))

(defrecord uniprotConnection [acc-list retype email]

  bios/biosequenceIO

  (bs-reader [this]
    (let [s (get-uniprot-stream (:acc-list this) (:retype this) (:email this))]
      (condp = (:retype this)
        :xml (init-uniprot-reader (io/reader s))
        :fasta (bios/init-fasta-reader (io/reader s) :iupacAminoAcids)))))

(defn init-uniprotxml-file
  "Initialises a uniprotXmlFile object."
  [path]
  (if (fs/exists? path)
    (->uniprotFile path)
    (throw (Throwable. (str "File not found: " path)))))

(defn init-uniprot-string
  [str]
  (->uniprotString str))

(defn init-uniprot-connection
  [accessions retype email]
  (if (#{:xml :fasta} retype)
    (let [l (if (coll? accessions) accessions (list accessions))]
      (->uniprotConnection l retype email))
    (throw (Throwable. (str retype " not allowed. "
                            "Only :xml and :fasta are allowed retype values.")))))

;; persistence

(defrecord uniprotStore [file])

(defrecord uniprotStoreDir [dir]

  bios/biosequenceStoreDir

  (load-store [this dbfile]
    (->uniprotStore dbfile)))

(defn index-uniprot-file
  "Indexes a uniprotXmlFile and returns a uniprotStore object."
  [uniprotfile]
  (let [st (ps/init-store (->uniprotStore
                           (ps/index-file-name (:file uniprotfile))))]
    (bios/index-biosequence-file uniprotfile st)))

(defn load-uniprot-store
  "Loads a uniprotStore."
  [dir]
  (bios/load-biosequence-store (->uniprotStoreDir dir)))

;; web search

(defn wget-uniprot-search
  "Returns a non-lazy list of uniprot accession numbers satisfying the search term. 
   The search term uses the same syntax as the uniprot web interface. For 
   example:
   - to get all Schistosoma mansoni proteins in the proteome reference set term
     would be:
     'organism:6183 AND keyword:1185'
   - get all Schistosoma mansoni proteins in the proteome set that are intrinsic
     to the membrane:
     'taxonomy:6183 AND keyword:1185 AND go:0031224'
   - get all reviewed human entries:
     'reviewed:yes AND organism:9606'
   And so on. Returns an empty list if no matches found. Uniprot
   requires an email so an email can be supplied using the email argument."
  ([term email] (wget-uniprot-search term email 0))
  ([term email offset]
     (let [r (remove #(= % "")
                     (-> (client/get
                          (str "http://www.uniprot.org/uniprot/?query="
                               term
                               "&format=list"
                               (str "&offset=" offset)
                               "&limit=1000")
                          {:client-params {"http.useragent"
                                           (str "clj-http " email)}})
                         (:body)
                         (string/split #"\n")))]
       (if (empty? r)
         nil
         (lazy-cat r (wget-uniprot-search term email (+ offset 1000)))))))

;; uniprot convienence functions

(defn organism
  "Returns organism information from Uniprot biosequence as xml elements.
   Includes organism name and taxonomy information."
  [uniprot]
  (zip/node (zf/xml1-> (zip/xml-zip (:src uniprot))
                       :organism)))

(defn prot-name
  "Returns the name of a uniprot as a string."
  [uniprot]
  (zf/xml1-> (zip/xml-zip (:src uniprot)) :name zf/text))

(defn amino-acids
  "Returns sequence information as xml elements. Information includes mass, 
   checksum, modified, version and amino-acids."
  [uniprot]
  (zip/node (zf/xml1-> (zip/xml-zip (:src uniprot))
                       :sequence)))

(defn nomenclature
  "Returns protein naming information as xml elements. Includes information on
   recommended, submitted, alternative, allergen, biotech, cdantigen names."
  [uniprot]
  (zip/node (zf/xml1-> (zip/xml-zip (:src uniprot))
                       :protein)))

(defn gene
  "Returns gene information as a list of xml elements. Includes type and name 
   information."
  [uniprot]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot)) :gene)))

(defn gene-location
  "Returns gene location information as a list of xml elements."
  [uniprot]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot))
                          :geneLocation)))

(defn citations
  "Returns citation information as a list of xml elements. Contains information
   country, last (page), date, pubmed, institute, name (journal), first 
   (page), title, city, scope, type, consortium, number, authors, source,
   editors, publisher, volume and db."
  [uniprot]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot))
                          :reference)))

(defn comment-value
  "'subcellular location', 'alternative products', 'interactions', 'mass spectrometry', "
  [uniprot comment]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot))
                          :comment
                          (zf/attr= :type comment))))

(defn comments
  "Change this to list all comment keys."
  [uniprot]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot))
                          :comment)))

(defn db-references
  "Returns all db-references as a list of xml elements."
  [uniprot]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot))
                          :dbReference)))

(defn go-terms
  [uniprot]
  (zf/xml-> (zip/xml-zip (:src uniprot))
            :dbReference
            (zf/attr= :type "GO")
            :property
            (zf/attr= :type "term")
            (zf/attr :value)))

(defn existence
  "Protein existence evidence as a list of xml elements."
  [uniprot]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot))
                          :proteinExistence)))

(defn keywords
  "Keywords as a list of xml elements."
  [uniprot]
  (map zip/node
       (zf/xml-> (zip/xml-zip (:src uniprot))
                 :keyword)))

(defn features
  "Features as a list of xml elements."
  [uniprot]
  (map zip/node (zf/xml-> (zip/xml-zip (:src uniprot))
                          :feature)))

(defn meta-data
  "Returns a map with uniprot meta-data. Keys - :dataset, :created, :modified
   and :version. All values are strings except :version, which is an integer."
  [uniprot]
  (let [s (zip/xml-zip (:src uniprot))]
    {:dataset (zf/attr s :dataset)
     :created (zf/attr s :created)
     :modified (zf/attr s :modified)
     :version (if-let [f (zf/attr s :version)]
                (Integer. ^String f))}))

;; utilities

(defn- uniprot-process-request
  [address params file]
  (let [p (client/post address params)]
    (letfn [(check [a c]
              (let [r (client/get a {:follow-redirects false
                                     :as :stream})]
                (cond
                 (nil? (get (:headers r) "retry-after"))
                 (if (= (:status r) 200)
                   r
                   (throw (Throwable. (str "Error in sequence retrieval: code"
                                           (:status r)))))
                 (> c 50)
                 (throw (Throwable. "Too many tries."))
                 :else
                 (recur
                  (do (Thread/sleep 3000)
                      a)
                  (+ 1 c)))))]
      (if (some #(= (:status p) %) '(302 303))
        (do
          (if (get (:headers p) "retry-after")
            (Thread/sleep (read-string (get (:headers p) "retry-after"))))
          (check (get (:headers p) "location") 0))
        (throw (Throwable. (str "Error in sequence retrieval" p)))))))

(defn- get-uniprot-stream
  "Returns a GZIPInputStream from Uniprot with the results of a batch fetch 
   command for the sequences in a collection of accessions."
  [accessions class email]
  (if (empty? accessions)
    nil
    (let [f (let [file (fs/temp-file "up-seq-")]
              (doseq [s accessions]
                (spit file (str s "\n") :append true))
              file)]
      (try
        (:body
         (uniprot-process-request
          "http://www.uniprot.org/batch/"
          {:client-params {"http.useragent"
                           (str "clj-http " email)}
           :multipart [{:name "file" :content f}
                       {:name "format" :content (name class)}]
           :follow-redirects false}
          f))
        (catch Exception e
          (throw e))
        (finally (fs/delete f))))))

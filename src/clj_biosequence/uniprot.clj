(ns clj-biosequence.uniprot
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clojure.java.io :as io]
            [clojure.pprint :as pp]
            [clj-biosequence.core :as bios]
            [clojure.string :as st]
            [clj-http.client :as client]
            [fs.core :as fs]
            [clojure.java.jdbc :as sql]))

(declare prot-names prot-name process-feature process-sequence process-cites init-uniprot-store meta-data amino-acids nomenclature organism uniprot-process-request uniprot-sequence-helper read-up-xml-from-stream org-scientific-name)

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
      (str " | "
           (zf/xml1-> nom
                      :recommendedName
                      :fullName
                      zf/text)
           " | "
           (zf/xml1-> nom
                      :alternativeName
                      :fullName
                      zf/text)
           " ["
           (org-scientific-name this)
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
    (throw (Throwable. "Action not defined for protein sequence."))))

;; file

(defn read-up-xml-from-stream
  [rdr]
  (letfn [(process [x]
            (if (empty? x)
              nil
              (lazy-seq (cons (->uniprotProtein (first x))
                              (process (rest x))))))]
    (let [xml (xml/parse rdr)]
      (if (= (:tag xml) :uniprot)
        (process (filter #(= (:tag %) :entry)
                         (:content xml)))
        (throw (Throwable. "Doesn't appear to be a uniprot XML file."))))))

(defrecord uniprotXmlFile [file])

(extend-protocol bios/biosequenceFile

  uniprotXmlFile
  
  (biosequence-seq-file [this rdr]
    (read-up-xml-from-stream rdr))

  (file-path [this]
    (:file this)))

(defn init-uniprotxml-file
  "Initialises a uniprotXmlFile object."
  [path]
  (->uniprotXmlFile path))

;; persistence

(defrecord uniprotStore [file])

(defn index-uniprotxml-file
  "Indexes a uniprotXmlFile and returns a uniprotStore object."
  ([uniprotfile] (index-uniprotxml-file uniprotfile false))
  ([uniprotfile memory]
     (let [st (init-uniprot-store uniprotfile memory)]
       (bios/with-connection-to-store [st]
         (bios/with-biosequences [l uniprotfile]
           (dorun (pmap #(bios/save-object %) l))))
       st)))

(defn load-uniprot-store
  "Loads a uniprotStore."
  [dir]
  (let [file (first (fs/glob (str dir "/" "*.h2.db")))
        db-file (second (re-find  #"(.+)\.h2.db" (fs/absolute-path file)))]
    (if (not (nil? db-file))
      (assoc (->uniprotStore db-file) :db (bios/make-db-connection db-file true))
      (throw (Throwable. "DB file not found!")))))

;; web

(defn get-uniprot-stream
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

(defmacro with-wget-uniprot-sequence
  "Takes a list of accessions and returns a lazy list of sequences corresponding
   to the accession numbers. The type of sequence entry is specified by 'retype'
   and can be either :xml for a uniprotProtein object or :fasta for a fastaSequence,
   no other values are allowed. Uniprot requires an email address so one should 
   be provided in the 'email' argument."
  [[handle accessions retype email] & code]
  `(if (not (some #(= ~retype %) '(:xml :fasta)))
     (throw (Throwable. (str ~retype
                             " not allowed. "
                             "Only :xml and :fasta are allowed retype values.")))
     (let [rdr# (get-uniprot-stream (if (coll? ~accessions)
                                      ~accessions
                                      (list ~accessions))
                                    ~retype ~email)]
       (with-open [str# (java.io.PushbackReader. (io/reader rdr#))]
         (let [~handle (condp = ~retype
                         :xml (read-up-xml-from-stream str#)
                         :fasta (bios/read-fasta-from-stream str# :protein))]
           (try
             ~@code
             (catch Exception e#
               (throw (Throwable. e#)))
             (finally
               (.close rdr#))))))))

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
   And so on. Returns an empty list if no matches found. Offset refers to the
   start of the retrieved results and may be of use if a download fails. Uniprot
   requires an email so an email can be supplied using the email argument."
  ([term] (wget-uniprot-search term "" 0))
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
                (Integer. f))}))

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
        (throw (Throwable. (str "Error in sequence retrieval"
                                p)))))))

(defn- init-uniprot-store
  [file memory]
  (if memory
    (bios/init-in-mem-store (->uniprotStore 'in-memory))
    (bios/init-store
     (->uniprotStore (bios/index-file-name (bios/file-path file))))))


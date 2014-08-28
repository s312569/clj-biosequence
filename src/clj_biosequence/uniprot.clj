(ns clj-biosequence.uniprot
  (:require [clojure.data.xml :refer [parse]]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :refer [reader]]
            [clojure.zip :refer [node xml-zip]]
            [clojure.string :refer [split]]
            [clj-biosequence.core :as bs]
            [fs.core :refer [file? temp-file delete absolute-path]]))

(declare prot-name meta-data recommended-name alternative-name get-uniprot-stream organism init-indexed-uniprot)

;; interval

(defrecord uniprotInterval [src]

  bs/biosequenceInterval

  (start [this]
    (Integer/parseInt
     (zf/xml1-> (xml-zip (:src this)) :begin (zf/attr :position))))

  (end [this]
    (Integer/parseInt
     (zf/xml1-> (xml-zip (:src this)) :end (zf/attr :position))))

  (comp? [this]
    false))

;; feature

(defrecord uniprotFeature [src]

  bs/biosequenceFeature

  (feature-type [this]
    (:type (:attrs (:src this))))

  (interval-seq [this]
    (map #(->uniprotInterval %)
         (filter #(= (:tag %) :location)
                 (:content (:src this))))))

;; citation

(defrecord uniprotCitation [src]

  bs/biosequenceCitation

  (ref-type [this]
    (zf/xml1-> (xml-zip (:src this))
               :citation (zf/attr :type)))

  (title [this]
    (zf/xml1-> (xml-zip (:src this))
               :citation :title zf/text))

  (journal [this]
    (zf/xml1-> (xml-zip (:src this))
               :citation (zf/attr :name)))

  (year [this]
    (Integer/parseInt (zf/xml1-> (xml-zip (:src this))
                                 :citation (zf/attr :date))))

  (volume [this]
    (Integer/parseInt (zf/xml1-> (xml-zip (:src this))
                                 :citation (zf/attr :volume))))

  (pstart [this]

    (Integer/parseInt (zf/xml1-> (xml-zip (:src this))
                                 :citation (zf/attr :first))))

  (pend [this]
    (Integer/parseInt (zf/xml1-> (xml-zip (:src this))
                                 :citation (zf/attr :last))))

  (authors [this]
    (zf/xml-> (xml-zip (:src this))
              :citation :authorList :person (zf/attr :name))))

;; protein

(defrecord uniprotProtein [src]

  bs/Biosequence
  
  (accessions [uniprot]
    (zf/xml-> (xml-zip (:src uniprot)) :accession zf/text))

  (accession
    [this]
    (first (bs/accessions this)))

  (def-line [this]
    (str (recommended-name this) " | "
         (alternative-name this) " ["
         (get (organism this) "scientific") "]"))

  (bs-seq [this]
    (bs/clean-sequence
     (zf/xml1-> (xml-zip (:src this)) :sequence zf/text) :iupacAminoAcids))

  (fasta-string [this]
    (let [db (condp = (:dataset (meta-data this))
               "Swiss-Prot" "sp"
               "TrEMBL" "tr"
               "bs")]
      (str ">" db "|" (bs/accession this) "|"
           (prot-name this)
           " " (bs/def-line this)
           \newline
           (bs/bioseq->string this)
           \newline)))

  (protein? [this] true)

  (alphabet [this]
    :iupacAminoAcids)

  (feature-seq [this]
    (map #(->uniprotFeature (node %))
         (zf/xml-> (xml-zip (:src this)) :feature))))

;; IO

(defrecord uniprotReader [strm]

  bs/biosequenceReader

  (bs/biosequence-seq [this]
    (let [xml (parse (:strm this))]
      (map (fn [x]
             (vec x) ;; realizing all laziness
             (->uniprotProtein x))
           (filter #(= (:tag %) :entry)
                   (:content xml)))))

  (parameters [this])

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defn- init-uniprot-reader
  "Initialises a uniprot reader."
  [strm]
  (->uniprotReader strm))

(defrecord uniprotFile [file]

  bs/biosequenceIO

  (bs-reader [this]
    (init-uniprot-reader (reader (:file this))))

  bs/biosequenceFile

  (bs-path [this]
    (absolute-path (:file this)))

  (index-file [this]
    (init-indexed-uniprot (bs/bs-path this)))

  (index-file [this ofile]
    (init-indexed-uniprot (absolute-path ofile))))

(defrecord uniprotString [str]

  bs/biosequenceIO

  (bs-reader [this]
    (init-uniprot-reader (java.io.BufferedReader.
                          (java.io.StringReader. (:str this))))))

(defrecord uniprotConnection [acc-list retype email]

  bs/biosequenceIO

  (bs-reader [this]
    (let [s (get-uniprot-stream (:acc-list this) (:retype this) (:email this))]
      (condp = (:retype this)
        :xml (init-uniprot-reader (reader s))
        :fasta (bs/init-fasta-reader (reader s) :iupacAminoAcids)))))

(defn init-uniprotxml-file
  "Initialises a uniprotXmlFile object for use with bs-reader."
  [path]
  (if (file? path)
    (->uniprotFile path)
    (throw (Throwable. (str "File not found: " path)))))

(defn init-uniprot-string
  "Initialises a uniprot string for use with bs-reader. String needs
  to be valid uniprot xml."
  [str]
  (->uniprotString str))

(defn init-uniprot-connection
  "Initialises a connection which can be used with bs-reader to open a
  lazy list of uniprotProteins from the Uniprot servers."
  [accessions retype email]
  (if (#{:xml :fasta} retype)
    (let [l (if (coll? accessions) accessions (list accessions))]
      (->uniprotConnection l retype email))
    (throw (Throwable. (str retype " not supported. "
                            "Only :xml and :fasta are allowed retype values.")))))

;; web search

(defn uniprot-search
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
  ([term email] (uniprot-search term email 0))
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
                         (split #"\n")))]
       (if (empty? r)
         nil
         (lazy-cat r (uniprot-search term email (+ offset 1000)))))))

;; uniprot convienence functions

(defn organism
  "Returns a hash-map of organism names keyed by name type, eg.
  scientific etc."
  [uniprot]
  (assoc (->> (map node (zf/xml-> (xml-zip (:src uniprot)) :organism :name))
              (map #(vector (:type (:attrs %)) (first (:content %))))
              (into {}))
    :taxid (Integer/parseInt
            (zf/xml1-> (xml-zip (:src uniprot)) :organism :dbReference
                       (zf/attr= :type "NCBI Taxonomy") (zf/attr :id)))))

(defn lineage
  "Returns a list of the taxon terms of an organism"
  [uniprot]
  (zf/xml-> (xml-zip (:src uniprot)) :organism :lineage :taxon zf/text))

(defn prot-name
  "Returns the name of a uniprot as a string."
  [uniprot]
  (zf/xml1-> (xml-zip (:src uniprot)) :name zf/text))

(defn sequence-info
  "Returns a list of sequence information, including \"length\",
  \"mass\", \"checksum\", \"modified\", and \"version\"."
  [uniprot]
  (let [a (:attrs (node (zf/xml1-> (xml-zip (:src uniprot))
                                   :sequence)))]
    (assoc a
      :length (Integer/parseInt (:length a))
      :mass (Float/parseFloat (:mass a))
      :version (Integer/parseInt (:version a)))))

(defn recommended-name
  "Returns the recommended name. If multiple returns the first."
  [uniprot]
  (zf/xml1-> (xml-zip (:src uniprot))
             :protein :recommendedName :fullName zf/text))

(defn alternative-name
  "Returns the alternative name. If multiple returns the first."
  [uniprot & kinds]
  (zf/xml1-> (xml-zip (:src uniprot))
             :protein :alternativeName :fullName zf/text))

(defn gene
  "Returns a list of maps comprised of gene names and types."
  [uniprot]
  (map #(merge (:attrs (node (zf/xml1-> % :name)))
               (hash-map :gene (zf/xml1-> % :name zf/text)))
       (zf/xml-> (xml-zip (:src uniprot)) :gene)))

(defn citations
  [uniprot]
  (->> (zf/xml-> (xml-zip (:src uniprot)) :reference)
       (map #(->uniprotCitation (node %)))))

(defn comments
  "Returns the XML tree representations of all comments in the
  sequence."
  [uniprot]
  (map node (zf/xml-> (xml-zip (:src uniprot)) :comment)))

(defn subcellular-location
  "Returns a list of hash-maps describing the sub-cellular location of
   the sequence. Keys are :text, which describes the location,
   and :status, which describes the evidence."
  [uniprot]
  (->> (zf/xml-> (xml-zip (:src uniprot))
                 :comment
                 (zf/attr= :type "subcellular location")
                 :subcellularLocation :location)
       (map node)
       (map #(assoc (:attrs %) :text (first (:content %))))))

(defn db-references
  "Returns all db-references as a list of xml elements."
  [uniprot]
  (map node (zf/xml-> (xml-zip (:src uniprot))
                          :dbReference)))

(defn go-terms
  "Returns a list of GO terms. Only returns the term itself and not
   any other information, for more information see `db-references`."
  [uniprot]
  (zf/xml-> (xml-zip (:src uniprot))
            :dbReference
            (zf/attr= :type "GO")
            :property
            (zf/attr= :type "term")
            (zf/attr :value)))

(defn bp-go-terms
  "Returns a list of process GO terms.Only returns the term itself and
   not any other information, for more information see
   `db-references`."
  [uprot]
  (->> (go-terms uprot)
       (filter #(= \P (first %)))
       (map (fn [[a b & rest]]
              (apply str rest)))))

(defn mf-go-terms
  "Returns a list of process GO terms.Only returns the term itself and
   not any other information, for more information see
   `db-references`."
  [uprot]
  (->> (go-terms uprot)
       (filter #(= \M (first %)))
       (map (fn [[a b & rest]]
              (apply str rest)))))

(defn cc-go-terms
  "Returns a list of process GO terms.Only returns the term itself and
   not any other information, for more information see
   `db-references`."
  [uprot]
  (->> (go-terms uprot)
       (filter #(= \C (first %)))
       (map (fn [[a b & rest]]
              (apply str rest)))))

(defn existence
  "Protein existence evidence as a list of strings."
  [uniprot]
  (zf/xml-> (xml-zip (:src uniprot))
            :proteinExistence
            (zf/attr :type)))

(defn keywords
  "Keywords as a list of maps."
  [uniprot]
  (map #(merge (:attrs (node (zf/xml1-> %)))
               (hash-map :text (zf/xml1-> % zf/text)))
       (zf/xml-> (xml-zip (:src uniprot))
                 :keyword)))

(defn meta-data
  "Returns a map with uniprot meta-data. Keys - :dataset, :created, :modified
   and :version. All values are strings except :version, which is an integer."
  [uniprot]
  (let [s (xml-zip (:src uniprot))]
    {:dataset (zf/attr s :dataset)
     :created (zf/attr s :created)
     :modified (zf/attr s :modified)
     :version (if-let [f (zf/attr s :version)]
                (Integer. ^String f))}))

;; utilities

(defn- uniprot-process-request
  [address params file]
  (let [p (bs/post-req address params)]
    (letfn [(check [a c]
              (let [r (bs/get-req a {:follow-redirects false
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
                  (do (Thread/sleep 10000)
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
    (let [f (let [file (temp-file "up-seq-")]
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
        (finally (delete f))))))

;; indexing

(defrecord indexedUniprotReader [index strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (bs/indexed-seq this map->uniprotProtein))

  (get-biosequence [this accession]
    (bs/get-object this accession map->uniprotProtein))

  java.io.Closeable

  (close [this]
    (bs/close-index-reader this)))

(defrecord indexedUniprotFile [index path]

  bs/biosequenceIO

  (bs-reader [this]
    (->indexedUniprotReader (:index this) (:alphabet this) (bs/open-index-reader (:path this))))
  
  bs/biosequenceFile

  (bs-path [this]
    (absolute-path (:path this)))

  (empty-instance [this path]
    (init-indexed-uniprot path)))

(defn init-indexed-uniprot [file]
  (->indexedUniprotFile {} file))

(defmethod print-method clj_biosequence.uniprot.indexedUniprotFile
  [this w]
  (bs/print-tagged-index this w))

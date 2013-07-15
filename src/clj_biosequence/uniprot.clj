(ns clj-biosequence.uniprot
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clojure.java.io :as io]
            [clojure.pprint :as pp]
            [clj-biosequence.core :as bios]
            [clj-http.client :as client]
            [fs.core :as fs]
            [clojure.java.jdbc :as sql]))

(declare prot-names prot-name process-feature process-sequence process-cites process-dbref text? init-uniprot-store meta-data amino-acids nomenclature organism-name uniprot-process-request uniprot-sequence-helper read-xml-from-stream)

;; protein

(defrecord uniprotProtein [accession src])

(extend-protocol bios/Biosequence
  
  uniprotProtein

  (accessions [uniprot]
    (zf/xml-> (zip/xml-zip (:src uniprot)) :accession zf/text))

  (accession
    [this]
    (first (bios/accessions this)))

  (def-line [this]
    (str (string/join " | " (map #(:fullname %) 
                                 (flatten (vals (nomenclature this)))))
         " ["
         (first (:scientific (organism-name this)))
         "]"))

  (sequence-string [this]
    (:amino-acids (amino-acids this)))

  (fasta-string [this]
    (let [db (condp = (:dataset (meta-data this))
               "Swiss-Prot" "sp"
               "TrEMBL" "tr")]
      (str ">" db "|"
           (string/join "|" (bios/accessions this))
           "|"
           (prot-name this)
           " " (bios/def-line this)
           \newline
           (:amino-acids (amino-acids this))
           \newline)))

  (protein? [this] true)
  
  (org-scientific-name [this]
    (first (:scientific (organism-name this))))

  (created [this]
    (:created (meta-data this)))
  
  (modified [this]
    (:modified (meta-data this)))
  
  (version [this]
    (:version (meta-data this)))

  (database [this]
    (:dataset (meta-data this)))

  (taxonomy [this]
    (zf/xml-> (zip/xml-zip (:src this)) :organism :lineage :taxon zf/text))
  
  (taxid [this]
    (:ncbi-taxid (organism-name this))))

;; uniprot

(defn organism-name
  "Returns a map with keys - :scientific, :common, :full, :synonym, :abbreviation 
   and :ncbi-taxid. Values are lists of names attributed to the organism, except 
   :ncbi-taxid which is an integer."
  [uniprot]
  (let [values '("scientific" "common" "full" "synonym" "abbreviation")
        s (zip/xml-zip (:src uniprot))]
    (assoc (zipmap (map #(keyword %) values)
                   (map #(zf/xml-> s
                                   :organism
                                   :name
                                   (zf/attr= :type %)
                                   zf/text) values))
      :ncbi-taxid (read-string (zf/xml1-> s
                                          :organism
                                          :dbReference
                                          (zf/attr= :type "NCBI Taxonomy")
                                          (zf/attr :id))))))

(defn prot-name
  "Returns the name of a uniprot as a string."
  [uniprot]
  (zf/xml1-> (zip/xml-zip (:src uniprot)) :name zf/text))

(defn amino-acids
  "Returns a map with the keys - :mass, :checksum, :modified, :version, :amino-acids.
   All are strings except :mass, which is a long, and :version, which is an integer."
  [uniprot]
  (zf/xml1-> (zip/xml-zip (:src uniprot))
             :sequence
             process-sequence))

(defn nomenclature
  "Returns a map of descriptions. Keys - :recommended, :alternate, :submitted, :allergen.
   Values of each is a list of maps with keys - :fullname, :shortname and :ecnumber. All 
   values are strings except for :recommended which is map with the keys :fullname, 
   :shortname and :ecnumber."
  [uniprot]
  (let [s (zip/xml-zip (:src uniprot))
        rec (zf/xml-> s :protein :recommendedName prot-names)
        sub (zf/xml-> s :protein :submittedName prot-names)
        alt (zf/xml-> s :protein :alternativeName prot-names)
        all (zf/xml-> s :protein :allergenName zf/text)
        bio (zf/xml-> s :protein :biotechName zf/text)
        cd (zf/xml-> s :protein :cdAntigenName zf/text)
        inn (zf/xml-> s :protein :innName zf/text)]
    (merge {}
           (if (not (empty? rec)) {:recommended rec})
           (if (not (empty? alt)) {:alternate alt})
           (if (not (empty? sub)) {:submitted sub})
           (if (not (empty? all)) {:allergen all})
           (if (not (empty? bio)) {:allergen bio})
           (if (not (empty? cd)) {:allergen cd})
           (if (not (empty? inn)) {:allergen inn}))))

(defn gene
  "Returns a map with gene information from a uniprot. Keys - :type and :name. Values 
   are strings."
  [uniprot]
  (let [item (zf/xml-> (zip/xml-zip (:src uniprot)) :gene :name)]
    (map #(assoc {}
            :type (zf/attr % :type)
            :name (zf/text %)) item)))

(defn gene-location
  "Returns a map of gene location information from a uniprot."
  [uniprot]
  (let [values '("apicoplast" "chloroplast" "organellar chromatophore" "cyanelle" "hydrogenosome"
                 "mitochondrion" "non photosynthetic plastid" "nucleomorph" "plasmid" "plastid")]
    (zipmap (map #(keyword (string/replace % " " "-"))  values)
            (map #(zf/xml-> (zip/xml-zip (:src uniprot))
                            :geneLocation
                            (zf/attr= :type %)
                            :name
                            zf/text) values))))

(defn citation
  "Returns a list of citations from a uniprot. Each citation is a map with the keys - 
   :country, :last (page), :date, :pubmed, :institute, :name (journal), :first (page),
   :title, :city, :scope, :type, :consortium, :number, :authors, :source, :editors, 
   :publisher, :volume, :db. All are strings, except for :scope, :editor, :strain and 
   :author, which are lists of strings. The value of :source is a map with the keys - 
   :tissue, :transposons, :plasmid and :strain which are all lists of strings."
  [uniprot]
  (zf/xml-> (zip/xml-zip (:src uniprot))
            :reference
            process-cites))

(defn subcellular-location
  "Returns a list of maps specifying the sub-cellualr location of a uniprot. Eahc map has
   the keys - :comments, :orientation, :topology, :location and :evidence. All of which are
   strings."
  [uniprot]
  (let [items (zf/xml-> (zip/xml-zip (:src uniprot))
                        :comment
                        (zf/attr= :type "subcellular location")
                        :subcellularLocation)]
    (map #(assoc {}
            :evidence (zf/xml1-> % :location (zf/attr :status))
            :location (zf/xml1-> % :location zf/text)
            :topology (zf/xml1-> % :topology zf/text)
            :orientation (zf/xml1-> % :orientation zf/text)
            :comments (zf/xml1-> % :text zf/text)) items)))

(defn alternative-products
  "Returns a list of maps describing alternative products of a uniprot. Each map has the keys
   - :event, :isoform-id, :isoform-name and :comments. All of which are lists of strings except
   :comments which is a string."
  [uniprot]
  (let [items (zf/xml-> (zip/xml-zip (:src uniprot))
                        :comment
                        (zf/attr= :type "alternative products"))]
    (map #(assoc {}
            :event (zf/xml-> % :event (zf/attr :type))
            :isoform-id (zf/xml-> % :isoform :id zf/text)
            :isoform-name (zf/xml-> % :isoform :name zf/text)
            :comments (zf/xml1-> % :text zf/text))
         items)))

(defn interactions
  "Returns a list of maps describing interactions of a uniprot. Each map has the keys - :id,
  :label, :interactors, :experiments, :comments. All of whcih are lists of strings, except 
  :experiments and :comments which are strings."
  [uniprot]
  (let [items (zf/xml-> (zip/xml-zip (:src uniprot))
                        :comment
                        (zf/attr= :type "interaction"))]
    (map #(assoc {}
            :id (zf/xml-> % :id zf/text)
            :label (zf/xml-> % :label zf/text)
            :interactors (zf/xml-> % :interactant (zf/attr :type))
            :experiments (zf/xml1-> % :experiments (zf/attr :type))
            :comments (zf/xml1-> % :text zf/text))
         items)))

(defn mass-spectroscopy
  "Returns a list of maps describing the mass spectroscopy of a uniprot. Each map has the keys
   - :mass, :error, :method, :begin, :end and :comments. "
  [uniprot]
  (let [items (zf/xml-> (zip/xml-zip (:src uniprot))
                        :comment
                        (zf/attr= :type "mass spectrometry"))]
    (map #(assoc {}
            :mass (read-string (zf/xml1-> % (zf/attr :mass)))
            :error (read-string (zf/xml1-> % (zf/attr :error)))
            :method (zf/xml1-> % (zf/attr :method))
            :begin (zf/xml1-> % :location :begin (zf/attr :position))
            :end (zf/xml1-> % :location :end (zf/attr :position))
            :comments (zf/xml1-> % :text zf/text))
         items)))

(defn comments
  "Returns a map of comments describing a uniprot. Keys - :similarity and :function both of
   which have strings as values."
  [uniprot]
  (let [s (zip/xml-zip (:src uniprot))]
    (zipmap
     (map #(keyword (string/replace % " " "-"))
          (zf/xml-> s
                    :comment
                    text?
                    (zf/attr :type)))
     (zf/xml-> s
               :comment
               :text
               zf/text))))

(defn db-references
  "Returns a list of maps describing database cross-references of a uniprot. Each map has
   the keys :data, which has a map as a value with a range of different keys corresponding
   to the xml file, :db, a string describing the database, :id, a string of the id of the 
   cross-ref."
  [uniprot]
  (zf/xml-> (zip/xml-zip (:src uniprot))
            :dbReference
            process-dbref))

(defn existence
  "Returns a list of strings describing the evidence for the existence of a uniprot."
  [uniprot]
  (zf/xml-> (zip/xml-zip (:src uniprot))
            :proteinExistence
            (zf/attr :type)))

(defn keywords
  "A map of uniprot keywords. Keys are the keyword id with a string value of the keyword."
  [uniprot]
  (let [s (zip/xml-zip (:src uniprot))] 
    (zipmap
     (zf/xml-> s
               :keyword
               (zf/attr :id)
               keyword)
     (zf/xml-> s
               :keyword
               zf/text))))

(defn features
  "A list of maps describing the features of a uniprot. Each map has the keys - :id, :description,
  :type, :begin, :end. All have string values except :begin and :end which are integers."
  [uniprot]
  (zf/xml-> (zip/xml-zip (:src uniprot))
            :feature
            process-feature))

;; file

(defn read-xml-from-stream
  [rdr]
  (letfn [(process [x]
            (if (empty? x)
              nil
              (lazy-seq (cons (->uniprotProtein (zf/xml1-> (zip/xml-zip (first x))
                                                           :accession
                                                           zf/text)
                                                (first x))
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
    (read-xml-from-stream rdr))

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
         (bios/with-biosequences-in-file [l uniprotfile]
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
  "Returns a GZIPInputStream from Uniprot with the results of a batch fetch command for
   the sequences in a collection of accessions."
  [accessions class email]
  (if (empty? accessions)
    nil
    (let [f (let [file (fs/temp-file "up-seq-")]
              (doseq [s accessions]
                (spit file (str s "\n") :append true))
              file)]
      (:body (uniprot-process-request "http://www.uniprot.org/batch/"
                                 {:client-params {"http.useragent"
                                                  (str "clj-http " email)}
                                  :multipart [{:name "file" :content f}
                                              {:name "format" :content (name class)}]
                                  :follow-redirects false}
                                 f)))))

(defmacro with-wget-uniprot-sequence
  "Takes a list of accessions and returns a lazy list of sequences corresponding to the 
   accession numbers. The type of sequence entry is specified by 'retype' and can be either
   :xml for a uniprotProtein object or :fasta for a fastaSequence, no other values are allowed.
   Uniprot requires an email address so one should be provided in the 'email' argument."
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
                         :xml (read-xml-from-stream str#)
                         :fasta (bios/read-fasta-from-stream str#))]
           (try
             ~@code
             (catch Exception e#
               (throw (Throwable. e#)))
             (finally
               (.close rdr#))))))))

(defn wget-uniprot-search
  "Returns a lazy list of uniprot accession numbers satisfying the search term. The search
   term uses the same syntax as the uniprot web interface. For example:
   - to get all Schistosoma mansoni proteins in the proteome reference set term would be:
     'organism:6183 AND keyword:1185'
   - get all Schistosoma mansoni proteins in the proteome set that are intrinsic to the
     membrane:
     'taxonomy:6183 AND keyword:1185 AND go:0031224'
   - get all reviewed human entries:
     'reviewed:yes AND organism:9606'
   And so on. Returns an empty list if no matches found. Offset refers to the start of the
   retrieved results and may be of use if a download fails. Uniprot requires an email so an email
   can be supplied using the email argument."
  ([term] (wget-uniprot-search term "" 0))
  ([term email] (wget-uniprot-search term email 0))
  ([term email offset]
     (let [r (remove #(= % "")
                     (string/split
                      (:body
                       (client/get
                        (str "http://www.uniprot.org/uniprot/?query="
                             term
                             "&format=list"
                             (str "&offset=" offset)
                             "&limit=1000")
                        {:client-params {"http.useragent"
                                         (str "clj-http " email)}}))
                      #"\n"))]
       (if (empty? r)
         nil
         (lazy-cat r (wget-uniprot-search term email (+ offset 1000)))))))

;; utilities

(defn- uniprot-process-request
  [address params file]
  (try
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
                    (do (Thread/sleep (* 10
                                         (read-string
                                          (get (:headers r) "retry-after"))))
                        a)
                    (+ 1 c)))))]
        (if (some #(= (:status p) %) '(302 303))
          (do
            (if (get (:headers p) "retry-after")
              (Thread/sleep (read-string (get (:headers p) "retry-after"))))
            (check (get (:headers p) "location") 0))
          (throw (Throwable. (str "Error in sequence retrieval"
                                  p))))))
    (catch Exception e
      (throw (Throwable. e)))
    (finally
      (fs/delete file))))

(defn- init-uniprot-store
  [file memory]
  (if memory
    (bios/init-in-mem-store (->uniprotStore 'in-memory))
    (bios/init-store
     (->uniprotStore (bios/index-file-name (bios/file-path file))))))

(defn- text?
  "Checks if an element has a text tag."
  [zipper]
  (some #(= :text (:tag %)) (zip/children zipper)))

(defn- process-dbref
  [zipper]
  (assoc {:db (zf/attr zipper :type)
          :id (zf/attr zipper :id)}
    :data (zipmap 
           (zf/xml-> zipper
                     :property
                     (zf/attr :type))
           (zf/xml-> zipper
                     :property
                     (zf/attr :value)))))

(defn- process-cites
  [zipper]
  {:title (zf/xml1-> zipper :citation :title zf/text)
   :type (zf/xml1-> zipper :citation (zf/attr :type))
   :date (zf/xml1-> zipper :citation (zf/attr :date))
   :name (zf/xml1-> zipper :citation (zf/attr :name))
   :db (zf/xml1-> zipper :citation (zf/attr :db))
   :publisher (zf/xml1-> zipper :citation (zf/attr :publisher))
   :city (zf/xml1-> zipper :citation (zf/attr :city))
   :number (zf/xml1-> zipper :citation (zf/attr :number))
   :institute (zf/xml1-> zipper :citation (zf/attr :institute))
   :country (zf/xml1-> zipper :citation (zf/attr :country))
   :volume (zf/xml1-> zipper :citation (zf/attr :volume))
   :first (zf/xml1-> zipper :citation (zf/attr :first))
   :last (zf/xml1-> zipper :citation (zf/attr :last))
   :source (let [values '(:strain :plasmid :transposon :tissue)]
             (zipmap values
                     (map #(zf/xml-> zipper
                                     :source
                                     %
                                     zf/text) values)))
   :authors (zf/xml-> zipper
                      :citation
                      :authorList 
                      :person 
                      (zf/attr :name))
   :editors (zf/xml-> zipper
                      :citation
                      :editorList 
                      :person 
                      (zf/attr :name))
   :consortium (zf/xml-> zipper
                      :citation
                      :consortium
                      :person 
                      (zf/attr :name))
   :pubmed (zf/xml1-> zipper
                      :citation
                      :dbReference 
                      (zf/attr= :type "PubMed")
                      (zf/attr :id))
   :scope (zf/xml-> zipper :scope zf/text)})

(defn- process-feature
  [zipper]
  {:id (zf/attr zipper :id)
   :description (zf/attr zipper :description)
   :type (zf/attr zipper :type)
   :begin (if-let [f (zf/xml1-> zipper :location :begin (zf/attr :position))]
            (Integer. f))
   :end (if-let [f (zf/xml1-> zipper :location :end (zf/attr :position))]
          (Integer. f))
   :position (if-let [f (zf/xml1-> zipper :location :position (zf/attr :position))]
               (Integer. f))})

(defn- process-sequence
  [zipper]
  {:mass (if-let [f (zf/attr zipper :mass)]
           (read-string f))
   :checksum (zf/attr zipper :checksum)
   :modified (zf/attr zipper :modified)
   :version (if-let [f (zf/attr zipper :version)]
              (Integer. f))
   :amino-acids (string/replace (zf/text zipper)
                                " "
                                "")})

(defn- prot-names 
  [zipper]
  {:fullname (zf/xml1-> zipper :fullName zf/text)
   :shortname (zf/xml1-> zipper :shortName zf/text)
   :ecnumber (zf/xml1-> zipper :ecNumber zf/text)})


(defn- meta-data
  "Returns a map with uniprot meta-data. Keys - :dataset, :created, :modified and :version. 
   All values are strings except :version, which is an integer."
  [uniprot]
  (let [s (zip/xml-zip (:src uniprot))]
    {:dataset (zf/attr s :dataset)
     :created (zf/attr s :created)
     :modified (zf/attr s :modified)
     :version (if-let [f (zf/attr s :version)]
                (Integer. f))}))
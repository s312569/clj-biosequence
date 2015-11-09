(ns clj-biosequence.uniprot
  (:require [clojure.data.xml :refer [parse]]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :refer [reader]]
            [clojure.zip :refer [node]]
            [clojure.string :refer [split]]
            [clj-biosequence.core :as bs]
            [fs.core :refer [file? temp-file delete absolute-path]]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; citation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord uniprotCitation [src])

(extend uniprotCitation
  bs/biosequenceCitation
  (assoc bs/default-biosequence-citation
    :title (fn [this] (bs/get-text this :citation :title))
    :journal (fn [this] (bs/get-one this :citation (zf/attr :name)))
    :year (fn [this] (Integer/parseInt
                      (bs/get-one this :citation (zf/attr :date))))
    :volume (fn [this] (bs/get-one this :citation (zf/attr :volume)))
    :pstart (fn [this] (bs/get-one this :citation (zf/attr :first)))
    :pend (fn [this] (bs/get-one this :citation (zf/attr :last)))
    :authors (fn [this]
               (bs/get-list this :citation :authorList
                            :person (zf/attr :name)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; interval
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord uniprotInterval [src])

(extend uniprotInterval
  bs/biosequenceInterval
  (assoc bs/default-biosequence-interval
    :start
    (fn [this]
      (if-let [s (bs/get-one this :begin (zf/attr :position))]
        (Integer/parseInt s) (bs/point this)))
    :end
    (fn [this]
      (if-let [e (bs/get-one this :end (zf/attr :position))]
        (Integer/parseInt e) (bs/point this)))
    :point
    (fn [this]
      (if-let [p (bs/get-one this :position (zf/attr :position))]
        (Integer/parseInt p)))
    :comp? (fn [this] false))
  bs/biosequenceStatus
  (assoc bs/default-biosequence-status
    :status
    (fn [this]
      (bs/get-one this :position (zf/attr :status))))
  bs/biosequenceEvidence
  (assoc bs/default-biosequence-evidence
    :evidence
    (fn [this]
      (split (bs/get-one this :position (zf/attr :evidence))
             #"\s+"))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; feature
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord uniprotFeature [src])

(extend uniprotFeature
  bs/biosequenceNameObject
  (assoc bs/default-biosequence-nameobject
    :obj-type
    (fn [this]
      (:type (:attrs (:src this))))
    :obj-description
    (fn [this]
      (:description (:attrs (:src this))))
    :obj-id
    (fn [this]
      (:id (:attrs (:src this)))))
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession
    (fn [this]
      (:id (:attrs (:src this)))))
  bs/biosequenceStatus
  (assoc bs/default-biosequence-status
    :status
    (fn [this]
      (:status (:attrs (:src this)))))
  bs/biosequenceDescription
  (assoc bs/default-biosequence-description
    :description
    (fn [this]
      (:description (:attrs (:src this)))))
  bs/biosequenceEvidence
  (assoc bs/default-biosequence-evidence
    :evidence
    (fn [this]
      (if-let [e (:evidence (:attrs (:src this)))]
        (split e #"\s+"))))
  bs/biosequenceCitations
  (assoc bs/default-biosequence-citations
    :citation-key
    (fn [this]
      (:ref (:attrs (:src this)))))
  bs/biosequenceIntervals
  (assoc bs/default-biosequence-intervals
    :intervals (fn [this]
                 (map #(->uniprotInterval %)
                      (filter #(= (:tag %) :location)
                              (:content (:src this))))))
  bs/Biosequence
  (assoc bs/default-biosequence-biosequence
    :bs-seq
    (fn [this]
      (vec (:sequence this))))
  bs/biosequenceVariant
  (assoc bs/default-biosequence-variant
    :original (fn [this] (bs/get-text this :original))
    :variant (fn [this] (bs/get-list this :variation
                                     zf/text))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; db ref
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- get-db-prop
  [this]
  (into {}
   (map #(vector (zf/xml1-> % (zf/attr :type))
                 (zf/xml1-> % (zf/attr :value)))
        (bs/get-list this :property))))

(defn- go-help
  [this f]
  (if-let [g ((bs/db-properties this) "term")]
    (f (split g #":"))))

(defrecord uniprotDbRef [src])

(extend uniprotDbRef
  bs/biosequenceDbRef
  (assoc bs/default-biosequence-dbref
    :database-name (fn [this] (:type (:attrs (:src this))))
    :object-id (fn [this] (:id (:attrs (:src this))))
    :db-properties (fn [this]
                     (get-db-prop this)))
  bs/biosequenceGoTerm
  (assoc bs/default-biosequence-goterm
    :go-id (fn [this] (bs/object-id this))
    :go-term (fn [this] (go-help this second))
    :go-component (fn [this] (go-help this first)))
  bs/biosequenceEvidence
  (assoc bs/default-biosequence-evidence
    :evidence
    (fn [this] ((bs/db-properties this) "evidence"))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; org ref
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord uniprotTaxRef [src])

(extend uniprotTaxRef
  bs/biosequenceTaxonomy
  (assoc bs/default-biosequence-tax
    :tax-name
    (fn [this] (bs/get-text this :name
                            (zf/attr= :type "scientific")))
    :common-name
    (fn [this] (bs/get-text this :name
                            (zf/attr= :type "common")))
    :lineage
    (fn [this] (let [l (bs/get-list this :lineage :taxon
                                    zf/text)]
                 (apply str (interpose ";" l)))))
  bs/biosequenceDbRefs
  (assoc bs/default-biosequence-dbrefs
    :get-db-refs
    (fn [this] (->> (bs/get-list this :dbReference)
                    (map #(->uniprotDbRef (node %))))))
  bs/biosequenceEvidence
  (assoc bs/default-biosequence-evidence
    :evidence
    (fn [this]
      (if-let [e (bs/get-text this (zf/attr :evidence))]
        (split e #"\s+")))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; genes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord uniprotGene [src])

(extend uniprotGene
  bs/biosequenceGene
  (assoc bs/default-biosequence-gene
    :locus (fn [this] (bs/get-text this :name
                                   (zf/attr= :type "ordered locus")))
    :orf (fn [this] (bs/get-text this :name
                                 (zf/attr= :type "ORF"))))
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession
    (fn [this]
      (bs/get-text this :name (zf/attr= :type "primary"))))
  bs/biosequenceSynonyms
  (assoc bs/default-biosequence-synonyms
    :synonyms
    (fn [this] (bs/get-list this :name
                            (zf/attr= :type "synonym")
                            zf/text))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; comment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- get-subcell
  [this tag]
  (let [ls (bs/get-list this :subcellularLocation)]
    (map #(vector (zf/xml1-> % tag zf/text)
                  (zf/xml1-> % tag (zf/attr :evidence)))
         ls)))

(defrecord uniprotComment [src])

(extend uniprotComment
  bs/biosequenceSubCellLoc
  (assoc bs/default-biosequence-subcell
    :subcell-location (fn [this] (get-subcell this :location))
    :subcell-topol (fn [this] (get-subcell this :topology))
    :subcell-orient (fn [this] (get-subcell this :orientation)))
  bs/biosequenceNameObject
  (assoc bs/default-biosequence-nameobject
    :obj-type (fn [this]
                (bs/get-one this (zf/attr :type)))
    :obj-text (fn [this]
                (bs/get-text this :text))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; protein
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- u-date
  [str]
  (bs/make-date str (bs/make-date-format "yyyy-MM-dd")))

(defn- get-short-full
  [this]
  (if this
    (remove nil?
            (vector (bs/get-text {:src (node this)} :fullName)
                    (bs/get-text {:src (node this)} :shortName)))))

(defn- all-names
  [this]
  (let [n (mapcat #(% this)
                  (list bs/alternate-names
                        bs/submitted-names
                        bs/allergen-names
                        bs/biotech-names
                        bs/cd-antigen-names
                        bs/innnames))]
    (if (seq n) n)))

(defrecord uniprotProtein [src alphabet])

(extend uniprotProtein
  bs/Biosequence
  (assoc bs/default-biosequence-biosequence
         :bs-seq (fn [this] (vec (:sequence this)))
         :protein? (fn [this] true)
         :alphabet (fn [this] (:alphabet this))
         :moltype (fn [this] "AA")
         :keywords (fn [this]
                     (map #(vector (zf/attr % :id)
                                   (zf/text %))
                          (bs/get-list this :keyword))))
  bs/biosequenceID
  (assoc bs/default-biosequence-id
         :accessions (fn [this]
                       (concat (bs/get-list this :accession zf/text)
                               (list (bs/get-text this :name))))
         :accession (fn [this] (first (bs/accessions this)))
         :version (fn [this] (:version (:attrs (:src this))))
         :creation-date (fn [this]
                          (u-date (:created (:attrs (:src this)))))
         :update-date (fn [this]
                        (u-date (:modified (:attrs (:src this)))))
         :dataset (fn [this] (:dataset (:attrs (:src this)))))
  bs/biosequenceName
  (assoc bs/default-biosequence-name
         :names
         (fn [this]
           (let [r (get-short-full (bs/get-one this :protein
                                               :recommendedName))]
             (or r (all-names this) (bs/accessions this))))
         :alternate-names
         (fn [this]
           (mapcat #(get-short-full %)
                   (bs/get-list this :protein :alternativeName)))
         :submitted-names
         (fn [this]
           (mapcat #(get-short-full %)
                   (bs/get-list this :protein :submittedName)))
         :allergen-names
         (fn [this]
           (mapcat #(get-short-full %)
                   (bs/get-list this :protein :allergenName)))
         :biotech-names
         (fn [this]
           (mapcat #(get-short-full %)
                   (bs/get-list this :protein :biotechName)))
         :cd-antigen-names
         (fn [this]
           (mapcat #(get-short-full %)
                   (bs/get-list this :protein :cdAntigenName)))
         :innnames
         (fn [this]
           (mapcat #(get-short-full %)
                   (bs/get-list this :protein :innName))))
  bs/biosequenceDescription
  (assoc bs/default-biosequence-description
         :description (fn [this]
                        (str (first (bs/names this)) " ["
                             (-> (bs/tax-refs this)
                                 first
                                 bs/tax-name) "]")))
  bs/biosequenceCitations
  (assoc bs/default-biosequence-citations
         :citations
         (fn [this]
           (->> (bs/get-list this :reference)
                (map #(->uniprotCitation (node %))))))
  bs/biosequenceFeatures
  (assoc bs/default-biosequence-features
         :feature-seq
         (fn [this]
           (->> (bs/get-list this :feature)
                (map #(->uniprotFeature (node %)))
                (map #(assoc % :sequence
                             (bs/get-one % :location (zf/attr :sequence))))
                (bs/clean-sequences (:alphabet this)))))
  bs/biosequenceTaxonomies
  (assoc bs/default-biosequence-taxonomies
         :tax-refs
         (fn [this] (->> (bs/get-list this :organism)
                         (map #(->uniprotTaxRef (node %))))))
  bs/biosequenceGenes
  (assoc bs/default-biosequence-genes
         :genes
         (fn [this] (->> (bs/get-list this :gene)
                         (map #(->uniprotGene (node %))))))
  bs/biosequenceComments
  (assoc bs/default-biosequence-comments
         :comments
         (fn [this]
           (->> (bs/get-list this :comment)
                (map #(->uniprotComment (node %))))))
  bs/biosequenceSubCellLocs
  (assoc bs/default-biosequence-subcells
         :subcell-locations
         (fn [this]
           (->> (bs/comments this)
                (filter #(= (bs/obj-type %) "subcellular location")))))
  bs/biosequenceDbRefs
  (assoc bs/default-biosequence-dbrefs
         :get-db-refs
         (fn [this]
           (->> (bs/get-list this :dbReference)
                (map #(->uniprotDbRef (node %))))))
  bs/biosequenceGoTerms
  (assoc bs/default-biosequence-goterms
         :gos (fn [this]
                (filter #(= (bs/database-name %) "GO")
                        (bs/get-db-refs this))))
  bs/biosequenceEvidence
  (assoc bs/default-biosequence-evidence
         :evidence
         (fn [this]
           (bs/get-list this :proteinExistence (zf/attr :type))))
  bs/biosequenceProtein
  (assoc bs/default-biosequence-protein
         :calc-mol-wt
         (fn [this]
           (Float/parseFloat
            (bs/get-one this :sequence (zf/attr :mass))))
         :ecs
         (fn [this]
           (mapcat #(bs/get-text {:src (node this)} :ecNumber)
                   (concat (bs/get-list this :protein :recommendedName)
                           (bs/get-list this :protein :alternativeName)
                           (bs/get-list this :protein :submittedName))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; IO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord uniprotReader [strm alphabet]
  bs/biosequenceReader
  (bs/biosequence-seq [this]
    (let [xml (parse (:strm this))]
      (->> (map (fn [x]
                  (vec x) ;; realizing all laziness
                  (->uniprotProtein x (:alphabet this)))
                (filter #(= (:tag %) :entry)
                        (:content xml)))
           (map #(assoc % :sequence
                        (bs/get-text % :sequence)))
           (bs/clean-sequences (:alphabet this)))))
  java.io.Closeable
  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord uniprotFile [file opts alphabet])

(extend uniprotFile
  bs/biosequenceIO
  {:bs-reader
   (fn [this]
     (->uniprotReader (apply bs/bioreader (bs/bs-path this)
                             (:opts this))
                      (:alphabet this)))}
  bs/biosequenceFile
  bs/default-biosequence-file)

(defn init-uniprot-file
  "Initialises a uniprotXmlFile object for use with bs-reader."
  [path & {:keys [alphabet opts] :or {alphabet :iupacAminoAcids
                                      opts ()}}]
  (if (file? path)
    (->uniprotFile path opts alphabet)
    (throw (Throwable. (str "File not found: " path)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; connection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- uniprot-process-request
  [address params file]
  (let [p (bs/post-req address params)]
    (letfn [(check [a c]
              (let [r (bs/get-req a {:follow-redirects false
                                     :as :stream})]
                (cond
                 (nil? (get (:headers r) "retry-after"))
                 (cond (= (:status r) 200)
                       r
                       (= (:status r) 302)
                       (recur (get (:headers r) "Location") 0)
                       :else
                       (throw (Throwable. (str "Error in sequence retrieval: code"
                                               (:status r)))))
                 (> c 50)
                 (throw (Throwable. "Too many tries."))
                 :else
                 (recur (do (Thread/sleep 10000) a) (inc c)))))]
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

(defrecord uniprotConnection [acc-list retype email alphabet]
  bs/biosequenceIO
  (bs-reader [this]
    (let [s (get-uniprot-stream (:acc-list this)
                                (:retype this) (:email this))]
      (condp = (:retype this)
        :xml (->uniprotReader (bs/bioreader s) (:alphabet this))
        :fasta (bs/init-fasta-reader (bs/bioreader s) (:alphabet this))))))

(defn init-uniprot-connection
  "Initialises a connection which can be used with bs-reader to open a
  lazy list of uniprotProteins from the Uniprot servers."
  [accessions retype email & {:keys [alphabet] :or {alphabet :iupacAminoAcids}}]
  {:pre [(#{:xml :fasta} retype)]}
  (let [l (if (coll? accessions) accessions (list accessions))]
    (->uniprotConnection l retype email alphabet)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; web search
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
                    (-> (bs/get-req
                         (str "http://www.uniprot.org/uniprot/?query="
                              term
                              "&format=list"
                              (str "&offset=" offset)
                              "&limit=1000")
                         {:client-params {"http.useragent"
                                          (str "clj-http " email)}})
                      (:body)
                      (split #"\n")))]
     (cond (< (count r) 1000) r
           (not (seq r)) nil
           :else
           (concat r (lazy-cat (uniprot-search term email (+ offset 1000))))))))

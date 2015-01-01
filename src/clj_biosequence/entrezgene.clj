(ns clj-biosequence.entrezgene
  (:require [clojure.data.xml :refer [parse]]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :refer [split]]
            [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.eutilities :as eu]
            [clj-biosequence.core :as bs]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; entrez interval
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defprotocol entrezComments
  (entrez-comments [this]
    "Returns entrezComment record.")
  (property-comments [this]
    "Returns property entrezComment records.")
  (product-comments [this]
    "Returns product entrezComment records.")
  (locus-comments [this]
    "Returns locus entrezComment records.")
  (ref-gene-comments [this]
    "Returns ref-gene entrezComment records.")
  (homology-comments [this]
    "Returns homology entrezComment records."))

(def default-entrez-comments
  {:entrez-comments bs/return-nil
   :property-comments bs/return-nil
   :product-comments bs/return-nil
   :locus-comments bs/return-nil
   :ref-gene-comments bs/return-nil
   :homology-comments bs/return-nil})

(defn- all-comments
  [this]
  {:pre [(satisfies? entrezComments this)]}
  [this]
  (concat (entrez-comments this)
          (property-comments this)
          (product-comments this)
          (locus-comments this)
          (ref-gene-comments this)
          (homology-comments this)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; entrez interval
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezInterval [src])

(extend entrezInterval
  bs/biosequenceInterval
  (assoc bs/default-biosequence-interval
    :start
    (fn [this]
      (Integer/parseInt (bs/get-text this :Seq-interval_from)))
    :end
    (fn [this]
      (Integer/parseInt (bs/get-text this :Seq-interval_to)))
    :comp?
    (fn [this]
      (if-let [s (bs/get-one this :Seq-interval_strand :Na-strand
                             (zf/attr :value))]
        (= s "minus")))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; entrez seq location
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezSeqLocation [src])

(extend entrezSeqLocation
  bs/biosequenceIntervals
  (assoc bs/default-biosequence-intervals
    :intervals
    (fn [this]
      (map #(->entrezInterval (zip/node %))
           (bs/get-list this :Seq-loc_int :Seq-interval)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; db tag
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezDbTag [src])

(extend entrezDbTag
  bs/biosequenceDbRef
  (assoc bs/default-biosequence-dbref
    :database-name (fn [this] (bs/get-text this :Dbtag_db))
    :object-id
    (fn [this]
      (or
       (bs/get-text this :Dbtag_tag :Object-id :Object-id_str)
       (bs/get-text this :Dbtag_tag :Object-id :Object-id_id)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; other src
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezOtherSource [src])

(extend entrezOtherSource
  bs/biosequenceUrl
  (assoc bs/default-biosequence-url
    :url (fn [this]
           (bs/get-text this :Other-source_url))
    :pre-text (fn [this]
                (bs/get-text this :Other-source_pre-text))
    :anchor (fn [this]
              (bs/get-text this :Other-source_anchor))
    :post-text (fn [this]
                 (bs/get-text this :Other-source_post-text)))
  bs/biosequenceDbRefs
  (assoc bs/default-biosequence-dbrefs
    :get-db-refs
    (fn [this]
      (if-let [d (bs/get-one this :Other-source_src :Dbtag)]
        (list (->entrezDbTag (zip/node d)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; xtra terms
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezExtraTerm [src])

(extend entrezExtraTerm
  bs/biosequenceNameObject
  (assoc bs/default-biosequence-nameobject
    :obj-name (fn [this] (bs/get-text this :Xtra-Terms_tag))
    :obj-value (fn [this] (bs/get-text this :Xtra-Terms_value))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; comments
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezGeneComment [src])

(extend entrezGeneComment
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession (fn [this]
                 (bs/get-text this :Gene-commentary_accession))
    :accessions (fn [this]
                  (remove nil? (list (bs/accession this))))
    :version (fn [this]
               (bs/get-text this :Gene-commentary_version)))
  entrezComments
  (assoc default-entrez-comments
    :entrez-comments
    (fn [this]
      (map #(->entrezGeneComment (zip/node %))
           (bs/get-list this :Gene-commentary_comment
                        :Gene-commentary)))
    :property-comments
    (fn [this]
      (map #(->entrezGeneComment (zip/node %))
           (bs/get-list this :Gene-commentary_properties
                        :Gene-commentary)))
    :product-comments
    (fn [this]
      (map #(->entrezGeneComment (zip/node %))
           (bs/get-list this :Gene-commentary_products
                        :Gene-commentary))))
  bs/biosequenceComments
  (assoc bs/default-biosequence-comments
    :comments (fn [this] (all-comments this)))
  bs/biosequenceNameObject
  (assoc bs/default-biosequence-nameobject
    :obj-type
    (fn [this]
      (Integer/parseInt
       (bs/get-text this :Gene-commentary_type)))
    :obj-value
    (fn [this]
      (bs/get-one this :Gene-commentary_type (zf/attr :value)))
    :obj-heading
    (fn [this]
      (bs/get-text this :Gene-commentary_heading))
    :obj-label
    (fn [this]
      (bs/get-text this :Gene-commentary_label))
    :obj-text
    (fn [this]
      (bs/get-text this :Gene-commentary_text))))

(defn xtra-terms
  "Returns a list of vectors with the tag and value pairs from an
  Xtra-Terms in an entrezGeneComment."
  [comment]
  (map #(let [x {:src (zip/node %)}]
          (vector (bs/get-text x :Xtra-Terms_tag)
                  (bs/get-text x :Xtra-Terms_value)))
       (bs/get-list comment :Gene-commentary_xtra-properties
                 :Xtra-Terms)))

(defn pmids
  "Returns a list of pubmed ids from an entrezGeneComment."
  [comment]
  (bs/get-list comment :Gene-commentary_refs :Pub :Pub_pmid
            :PubMedId zf/text))

(defn other-sources
  "Returns a list of other entrezOtherSource records from an
  entrezGeneComment."
  [comment]
  (map #(->entrezOtherSource (zip/node %))
       (bs/get-list comment :Gene-commentary_source :Other-source)))

(defn seq-locations
  "Returns a list of entrezSeqLocation records from an
  entrezGeneComment."
  [comment]
  (map #(->entrezSeqLocation (zip/node %))
       (concat (bs/get-list comment :Gene-commentary_seqs
                            :Seq-loc)
               (bs/get-list comment
                            :Gene-commentary_genomic-coords
                            :Seq-loc))))

(defn walk-comments
  "Recursively steps through levels of comments. For each level a
  function and a value can be supplied to filter comments. A final
  argument can be supplied to extract a field from the final set of
  comments."
  [start & args]
  (letfn [(f [s [f v]]
            (if (satisfies? entrezComments s)
              (->> (filter #(= (f %) v)
                           (bs/comments s))
                   (remove nil?))
              (->> (mapcat #(bs/comments %) s)
                   (filter #(= (f %) v))
                   (remove nil?))))]
    (loop [s start a args]
      (if (and (seq s) (seq a))
        (if (= (count a) 1)
          (remove nil? (map #((first a) %) s))
          (recur (f s (take 2 a)) (drop 2 a)))
        s))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; entrez gene source
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezGeneSource [src])

(defn source-name
  [genesource]
  "Returns the name of a entrezGeneSource."
  (bs/get-text genesource :Gene-source_src))

(defn source-id
  [genesource]
  "Returns the id of a entrezGeneSource."
  (bs/get-text genesource :Gene-source_src-int))

(defn string1
  "Returns the string1 entry from a entrezGeneSource."
  [genesource]
  (bs/get-text genesource :Gene-source_src-str1))

(defn string2
  "Returns the string1 entry from a entrezGeneSource."
  [genesource]
  (bs/get-text genesource :Gene-source_src-str2))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; entrez map
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezMap [src])

(defn display-string
  "Returns list of display strings from an entrezMap."
  [entrezmap]
  (bs/get-one entrezmap :Maps_display-str zf/text))

(defn map-proxy
  "Returns list of urls to non mapviewer mapviewing resource from an
  entrezMap."
  [entrezmap]
  (bs/get-one entrezmap :Maps_method_map-type :Maps_proxy zf/text))

(defn map-method-type
  "Returns list of units used in display-str to query mapviewer from
  an entrezMap."
  [entrezmap]
  (bs/get-one entrezmap :Maps_method :Maps_method_map-type
           (zf/attr :value)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; gene track
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- e-date
  [y m d]
  (let [f (bs/make-date-format "yyyy-MM-dd")]
    (bs/make-date (apply str (interpose "-" (list y m d))) f)))

(defrecord entrezGeneTrack [src])

(extend entrezGeneTrack
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession (fn [this]
                 (Integer/parseInt
                  (bs/get-text this :Gene-track_geneid)))
    :accessions (fn [this]
                  (remove nil?
                          (list (bs/accession this))))
    :creation-date
    (fn [this]
      (if-let [cd (bs/get-one this :Gene-track_create-date
                              :Date :Date_std :Date-std)]
        (e-date (zf/xml1-> cd :Date-std_year zf/text)
                (zf/xml1-> cd :Date-std_month zf/text)
                (zf/xml1-> cd :Date-std_day zf/text))))
    :update-date
    (fn [this]
      (if-let [cd (bs/get-one this :Gene-track_update-date
                              :Date :Date_std :Date-std)]
        (e-date (zf/xml1-> cd :Date-std_year zf/text)
                (zf/xml1-> cd :Date-std_month zf/text)
                (zf/xml1-> cd :Date-std_day zf/text))))
    :discontinue-date
    (fn [this]
      (if-let [cd (bs/get-one this :Gene-track_discontinue-date
                              :Date :Date_std :Date-std)]
        (e-date (zf/xml1-> cd :Date-std_year zf/text)
                (zf/xml1-> cd :Date-std_month zf/text)
                (zf/xml1-> cd :Date-std_day zf/text)))))
  bs/biosequenceStatus
  (assoc bs/default-biosequence-status
  :status
  (fn [this]
    (if-let [s (bs/get-one this :Gene-track_status)]
      (vector (zf/xml1-> s (zf/attr :value))
              (Integer/parseInt (zf/xml1-> s zf/text)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; org name
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezOrgName [src])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; org-ref
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn org-name
  "Returns a org-name record from an org ref."
  [org]
  (->entrezOrgName
   (zip/node (bs/get-one org :Org-ref_orgname :OrgName))))

(defrecord entrezOrgRef [src])

(extend entrezOrgRef
  bs/biosequenceTaxonomy
  (assoc bs/default-biosequence-tax
    :tax-name (fn [this]
                (bs/get-text this :Org-ref_taxname))
    :common-name (fn [this]
                   (bs/get-text this :Org-ref_common))
    :mods (fn [this]
            (bs/get-list this :Org-ref_mod :Org-ref_mod_E zf/text))
    :lineage (fn [this]
               (bs/get-text (org-name this) :OrgName_lineage)))
  bs/biosequenceDbRefs
  (assoc bs/default-biosequence-dbrefs
    :get-db-refs
    (fn [this]
      (map #(->entrezDbTag (zip/node %))
           (bs/get-list this :Org-ref_db :Dbtag))))
  bs/biosequenceSynonyms
  (assoc bs/default-biosequence-synonyms
    :synonyms
    (fn [this]
      (bs/get-list this :Org-ref_syn :Org-ref_syn_E zf/text))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; subsource
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezSubSource [src])

(defn subsource-type
  "Returns a vector of the value and integer code of the subtype of a
  subsource."
  [sub]
  (let [s (bs/get-one sub :SubSource_subtype)]
    (vector (zf/xml1-> s (zf/attr :value))
            (zf/xml1-> s zf/text))))

(defn subsource-name
  "Returns the name of a subsource."
  [sub]
  (bs/get-text sub :SubSource_name))

(defn subsource-attribute
  "Returns the attribute from a subsource."
  [sub]
  (bs/get-text sub :SubSource_attrib))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; pcr primers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezPcrPrimers [src])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; biosource
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezBiosource [src])

(defn genome
  "Returns a vector of the value and integer code for the genome from
  an entrezBiosource."
  [biosource]
  (let [s (bs/get-one biosource :BioSource_genome)]
    (vector (zf/xml1-> s (zf/attr :value))
            (Integer/parseInt (zf/xml1-> s zf/text)))))

(defn origin
  "Returns a vector of the value and integer code for the origin from
  an entrezBiosource."
  [biosource]
  (let [s (bs/get-one biosource :BioSource_origin)]
    (vector (zf/xml1-> s (zf/attr :value))
            (Integer/parseInt (zf/xml1-> s zf/text)))))

(defn org-ref
  "Returns an org-ref from an entrezBiosource."
  [biosource]
  (if-let [d (bs/get-one biosource :BioSource_org :Org-ref)]
    (->entrezOrgRef (zip/node d))))

(defn subsource
  "Returns a list of subsource records from an entrezBiosource."
  [biosource]
  (map #(->entrezSubSource (zip/node %))
       (bs/get-list biosource :BioSource_subtype :SubSource)))

(defn pcr-primers
  "Returns a pcr primers record from an entrezBiosource."
  [biosource]
  (if-let [d (bs/get-one biosource :BioSource_pcr-primers
                         :PCRReactionSet)]
    (->entrezOrgRef (zip/node d))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; protein
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezProtein [src])

(extend entrezProtein
  bs/biosequenceProtein
  (assoc bs/default-biosequence-protein
    :ecs (fn [this]
          (bs/get-list this :Prot-ref_ec :Prot-ref_ec_E zf/text))
    :activities
    (fn [this]
      (bs/get-list this :Prot-ref_activity :Prot-ref_activity_E
                   zf/text)))
  bs/biosequenceSynonyms
  (assoc bs/default-biosequence-synonyms
    :synonyms
    (fn [this]
      (bs/get-list this :Prot-ref_name :Prot-ref_name_E zf/text)))
  bs/biosequenceDescription
  (assoc bs/default-biosequence-description
    :description (fn [this]
                   (bs/get-text this :Prot-ref_desc)))
  bs/biosequenceDbRefs
  (assoc bs/default-biosequence-dbrefs
    :get-db-refs
    (fn [this]
      (map #(->entrezDbTag (zip/node %))
           (bs/get-list this :Prot-ref_db :Dbtag)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; gene
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn entrez-gene-track
  "Returns an entrezGeneTrack record from a entrezgene."
  [entrezgene]
  (->entrezGeneTrack
   (zip/node (bs/get-one entrezgene :Entrezgene_track-info :Gene-track))))

(defn entrez-biosource
  "Returns a entrezgeneBiosource from an entrezgene."
  [entrezgene]
  (->entrezBiosource
   (zip/node (bs/get-one entrezgene :Entrezgene_source :BioSource))))

(defn entrez-location
  "Returns a list of location maps from an entrezgene"
  [entrezgene]
  (map #(->entrezMap (zip/node %))
       (bs/get-list entrezgene :Entrezgene_location :Maps)))

(defn entrez-gene-source
  "Returns a gene source record from an entrezgene."
  [entrezgene]
  (->entrezGeneSource
   (zip/node
    (bs/get-one entrezgene :Entrezgene_gene-source :Gene-source))))

(defn entrez-unique-keys
  "Returns a list of entrezDbtags representing unique keys."
  [entrezgene]
  (map #(->entrezDbTag (zip/node %))
       (bs/get-list entrezgene :Entrezgene_unique-keys :Dbtag)))

(defn entrez-non-unique-keys
  "Returns a list of entrezDbtags representing non-unique keys."
  [entrezgene]
  (map #(->entrezDbTag (zip/node %))
       (bs/get-list entrezgene :Entrezgene_non-unique-keys :Dbtag)))

(defrecord entrezGene [src])

(extend entrezGene
  bs/biosequenceGene
  (assoc bs/default-biosequence-gene
    :locus
    (fn [this]
      (bs/get-text this :Entrezgene_gene :Gene-ref
                   :Gene-ref_locus))
    :map-location
    (fn [this]
      (bs/get-text this :Entrezgene_gene :Gene-ref
                   :Gene-ref_maploc))
    :locus-tag
    (fn [this]
      (bs/get-text this :Entrezgene_gene :Gene-ref
                   :Gene-ref_locus-tag)))
  bs/biosequenceSynonyms
  (assoc bs/default-biosequence-synonyms
    :synonyms
    (fn [this]
      (bs/get-list this :Entrezgene_gene :Gene-ref
                   :Gene-ref_syn :Gene-ref_syn_E zf/text)))
  bs/biosequenceDescription
  (assoc bs/default-biosequence-description
    :description (fn [this]
                   (bs/get-text this :Entrezgene_gene :Gene-ref
                                :Gene-ref_desc)))
  bs/biosequenceDbRefs
  (assoc bs/default-biosequence-dbrefs
    :get-db-refs (fn [this]
                   (map #(->entrezDbTag (zip/node %))
                        (bs/get-list this :Entrezgene_gene :Gene-ref
                                     :Gene-ref_db :Dbtag))))
  entrezComments
  (assoc default-entrez-comments
    :entrez-comments
    (fn [this]
      (map #(->entrezGeneComment (zip/node %))
           (bs/get-list this :Entrezgene_comments :Gene-commentary)))
    :property-comments
    (fn [this]
      (map #(->entrezGeneComment (zip/node %))
           (bs/get-list this :Entrezgene_properties
                        :Gene-commentary)))
    :homology-comments
    (fn [this]
      (map #(->entrezGeneComment (zip/node %))
           (bs/get-list this :Entrezgene_homology
                        :Gene-commentary)))
    :locus-comments
    (fn [this]
      (map #(->entrezGeneComment (zip/node %))
           (bs/get-list this :Entrezgene_locus
                        :Gene-commentary))))
  bs/biosequenceComments
  (assoc bs/default-biosequence-comments
    :comments (fn [this] (all-comments this)))
  bs/biosequenceTranslation
  (assoc bs/default-biosequence-translation
    :trans-table (fn [this]
                   (bs/get-text (-> (bs/tax-refs this)
                                    first org-name)
                                :OrgName_gcode)))
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession (fn [this] (bs/accession (entrez-gene-track this)))
    :accessions (fn [this] (remove nil? (list (bs/accession this))))
    :version (fn [this] (bs/version (entrez-gene-track this)))
    :creation-date (fn [this]
                     (bs/creation-date (entrez-gene-track this)))
    :update-date (fn [this]
                   (bs/update-date (entrez-gene-track this)))
    :discontinue-date (fn [this]
                        (bs/discontinue-date
                         (entrez-gene-track this))))
  bs/biosequenceStatus
  (assoc bs/default-biosequence-status
    :status (fn [this] (bs/status (entrez-gene-track this))))
  bs/biosequenceSummary
  (assoc bs/default-biosequence-summary
    :summary (fn [this]
               (bs/get-text this :Entrezgene_summary)))
  bs/biosequenceTaxonomies
  (assoc bs/default-biosequence-taxonomies
    :tax-refs
    (fn [this]
      (list (org-ref (entrez-biosource this)))))
  bs/biosequenceProteins
  (assoc bs/default-biosequence-proteins
    :proteins
    (fn [this]
      (list (->entrezProtein
             (zip/node (bs/get-one this :Entrezgene_prot
                                   :Prot-ref))))))
  bs/Biosequence
  (assoc bs/default-biosequence-biosequence
    :alphabet (fn [_] :iupacNucleicAcids)
    :moltype (fn [this]
               (bs/get-text this :Entrezgene_type))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; reader
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezGeneReader [strm]
  bs/biosequenceReader
  (biosequence-seq [this]
    (let [xml (parse (:strm this) :support-dtd false)]
      (map (fn [x]
             (vec x) ;; realising all the laziness
             (->entrezGene x))
           (filter #(= (:tag %) :Entrezgene)
                   (:content xml)))))
  java.io.Closeable
  (close [this]
    (.close (:strm this))))

(defn init-gene-reader
  [strm]
  (->entrezGeneReader strm))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezgeneFile [file encoding])

(extend entrezgeneFile
  bs/biosequenceIO
  {:bs-reader (fn [this]
                (init-gene-reader (io/reader (:file this))))}
  bs/biosequenceFile
  bs/default-biosequence-file)

(defn init-entrezgene-file
  [file & {:keys [encoding] :or {encoding "UTF-8"}}]
  (->entrezgeneFile file encoding))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; connection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord entrezGeneConnection [acc-list]
  bs/biosequenceIO
  (bs-reader [this]
    (init-gene-reader
     (io/reader
      (eu/e-fetch (:acc-list this) "gene" nil "xml")))))

(defn init-entrezgene-connection
  [acc-list]
  (->entrezGeneConnection acc-list))

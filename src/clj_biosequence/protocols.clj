(in-ns 'clj-biosequence.core)

(def return-nil (fn [_] nil))

(defprotocol fastaReduce
  (fasta-reduce [this func fold]
    "Applies a function to sequence data streamed line-by-line and
    reduces the results using the supplied `fold` function. Uses the
    core reducers library so the fold function needs to have an
    'identity' value that is returned when the function is called with
    no arguments."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; readers and files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defprotocol biosequenceIO
  (bs-reader [this]
    "Returns a reader for a file containing biosequences. Use with
    `with-open'"))

(defprotocol biosequenceParameters
  (parameters [this]
    "Returns parameters from a reader."))

(defprotocol biosequenceReader
  (biosequence-seq [this]
    "Returns a lazy sequence of biosequence objects.")
  (get-biosequence [this accession]
    "Returns the biosequence object with the corresponding
    accession."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; annotations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defprotocol biosequenceGoTerms
  (gos [this]
    "Returns go term records."))

(def ^{:doc "Default implementation of biosequenceGoTerms protocol."}
  default-biosequence-goterms
  {:gos return-nil})

(defprotocol biosequenceGoTerm
  (go-id [this]
    "The GO id.")
  (go-term [this]
    "The GO term.")
  (go-component [this]
    "The GO component, molecular function etc."))

(def ^{:doc "Default implementation of biosequenceGoTerm protocol."}
  default-biosequence-goterm
  {:go-id return-nil
   :go-term return-nil
   :go-component return-nil})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; others
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defprotocol biosequenceDbRefs
  (get-db-refs [this]
    "Returns db ref records."))

(def ^{:doc "Default implementation of biosequenceDbRefs protocol."}
  default-biosequence-dbrefs
  {:get-db-refs return-nil})

(defprotocol biosequenceDbRef
  (database-name [this]
    "Returns the name of the database.")
  (object-id [this]
    "Returns the ID of an database object.")
  (db-properties [this]
    "Returns properties of the reference."))

(def ^{:doc "Default implementation of biosequenceDbRef protocol."}
  default-biosequence-dbref
  {:database-name return-nil
   :object-id return-nil
   :db-properties return-nil})

(defprotocol biosequenceNotes
  (notes [this]
    "Returns notes."))

(def ^{:doc "Default implementation of biosequenceNotes protocol."}
  default-biosequence-notes
  {:notes return-nil})

(defprotocol biosequenceUrl
  (url [this]
    "Returns the url.")
  (pre-text [this]
    "Text before anchor.")
  (anchor [this]
    "Text to show as highlight")
  (post-text [this]
    "Text after anchor"))

(def ^{:doc "Default implementation of biosequenceUrl protocol."}
  default-biosequence-url
  {:url return-nil
   :pre-text return-nil
   :anchor return-nil
   :post-text return-nil})

(defprotocol biosequenceDescription
  (description [this]
    "Returns the description of a biosequence object."))

(def ^{:doc "Default implementation of biosequenceDescription protocol."}
  default-biosequence-description
  {:description return-nil})

(defprotocol biosequenceSynonyms
  (synonyms [this]
    "Returns a list of synonyms."))

(def ^{:doc "Default implementation of biosequenceSynonyms protocol."}
  default-biosequence-synonyms
  {:synonyms return-nil})

(defprotocol biosequenceStatus
  (status [this]
    "Status of a biosequence."))

(def ^{:doc "Default implementation of biosequenceStatus protocol."}
  default-biosequence-status
  {:status return-nil})

(defprotocol biosequenceID
  (accession [this]
    "Returns the accession of a biosequence object.")
  (accessions [this]
    "Returns a list of accessions for a biosequence object.")
  (dataset [this]
    "Returns the dataset.")
  (version [this]
    "Returns the version of the accession nil if none.")
  (creation-date [this]
    "Returns a java date object.")
  (update-date [this]
    "Returns a java date object.")
  (discontinue-date [this]
    "Returns a date object."))

(def ^{:doc "Default implementation of biosequenceID protocol."}
  default-biosequence-id
  {:accession return-nil
   :accessions return-nil
   :dataset return-nil
   :version return-nil
   :creation-date return-nil
   :update-date return-nil
   :discontinue-date return-nil})

(defprotocol biosequenceTaxonomies
  (tax-refs [src]
    "Returns taxonomy records."))

(def ^{:doc "Default implementation of biosequenceTaxonomies protocol."}
  default-biosequence-taxonomies
  {:tax-refs return-nil})

(defprotocol biosequenceTaxonomy
  (tax-name [this]
    "Returns the taxonomic name.")
  (common-name [this]
    "Returns the common name.")
  (mods [this]
    "Returns a list of modifications to taxonomy.")
  (lineage [this]
    "Returns a lineage string."))

(def ^{:doc "Default implementation of biosequenceTaxonomy protocol."}
  default-biosequence-tax
  {:tax-name (fn [_] nil)
   :common-name (fn [_] nil)
   :mods (fn [_] nil)
   :lineage (fn [_] nil)})

(defprotocol biosequenceGenes
  (genes [this]
    "Returns sub-seq gene records."))

(def ^{:doc "Default implementation of biosequenceGenes protocol."}
  default-biosequence-genes
  {:genes return-nil})

(defprotocol biosequenceGene
  (locus [this]
    "Returns the gene locus.")
  (map-location [this]
    "Returns the map location.")
  (locus-tag [this]
    "Returns a locus tag.")
  (products [this]
    "The products of a gene.")
  (orf [this]
    "ORF associated with the gene."))

(def ^{:doc "Default implementation of biosequenceGene protocol."}
  default-biosequence-gene
  {:locus return-nil
   :orf return-nil
   :map-location return-nil
   :locus-tag return-nil
   :products return-nil})

(defprotocol biosequenceSummary
  (summary [this]
    "Returns the summary of a sequence."))

(def ^{:doc "Default implementation of biosequenceSummary protocol."}
  default-biosequence-summary
  {:summary return-nil})

(defprotocol biosequenceVariant
  (original [this]
    "Returns the original.")
  (variant [this]
    "Returns the variant."))

(def ^{:doc "Default implementation of biosequenceVariant protocol."}
  default-biosequence-variant
  {:original return-nil
   :variant return-nil})

(defprotocol biosequenceEvidence
  (evidence [this]
    "Returns evidence records."))

(def default-biosequence-evidence
  ^{:doc "Default implementation of biosequenceEvidence protocol."}
  {:evidence return-nil})

(defprotocol biosequenceProteins
  (proteins [this]
    "Returns protein sub-seq records."))

(def ^{:doc "Default implementation of biosequenceProteins protocol."}
  default-biosequence-proteins
  {:proteins return-nil})

(defprotocol biosequenceProtein
  (ecs [this]
    "Returns list of E.C numbers.")
  (activities [this]
    "Returns a lit of activities.")
  (processed [this]
    "Processing of the protein.")
  (calc-mol-wt [this]
    "The calculated molecular weight."))

(def ^{:doc "Default implementation of biosequenceProtein protocol."}
  default-biosequence-protein
  {:ec return-nil
   :protein-ids return-nil
   :activities return-nil
   :processed return-nil})

(defprotocol biosequenceNameObject
  (obj-name [this])
  (obj-id [this])
  (obj-description [this])
  (obj-type [this])
  (obj-value [this])
  (obj-label [this])
  (obj-heading [this])
  (obj-text [this]))

(def ^{:doc "Default implementation of biosequenceNameObject protocol."}
  default-biosequence-nameobject
  {:obj-name return-nil
   :obj-id return-nil
   :obj-description return-nil
   :obj-type return-nil
   :obj-value return-nil
   :obj-label return-nil
   :obj-heading return-nil
   :obj-text return-nil})

(defprotocol biosequenceComments
  (comments [this]
    "Returns comments.")
  (filter-comments [this value]
    "Filters comments based on the return value of obj-type."))

(def ^{:doc "Default implementation of biosequenceComments protocol."}
  default-biosequence-comments
  {:comments return-nil
   :filter-comments (fn [this value]
                      (filter #(= value (obj-type %))
                              (comments this)))})

(defprotocol biosequenceName
  (names [this]
    "Returns the default names of a record.")
  (alternate-names [this]
    "Returns the alternate names.")
  (allergen-names [this]
    "Returns the allergen names.")
  (submitted-names [this]
    "Returns the submitted names.")
  (biotech-names [this]
    "Returns the biotech names.")
  (cd-antigen-names [this]
    "Returns the cd names.")
  (innnames [this]
    "Returns the innname."))

(def ^{:doc "Default implementation of biosequenceName protocol."}
  default-biosequence-name
  {:names return-nil
   :alternate-names return-nil
   :submitted-names return-nil
   :allergen-names return-nil
   :biotech-names return-nil
   :cd-antigen-names return-nil
   :innnames return-nil})

(defprotocol Biosequence
  (bs-seq [this]
    "Returns the sequence of a biosequence as a vector.")
  (protein? [this]
    "Returns true if a protein and false otherwise.")
  (alphabet [this]
    "Returns the alphabet of a biosequence.")
  (moltype [this]
    "Returns the moltype of a biosequence.")
  (keywords [this]
    "Returns a collection of keywords."))

(def ^{:doc "Default implementation of Biosequence protocol."}
  default-biosequence-biosequence
  {:bs-seq return-nil
   :protein? (fn [_] false)
   :alphabet return-nil
   :moltype return-nil
   :keywords return-nil})

(defprotocol biosequenceCitations
  (citations [this]
    "Returns a collection of references in a sequence record.")
  (citation-key [this]
    "Returns a citation key from a record."))

(def ^{:doc "Default implementation of biosequenceCitations protocol."}
  default-biosequence-citations
  {:citations return-nil
   :citation-key return-nil})

(defprotocol biosequenceCitation
  (title [this]
    "Returns the title of a citation object.")
  (journal [this]
    "Returns the journal of a citation object.")
  (year [this]
    "Returns the year of a citation object.")
  (volume [this]
    "Returns the volume of a citation object.")
  (pstart [this]
    "Returns the start page of a citation object.")
  (pend [this]
    "Returns the end page of a citation object.")
  (authors [this]
    "Returns the authors from a citation object.")
  (pubmed [this]
    "Returns the pubmed id of a reference if there is one.")
  (crossrefs [this]
    "Returns crossrefs - DOI, ISBN etc")
  (abstract [this]
    "Returns the abstract"))

(def ^{:doc "Default implementation of biosequenceCitation protocol."}
  default-biosequence-citation
  {:title return-nil
   :journal return-nil
   :year return-nil
   :volume return-nil
   :pstart return-nil
   :pend return-nil
   :authors return-nil
   :pubmed return-nil
   :crossrefs return-nil
   :notes return-nil
   :abstract return-nil})

(defprotocol biosequenceTranslation
  (frame [this]
    "Returns the frame a sequence should be translated in.")
  (trans-table [this]
    "Returns the translation table code to be used.")
  (translation [this]
    "Returns the translation of a sequence.")
  (codon-start [this]
    "The start codon."))

(def ^{:doc "Default implementation of biosequenceTranslation protocol."}
  default-biosequence-translation
  {:frame return-nil
   :codon-start return-nil
   :trans-table return-nil
   :translation return-nil})

(defprotocol biosequenceFile
  (bs-path [this]
    "Returns the path of the file as string."))

(def ^{:doc "Default implementation of biosequenceFile protocol."}
  default-biosequence-file
  {:bs-path (fn [this] (absolute-path (:file this)))})

(defprotocol biosequenceFeatures
  (feature-seq [this]
    "Returns a lazy list of features in a sequence.")
  (filter-features [this name]
    "Returns a list of features that return 'name' when called using
    bs/bs-name from biosequenceName protocol."))

(def ^{:doc "Default implementation of biosequenceFeatures protocol."}
  default-biosequence-features
  {:feature-seq return-nil
   :filter-features return-nil})

(defprotocol biosequenceFeature
  (operator [this]
    "Returns an operator for dealing with intervals."))

(def ^{:doc "Default implementation of biosequenceFeature protocol."}
  default-biosequence-feature
  {:operator return-nil})

(defprotocol biosequenceIntervals
  (intervals [this]
    "Returns a list of intervals."))

(def ^{:doc "Default implementation of biosequenceIntervals protocol."}
  default-biosequence-intervals
  {:intervals return-nil})

(defprotocol biosequenceInterval
  (start [this]
    "Returns the start position of an interval as an integer.")
  (end [this]
    "Returns the end position of an interval as an integer.")
  (point [this]
    "Returns a point interval.")
  (comp? [this]
    "Is the interval complementary to the biosequence
    sequence. Boolean"))

(def ^{:doc "Default implementation of biosequenceInterval protocol."}
  default-biosequence-interval
  {:start return-nil
   :end return-nil
   :point return-nil
   :comp? (fn [_]
            (throw (Throwable. "No strand information.")))})

(defprotocol biosequenceSubCellLocs
  (subcell-locations [this]))

(def ^{:doc "Default implementation of biosequenceSubCellLocs protocol."}
  default-biosequence-subcells
  {:subcell-locations return-nil})

(defprotocol biosequenceSubCellLoc
  (subcell-location [this])
  (subcell-topol [this])
  (subcell-orient [this]))

(def
  ^{:doc "Default implementation of biosequenceSubCellLoc protocol."}
  default-biosequence-subcell
  {:subcell-location return-nil
   :subcell-topol return-nil
   :subcell-orient return-nil})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; nil
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(extend nil
  biosequenceDbRefs default-biosequence-dbrefs
  biosequenceDbRef default-biosequence-dbref
  biosequenceComments default-biosequence-comments
  biosequenceUrl default-biosequence-url
  biosequenceDescription default-biosequence-description
  biosequenceSynonyms default-biosequence-synonyms
  biosequenceID default-biosequence-id
  biosequenceTaxonomy default-biosequence-tax
  biosequenceGenes default-biosequence-genes
  biosequenceGene default-biosequence-gene
  biosequenceSummary default-biosequence-summary
  biosequenceProteins default-biosequence-proteins
  biosequenceProtein default-biosequence-protein
  Biosequence default-biosequence-biosequence
  biosequenceCitations default-biosequence-citations
  biosequenceCitation default-biosequence-citation
  biosequenceTranslation default-biosequence-translation
  biosequenceFeatures default-biosequence-features
  biosequenceFeature default-biosequence-feature
  biosequenceIntervals default-biosequence-interval
  biosequenceInterval default-biosequence-interval
  biosequenceName default-biosequence-name
  biosequenceNotes default-biosequence-notes
  biosequenceVariant default-biosequence-variant
  biosequenceEvidence default-biosequence-evidence
  biosequenceStatus default-biosequence-status
  biosequenceTaxonomies default-biosequence-taxonomies
  biosequenceNameObject default-biosequence-nameobject
  biosequenceSubCellLocs default-biosequence-subcells
  biosequenceSubCellLoc default-biosequence-subcell
  biosequenceGoTerms default-biosequence-goterms
  biosequenceGoTerm default-biosequence-goterm)

(ns clj-biosequence.genbank
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :refer [split]]
            [clj-http.client :as client]
            [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.core :as bs]))

(declare genbank-search-helper genbank-sequence-helper moltype get-genbank-stream check-db check-rt)

; interval

(defrecord genbankInterval [src]

  bs/biosequenceInterval

  (start [gb-interval]
    (let [r (or (zf/xml1-> (zip/xml-zip (:src gb-interval))
                           :GBInterval_from
                           zf/text)
                (zf/xml1-> (zip/xml-zip (:src gb-interval))
                           :GBInterval_point
                           zf/text))]
      (if r
        (Integer/parseInt r)
        (throw (Throwable. "No start value in interval!")))))

  (end [gb-interval]
    (let [r (or (zf/xml1-> (zip/xml-zip (:src gb-interval))
                           :GBInterval_to
                           zf/text)
                (zf/xml1-> (zip/xml-zip (:src gb-interval))
                           :GBInterval_point
                           zf/text))]
      (if r
        (Integer/parseInt r)
        (throw (Throwable. "No end value in interval!")))))

  (comp? [gb-interval]
    (let [c (zf/xml1-> (zip/xml-zip (:src gb-interval))
                       :GBInterval_iscomp
                       (zf/attr :value))]
      (if (= "true" c) true false))))

;; qualifier

(defrecord genbankQualifier [src])

(defn qualifier-name
  [qual]
  (zf/xml1-> (zip/xml-zip (:src qual)) :GBQualifier_name zf/text))

(defn qualifier-value
  [qual]
  (zf/xml1-> (zip/xml-zip (:src qual)) :GBQualifier_value zf/text))

; feature

(defrecord genbankFeature [src]

  bs/biosequenceFeature

  (feature-type [feature]
    (zf/xml1-> (zip/xml-zip (:src feature))
               :GBFeature_key
               zf/text))

  (interval-seq [gb-feature]
    (map #(->genbankInterval (zip/node %))
         (zf/xml-> (zip/xml-zip (:src gb-feature))
                   :GBFeature_intervals
                   :GBInterval))))

(defn qualifier-seq
  [feature]
  (->> (:content (some #(if (= (:tag %) :GBFeature_quals) %)
                       (:content (:src feature))))
       (map #(->genbankQualifier %))))

(defn feature-operator
  [feat]
  (zf/xml1-> (zip/xml-zip (:src feat)) :GBFeature_operator zf/text))

(defn feature-location
  "Returns the value corresponding to the 'GBFeature_location' element of a 
   genbank feature element."
  [feature]
  (zf/xml1-> (zip/xml-zip (:src feature))
             :GBFeature_location
             zf/text))
; sequence

(defrecord genbankSequence [src]

  bs/Biosequence

  (feature-seq [this]
    (map #(->genbankFeature %)
         (:content (some #(if (= (:tag %) :GBSeq_feature-table)
                            %) (:content (:src this))))))

  (accession [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :GBSeq_primary-accession
               zf/text))

  (accessions [this]
    (zf/xml-> (zip/xml-zip (:src this))
              :GBSeq_other-seqids
              :GBSeqid
              zf/text))

  (def-line [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :GBSeq_definition
               zf/text))

  (protein? [this]
    (= (bs/alphabet this) :iupacAminoAcids))

  (bs-seq [this]
    (bs/clean-sequence
     (zf/xml1-> (zip/xml-zip (:src this))
                :GBSeq_sequence
                zf/text)
     (bs/alphabet this)))

  (fasta-string [this]
    (str ">gb|" (bs/accession this) "|"
         (apply str (interpose "|" (bs/accessions this)))
         "| " (bs/def-line this) \newline
         (bs/bioseq->string this) \newline))

  (alphabet [this]
    (cond (#{"genomic" "precursor RNA" "mRNA" "rRNA" "tRNA" "snRNA" "scRNA"
             "other-genetic" "DNA" "cRNA" "snoRNA" "transcribed RNA"} (moltype this))
          :iupacNucleicAcids
          (#{"AA"} (moltype this))
          :iupacAminoAcids
          :else
          (throw (Throwable. (str "Unknown moltype: " (moltype this)))))))

; IO

(defrecord genbankReader [strm]

  bs/Biosequence

  (feature-seq [this]
    (let [xml (xml/parse (:strm this) :support-dtd false)]
      (->> (:content (first (filter #(= (:tag %) :GBSeq_feature-table)
                                    (:content (first (:content xml))))))
           (filter #(= (:tag %) :GBFeature))
           (map #(->genbankFeature %)))))

  bs/biosequenceReader

  (biosequence-seq [this]
    (let [xml (xml/parse (:strm this) :support-dtd false)]
      (map (fn [x]
             (vec x) ;; realising all the laziness
             (->genbankSequence x))
           (filter #(= (:tag %) :GBSeq)
                   (:content xml)))))

  (parameters [this]
    ())

  java.io.Closeable

  (close [this]
    (.close (:strm this))))

(defn init-genbank-reader
  [strm]
  (->genbankReader strm))

(defrecord genbankFile [file]
  
  bs/biosequenceIO

  (bs-reader [this]
    (init-genbank-reader (io/reader (:file this))))

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this))))

(defrecord genbankString [str]

  bs/biosequenceIO

  (bs-reader [this]
    (init-genbank-reader (io/reader (:str this)))))

(defrecord genbankConnection [acc-list db retype]

  bs/biosequenceIO

  (bs-reader [this]
    (let [s (get-genbank-stream (:acc-list this)
                                (:db this)
                                (:retype this))]
      (condp = (:retype this)
        :xml (init-genbank-reader (io/reader s))
        :fasta (bs/init-fasta-reader (io/reader s)
                                     (cond (#{:nucest :nuccore :nucgss :popset} db)
                                           :iupacNucleicAcids
                                           (= :protein db)
                                           :iupacAminoAcids))))))

(defn init-genbank-file
  [file]
  (if (fs/file? file)
    (->genbankFile file)
    (throw (Throwable. (str "File not found: " file)))))

(defn init-genbank-string
  [str]
  (->genbankString str))

(defn init-genbank-connection
  [accessions db retype]
  (and (check-db db)
       (check-rt retype)
       (->genbankConnection accessions db retype)))

;; web

(defn genbank-search
  "Returns a non-lazy list of result ids from NCBI for a particular search term and
   database. Search term syntax is the same as is used at the NCBI. For example,
   to retreive all Schistosoma mansoni protein sequences the 'db' argument would
   be 'protein' and the search term 'txid6183[Organism:noexp]'. Database arguments
   are restricted to the following keyword arguments and any other value will cause
   an error:
   :protein    - a collection of sequences from several sources, including 
                 translations from annotated coding regions in GenBank, RefSeq 
                 and TPA, as well as records from SwissProt, PIR, PRF, and PDB.
   :nucest     - a collection of short single-read transcript sequences from GenBank.
   :nuccore    - a collection of sequences from several sources, including GenBank, 
                 RefSeq, TPA and PDB.
   :nucgss     - a collection of unannotated short single-read primarily genomic
                 sequences from GenBank including random survey sequences clone-end
                 sequences and exon-trapped sequences.
   :popset     - a set of DNA sequences that have been collected to analyse the 
                 evolutionary relatedness of a population. The population could 
                 originate from different members of the same species, or from 
                 organisms from different species. 

   Note that when retrieving popset entries multiple sequences are returned from a
   single accession number. Returns an empty list if no matches found."
  ([term db] (genbank-search term db 0 nil))
  ([term db restart key]
     (if (check-db db)
         (let [r (genbank-search-helper term db restart key)
               k (zf/xml1-> (zip/xml-zip r) :WebEnv zf/text)
               c (Integer/parseInt (zf/xml1-> (zip/xml-zip r) :Count zf/text))]
           (if (> restart c)
             nil
             (lazy-cat (zf/xml-> (zip/xml-zip r) :IdList :Id zf/text)
                       (genbank-search term db (+ restart 1000) k))))
         (throw (Throwable. (str "'" db "' " "not allowed. Only :protein, :nucleotide, :nucest, :nuccore, :nucgss and :popset are acceptable database arguments. See the documentation for 'genbank-search' for an explanation of the different databases."))))))

;; convenience functions

(defn created
  "Returns the creation date of the sequence. Content of the GBSeq_create-date tag."
  [this]
  (zf/xml1-> (zip/xml-zip (:src this))
             :GBSeq_create-date zf/text))

(defn modified
  "Returns modification date. Content of the GBSeq_update-date tag."
  [this]
  (zf/xml1-> (zip/xml-zip (:src this))
             :GBSeq_update-date zf/text))

(defn version
  "Version number. Content of GBSeq_accession-version tag."
  [this]
  (Integer. (second (re-find #".+\.(\d+)"
                             (zf/xml1-> (zip/xml-zip (:src this))
                                        :GBSeq_accession-version zf/text)))))

(defn taxid
  "NCBI taxonomic id of the organism from which the sequence is derived."
  [this]
  (Integer/parseInt
   (second (split (->> (filter #(= "source" (bs/feature-type %)) (bs/feature-seq this))
                       first
                       (qualifier-seq)
                       (filter #(= "db_xref" (qualifier-name %)))
                       first
                       qualifier-value)
                  #":"))))

(defn taxonomy
  "Returns a list of the taxonomy of the organism from which the
  sequence is derived."
  [this]
  (seq (split (zf/xml1-> (zip/xml-zip (:src this))
                                :GBSeq_taxonomy zf/text)
                     #";")))

(defn org-scientific-name
  "Scientific name of the organism from which the sequence is derived."
  [this]
  (zf/xml1-> (zip/xml-zip (:src this))
             :GBSeq_organism
             zf/text))

(defn moltype
  "The type of the sequence. Content of the GBSeq_moltype tag."
  [gbseq]
  (zf/xml1-> (zip/xml-zip (:src gbseq))
             :GBSeq_moltype
             zf/text))

(defn gb-locus
  "Returns the locus. Contents of GBSeq_locus tag."
  [gbseq]
  (zf/xml1-> (zip/xml-zip (:src gbseq)) :GBSeq_locus zf/text))

;private

(defn- check-db
  [db]
  (if (not (#{:nucest :nuccore :nucgss :popset :protein} db))
    (throw (Throwable. (str "DB not supported: "
                            db
                            ". Only :protein, :nucest, :nuccore, :nucgss and :popset are supported.")))
    true))

(defn- check-rt
  [rt]
  (if (not (#{:xml :fasta} rt))
    (throw (Throwable. (str rt " not supported. "
                            "Only :xml and :fasta are allowed retype values.")))
    true))

(defn- get-genbank-stream
  [a-list db rettype]
  (if (empty? a-list)
    nil
    (let [r (client/post "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                         {:query-params
                          {:db (name db)
                           :id (apply str (interpose "," a-list))
                           :rettype (condp = rettype
                                      :xml "gb"
                                      :fasta "fasta")
                           :retmode (condp = rettype
                                      :xml "xml"
                                      :fasta "text")}
                          :as :stream})]
      (:body r))))

(defn- genbank-search-helper
  ([term db retstart] (genbank-search-helper term db retstart nil))
  ([term db retstart key]
   (xml/parse-str
    (:body
     (client/get
      (str "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db="
           (java.net.URLEncoder/encode (name db))
           "&term="
           (java.net.URLEncoder/encode term)
           "&retmax=" 1000
           "&retstart=" retstart
           (if key (str "&WebEnv=" key))
           "&usehistory=y"))))))

(defn- extract-features-genbank
  [rdr feature]
  "Returns a lazy list of features, corresponding to 'feature' from a file handle opened on a Genbank xml file."
  (filter #(and (= (:tag %) :GBFeature)
                (= (zf/xml1-> (zip/xml-zip %) :GBFeature_key zf/text) feature))
          (:content (some #(if (= (:tag %) :GBSeq_feature-table)
                             %)
                          (:content (first (:content (xml/parse rdr))))))))

(defn- gb-sequence
  [rdr]
  (:content (some #(if (= (:tag %) :GBSeq_sequence)
                     %)
                  (:content (first (:content (xml/parse rdr)))))))

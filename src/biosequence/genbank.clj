(ns biosequence.genbank
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clj-http.client :as client]
            [clojure.java.io :as io]
            [clojure.pprint :as pp]
            [fs.core :as fs]
            [clojure.java.jdbc :as sql]
            [biosequence.core :as bs]))

(declare qualifier-extract init-genbank-store feature-seq genbank-search-helper genbank-sequence-helper)

; interval

(defrecord genbankInterval [src])

(defn start
  "Start index of a genbankInterval."
  [gb-interval]
  (let [r (or (zf/xml1-> (zip/xml-zip (:src gb-interval))
                         :GBInterval_from
                         zf/text)
              (zf/xml1-> (zip/xml-zip (:src gb-interval))
                         :GBInterval_point
                         zf/text))]
    (if r
      (read-string r)
      (throw (Throwable. "No start value in interval!")))))

(defn end
  "End index of a genbankInterval."
  [gb-interval]
  (let [r (or (zf/xml1-> (zip/xml-zip (:src gb-interval))
                         :GBInterval_to
                         zf/text)
              (zf/xml1-> (zip/xml-zip (:src gb-interval))
                         :GBInterval_point
                         zf/text))]
    (if r
      (read-string r)
      (throw (Throwable. "No end value in interval!")))))

(defn comp?
  "Is a genbankInterval complementary?"
  [gb-interval]
  (let [c (zf/xml1-> (zip/xml-zip (:src gb-interval))
                     :GBInterval_iscomp
                     (zf/attr :value))]
    (if (= "true" c) true false)))

(defn get-interval-sequence
  "Returns a fastaSequence object containing the sequence specified in the genbankInterval from a genbankSequence. Designed for getting fastasequences by applying  genbankInterval to the sequence entry it originates from."
  [gb-interval gb-sequence]
  (let [dna (bs/sequence-string gb-sequence) 
        start (start gb-interval)
        end (end gb-interval)]
    (bs/init-fasta-sequence
     (bs/accession gb-sequence)
     (str (bs/def-line gb-sequence) " [" start "-" end "]")
     (if (bs/protein? gb-sequence) :protein :nucleotide)
     (if (false? (comp? gb-interval))
       (subs dna (- start 1) end)
       (bs/revcom-dna-string (subs dna (- end 1) start))))))

; feature

(defrecord genbankFeature [src])

(defn feature-type
  "Returns the feature key. For example: protein, Region, Site, CDS etc."
  [feature]
  (zf/xml1-> (zip/xml-zip (:src feature))
             :GBFeature_key
             zf/text))

(defn qualifier-extract
  "Takes a genbankFeature and returns the value specified by 'element'. For instance to retrieve the organism name from a 'source' feature arguments would be the feature object and 'organism'."
  [feature element]
  (let [quals (:content (some #(if (= (:tag %) :GBFeature_quals)
                                 %)
                              (:content (:src feature))))]
    (some #(let [z (zip/xml-zip %)]
             (if (= (zf/xml1-> z :GBQualifier_name zf/text)
                    element)
               (zf/xml1-> z :GBQualifier_value zf/text)))
          quals)))

(defn interval-seq
  "Returns a non-lazy list of intervals in a feature. It also calculates a frame for each interval so that individual translations of a DNA interval provide the correct protein sequence. This value is accessed through the :frame keyword for each interval. These values will have no meaning if the sequence is a protein."
  [gb-feature]
  (let [ints (map #(->genbankInterval %)
                  (:content (zip/node (zf/xml1-> (zip/xml-zip (:src gb-feature))
                                                 :GBFeature_intervals))))]
    (loop [i ints
           f 1
           acc nil]
      (if (empty? i)
        (reverse acc)
        (let [s (start (first i))
              e (end (first i))
              aint (- (+ 1 (if-not (comp? (first i))
                            (- e s)
                            (- s e)))
                      f)
              frame (cond (< aint 0)
                          (* -1 aint)
                          (< aint 3)
                          (- 3 aint)
                          :else
                          (- 3 (mod aint 3)))]
          (recur (rest i) frame
                 (cons (assoc (first i) :frame f) acc)))))))

(defn get-feature-sequence
  "Returns a fastaSequence object containing the sequence specified in a genbankFeature object from a genbankSequence object. Designed for applying intervals to the sequence entry they originate from."
  [gb-feat gbseq]
  (let [intervals (interval-seq gb-feat)]
    (bs/init-fasta-sequence
     (bs/accession gbseq)
     (str (bs/def-line gbseq) " - Feature: " (feature-type gb-feat)
          " - [" (start (first intervals)) "-" (end (last intervals)) "]")
     (if (bs/protein? gbseq) :protein :nucleotide)
     (apply str
            (map #(if (comp? %)
                    (bs/revcom-dna-string
                     (subs (bs/sequence-string gbseq)
                           (- (end %) 1)
                           (start %)))
                    (subs (bs/sequence-string gbseq)
                          (- (start %) 1)
                          (end %))) intervals)))))

; sequence

(defrecord genbankSequence [accession src])

(extend-protocol bs/Biosequence

  genbankSequence

  (accession [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :GBSeq_primary-accession
               zf/text))

  (accessions [this]
    (zf/xml-> (zip/xml-zip (:src this))
              :GBSeq_other-seqids
              :GBSeqid
              zf/text))

  (created [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :GBSeq_create-date zf/text))

  (modified [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :GBSeq_update-date zf/text))

  (version [this]
    (Integer. (second (re-find #".+\.(\d+)"
                               (zf/xml1-> (zip/xml-zip (:src this))
                                          :GBSeq_accession-version zf/text)))))

  (taxid [this]
    (Integer. (second
               (string/split
                (qualifier-extract (first (filter #(= (feature-type %) "source")
                                                  (feature-seq this)))
                                   "db_xref")
                #":"))))

  (taxonomy [this]
    (seq (string/split (zf/xml1-> (zip/xml-zip (:src this))
                                  :GBSeq_taxonomy zf/text)
                       #";")))

  (org-scientific-name [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :GBSeq_organism
               zf/text))

  (def-line [this]
    (zf/xml1-> (zip/xml-zip (:src this))
              :GBSeq_definition
              zf/text))

  (protein? [this]
    (= "AA"
       (zf/xml1-> (zip/xml-zip (:src this))
                  :GBSeq_moltype
                  zf/text)))

  (sequence-string [this]
    (string/upper-case (apply str (remove #{\space\newline}
                                          (zf/xml1-> (zip/xml-zip (:src this))
                                                     :GBSeq_sequence
                                                     zf/text)))))

  (fasta-string [this]
    (str ">gb|" (bs/accession this) "|"
         (apply str (bs/accessions this))
         "| " (bs/def-line this) \newline
         (bs/sequence-string this)))

  (fasta-to-stream [this wrt]
    (pp/cl-format wrt "~A~%" (bs/fasta-string this))))

(defn moltype
  [gbseq]
  (zf/xml1-> (zip/xml-zip (:src gbseq))
             :GBSeq_moltype
             zf/text))

(defn gb-locus
  "Returns the locus of a genbankSequence."
  [gbseq]
  (zf/xml1-> (zip/xml-zip (:src gbseq)) :GBSeq_locus zf/text))

(defn feature-seq
  "Returns a lazy list of features from a genbankSequence."
  [gbseq]
  (map #(->genbankFeature %)
       (:content (some #(if (= (:tag %) :GBSeq_feature-table)
                          %) (:content (:src gbseq))))))

; file

(defrecord genbankFile [file])

(extend-protocol bs/biosequenceFile

  genbankFile

  (biosequence-seq-file [this rdr]
    (let [xml (xml/parse rdr)]
      (if (= (:tag xml) :GBSet)
        (map #(let [a (zf/xml1-> (zip/xml-zip %)
                                 :GBSeq_primary-accession
                                 zf/text)]
                (->genbankSequence a %))
             (filter #(= (:tag %) :GBSeq)
                     (:content xml)))
        (throw (Throwable. "Doesn't appear to be a genbank XML file.")))))
  
  (file-path [this]
    (:file this)))

(defn init-genbank-file
  [file]
  (->genbankFile file))

(defmacro with-genbank-features-in-file
  "For use with very large sequences such as genomes."
  [[gbseq handle feature] & code]
  `(let [fl# (extract-features-genbank rdr# ~feature)
           ~handle (map #(->genbankFeature %) fl#)]
     ~@code))

;; persistance

(defrecord genbankStore [file])

(defn index-genbank-file
  "Indexes a genbankFile object and returns a uniprotStore object. Duplicate primary accessions in the file will cause a 'Unique index or primary key violation'."
  ([genbankfile] (index-genbank-file genbankfile false))
  ([genbankfile memory]
     (let [st (init-genbank-store genbankfile memory)]
       (bs/with-connection-to-store [st]
         (bs/with-biosequences-in-file [l genbankfile]
           (dorun (map #(bs/save-object %) l))))
       st)))

(defn load-genbank-store
  [dir]
  (let [file (first (fs/glob (str dir "/" "*.h2.db")))
        db-file (second (re-find  #"(.+)\.h2.db" (fs/absolute-path file)))]
    (if (not (nil? db-file))
      (assoc (->genbankStore db-file) :db (bs/make-db-connection db-file true))
      (throw (Throwable. "DB file not found!")))))

;; web

(defn wget-genbank-sequence
  "Returns a 'semi-lazy' (1000 sequences are fetched at a time) list of genbankSequences corresponding to the 'accessions' list argument. Accepts any identification that genbank accepts but will return an exception (status code 400) if the id is not recognised. Can return genbankSequences or fastaSequence objects depending whether the 'class' argument is :xml or :fasta respectively. The 'db' argument has the same format and allowed values as the wget-genbank-search function."
  [accessions db rettype]
  (if (some #(= db %) '(:protein :nucleotide :nucest :nuccore :nucgss :popset))
    (if-not (some #(= rettype %)
                   '(:xml :fasta))
      (throw (Throwable. (str rettype
                              " not allowed. "
                              "Only :xml and :fasta are allowed rettype values.")))
      (genbank-sequence-helper (partition-all 1000 (if (coll? accessions)
                                                     accessions
                                                     (list accessions)))
                               db rettype))
    (throw (Throwable. (str "'" db "' " "not allowed. Only :protein, :nucleotide, :nucest, :nuccore, :nucgss and :popset are acceptable database arguments. See the documentation for 'wget-genbank-search' for an explanation of the different databases.")))))

(defn wget-genbank-search
  "Returns a 'semi-lazy' (in that it fetches 1000 results at a time) list of result ids from NCBI for a particular search term and database. Search term syntax is the same as is used at the NCBI. For example, to retreive all Schistosoma mansoni protein sequences the 'db' argument would be 'protein' and the search term 'txid6183[Organism:noexp]'. Database arguments are restricted to the following keyword arguments and any other value will cause an error:
   :protein    - a collection of sequences from several sources, including translations 
                 from annotated coding regions in GenBank, RefSeq and TPA, as well as 
                 records from SwissProt, PIR, PRF, and PDB.
   :nucleotide - a collection of sequences from several sources, including GenBank, 
                 RefSeq, TPA and PDB.
   :nucest     - a collection of short single-read transcript sequences from GenBank.
   :nuccore    - a collection of sequences from several sources, including GenBank, 
                 RefSeq, TPA and PDB (same as :nucleotide).
   :nucgss     - a collection of unannotated short single-read primarily genomic sequences
                 from GenBank including random survey sequences clone-end sequences and 
                 exon-trapped sequences.
   :popset     - a set of DNA sequences that have been collected to analyse the evolutionary
                 relatedness of a population. The population could originate from different
                 members of the same species, or from organisms from different species. 

Note that when retrieving popset entries multiple sequences are returned from a single accession number. Returns an empty list if no matches found."
  ([term db] (wget-genbank-search term db 0 nil))
  ([term db restart key]
     (if (some #(= db %) '(:protein :nucleotide :nucest :nuccore :nucgss :popset))
         (let [r (genbank-search-helper term db restart key)
               k (zf/xml1-> (zip/xml-zip r) :WebEnv zf/text)
               c (read-string (zf/xml1-> (zip/xml-zip r) :Count zf/text))]
           (if (> restart c)
             nil
             (lazy-cat (zf/xml-> (zip/xml-zip r)
                                 :IdList
                                 :Id
                                 zf/text)
                       (wget-genbank-search term db (+ restart 1000) k))))
         (throw (Throwable. (str "'" db "' " "not allowed. Only :protein, :nucleotide, :nucest, :nuccore, :nucgss and :popset are acceptable database arguments. See the documentation for 'wget-genbank-search' for an explanation of the different databases."))))))

;private

(defn- genbank-sequence-helper
  [a-list db class]
  (if (empty? a-list)
    nil
    (let [r (client/post "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                         {:query-params
                          {:db (name db)
                           :id (apply str (interpose "," (first a-list)))
                           :rettype (condp = class
                                      :xml "gb"
                                      :fasta "fasta")
                           :retmode (condp = class
                                      :xml "xml"
                                      :fasta "text")}})
          s (condp = class
              :xml (filter #(= (:tag %) :GBSeq)
                           (:content
                            (xml/parse-str (:body r))))
              :fasta (bs/fasta-seq-string (:body r) (if (= db :protein)
                                                      :protein
                                                      :nucleotide)))]
      (lazy-cat (map #(condp = class
                        :xml (->genbankSequence (zf/xml1-> (zip/xml-zip %)
                                                           :GBSeq_primary-accession
                                                           zf/text)
                                                %)
                        :fasta %) 
                     s)
                (genbank-sequence-helper (rest a-list) db class)))))

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

(defn- init-genbank-store
  [file memory]
  (if memory
    (bs/init-in-mem-store (->genbankStore (bs/file-path file)))
    (bs/init-store
     (->genbankStore (bs/index-file-name (bs/file-path file))))))

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




(ns clj-biosequence.core
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-http.client :as client]
            [clojure.string :as string]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.persistence :as ps]))

(declare init-fasta-store init-fasta-sequence translate)

(defprotocol biosequenceIO
  (bs-reader [this]))

(defprotocol biosequenceReader
  (biosequence-seq [this]))

(defprotocol Biosequence
  (accession [this]
    "Returns the accession of a biosequence.")
  (accessions [this]
    "Returns a list of strings describing the accessions of a biosequence object.")
  (def-line [this]
    "Returns a description for a biosequence object.")
  (bs-seq [this]
    "Returns the sequence of a biosequence as a vector.")
  (fasta-string [this]
    "Returns the biosequence as a string in fasta format.")
  (protein? [this]
    "Returns true if a protein and false otherwise.")
  (alphabet [this]
    "Returns the alphabet of a biosequence.")
  (reverse-seq [this]
    "Returns a new fastaBiosequence with the reversed sequence of the original.")
  (reverse-comp [this]
    "Returns a new fastaBiosequence with the reverse complement of the original."))

(defprotocol biosequenceStoreDir
  (load-store [this dbfile]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; store

(defn save-biosequence
  "Saves an object to a store."
  [obj]
  (ps/save-object (accession obj) obj))

(defmacro with-biosequences-in-store
  "Provides a handle to a lazy list of biosequences in an object store."
  [[handle store] & code]
  `(ps/with-objects [~handle ~store]
     ~@code))

(defn update-biosequence 
  "Updates an object in the store with a current connection."
  [obj]
  (ps/update-object (accession obj) obj))

(defn update-biosequence-by-accession
  "Takes an accession number and key value pairs. If the biosequence exists in
   the current store it will be updated with the key value pairs and saved. 
   Throws an exception if a corresponding object is not found in the store."
  [store accession & args]
  (ps/update-object-by-id store accession args))

(defn get-biosequence [accession]
  "Returns a biosequence object from a store implementing the biosequence"
  (ps/get-object accession))

(defn index-biosequence-file
  [file store]
  (ps/with-store [store]
    (with-open [^java.io.BufferedReader rdr (bs-reader file)]
      (dorun
       (pmap #(save-biosequence %)
             (biosequence-seq rdr))))
    store))

(defmacro with-biosequences
  [[handle store] & code]
  `(ps/with-objects [~handle ~store]
     ~@code))

(defmacro with-biosequence-store
  [[store] & code]
  `(ps/with-store [store]
     ~@code))

(defn load-biosequence-store
  "Loads a fastaStore."
  [dir]
  (let [file (first (fs/glob (str (:dir dir) "/" "*.h2.db")))
        db-file (second (re-find  #"(.+)\.h2.db" (fs/absolute-path file)))
        d (load-store dir db-file)]
    (if (not (nil? db-file))
      (assoc d :db
             (ps/make-db-connection db-file false))
      (throw (Throwable. "DB file not found!")))))

; other

(defn bioseq->string
  [bs]
  (apply str (bs-seq bs)))

(defn residue-frequencies
  [bs]
  (frequencies (bs-seq bs)))

(defn sub-bioseq
  "Returns a new fasta sequence object with the sequence corresponding to
   'beg' (inclusive) and 'end' (exclusive) of 'bs'. If no 'end' argument 
   returns from 'start' to the end of the sequence. Indexes start at zero."
  ([bs beg] (sub-bioseq bs beg nil))
  ([bs beg end]
     (init-fasta-sequence (accession bs)
                          (str (def-line bs)
                               "[" beg " - "
                               (if end end "End") "]")
                          (alphabet bs)
                          (if end
                            (subvec (bs-seq bs) beg end)
                            (subvec (bs-seq bs) beg)))))

(defn partition-bioseq
  "Partitions a sequence into a lazy list of lists of 'n' size. Default
   partition size is 3."
  ([bs] (partition-bioseq bs 3))
  ([bs n]
     (partition-all n (bs-seq bs))))

(defn concat-bioseqs
  [s1 s2]
  (if (= (alphabet s1) (alphabet s2))
    (init-fasta-sequence (str (accession s1) "/" (accession s2))
                         (str (def-line s1) "/" (def-line s2))
                         (alphabet s1)
                         (vec (concat (bs-seq s1) (bs-seq s2))))
    (throw (Throwable. "Incompatible alphabets for concatenation of biosequence."))))

(defn translate
  "Returns a fastaSequence object corresponding to the protein translation 
   of the sequence in the specified frame."
  ([bs frame] (translate bs frame (ala/codon-tables 1)))
  ([bs frame table]
     (cond (protein? bs)
           (throw (Throwable. "Can't translate a protein sequence!"))
           (not (#{1 2 3 4 5 6 -1 -2 -3} frame))
           (throw (Throwable. "Invalid frame."))
           :else
           (init-fasta-sequence (str (accession bs) "-"  frame)
                                (str (def-line bs) " - Translated frame: " frame)
                                :iupacAminoAcids
                                (let [v (cond (#{1 2 3} frame)
                                              (sub-bioseq bs (- frame 1))
                                              (#{-1 -2 -3} frame)
                                              (-> (reverse-comp bs)
                                                  (sub-bioseq (- (* -1 frame) 1)))
                                              (#{4 5 6} frame)
                                              (-> (reverse-comp bs)
                                                  (sub-bioseq ( - frame 4))))]
                                  (vec (map #(ala/codon->aa % table)
                                            (partition-bioseq v 3))))))))

(defn six-frame-translation
  "Returns a lazy list of fastaSequence objects representing translations of
   a nucleotide biosequence object in six frames."
  ([nucleotide] (six-frame-translation nucleotide (ala/codon-tables 1)))
  ([nucleotide table]
     (map #(translate nucleotide % table)
          '(1 2 3 -1 -2 -3))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; an implementation of biosequence for fasta sequences
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastaSequence [accession description alphabet sequence]

  Biosequence
  
  (accession [this]
    (:accession this))

  (accessions [this]
    (list (:accession this)))

  (bs-seq [this]
    (:sequence this))

  (def-line [this]
    (:description this))

  (protein? [this]
    (if (= :iupacAminoAcids (:alphabet this))
      true
      false))

  (reverse-comp [this]
    (if (protein? this)
      (throw (Throwable. "Can't reverse/complement a protein sequence."))
      (init-fasta-sequence (accession this)
                           (str (def-line this) " - reverse-comp")
                           (alphabet this)
                           (ala/revcom (bs-seq this)))))

  (reverse-seq [this]
    (init-fasta-sequence (accession this)
                         (str (def-line this) " - reverse")
                         (alphabet this)
                         (vec (reverse (bs-seq this)))))

  (fasta-string [this]
    (if (:description this)
      (str ">" (:accession this) " " (def-line this) "\n" (bioseq->string this) "\n")))

  (alphabet [this]
    (:alphabet this)))

(defn init-fasta-sequence
  "Returns a fastaSequence object with the specified information. Alphabet can be
   one of :iupacNucleicAcids or :iupacAminoAcids."
  [accession description alphabet sequence]
  (->fastaSequence accession description alphabet sequence))

;; IO

(defrecord fastaReader [strm alphabet]

  biosequenceReader
  
  (biosequence-seq [this]
    (map (fn [[d s]]
           (let [seqs (apply str s)]
             (cond (not (re-find #"^>" (first d)))
                   (throw (Throwable. (str "Data corrupted at " (first d))))
                   (> (count d) 1)
                   (throw (Throwable. (str "No seqeunce for entry " (first d))))
                   :else
                   (init-fasta-sequence (second (re-find #"^>([^\s]+)" (first d)))
                                        (second (re-find #">[^\s]+\s+(.+)" (first d)))
                                        (:alphabet this)
                                        (->> seqs
                                             vec
                                             (remove (complement
                                                      (-> (ala/alphabet (:alphabet this))
                                                          keys
                                                          set)))
                                             vec)))))
         (partition 2 (partition-by #(re-find #"^>" %) (line-seq (:strm this))))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defn init-fasta-reader
  [strm alphabet]
  (->fastaReader strm alphabet))

(defrecord fastaFile [file alphabet]

  biosequenceIO

  (bs-reader [this]
    (init-fasta-reader (io/reader (:file this))
                       (:alphabet this))))

(defrecord fastaString [str alphabet]

  biosequenceIO

  (bs-reader [this]
    (init-fasta-reader (io/reader (:str this))
                       (alphabet this))))

(defn init-fasta-file
  "Initialises fasta protein file. Accession numbers and description are 
   processed by splitting the string on the first space, the accession 
   being the first value and description the second."
  [path alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (if (fs/exists? path)
      (->fastaFile path alphabet)
      (throw (Throwable. (str "File not found: " path))))))

(defn init-fasta-string
  "Initialises a fasta string."
  [str alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (->fastaString str alphabet)))

;; persistence

(defrecord fastaStore [file])

(defn index-fasta-file
  "Indexes a fastaFile object and returns a fastaStore object."
  [fastafile]
  (let [st (ps/init-store (->fastaStore 
                           (ps/index-file-name (:file fastafile))))]
    (index-biosequence-file fastafile st)))

(defrecord fastaStoreDir [dir]

  biosequenceStoreDir

  (load-store [this dbfile]
    (->fastaStore dbfile)))

(defn load-fasta-store
  "Loads a fastaStore."
  [dir]
  (load-biosequence-store (->fastaStoreDir dir)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; id mapping
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn id-convert
  "Takes either a single accession or a list of accessions and returns a hash-map
   mapping the accession numbers to the corresponding identification number in the 
   specified 'to' database. 'From' database also needs to be specified. If not 
   found returns an empty hash-map. Uses the Uniprot id mapping utility and a list 
   of supported databases is supplied at http://www.uniprot.org/faq/28#id_mapping_examples.
   Some common mappings include:
   DB Name                  Abbreviation     Direction
   UniProtKB AC/ID	    ACC+ID	     from
   UniProtKB AC	            ACC              to
   EMBL/GenBank/DDBJ	    EMBL_ID	     both
   EMBL/GenBank/DDBJ CDS    EMBL	     both
   Entrez Gene (GeneID)     P_ENTREZGENEID   both
   GI number	            P_GI	     both
   RefSeq Protein	    P_REFSEQ_AC	     both
   RefSeq Nucleotide	    REFSEQ_NT_ID     both
   WormBase	            WORMBASE_ID	     both

   There is a 100,000 limit on accessions in a single query imposed by Uniprot."
  [ids from to email]
  (let [i (if (seq? ids)
            ids (list ids))]
    (if (<= (count ids) 100000)
      (let [param {:from from :to to :format "tab" 
                   :query (apply str (doall (interpose "," (if (list? ids)
                                                             ids
                                                             (list ids)))))}
            address "http://www.uniprot.org/mapping/"
            r (client/post address 
                           {:client-params
                            {"http.useragent" (str "clj-http " email)}
                            :follow-redirects true
                            :force-redirects true
                            :form-params param})]
        (if-not (= 200 (:status r))
          (throw (Throwable. (str "Error in mapping request: " (:body r))))
          (into {} (map #(string/split % #"\t") (rest (string/split (:body r) #"\n"))))))
      (throw (Throwable. "No more than 100,000 mappings per query allowed.")))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn now
  []
  (.getTime (java.util.Date.)))

(defn time-stamped-file
  ([base] (time-stamped-file base nil))
  ([base ext]
     (let [nf (if ext
                (fs/file (str base "-" (now) "." ext))
                (fs/file (str base "-" (now))))]
       (if-not (fs/exists? nf)
         nf
         (time-stamped-file base ext)))))


(ns clj-biosequence.core
  (:require [clojure.java.io :as io]
            [clojure.java.jdbc :as sql]
            [fs.core :as fs]
            [clojure.pprint :as pp]
            [clj-http.client :as client]
            [clojure.string :as string]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.persistence :as ps]))

(declare read-seq pb-read-line init-fasta-store init-fasta-sequence translate translate-string adjust-dna-frame map-frame)

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
    "Returns true if a protein and false otherwise."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; store

(defn save-biosequence
  "Saves an object to a store."
  [store obj]
  (ps/save-object store (accession obj) obj))

(defmacro with-biosequences-in-store
  "Provides a handle to a lazy list of biosequences in an object store."
  [[handle store] & code]
  `(ps/with-objects [~handle ~store]
     ~@code))

(defn update-biosequence 
  "Updates an object in the store with a current connection."
  [store obj]
  (ps/update-object store (accession obj) obj))

(defn update-biosequence-by-accession
  "Takes an accession number and key value pairs. If the biosequence exists in
   the current store it will be updated with the key value pairs and saved. 
   Throws an exception if a corresponding object is not found in the store."
  [store accession & args]
  (ps/update-object-by-id store accession args))

(defn get-biosequence [store accession]
  "Returns a biosequence object from a store implementing the biosequence"
  (ps/get-object store accession))

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
                          (:alphabet bs)
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
  (if (= (:alphabet s1) (:alphabet s2))
    (init-fasta-sequence (str (accession s1) "/" (accession s2))
                         (str (def-line s1) "/" (def-line s2))
                         (:alphabet s1)
                         (vec (concat (bs-seq s1) (bs-seq s2))))
    (throw (Throwable. "Incompatible alphabets for concatenation of biosequence."))))

(defn revcom-bioseq
  [bs]
  (if (protein? bs)
    (throw (Throwable. "Can't reverse/complement a protein sequence."))
    (init-fasta-sequence (accession bs)
                         (str (def-line bs) " - reverse-comp")
                         (:alphabet bs)
                         (ala/revcom (bs-seq bs)))))

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
                                              (-> (revcom-bioseq bs)
                                                  (sub-bioseq (- (* -1 frame) 1)))
                                              (#{4 5 6} frame)
                                              (-> (revcom-bioseq bs)
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
  
  (fasta-string [this]
    (if (:description this) 
      (pp/cl-format nil ">~A ~A~%~A~%"
                    (:accession this)
                    (:description this)
                    (apply str (:sequence this)))
      (pp/cl-format nil ">~A~%~A~%"
                    (:accession this)
                    (apply str (:sequence this))))))

(defn init-fasta-sequence
  "Returns a fastaSequence object with the specified information. Alphabet can be
   one of :iupacNucleicAcids or :iupacAminoAcids."
  [accession description alphabet sequence]
  (->fastaSequence accession description alphabet sequence))


;; reader

(defrecord fastaReader [strm alphabet]

  biosequenceReader
  
  (biosequence-seq [this]
    (letfn [(process [x]
              (if (empty? x)
                nil
                (let [s (first x)]
                  (lazy-seq (cons (fastaSequence. (first s)
                                                   (second s)
                                                   alphabet
                                                   (vec (nth s 2)))
                                  (process (rest x)))))))]
      (process (take-while (complement nil?)
                           (repeatedly #(read-seq (:strm this)))))))

  java.io.Closeable

  (close [this]
    (.close (:strm this))))

;; files and strings

(defrecord fastaFile [file alphabet]

  biosequenceIO

  (bs-reader [this]
    (->fastaReader
     (java.io.PushbackReader.
      (java.io.BufferedReader.
       (java.io.FileReader. (:file this))))
     (:alphabet this))))

(defrecord fastaString [str alphabet]

  biosequenceIO

  (bs-reader [this]
    (->fastaReader
     (java.io.PushbackReader.
      (java.io.StringReader. (:str this)))
     (:alphabet this))))

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
    (with-open [rdr (io/reader fastafile)]
      (let [seqs (biosequence-seq rdr)]
        (doseq [s seqs]
          (save-biosequence st s))))
    st))

(defn load-fasta-store
  "Loads a fastaStore."
  [dir]
  (let [file (first (fs/glob (str dir "/" "*.h2.db")))
        db-file (second (re-find  #"(.+)\.h2.db" (fs/absolute-path file)))]
    (if (not (nil? db-file))
      (assoc (->fastaStore db-file) :db
             (ps/make-db-connection db-file false))
      (throw (Throwable. "DB file not found!")))))

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

(defn- read-seq
  [^java.io.PushbackReader strm]
  (let [c (.read strm)]
    (cond (= c -1) nil
          (= (char c) \>)
          (let [d (pb-read-line strm)]
            (vector (get (re-find #"^([^\s]+)" d) 1)
                    (get (re-find #"^[^\s]+\s+(.+)" d) 1)
                    (loop [e (.read strm)
                           acc []]
                      (cond (= e -1) (apply str acc)
                            (= (char e) \>) (do
                                              (.unread strm e)
                                              (apply str acc))
                            :else
                            (recur (.read strm) (if (= (char e) \newline)
                                                  acc
                                                  (conj acc (char e))))))))
          :else
          (do (println (char c))
              (throw (Throwable. "Format error in fasta file."))))))

(defn- pb-read-line
  [^java.io.PushbackReader strm]
  (loop [c (char (.read strm))
         acc []]
    (if (= c \newline)
      (apply str acc)
      (recur (char (.read strm)) (conj acc c)))))

;; translation










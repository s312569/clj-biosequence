(ns clj-biosequence.core
  (:require [clojure.java.io :as io]
            [clojure.java.jdbc :as sql]
            [fs.core :as fs]
            [clojure.pprint :as pp]
            [clj-http.client :as client]
            [clojure.string :as string]))

(declare codon-tables declob with-connection-to-store read-seq pb-read-line init-fasta-store translate translate-string adjust-dna-frame map-frame db-parsing)

; macros

(defmacro with-connection-to-store
  "Provides a connection to a biosequence store."
  [[store] & code]
  `(sql/with-connection (:db ~store)
     (sql/transaction
      ~@code)))

(defmacro with-biosequences-in-file
  "Returns a handle to a lazy list of biosequences in a biosequence file object."
  [[handle file] & body]
  `(with-open [rdr# (java.io.PushbackReader. (io/reader (:file ~file)))]
     (let [~handle (biosequence-seq-file ~file rdr#)]
       ~@body)))

(defmacro with-biosequences
  "Provides a handle to a lazy list of biosequences in a biosequence store object."
  [[handle store] & body]
  `(with-connection-to-store [~store]
     (sql/with-query-results res#
       ["select * from sequence"]
       {:fetch-size 10 :concurrency :read-only :result-type :forward-only}
       (let [~handle (map #(assoc (declob (:src %))
                             :db (:db ~store))
                          res#)]
         ~@body))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; biosequence protocol
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; file

(defprotocol biosequenceFile
  (biosequence-seq-file [this rdr]
    "Returns a lazy sequence of biosequences from a biosequence file object.")
  (file-path [this]
    "Returns the path of a biosequence file."))

; biosequence

(defprotocol Biosequence
  (accession [this]
    "Returns the accession of a biosequence.")
  (accessions [this]
    "Returns a list of strings describing the accessions of a uniprot. 
     Ordered from most recent to oldest.")
  (def-line [this]
    "Returns the description of a biosequence.")
  (sequence-string [this]
    "Returns the sequence of a biosequence as a string.")
  (fasta-string [this]
    "Returns the biosequence as a string in fasta format.")
  (protein? [this]
    "Returns true if a protein and false otherwise.")
  (org-scientific-name [this]
    "Returns the scientific name of the organism possessing the sequence as a string.")
  (created [this]
    "Date sequence created.")
  (modified [this]
    "Date sequence modified.")
  (version [this]
    "Version of the sequence (Integer).")
  (database [this]
    "Returns the database that contains the sequence.")
  (taxonomy [this]
    "Returns the lineage from a uniprot as a list of strings. Strings in order from kingdom to species.")
  (taxid [this]
    "Returns the taxid of a sequence (Integer)."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-biosequence [accession]
  "Returns a biosequence object from a store implementing the biosequence"
  (sql/with-query-results row
    ["select * from sequence where id=?" accession]
    {:fetch-size 1}
    (if-not (empty? row)
      (declob (:src (first row))))))

(defn translate-biosequence
  "Returns a fastaSequence object corresponding to the protein translation 
   of the sequence in the specified frame."
  ([nucleotide frame] 
     (translate-biosequence nucleotide frame (codon-tables :standard)))
  ([nucleotide frame table]
     (translate nucleotide frame table)))

(defn six-frame-translation
  "Returns a lazy list of fastaSequence objects representing translations of
   a nucleotide biosequence object in six frames."
  ([nucleotide] (six-frame-translation nucleotide (codon-tables :standard)))
  ([nucleotide table]
     (map #(translate nucleotide % table)
          '(1 2 3 -1 -2 -3))))

;; persistance

(defn declob [^java.sql.Clob clob]
  "Turn a Derby 10.6.1.0 EmbedClob into a String"
  (binding [*read-eval* false]
    (read-string 
     (with-open [rdr (java.io.BufferedReader. (.getCharacterStream clob))]
       (apply str (line-seq rdr))))))

(defn make-db-connection 
  "Takes a file path and returns a specification for a database based on that file.
   Exists is a boolean and specifies whether an error is signalled if file already exists."
  [file exists]
  {
   :classname   "org.h2.Driver"
   :subprotocol "h2:file"
   :subname     (str file ";IFEXISTS=" exists)
   :user        "sa"
   :password    ""
   })

(defn make-mem-db-connection
  "Returns a in-memory database specification."
  []
  {:classname   "org.h2.Driver"
   :subprotocol "h2:mem:test"
   :subname     ";DB_CLOSE_DELAY=-1"
   :user        "sa"
   :password    ""
   })

(defn init-store
  "Initialises a permanent store."
  [store]
  (let [db (assoc store :db
                  (make-db-connection (:file store) false))]
    (try
      (do
        (sql/with-connection (:db db)
          (sql/create-table :sequence
                            [:id "varchar(255)" "PRIMARY KEY" "NOT NULL"]
                            [:src :clob])
          (sql/create-table :meta
                            [:id "varchar(50)" "PRIMARY KEY" "NOT NULL"]
                            [:src :clob]))
        db)
      (catch Exception e
        (fs/delete-dir (fs/parent (:file store)))))))

(defn init-in-mem-store
  "Initialises an in-memory permanent store"
  [store]
  (let [db (assoc store :db
                  (make-mem-db-connection))]
    (sql/with-connection (:db db)
      (sql/transaction
       (println db)
       (sql/create-table :sequence
                         [:id "varchar(255)" "PRIMARY KEY" "NOT NULL"]
                         [:src :clob])))
    db))

(defn save-object [obj]
  "Saves an object to the currently opened store."
  (sql/insert-record :sequence
                     {:id (accession obj)
                      :src (pr-str obj)}))

(defn update-object [obj]
  "Updates an object in the store with a current connection."
  (if (get-biosequence (accession obj))
    (do (sql/update-values :sequence 
                           ["id=?" (accession obj)]
                           {:src (pr-str obj)})
        obj)
    (throw (Throwable. "Object not currently in database."))))

(defn index-file-name
  "Returns a path suitable for a persistent store."
  [file]
  (let [dir (if (string? file)
              (str file ".ind")
              (str (:file file) ".ind"))]
    (if (.exists (io/file dir))
      (throw (Throwable. (str "Index directory already exists: " dir)))
      (str dir "/" (fs/base-name (if (string? file)
                                   file
                                   (:file file)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; an implementation of biosequence for fasta sequences
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastaSequence [accession description type sequence])

(defn init-fasta-sequence
  "Returns a fastaSequence object with the specified information."
  [accession description type sequence]
  (->fastaSequence accession description type
                   (string/replace sequence #"[^a-z|A-Z|*]" "")))

(extend-protocol Biosequence
  
  fastaSequence

  (accession [this]
    (:accession this))

  (accessions [this]
    (list (:accession this)))
  
  (sequence-string [this]
    (:sequence this))
  
  (def-line [this]
    (:description this))

  (protein? [this]
    (if (= :protein (:type this))
      true))
  
  (fasta-string [this]
    (if (:description this)
      (pp/cl-format nil ">~A ~A~%~A~%"
                    (:accession this)
                    (:description this)
                    (:sequence this))
      (pp/cl-format nil ">~A~%~A~%"
                    (:accession this)
                    (:sequence this))))

  (org-scientific-name [this]
    nil)

  (created [this]
    nil)

  (modified [this]
    nil)

  (version [this]
    nil)

  (database [this]
    nil)

  (taxonomy [this]
    nil)
  
  (taxid [this]
    nil))

;; files

(defn read-fasta-from-stream
  [rdr type]
  (letfn [(process [x]
            (if (empty? x)
              nil
              (let [s (first x)]
                (lazy-seq (cons (->fastaSequence (first s) 
                                             (second s)
                                             type
                                             (nth s 2))
                                (process (rest x)))))))]
    (process (take-while (complement nil?)
                         (repeatedly #(read-seq rdr))))))

(defrecord fastaFile [file type])

(extend-protocol biosequenceFile
  
  fastaFile

  (file-path [this]
    (:file this))

  (biosequence-seq-file [this rdr]
    (read-fasta-from-stream rdr (:type this))))

(defn init-fasta-file
  "Initialises fasta protein file. Accession numbers and description are processed by splitting 
   the string on the first space, the accession being the first value and description the second."
  [path type]
  (if-not (or (= :protein type) (= :nucleotide type))
    (throw (Throwable. "Fasta file type can be :protein or :nucleotide only."))
    (if (fs/exists? path)
      (->fastaFile path type)
      (throw (Throwable. (str "File not found: " path))))))

(defn fasta-seq-string
  "Returns a non-lazy list of fastaSequence objects from a string."
  [string type]
  (with-in-str string
    (with-open [rdr (java.io.PushbackReader. *in*)]
      (doall (read-fasta-from-stream rdr type)))))

;; persistence

(defrecord fastaStore [file type])

(defn index-fasta-file
  "Indexes a fastaFile object and returns a fastaStore object."
  ([fastafile] (index-fasta-file fastafile false))
  ([fastafile memory]
     (let [st (init-fasta-store fastafile memory)]
       (with-connection-to-store [st]
         (with-biosequences-in-file [l fastafile]
           (dorun (pmap #(future (save-object %)) l))))
       st)))

(defn load-fasta-store
  "Loads a fastaStore."
  [dir type]
  (let [file (first (fs/glob (str dir "/" "*.h2.db")))
        db-file (second (re-find  #"(.+)\.h2.db" (fs/absolute-path file)))]
    (if (not (nil? db-file))
      (assoc :db (->fastaStore db-file type)
             (make-db-connection file false))
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
   WormBase	            WORMBASE_ID	     both"
  [ids from to email]
  (let [param {:from from :to to :format "tab" 
               :query (apply str (doall (interpose "," (if (list? ids)
                                                         ids
                                                         (list ids)))))}
        address "http://www.uniprot.org/mapping/"
        r (client/post address 
                       {:client-params {"http.useragent" (str "clj-http " email)}
                        :follow-redirects true
                        :force-redirects true
                        :form-params param})]
    (if-not (= 200 (:status r))
      (throw (Throwable. (str "Error in mapping request: " (:body r))))
      (into {} (map #(string/split % #"\t") (rest (string/split (:body r) #"\n")))))))

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

(defn- db-parsing
  [st k]
  (with-connection-to-store [st]
    (sql/with-query-results row
      ["select * from parsing"]
      {:fetch-size 1}
      (if-not (empty? row)
        (read-string (k (first row)))))))

(defn- init-fasta-store
  "Initialises a fasta store."
  [fastafile memory]
  (if memory
    (init-in-mem-store (assoc 
                           (->fastaStore 
                            'in-memory 
                            (:type fastafile))))
    (init-store (->fastaStore 
                 (index-file-name (file-path fastafile)) 
                 (:type fastafile)))))

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

(defn- translate
  [seq frame table]
  (if (protein? seq)
    (throw (Throwable. "Can't translate a protein sequence!"))
    (let [trans (translate-string
                 (adjust-dna-frame 
                  (sequence-string seq) frame)
                 table)]
      (->fastaSequence (str (accession seq) "-"  (map-frame frame))
                       (str (def-line seq) " - Translated frame: " frame)
                       :protein
                       trans))))

(defn codon-tables
  "Returns a codon table suitable for use as an argument in the functions 
  'translate-biosequence' and 'six-frame-translation'. Takes a keyword argument
   to denote different codon tables. So far only provides the standard (:standard)
   codon table but further are planned."
  [key]
  (key {:standard
        '((\T
           (\T (\T \C "F") (\A \G "L"))
           (\C "S")
           (\A (\T \C "Y") (\A \G "*"))
           (\G (\T \C "C") (\A "*") (\G "W")))
          (\C
           (\T "L")
           (\C "P")
           (\A (\T \C "H") (\A \G "Q"))
           (\G "R"))
          (\A
           (\T (\T \C \A "I") (\G "M"))
           (\C "T")
           (\A (\T \C "N") (\A \G "K"))
           (\G (\T \C "S") (\A \G "R")))
          (\G
           (\T "V")
           (\C "A")
           (\A (\T \C "D") (\A \G "E"))
           (\G "G")))}))

(defn- get-amino-acid
  [lst table]
  (loop [l lst
         t table]
    (if (not (list? (last t)))
      (if (empty? t) "X" (last t))
      (recur (rest l)
             (remove #(nil? %)
                     (mapcat #(if (list? %)
                                (if (some #{(first l)} %) (rest %))
                                (if (= (first l) %) (rest %)))
                             t))))))

(defn- translate-string
  [string table]
  (let [s (string/upper-case (apply str (remove #{\space\newline} string)))]
    (apply str (map #(get-amino-acid % table) (partition-all 3 s)))))

(defn revcom-dna-string
  "Provides the reverse, complement of a string representing a DNA sequence. 
   Converts all 'U's to 'T's and upcases the string."
  [string]
  (apply str (map #(condp = %
                     \A \T
                     \T \A
                     \G \C
                     \C \G
                     \X) (reverse (string/replace (string/upper-case string)
                                                  #"U" "T")))))

(defn- map-frame
  [frame]
  (if (> frame 0)
    frame
    (condp = frame
      -1 4
      -2 5
      -3 6
      :else (throw (Throwable. "Frame out of bounds")))))

(defn- adjust-dna-frame
  [string frame]
  (let [s (string/upper-case string)]
    (cond
     (= frame 1) s
     (= frame -1) (revcom-dna-string s)
     (> frame 0) (subs s (- frame 1))
     (< frame 0) (subs (revcom-dna-string s) (- (* -1 frame) 1))
     :else (throw (Throwable. "Frame out of bounds")))))


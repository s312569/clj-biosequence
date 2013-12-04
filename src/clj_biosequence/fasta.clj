(in-ns 'clj-biosequence.core)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; an implementation of biosequence for fasta sequences
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastaSequence [acc description alphabet sequence]

  Biosequence
  
  (accession [this]
    (:acc this))
  
  (accessions [this]
    (list (:acc this)))
  
  (bs-seq [this]
    (:sequence this))
  
  (def-line [this]
    (:description this))
  
  (protein? [this]
    (if (= :iupacAminoAcids (:alphabet this))
      true
      false))
  
  (fasta-string [this]
    (str ">" (accession this) " " (def-line this) "\n" (bioseq->string this) "\n"))
  
  (alphabet [this]
    (:alphabet this)))

(defn init-fasta-sequence
  "Returns a new fastaSequence. Currently :iupacNucleicAcids
  and :iupacAminoAcids are supported alphabets."
  [accession description alphabet sequence]
  (->fastaSequence accession description alphabet (clean-sequence sequence alphabet)))

;; IO

(defrecord fastaReader [strm alphabet]

  biosequenceReader
  
  (biosequence-seq [this]
    (map (fn [[d s]]
           (let [seqs (apply str s)]
             (cond (not (re-find #"^>" (first d)))
                   (throw (Throwable. (str "Data corrupted at " (first d))))
                   (> (count d) 1)
                   (throw (Throwable. (str "No sequence for entry " (first d))))
                   :else
                   (init-fasta-sequence (second (re-find #"^>([^\s]+)" (first d)))
                                        (second (re-find #">[^\s]+\s+(.+)" (first d)))
                                        (:alphabet this)
                                        (clean-sequence seqs (:alphabet this))))))
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
                       (:alphabet this)))

  biosequenceFile

  (bs-path [this]
    (:file this)))

(defrecord fastaString [str alphabet]

  biosequenceIO

  (bs-reader [this]
    (init-fasta-reader (java.io.BufferedReader. (java.io.StringReader. (:str this)))
                       (:alphabet this))))

(defn init-fasta-file
  "Initialises fasta protein file. Accession numbers and description are 
   processed by splitting the string on the first space, the accession 
   being the first value and description the second."
  [path alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (if (fs/file? path)
      (->fastaFile path alphabet)
      (throw (Throwable. (str "File not found: " path))))))

(defn init-fasta-string
  "Initialises a fasta string. Accession numbers and description are
   processed by splitting the string on the first space, the accession
   being the first value and description the second."
  [str alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (->fastaString str alphabet)))

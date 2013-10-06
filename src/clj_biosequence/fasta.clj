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
  
  (reverse-comp [this]
    (if (protein? this)
      (throw (Throwable. "Can't reverse/complement a protein sequence."))
      (init-fasta-sequence (accession this)
                           (str (def-line this) " - Reverse-comp")
                           (alphabet this)
                           (ala/revcom (bs-seq this)))))
  
  (reverse-seq [this]
    (init-fasta-sequence (accession this)
                         (str (def-line this) " - Reversed")
                         (alphabet this) 
                         (vec (reverse (bs-seq this)))))
  
  (fasta-string [this]
    (if (:description this)
      (str ">" (accession this) " " (def-line this) "\n" (bioseq->string this) "\n")))
  
  (alphabet [this]
    (:alphabet this)))

(defmethod print-method clj_biosequence.core.fastaSequence
  [this ^java.io.Writer w]
  (print-tagged this w))

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
                       (:alphabet this))))

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
  "Initialises a fasta string."
  [str alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (->fastaString str alphabet)))

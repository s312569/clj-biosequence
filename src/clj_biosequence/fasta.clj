(in-ns 'clj-biosequence.core)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; sequence
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastaSequence [acc description alphabet sequence])

(extend fastaSequence
  biosequenceID
  (assoc default-biosequence-id
    :accession (fn [this] (:acc this))
    :accessions (fn [this] (list (accession this))))
  biosequenceDescription
  (assoc default-biosequence-description
    :description (fn [this] (:description this)))
  Biosequence
  (assoc default-biosequence-biosequence
    :bs-seq (fn [this] (:sequence this))
    :protein? (fn [this]
                (ala/alphabet-is-protein (:alphabet this)))
    :alphabet (fn [this] (:alphabet this))
    :moltype
    (fn [this] (if (protein? this) "AA" "Nucleic acid"))))

(defn init-fasta-sequence
  "Returns a new fastaSequence. Currently :iupacNucleicAcids
  and :iupacAminoAcids are supported alphabets."
  [accession description alphabet sequence]
  (->fastaSequence accession description alphabet
                   (clean-sequence sequence alphabet)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; IO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- parse-fasta
  [this]
  (let [l (if (:junk this)
            (drop-while #(not (= \> (first %))) (line-seq (:strm this)))
            (line-seq (:strm this)))]
    (map (fn [[d s]]
           (let [seqs (apply str (map trim s))]
             (cond (not (re-find #"^>" (first d)))
                   (throw
                    (Throwable. (str "Data corrupted at "
                                     (first d))))
                   (> (count d) 1)
                   (throw
                    (Throwable. (str "No sequence for entry "
                                     (first d))))
                   :else
                   (init-fasta-sequence
                    (second (re-find #"^>([^\s]+)" (first d)))
                    (second (re-find #">[^\s]+\s+(.+)" (first d)))
                    (:alphabet this)
                    seqs))))
         (partition 2 (partition-by #(re-find #"^>" %) l)))))

(defrecord fastaReader [strm alphabet]
  biosequenceReader
  (biosequence-seq [this] (parse-fasta this))
  java.io.Closeable
  (close [this] (.close ^java.io.BufferedReader (:strm this))))

(defn init-fasta-reader
  [strm alphabet & {:keys [junk] :or {junk false}}]
  (assoc (->fastaReader strm alphabet)
         :junk junk))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fasta files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastaFile [file alphabet opts])

(extend fastaFile
  biosequenceIO
  {:bs-reader
   (fn [this]
     (init-fasta-reader (apply bioreader (bs-path this) (:opts this))
                        (:alphabet this)
                        :junk (apply hash-map (:opts this))))}
  biosequenceFile
  default-biosequence-file
  fastaReduce
  {:fasta-reduce
   (fn [this func fold]
     (->> (iot/seq (:file this))
          (r/filter #(not (= \> (first %))))
          (r/map #(clean-sequence % (:alphabet this)))
          (r/map func)
          (r/fold fold)))})

(defn init-fasta-file
  "Initialises fasta protein file. Accession numbers and description
  are processed by splitting the string on the first space, the
  accession being the first value and description the second. Encoding
  can be specified using the :encoding keyword, defaults to UTF-8."
  [path alphabet & opts]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (if (file? path)
      (->fastaFile path alphabet opts)
      (throw (Throwable. (str "File not found: " path))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; fasta string
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord fastaString [str alphabet]
  biosequenceIO
  (bs-reader [this]
    (->fastaReader (java.io.BufferedReader.
                    (java.io.StringReader. (:str this)))
                   (:alphabet this))))

(defn init-fasta-string
  "Initialises a fasta string. Accession numbers and description are
   processed by splitting the string on the first space, the accession
   being the first value and description the second."
  [str alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (->fastaString str alphabet)))

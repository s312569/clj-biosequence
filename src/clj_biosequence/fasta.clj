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
    (ala/alphabet-is-protein (:alphabet this)))
  
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
    (let [l (line-seq (:strm this))]
      (map (fn [[d s]]
             (let [seqs (apply str (map trim s))]
               (cond (not (re-find #"^>" (first d)))
                     (throw (Throwable. (str "Data corrupted at " (first d))))
                     (> (count d) 1)
                     (throw (Throwable. (str "No sequence for entry " (first d))))
                     :else
                     (init-fasta-sequence (second (re-find #"^>([^\s]+)" (first d)))
                                          (second (re-find #">[^\s]+\s+(.+)" (first d)))
                                          (:alphabet this)
                                          seqs))))
           (partition 2 (partition-by #(re-find #"^>" %) l)))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defn init-fasta-reader
  [strm alphabet]
  (->fastaReader strm alphabet))

;; fasta files

(defprotocol fastaReduce
  (fasta-reduce [this func fold]
    "Applies a function to sequence data streamed line-by-line and
    reduces the results using the supplied `fold` function. Uses the
    core reducers library so the fold function needs to have an
    'identity' value that is returned when the function is called with
    no arguments."))

(defrecord fastaFile [file alphabet]

  biosequenceIO

  (bs-reader [this]
    (init-fasta-reader (reader (:file this))
                       (:alphabet this)))

  biosequenceFile

  (bs-path [this]
    (absolute-path (:file this)))

  (index-file [this]
    (init-indexed-fasta (bs-path this) (:alphabet this)))

  (index-file [this ofile]
    (init-indexed-fasta (absolute-path ofile) (:alphabet this)))

  fastaReduce

  (fasta-reduce
    [this func fold]
    (->> (iot/seq (:file this))
         (r/filter #(not (= \> (first %))))
         (r/map #(clean-sequence % (:alphabet this)))
         (r/map func)
         (r/fold fold))))

(defn init-fasta-file
  "Initialises fasta protein file. Accession numbers and description are 
   processed by splitting the string on the first space, the accession 
   being the first value and description the second."
  [path alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (if (file? path)
      (->fastaFile path alphabet)
      (throw (Throwable. (str "File not found: " path))))))

;; strings

(defrecord fastaString [str alphabet]

  biosequenceIO

  (bs-reader [this]
    (init-fasta-reader (java.io.BufferedReader. (java.io.StringReader. (:str this)))
                       (:alphabet this))))

(defn init-fasta-string
  "Initialises a fasta string. Accession numbers and description are
   processed by splitting the string on the first space, the accession
   being the first value and description the second."
  [str alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet keyword. Currently :iupacNucleicAcids :iupacAminoAcids are allowed."))
    (->fastaString str alphabet)))

;; indexed files

(defrecord indexedFastaFile [index path alphabet]

  biosequenceFile

  (bs-path [this]
    (absolute-path (:path this)))

  indexFileIO

  (bs-writer [this]
    (init-index-writer this))

  biosequenceReader

  (biosequence-seq [this]
    (map (fn [[o l]]
           (map->fastaSequence (read-one o l (str (bs-path this) ".bin"))))
         (vals (:index this))))

  (get-biosequence [this accession]
    (let [[o l] (get (:index this) accession)]
      (if o
        (map->fastaSequence (read-one o l (str (bs-path this) ".bin")))))))

(defn init-indexed-fasta [file alphabet]
  (->indexedFastaFile {} file alphabet))

(defmethod print-method clj_biosequence.core.indexedFastaFile
  [this w]
  (print-tagged this w))

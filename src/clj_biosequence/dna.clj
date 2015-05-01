(in-ns 'clj-biosequence.core)

(defn reverse-comp
  "Returns a new fastaSequence with the reverse complement sequence."
  [this]
  {:pre [(not (protein? this))]}
  (init-fasta-sequence (accession this)
                       (str (description this) " - Reverse-comp")
                       (alphabet this)
                       (apply str (ala/revcom (bs-seq this)))))

(defn reverse-seq
  "Returns a new fastaSequence with the reverse sequence."
  [this]
  {:pre [(not (protein? this))]}
  (init-fasta-sequence (accession this)
                       (str (description this) " - Reversed")
                       (alphabet this) 
                       (apply str (reverse (bs-seq this)))))

(defn translate
  "Returns a fastaSequence sequence representing the translation of
  the specified biosequence in the specified frame."
  [bs frame & {:keys [table id-alter] :or {table (ala/codon-tables 1)
                                           id-alter true}}]
  {:pre [(not (protein? bs))]}
  (let [f ({1 1 2 2 3 3 4 4 5 5 6 6 -1 4 -2 5 -3 6} frame)]
    (if (not (#{1 2 3 4 5 6} f))
      (throw (IllegalArgumentException.
              (str "Invalid frame: " frame))))
    (init-fasta-sequence
     (if id-alter (str (accession bs) "-" f) (accession bs))
     (str (description bs) " - Translated frame: " f)
     :iupacAminoAcids
     (let [v (cond (#{1 2 3} f)
                   (sub-bioseq bs (- f 1))
                   (#{4 5 6} f)
                   (-> (reverse-comp bs)
                       (sub-bioseq ( - f 4))))]
       (apply str (map #(ala/codon->aa % table)
                       (partition-all 3 (bs-seq v))))))))

(defn six-frame-translation
  "Returns a lazy list of fastaSequence objects representing
  translations of a nucleotide biosequence object in six frames."
  ([nucleotide] (six-frame-translation nucleotide (ala/codon-tables 1)))
  ([nucleotide table]
   (map #(translate nucleotide % :table table)
        '(1 2 3 -1 -2 -3))))

(defn n50
  "Takes anything that can have `biosequence-seq' called on it and
  returns the N50 of the sequences therein."
  [reader]
  (let [sa (sort > (pmap #(count (bs-seq %))
                         (biosequence-seq reader)))
        t (/ (reduce + sa) 2)
        n50 (atom 0)]
    (loop [l sa]
      (if (seq l)
        (do (swap! n50 + (first l))
            (if (>= @n50 t)
              (first l)
              (recur (rest l))))
        (throw (Throwable. "N50 calculation failed!"))))))

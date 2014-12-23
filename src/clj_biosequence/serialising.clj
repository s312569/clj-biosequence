(in-ns 'clj-biosequence.core)

(defrecord biosequenceIndexReader [index strm]

  biosequenceReader

  (biosequence-seq [this]
    (ind/object-seq (:strm this) (vals (:index this))))

  (get-biosequence [this accession]
    (if-let [i ((:index reader) accession)]
      (first (ind/object-seq (:strm this) (list i)))))

  java.io.Closeable

  (close [this]
    (ind/close-index-reader (:strm this))))

(defrecord biosequenceIndex [index path]

  biosequenceIO

  (bs-reader [this]
    (->biosequenceIndexReader (:index this)
                              (ind/index-reader (bs-path this))))

  biosequenceFile

  (bs-path [this]
    (absolute-path (:path this))))

(defn index-biosequence-file
  [file & {:keys [func] :or {func accession}}]
  (try
    (let [index (->biosequenceIndex {} (bs-path file))
          i (with-open [w (ind/index-writer (bs-path index))]
              (assoc index :index
                     (with-open [r (bs-reader file)]
                       (ind/index-objects w (biosequence-seq r) func))))]
      (ind/save-index (bs-path file) (:index i))
      i)
    (catch Exception e
      (ind/delete-index (bs-path file))
      (println (str "Exception: " (.getMessage e))))))

(defn load-biosequence-index
  [path]
  (->biosequenceIndex (ind/load-indexed-file path) (absolute-path path)))

(defn delete-indexed-biosequence
  [index-file]
  (ind/delete-index (bs-path index-file)))

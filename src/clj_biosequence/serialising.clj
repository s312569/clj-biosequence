(in-ns 'clj-biosequence.core)

(defrecord biosequenceIndexReader [index strm]

  biosequenceReader

  (biosequence-seq [this]
    (ind/object-seq (:strm this) (vals (dissoc (:index this) :parameters))))

  (get-biosequence [this accession]
    (if-let [i ((:index this) accession)]
      (first (ind/object-seq (:strm this) (list i)))))

  biosequenceParameters

  (parameters [this]
    (get-biosequence this :parameters))

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

(defn- parameters?
  [r]
  (if (satisfies? biosequenceParameters r)
    (parameters r)))

(defn p-key
  [_]
  :parameters)

(defn index-biosequence-file
  [file & {:keys [func] :or {func accession}}]
  (try
    (let [index (->biosequenceIndex {} (bs-path file))
          i (with-open [w (ind/index-writer (bs-path index))
                        r (bs-reader file)]
              (assoc index :index
                     (merge (ind/index-objects w (biosequence-seq r) func)
                            (if-let [p (parameters? r)]
                              (ind/index-objects w (list p) p-key)))))]
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

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

(defn is-indexed?
  [file]
  (exists? (str (bs-path file) ".idx")))

(defn load-biosequence-index
  [path]
  (->biosequenceIndex (ind/load-indexed-file path) (absolute-path path)))

(defn index-biosequence-collection
  [collection file & {:keys [func] :or {func accession}}]
  (try
    (let [index (->biosequenceIndex {} file)
          i (with-open [w (ind/index-writer (bs-path index))]
              (assoc index :index
                     (ind/index-objects w collection func)))]
      (ind/save-index file (:index i))
      i)
    (catch Exception e
      (ind/delete-index file)
      (println (str "Exception: " (.getMessage e))))))

(defn index-biosequence-file
  [file & {:keys [func] :or {func accession}}]
  (if (is-indexed? file)
    (load-biosequence-index (bs-path file))
    (with-open [r (bs-reader file)]
      (-> (if-let [p (parameters? r)]
            (cons p (biosequence-seq r))
            (biosequence-seq r))
          (index-biosequence-collection (bs-path file))))))

(defn delete-indexed-biosequence
  [index-file]
  (ind/delete-index (bs-path index-file)))

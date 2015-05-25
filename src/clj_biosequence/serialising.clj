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

(defn init-biosequence-index
  [path]
  (->biosequenceIndex {} path))

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
  [collection file & {:keys [func index] :or {func accession index nil}}]
  (try
    (let [index (if index index
                    (->biosequenceIndex {} file))
          i (with-open [w (ind/index-writer (bs-path index))]
              (assoc index :index
                     (merge (:index index)
                            (ind/index-objects w collection func))))]
      (ind/save-index (if index (bs-path index) file) (:index i))
      i)
    (catch Exception e
      (ind/delete-index file)
      (println (str "Exception: " (.getMessage e))))))

(defn index-biosequence-file
  [file & {:keys [func index] :or {func accession index nil}}]
  (if (is-indexed? file)
    (load-biosequence-index (bs-path file))
    (with-open [r (bs-reader file)]
      (-> (if-let [p (parameters? r)]
            (cons p (biosequence-seq r))
            (biosequence-seq r))
          (index-biosequence-collection (bs-path file)
                                        :index (if index index))))))

(defn index-combine-files
  "Takes a collection of biosequence files and generates an index
  combining all the files. If files have parameters only the parameter
  object from the last file will be retained."
  [files outpath]
  (loop [f files
         i (init-biosequence-index outpath)]
    (if-not (seq f)
      i
      (recur (rest f) (index-biosequence-file (first f) :index i)))))

(defn delete-indexed-biosequence
  [index-file]
  (ind/delete-index (bs-path index-file)))

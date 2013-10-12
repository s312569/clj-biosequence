(in-ns 'clj-biosequence.core)

(defrecord storeReader [cursor]

  biosequenceReader

  (biosequence-seq [this]
    (ps/record-seq (:cursor this)))

  java.io.Closeable

  (close [this]
    nil))

(defrecord biosequenceStore [name]

  biosequenceIO

  (bs-reader [this]
    (->storeReader (ps/find-all (:name this)))))

(defn update-biosequence
  [bs s]
  (ps/update-record (bs-save bs) (:name s)))

(defn save-biosequences
  [lst s]
  (ps/save-records (pmap bs-save lst) (:name s)))

(defn get-biosequence
  [a s]
  (ps/get-record a (:name s)))

(defn index-biosequence-file
  [file name]
  (let [s (->biosequenceStore name)]
    (do
      (with-open [r (bs-reader file)]
        (save-biosequences (biosequence-seq r) s))
      s)))

(defn bs-collections
  []
  (ps/get-collections))

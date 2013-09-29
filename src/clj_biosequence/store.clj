(in-ns 'clj-biosequence.core)

(defrecord storeReader [conn stm]

  biosequenceReader

  (biosequence-seq [this]
    (map #(ps/declob (:src %))
         (ps/store-result-map (:stm this) "select * from object")))

  java.io.Closeable

  (close [this]
    (.close (:stm this))
    (.close (:conn this))))

(defrecord biosequenceStore [path]

  biosequenceIO

  (bs-reader [this]
    (let [c (ps/store-connection (:db this))]
      (->storeReader c
                     (ps/store-statement c)))))

(defn new-store
  [path]
  (ps/init-store (->biosequenceStore
                  (ps/index-file-name path))))

(defn load-store
  "Loads a fastaStore."
  [dir]
  (let [file (first (fs/glob (str dir "/" "*.h2.db")))
        db-file (second (re-find  #"(.+)\.h2.db" (fs/absolute-path file)))]
    (if (not (nil? db-file))
      (->biosequenceStore
       (ps/make-db-connection db-file false))
      (throw (Throwable. "DB file not found!")))))

(defn save-biosequence
  "Saves an object to a store."
  [obj]
  (ps/save-object (accession obj) obj))

(defn update-biosequence 
  "Updates an object in the store with a current connection."
  [obj]
  (ps/update-object (accession obj) obj))

(defn update-biosequence-by-accession
  "Takes an accession number and key value pairs. If the biosequence exists in
   the current store it will be updated with the key value pairs and saved. 
   Throws an exception if a corresponding object is not found in the store."
  [store accession & args]
  (ps/update-object-by-id store accession args))

(defn get-biosequence [accession]
  "Returns a biosequence object from a store implementing the biosequence"
  (ps/get-object accession))

(defmacro with-biosequence-store
  [[st] & code]
  `(ps/with-store [~st]
     ~@code))

(defn index-biosequence-file
  [file store]
  (with-biosequence-store [store]
    (with-open [rdr (bs-reader file)]
      (dorun
       (pmap #(save-biosequence %)
             (biosequence-seq rdr))))
    store))

(defmacro with-biosequences
  [[handle st] & code]
  `(ps/with-objects [~handle ~st]
     ~@code))

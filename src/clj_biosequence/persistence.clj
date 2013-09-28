(ns clj-biosequence.persistence
  (:require [clojure.java.jdbc :as sql]
            [fs.core :as fs]
            [clojure.data.xml :as xml]
            [clojure.java.io :as io]
            [clojure.edn :as ed]))

;; persistance

(defmacro with-store
  "Provides a connection to an object store."
  [[store] & code]
  `(sql/with-connection (:db ~store)
     (sql/transaction
      ~@code)))

(defmacro with-objects
  "Provides a handle to a lazy list of all objects in a store."
  [[handle store] & code]
  `(with-store [~store]
     (sql/with-query-results res#
       ["select * from object"]
       {:fetch-size 10 :concurrency :read-only :result-type :forward-only}
       (let [~handle (map #(declob (:src %))
                          res#)]
         ~@code))))

(defn my-tag->factory
  "Returns the map-style record factory for the `tag` symbol.  Returns nil if `tag` does not
  refer to a record. Blatantly copied from miner.tagged."
  [tag]
  (when (namespace tag)
    (resolve (symbol (str (namespace tag) "/map->" (name tag))))))

(defn default-reader
  [tag value]
  (if-let [factory (and (map? value)
                     (my-tag->factory tag))]
    (factory value)
    (throw (Throwable. (str "Record not supported: " tag)))))

(defn declob [^java.sql.Clob clob]
  "Turn a Derby 10.6.1.0 EmbedClob into a String"
  (ed/read-string
   {:default #'default-reader
    :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
   (with-open [rdr (java.io.BufferedReader. (.getCharacterStream clob))]
     (apply str (line-seq rdr)))))

(defn make-db-connection 
  "Takes a file path and returns a specification for a database based on 
   that file. Exists is a boolean and specifies whether an error is 
   signalled if file already exists."
  [file exists]
  {
   :classname   "org.h2.Driver"
   :subprotocol "h2:file"
   :subname     (str file ";IFEXISTS=" exists)
   :user        "sa"
   :password    ""
   })

(defn init-store
  "Initialises a permanent store."
  [store]
  (let [db (assoc store :db
                  (make-db-connection (:file store) false))]
    (try
      (do
        (sql/with-connection (:db db)
          (sql/create-table :object
                            [:id "varchar(255)" "PRIMARY KEY" "NOT NULL"]
                            [:src :clob])
          (sql/create-table :meta
                            [:id "varchar(50)" "PRIMARY KEY" "NOT NULL"]
                            [:src :clob]))
        db)
      (catch Exception e
        (fs/delete-dir (fs/parent (:file store)))))))

(defn get-object [id]
  "Returns an object from the currently open store."
  (sql/with-query-results row
    ["select * from object where id=?" id]
    {:fetch-size 1}
    (if-not (empty? row)
      (declob (:src (first row))))))

(defn save-object [id obj]
  "Saves an object to the currently opened store."
  (sql/insert-record :object
                     {:id id
                      :src (pr-str obj)}))

(defn update-object [id obj]
  "Updates an object in the store with a current connection. Returns
   the object."
  (sql/update-values :object
                     ["id=?" id]
                     {:src (pr-str obj)})
  obj)

(defn update-object-by-id
  "Takes an accession number and key value pairs. If the object exists in
   the current store the object with the corresponding accession number
   will be updated with the key value pairs and saved. Throws an exception
   if a corresponding object is not found in the store."
  [store id & args]
  (if-let [obj (get-object id)]
    (update-object (apply assoc obj args))
    (throw (Throwable. (str "No object found with id: " id)))))

(defn index-file-name
  "Returns a path suitable for a persistent store."
  [file]
  (let [dir (if (string? file)
              (str file ".ind")
              (str (:file file) ".ind"))]
    (if (.exists (io/file dir))
      (throw (Throwable. (str "Index directory already exists: " dir)))
      (str dir "/" (fs/base-name (if (string? file)
                                   file
                                   (:file file)))))))

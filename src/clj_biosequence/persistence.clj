(ns clj-biosequence.persistence
  (:require [clojure.java.jdbc :as sql]
            [fs.core :as fs]
            [clojure.data.xml :as xml]
            [clojure.java.io :as io]
            [clojure.edn :as ed]))

(defn store-connection
  [spec]
  (sql/get-connection spec))

(defn store-statement
  [conn]
  (.createStatement conn))

(defn store-result-map
  [s q]
  (try
    (sql/result-set-seq (.executeQuery s q))
    (catch java.sql.SQLException e (str "SQL error: " (.getMessage e)))))

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
                  (make-db-connection (:path store) false))]
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
        (fs/delete-dir (fs/parent (:path store)))))))

(defn get-object [a s]
  "Returns an object from the currently open store. Returns nil if not found."
  (let [r (-> (sql/query (:db s)
                         ["select * from object where id=?" a]
                         :row-fn (fn [x] (declob (second x)))
                         :as-arrays? true))]
    (second r)))

(defn save-records
  [l s]
  (sql/db-transaction
   [tcon (:db s)]
   (dorun (pmap #(sql/insert! tcon :object %) l))))

(defn update-records
  "Updates an object in the store with a current connection."
  [l s]
  (letfn [(ud [x]
            (sql/update! (:db s) :object x ["id=?" (:id x)]))]
    (sql/db-transaction
     [tcon (:db s)]
     (dorun (map ud l)))))

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

(ns clj-biosequence.persistence
  (:require [fs.core :as fs]
            [clojure.data.xml :as xml]
            [clojure.java.io :as io]
            [monger.core :as mg]
            [monger.collection :as mc]
            [monger.conversion :as con]
            [monger.db :as mdb]
            [clojure.edn :as ed])
  (:import [com.mongodb MongoOptions ServerAddress]
           [org.bson.types ObjectId]))

(declare bs-read)

(mg/connect!)
(mg/use-db! "clj-biosequence")

(defn my-tag->factory
  "Returns the map-style record factory for the `tag` symbol.  Returns nil if `tag` does not
  refer to a record."
  [tag]
  (when (namespace tag)
    (resolve (symbol (str (namespace tag) "/map->" (name tag))))))

(defn default-reader
  [tag value]
  (if-let [factory (and (map? value)
                     (my-tag->factory tag))]
    (factory value)
    (throw (Throwable. (str "Record not supported: " tag)))))

(defn find-all
  [c]
  (mc/find c))

(defn record-seq
  [c]
  (pmap #(let [r (con/from-db-object % true)]
           (merge (bs-read (:src r)) (dissoc r :src))) (seq c)))

(defn bs-read
  [s]
  (ed/read-string
   {:default #'default-reader
    :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
   s))

(defn update-record
  [m c]
  (mc/update-by-id c (:_id m) m))

(defn save-records
  [l c]
  (mc/insert-batch c (doall (pmap #(assoc % :_id (ObjectId.)) l)))
  (mc/ensure-index c {"acc" 1} {:unique true :name "unique_acc"}))

(defn get-record
  [a c]
  (let [r (first (mc/find-maps "sequences" {:acc a}))]
    (merge (bs-read (:src r)) (dissoc r :src))))

(defn get-collections
  []
  (mdb/get-collection-names))

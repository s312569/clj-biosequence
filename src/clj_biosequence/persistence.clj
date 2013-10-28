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

(mg/connect!)

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
  [d c & {:keys [query] :or {query nil}}]
  (mg/use-db! d)
  (if query
    (mc/find c query)
    (mc/find c)))

(defn bs-read
  [s]
  (ed/read-string
   {:default #'default-reader
    :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
   s))

(defn record-seq
  [c]
  (map #(let [r (con/from-db-object % true)]
           (merge (bs-read (:src r)) (dissoc r :src))) (seq c)))

(defn update-record
  [m d c]
  (mg/use-db! d)
  (mc/update-by-id c (:_id m) m))

(defn save-records
  [l d c & {:keys [indexes] :or {indexes '(:acc)}}]
  (mg/use-db! d)
  (mc/insert-batch c (pmap #(assoc % :_id (ObjectId.)) l))
  (mc/ensure-index c (apply array-map (interleave indexes (iterate identity 1)))
                   {:unique true :name "unique_acc"}))

(defn get-record
  [h d c]
  (mg/use-db! d)
  (let [r (first (mc/find-maps c h))]
    (merge (bs-read (:src r)) (dissoc r :src))))

(defn get-collections
  [d]
  (mg/use-db! d)
  (mdb/get-collection-names))

(defn drop-collection
  [d n]
  (mg/use-db! d)
  (mc/drop n))

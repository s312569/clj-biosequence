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

;; printing objects

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

(defn bs-read
  [s]
  (ed/read-string
   {:default #'default-reader
    :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
   s))

;; db interactions

(defn find-all
  [d c & {:keys [query] :or {query nil}}]
  (mg/use-db! d)
  (if query
    (mc/find c query)
    (mc/find c)))

(defn record-seq
  "Returns a lazy sequence of all records in the collection `c`."
  [c]
  (map #(let [r (con/from-db-object % true)]
          (assoc (bs-read (:src r)) :_id (:_id r))) (seq c)))

(defn update-record
  "Updates record `m` in database `d` and collection `c`."
  [m d c]
  (mg/use-db! d)
  (let [w  (mc/update-by-id c (:_id m) m)]))

(defn save-records
  "Saves a list of records `l` into database `d` and collection `c`.
  `indexes` keyword allows for the specification of slots that are to
  be used as indexes."
  [l d c]
  (mg/use-db! d)
  (mc/insert-batch c (pmap #(assoc % :_id (ObjectId.)) l))
  (mc/ensure-index c (array-map :acc 1)
                   {:unique true :name "unique_acc"}))

(defn get-record
  "Returns a record corresponding to the query hash, `h`, from
   database, `d`, and collection, `c`."
  [h d c]
  (mg/use-db! d)
  (let [r (first (mc/find-maps c h))]
    (assoc (bs-read (:src r)) :_id (:_id r))))

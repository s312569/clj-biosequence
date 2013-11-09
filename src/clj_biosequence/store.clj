(ns clj-biosequence.store
  (:require [clj-biosequence.write :as wr]
            [monger.core :as mg]
            [monger.collection :as mc]
            [monger.conversion :as con]
            [monger.db :as mdb]
            [monger.util :as mu])
  (:import [com.mongodb MongoOptions ServerAddress WriteConcern]
           [org.bson.types ObjectId]))

(declare save-list store-read)

;; interface

(defprotocol storeCollectionIO
  (mongo-save-file [this project name]
    "Saves a biosequenceFile to a mongoDB for random access."))

(defprotocol storeCollectionAccess
  (collection-seq [this]
    "Returns a lazy list of entries in a biosequenceCollection."))

;; project

(defrecord mongoProject [name])

(defrecord biosequenceCollection [name pname type]

  storeCollectionAccess

  (collection-seq [this]
    (map store-read
         (mc/find-maps "sequences" {:batch_id this}))))

(defmethod print-method clj_biosequence.store.biosequenceCollection
  [this ^java.io.Writer w]
  (wr/print-tagged this w))

(defn init-biosequence-collection
  [name pname type]
  (->biosequenceCollection name pname type))

;; functions

(defn mongo-connect
  "Connects to the default mongo database server."
  []
  (mg/connect!))

(defn mongo-disconnect
  []
  "Disconnects from the default mongo database server."
  (mg/disconnect!))

(defn init-project
  "Returns a mongoProject for storing sequence collections. Used for
  accessing existing projects or initialising new ones."
  [name]
  (when (not (mc/exists? "sequence"))
    (mg/use-db! "clj-projects")
    (mc/create "sequence" {})
    (mc/ensure-index "sequences"
                     (array-map :acc 1 :batch_id 1 :_id 1)
                     {:unique true})
    (mc/ensure-index "sequences"
                     (array-map :pname 1 :iname 1 :coll 1 :_id 1)
                     {:sparse true}))
  (->mongoProject name))

(defn list-projects
  "Returns a set of projects on the server."
  []
  (mg/use-db! "clj-projects")
  (distinct (map :pname (mc/find-maps "sequences" {:coll "t"} [:pname]))))

(defn drop-project
  "Takes a mongoProject and drops it from the database."
  [project]
  (mg/use-db! "clj-projects")
  (mc/remove "sequences" {:pname (:name project)}))

(defn list-collections
  "Takes a mongoProject and returns a hash-map of collection names and
  types in the project."
  [project]
  (mg/use-db! "clj-projects")
  (->> (get-collections project)
       (map (fn [{n :name t :type}] {n t}))))

(defn get-collections
  "Returns a list of collections in a project."
  ([project] (get-collections project nil))
  ([project collection]
     (mg/use-db! "clj-projects")
     (let [c (if collection {:iname collection} {})]
       (->> (mc/find-maps "sequences"
                          (merge {:pname (:name project) :coll "t"} c)
                          [:src])
            (map #(wr/bs-read (:src %)))))))

(defn drop-collection
  "takes a collection object and drops it from the database."
  [collection]
  (mg/use-db! "clj-projects")
  (mc/remove "sequences"
             {:batch_id (-> (get-collections (init-project (:pname collection))
                                             (:name collection))
                            first
                            :batch_id)}))

(defn get-record
  ([collection value] (get-record collection :acc value))
  ([collection key value & kv]
     (mg/use-db! "clj-projects")
     (map store-read (mc/find-maps "sequences" (merge {key value
                                                       :batch_id (:batch_id collection)}
                                                      (apply hash-map kv))))))

(defn save-list
  "Takes a list of hash-maps for insertion into a mongoDB and a
  collection object and inserts all members of the list."
  [l i]
  (mg/use-db! "clj-projects")
  (let [u (mu/random-uuid)]
    (try
      (mc/insert-batch "sequences"
                       (cons
                        {:_id (ObjectId.) :acc (mu/random-uuid)
                         :src (pr-str (assoc i :batch_id u))
                         :pname (:pname i) :iname (:name i)
                         :coll "t" :batch_id u}
                        (pmap #(assoc % :_id (ObjectId.) :batch_id u
                                      :pname (:pname i) :iname (:name i))
                              l))
                       WriteConcern/JOURNAL_SAFE)
      (catch Exception e
        (mc/remove "sequences" {:batch_id u})
        (throw e)))
    (assoc i :batch_id u)))

;; utilities

(defn- store-read
  [h]
  (if-let [o (wr/bs-read (:src h))]
    (if (string? o)
      o
      (merge o (dissoc h :src)))))


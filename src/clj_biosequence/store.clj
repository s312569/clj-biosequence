(ns clj-biosequence.store
  (:require [clj-biosequence.write :as wr]
            [monger.core :as mg]
            [monger.collection :as mc]
            [monger.conversion :as con]
            [monger.db :as mdb]
            [monger.util :as mu])
  (:import [com.mongodb MongoOptions ServerAddress WriteConcern]
           [org.bson.types ObjectId]))

(declare save-list store-read get-record)

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
    (get-record this :element "sequence")))

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
  (mg/use-db! "clj-projects")
  (if (not (mc/exists? "sequences"))
    (let [p (mc/insert-and-return "sequences"
                                  {:pname name :project "t" :started (new java.util.Date)})]
      (assoc (->mongoProject name) :started (:started p)))
    (let [p (first (mc/find-maps "sequences" {:pname name :project "t"}))]
      (if p
        (assoc (->mongoProject name) :started (:started p))
        (assoc (->mongoProject name) :started
               (mc/insert-and-return "sequences"
                                     {:pname name :project "t" :started (new java.util.Date)}))))))

(defn list-projects
  "Returns a set of projects on the server."
  []
  (mg/use-db! "clj-projects")
  (distinct (map :pname (mc/find-maps "sequences" {:project "t"} [:pname]))))

(defn drop-project
  "Takes a mongoProject and drops it from the database."
  [project]
  (mg/use-db! "clj-projects")
  (mc/remove "sequences" {:pname (:name project)}))

;; collection functions

(defn get-collection
  "Returns a list of collections in a project."
  [project collection]
  (mg/use-db! "clj-projects")
  (-> (mc/find-maps "sequences"
                    {:pname (:name project) :coll "t" :name collection}
                    [:src])
      first
      :src
      wr/bs-read))

(defn list-collections
  "Takes a mongoProject and returns a hash-map of collection names and
  types in the project."
  [project]
  (mg/use-db! "clj-projects")
  (into {} (for [x (mc/find-maps "sequences"
                                 {:pname (:name project) :coll "t"}
                                 [:name :type])]
             (vector (:name x) (:type x)))))

(defn drop-collection
  "takes a collection object and drops it from the database."
  [collection]
  (mg/use-db! "clj-projects")
  (mc/remove "sequences"
             {:pname (:pname collection) :name (:name collection)}))

;; record functions

(defn get-record
  ([collection value] (get-record collection :acc value))
  ([collection key value & kv]
     (mg/use-db! "clj-projects")
     (store-read (first (mc/find-maps "sequences" (merge {key value
                                                          :pname (:pname collection)
                                                          :name (:name collection)}
                                                         (apply hash-map kv)))))))

;; saving to store

(defn save-list
  "Takes a list of hash-maps for insertion into a mongoDB and a
  collection object and inserts all members of the list."
  [l i]
  (mg/use-db! "clj-projects")
  (let [u (mu/random-uuid)
        ci (ObjectId.)]
    (try
      (do (mc/ensure-index "sequences"
                           (array-map :acc 1 :pname 1 :name 1 :_id 1))
          (mc/ensure-index "sequences"
                           (array-map :acc 1 :pname 1 :name 1)
                           {:unique true})
          (dorun (pmap #(mc/insert-batch "sequences"
                                         %
                                         WriteConcern/JOURNAL_SAFE)
                       (partition-all 1000
                                      (cons
                                       (merge i
                                              {:_id ci :acc ci
                                               :src (pr-str (assoc i :batch_id u))
                                               :coll "t" :batch_id u})
                                       (pmap #(merge {:_id (ObjectId.) :batch_id u
                                                      :pname (:pname i) :name (:name i)
                                                      :element "sequence"} %)
                                             l))))))
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


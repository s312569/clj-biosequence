(ns clj-biosequence.store
  (:require [clj-biosequence.core :as bs]
            [fs.core :as fs]
            [clojure.data.xml :as xml]
            [clojure.java.io :as io]
            [monger.core :as mg]
            [monger.collection :as mc]
            [monger.conversion :as con]
            [monger.db :as mdb]
            [monger.util :as mu])
  (:import [clj_biosequence.core fastaSequence fastaFile]
           [com.mongodb MongoOptions ServerAddress WriteConcern]
           [org.bson.types ObjectId]))

(declare save-list store-read biosequence-save-helper collection-seq-helper)

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
    (collection-seq-helper this)))

(defmethod print-method clj_biosequence.store.biosequenceCollection
  [this ^java.io.Writer w]
  (bs/print-tagged this w))

;; extend persistence

(extend-protocol storeCollectionIO
  fastaFile
  (mongo-save-file [this project name]
    (biosequence-save-helper this project name "biosequence/fasta")))

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
  (->mongoProject name))

(defn list-projects
  "Returns a set of projects on the server."
  []
  (mg/use-db! "clj-projects")
  (set (mc/distinct "sequences" :pname)))

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
  (map #(hash-map % (:type (mc/find-one-as-map "sequences" {:iname %} [:type])))
       (mc/distinct "sequences" :iname)))

(defn get-collection
  "Returns an collection object."
  [name project]
  (mg/use-db! "clj-projects")
  (bs/bs-read (:i (mc/find-one-as-map "sequences"
                                      {:pname (:name project) :iname name}))))

(defn drop-collection
  "takes a collection object and drops it from the database."
  [collection]
  (mc/remove "sequences" {:iname (:name collection) :pname (:pname collection)}))

(defn save-list
  "Takes a list of hash-maps for insertion into a mongoDB and a
  collection object and inserts all members of the list."
  [l i]
  (let [u (mu/random-uuid)]
    (try
      (do (mc/ensure-index "sequences"
                           (array-map :acc 1 :pname -1 :iname -1)
                           {:unique true :sparse true})
          (mc/insert-batch "sequences"
                           (pmap #(assoc % :_id (ObjectId.) :type (:type i)
                                         :pname (:pname i) :iname (:name i)
                                         :i (pr-str i) :batch_id u)
                                 l)
                           WriteConcern/JOURNAL_SAFE))
      (catch Exception e
        (mc/remove "sequences" {:batch_id u})
        (throw e)))))

;; utilities

(defn- store-read
  [h]
  (if-let [o (bs/bs-read (:src h))]
    (merge o (dissoc h :src))))

(defn- biosequence-save-helper
  [this project name type]
  (let [i (->biosequenceCollection name (:name project) type)]
    (with-open [r (bs/bs-reader this)]
      (save-list (pmap #(hash-map :acc (bs/accession %)
                                  :src (pr-str %))
                       (bs/biosequence-seq r))
                 i))
    i))

(defn collection-seq-helper
  [collection]
  (map store-read
       (mc/find-maps "sequences" {:pname (:pname collection) :iname (:name collection)})))


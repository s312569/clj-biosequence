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

(declare save-list store-read biosequence-save-helper)

;; interface

(defprotocol storeIndexIO
  (mongo-save-file [this project name]))

;; project

(defrecord mongoProject [name])

(defrecord biosequenceIndex [name pname type])

(defmethod print-method clj_biosequence.store.biosequenceIndex
  [this ^java.io.Writer w]
  (bs/print-tagged this w))

;; extend persistence

(extend-protocol storeIndexIO
  fastaFile
  (mongo-save-file [this project name]
    (biosequence-save-helper this project name "biosequence/fasta")))

;; functions

(defn mongo-connect
  []
  (mg/connect!))

(defn init-project
  [name]
  (->mongoProject name))

(defn list-projects
  []
  (mg/use-db! "clj-projects")
  (set (mc/distinct "sequences" :pname)))

(defn drop-project
  [project]
  (mg/use-db! "clj-projects")
  (mc/remove "sequences" {:pname (:name project)}))

(defn list-indexes
  [project]
  (mg/use-db! "clj-projects")
  (map #(hash-map % (:type (mc/find-one-as-map "sequences" {:iname %} [:type])))
       (mc/distinct "sequences" :iname)))

(defn get-index
  [name project]
  (mg/use-db! "clj-projects")
  (bs/bs-read (:i (mc/find-one-as-map "sequences"
                                      {:pname (:name project) :iname name}))))

(defn drop-index
  [index]
  (mc/remove "sequences" {:iname (:name index) :pname (:pname index)}))

(defn index-seq
  [index]
  (map store-read
       (mc/find-maps "sequences" {:pname (:pname index) :iname (:name index)})))

(defn index-file
  [file alphabet project name]
  (index-source (bs/init-fasta-file file alphabet) project name))

(defn save-list
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
  (let [i (->biosequenceIndex name (:name project) type)]
    (with-open [r (bs/bs-reader this)]
      (save-list (pmap #(hash-map :acc (bs/accession %)
                                  :src (pr-str %))
                       (bs/biosequence-seq r))
                 i))
    i))

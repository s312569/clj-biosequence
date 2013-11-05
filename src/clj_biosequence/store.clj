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

;; interface

(defprotocol storeItemIO
  (save-rep [this]))

(defprotocol storeIndexIO
  (save-index [this project name]))

;; project

(defrecord mongoProject [name])

(defrecord biosequenceIndex [name pname type])

(defmethod print-method clj_biosequence.store.biosequenceIndex
  [this ^java.io.Writer w]
  (bs/print-tagged this w))

;; extend persistence

(extend-protocol storeItemIO

 fastaSequence

  (save-rep [this]
    (hash-map :acc (bs/accession this) :src (pr-str this) :_id (:_id this))))

(extend-protocol storeIndexIO

  fastaFile

  (save-index [this project name]
    (->biosequenceIndex name (:name project) "biosequence/fasta")))

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
  (bs/bs-read {:src (:i (mc/find-one-as-map "sequences"
                                            {:pname (:name project) :iname name}))}))

(defn drop-index
  [index]
  (mc/remove "sequences" {:iname (:name index) :pname (:pname index)}))

(defn index-source
  [src project name]
  (mg/use-db! "clj-projects")
  (let [h (save-index src project name)
        u (mu/random-uuid)]
    (with-open [r (bs/bs-reader src)]
      (try
        (do (mc/ensure-index "sequences"
                             (array-map :acc 1 :pname -1 :iname -1)
                             {:unique true :sparse true})
            (mc/insert-batch "sequences"
                             (pmap #(assoc (save-rep %)
                                      :_id (ObjectId.) :type (:type h)
                                      :pname (:pname h) :iname (:name h)
                                      :i (pr-str h)
                                      :batch_id u)
                                   (bs/biosequence-seq r))
                             WriteConcern/JOURNAL_SAFE))
        (catch Exception e
          (mc/remove "sequences" {:batch_id u})
          (throw e))))
    h))

(defn index-seq
  [index]
  (map bs/bs-read
       (mc/find-maps "sequences" {:pname (:pname index) :iname (:name index)})))

(defn index-file
  [file alphabet project name]
  (index-source (bs/init-fasta-file file alphabet) project name))

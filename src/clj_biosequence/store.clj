(ns clj-biosequence.store
  (:require [monger.core :as mg]
            [monger.collection :as mc]
            [monger.conversion :as con]
            [monger.result :as re]
            [monger.db :as mdb]
            [monger.util :as mu]
            [clj-biosequence.core :as bs]
            [fs.core :refer [base-name]])
  (:import [com.mongodb MongoOptions ServerAddress WriteConcern]
           [org.bson.types ObjectId]
           [clj_biosequence.core fastaSequence]))

(declare store-read all-bs-collections)

;; interface

(defprotocol mongoBSRecordIO
  (mongo-bs-save [this pname cname]))

(defprotocol mongoBSCollectionIO
  (save-bs-collection [this coll])
  (run-type [this])
  (drop-project [this])
  (drop-collection [this]))

;; connect

(defn mongo-connect
  []
  (mg/connect!))

;; run

(defrecord mongoBSCollection [pname cname]

  mongoBSCollectionIO

  (save-bs-collection [this coll]
    (if (or (nil? (:cname (first coll)))
            (= (:cname this) (:cname (first coll)))
            (not (some #(= % this) (all-bs-collections))))
      (do (mg/use-db! "clj-biosequence")
          (mc/ensure-index "sequences" (array-map :pname 1 :cname 1 :acc -1)
                           {:unique true})
          (let [b (mu/random-uuid)]
            (try
              (do (dorun (->> (map #(mongo-bs-save % (:pname this) (:cname this)) coll)
                              (pmap #(mc/save "sequences" (assoc % :batch b)
                                              WriteConcern/SAFE))))
                  this)
              (catch Exception e
                (mc/remove "sequences" {:batch b})
                (throw e)))))
      (throw (Throwable. "Illegal operation"))))

  (run-type [this]
    (mc/distinct "sequences" :type {:pname (:pname this)
                                    :cname (:ename this)}))

  (drop-project [this]
    (mg/use-db! "clj-biosequence")
    (mc/remove "sequences" {:pname (:pname this)}))

  (drop-collection [this]
    (mg/use-db! "clj-biosequence")
    (mc/remove "sequences"
             {:pname (:pname this)
              :cname (:cname this)}))

  bs/biosequenceReader

  (biosequence-seq [this]
    (mg/use-db! "clj-biosequence")
    (map store-read (mc/find-maps "sequences" {:pname (:pname this) :cname (:cname this)
                                               :element "sequence"}))))

(defn init-mongo-collection
  [pname cname]
  (->mongoBSCollection pname cname))

;; MS interface

(defn save-biosequence-file
  [file pname]
  (let [run (init-mongo-collection pname (base-name (bs/bs-path file)))]
    (with-open [r (bs/bs-reader file)]
      (save-bs-collection run (bs/biosequence-seq r)))))

(defn all-biosequence-collections
  []
  (mg/use-db! "clj-biosequence")
  (map (fn [{:keys [pname cname]}]
         (init-mongo-collection pname cname))
       (distinct (map #(dissoc % :_id)
                      (mc/find-maps "sequences" {:cname {"$regex" ".*"}}
                                    [:pname :cname])))))

(defn biosequence-collections
  [pname cname]
  (mg/use-db! "clj-biosequence")
  (->> (all-biosequence-collections)
       (filter #(= (:pname %) pname))))

(defn biosequence-projects
  []
  (mg/use-db! "clj-biosequence")
  (->> (all-biosequence-collections)
       (map :pname)
       distinct))

;; utilities

(defn store-read
  [h]
  (if-let [o (bs/bs-thaw (:src h))]
    (if (string? o)
      o
      (merge o (dissoc h :src)))))

(defn unique-uuid
  []
  (mu/random-uuid))

;; extending fasta

(extend-protocol mongoBSRecordIO

  fastaSequence

  (mongo-bs-save [this pname cname]
    (let [s (hash-map :acc (bs/accession this) :element "sequence"
                      :pname pname :cname cname
                      :type "biosequence/fasta" :src (bs/bs-freeze this))]
      (if (:_id this)
        (assoc s :_id (:_id this))
        s))))

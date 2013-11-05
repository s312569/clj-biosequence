(ns clj-biosequence.store
  (:require [fs.core :as fs]
            [clojure.data.xml :as xml]
            [clojure.java.io :as io]
            [monger.core :as mg]
            [monger.collection :as mc]
            [monger.conversion :as con]
            [monger.db :as mdb])
  (:import [com.mongodb MongoOptions ServerAddress]
           [org.bson.types ObjectId]))

(declare prep-obj)

;; db interactions

(defn connect
  []
  (mg/connect!))

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
  "Saves a list of records `l` into database `d` and collection `c`."
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

;; housekeeping

(defn register-project
  "Registers a project as a clj-biosequence project."
  [name]
  (mg/use-db! "clj-projects")
  (let [r (mc/insert-and-return "projects" {:name name})]
    (mc/ensure-index "clj-projects" (array-map :name 1)
                     {:unique true})
    name))

(defn register-index
  "Registers an index in a project."
  [project name type]
  )

(defn list-projects
  "Returns a set of clj-biosequence projects on a server."
  []
  (mg/use-db! "clj-projects")
  (set (map :name (mc/find-maps "projects"))))

(defn list-indexes
  "Returns a set of indexes in a project."
  [project]
  (mg/use-db! project)
  (mc/find-maps "indexes"))

;; project

(defrecord biosequenceProject [name])

(defn init-project
  "Returns a new project with specified name."
  [name]
  (->biosequenceProject (ps/register-project name)))

(defn init-index
  [db collection type]
  (mg/use-db! project)
  (let [r (mc/insert-and-return "indexes" {:name collection :type type})]
    (mc/ensure-index "indexes" (array-map :name 1)
                     {:unique true})
    r))

(defn get-index
  "Reteieves a biosequenceIndex from a biosequenceProject."
  [bp name]
  (let [i (first (filter (fn [{n :name t :type}]
                           (and (= name n) (= t "biosequence")))
                         (ps/list-indexes (:name bp))))]
    (if i
      (->biosequenceIndex (:name i) (:name bp) (:type i)))))

(defn list-all-projects
  "Lists all projects."
  []
  (ps/list-projects))



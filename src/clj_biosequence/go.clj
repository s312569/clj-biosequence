(ns clj-biosequence.go
  (:require [clojure.java.io :refer [writer reader]]
            [fs.core :refer [file? absolute-path]]
            [clojure.string :refer [trim split]]
            [clojure.walk :refer [keywordize-keys]]
            [clojure.zip :as zip]
            [clj-biosequence.indexing :as ind]
            [clj-biosequence.core :as bs]))

(declare init-indexed-obo)

;; term

(defn stanza-type
  [s]
  (:type s))

(defn term-name
  [t]
  (:name (:info t)))

(defn term-namespace
  [t]
  (:namespace (:info t)))

(defn alternate-ids
  [t]
  (:alt_id (:info t)))

(defn is-anonymous
  [t]
  (:is_anonymous (:info t)))

(defn term-comment
  [t]
  (:comment (:info t)))

(defn subsets
  [t]
  (:subset (:info t)))

(defn synonyms
  [t]
  (:synonym (:info t)))

(defn xrefs
  [t]
  (:xref (:info t)))

(defn is-a-relations
  [t]
  (map #(trim (first (split % #"!")))
       (:is_a (:info t))))

(defn intersection-of-relations
  [t]
  (:intersection_of (:info t)))

(defn union-of-relations
  [t]
  (:union_of (:info t)))

(defn disjoint-from-relations
  [t]
  (:disjoint_from (:info t)))

(defn obsolete?
  [t]
  (:is_obsolete (:info t)))

(defn replace-by-terms
  [t]
  (:replaced_by (:info t)))

(defn consider-terms
  [t]
  (:consider (:info t)))

(defn created-by
  [t]
  (:created_by (:info t)))

(defn created-date
  [t]
  (:creation_date (:info t)))

(defn term-seq
  [obo-reader]
  (filter #(= "Term" (stanza-type %))
          (bs/biosequence-seq obo-reader)))

(defrecord goStanza [id type info]

  bs/Biosequence

  (accession [this]
    (:id this))

  (accessions [this]
    (cons (:id this) (alternate-ids this)))

  (def-line [this]
    (:def (:info this))))

;; reader

(defn- tokenise
  [m]
  (let [l (:remaining m)
        r (rest (->> (drop-while #(not (re-find #"^\[.+\]" %)) l)
                     (drop-while #(not (= "" (trim %))))))
        y (->> (drop-while #(not (re-find #"^\[.+\]" %)) l)
               (take-while #(not (= "" (trim %))))
               (filter #(not (re-find #"^!" %))))]
    (if (seq y)
      {:yield y :remaining r}
      {:end true})))

(defn- get-id [l]
  (-> (some #(re-find #"^id:.+" %) l)
      (split #":" 2)
      second
      trim))

(defn- get-info [l]
  (keywordize-keys
   (apply merge-with (fn [x y] (cond (and (string? x) (string? y))
                                     (vector x y)
                                     (string? x)
                                     (conj y x)
                                     (string? y)
                                     (conj x y)
                                     :else
                                     (vec (concat x y))))
          (->> (remove #(or (re-find #"^\[" %) (re-find #"^id:" %)) l)
               (map #(->> (split % #":" 2)
                          (map trim)
                          vec))
               (map #(apply assoc {} %))))))

(defn- compliance
  [s]
  (condp = (stanza-type s)
    "Term" (let [i (:info s)
                 r (remove #(not (string? (get i %)))
                           [:subset :synonym :xref :is_a :intersection_of :union_of
                            :disjoint_from :replaced_by :consider :relationship])
                 cs (if (seq r)
                      (assoc s :info
                             (apply assoc i (mapcat #(vector % (vector (get i %))) r)))
                      s)]
             (if-not (every? #(or (string? %) (nil? %))
                             (vector (:id cs) (:def (:info cs)) (:comment (:info cs))))
               (throw (Throwable. (str "Duplicate tags in term "  (:id s)))))
             (if-not (every? #(= 0 (mod (count (get i %)) 2))
                             (map #(% (:info cs)) [:intersection_of :union_of]))
               (throw (Throwable. (str "Odd number of tags error in term " (:id s)))))
             cs)
    s))

(defn- parse-obo [l]
  (->> {:remaining l}
       (iterate tokenise)
       (take-while #(not (contains? % :end)))
       (map :yield)
       (filter #(not (nil? %)))
       (map #(->goStanza (get-id %) (trim (second (re-find #"^\[(.+)\]" (first %))))
                         (get-info %)))
       (map compliance)))

(defrecord goOBOReader [strm path]

  bs/biosequenceReader

  (biosequence-seq [this]
    (let [l (line-seq (:strm this))]
      (parse-obo l)))

  bs/biosequenceFile

  (bs-path [this]
    (absolute-path (:path this)))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defn- init-obo-reader [strm path]
  (->goOBOReader strm path))

;; file

(defrecord goOBOFile [path]

  bs/biosequenceIO

  (bs-reader [this]
    (init-obo-reader (reader (bs/bs-path this))
                     (bs/bs-path this)))

  bs/biosequenceFile

  (bs-path [this]
    (absolute-path (:path this)))

  (index-file [this]
    (init-indexed-obo (bs/bs-path this)))

  (index-file [this ofile]
    (init-indexed-obo (absolute-path ofile))))

(defn init-obo-file [path]
  (->goOBOFile path))

;; indexed file

(defrecord indexedOBOReader [index strm path]

  bs/biosequenceReader

  (biosequence-seq [this]
    (bs/indexed-seq this map->goStanza))

  (get-biosequence [this accession]
    (bs/get-object this accession map->goStanza))

  bs/biosequenceFile

  (bs-path [this]
    (absolute-path (:path this)))

  java.io.Closeable

  (close [this]
    (bs/close-index-reader this)))

(defrecord indexedOBOFile [index path]

  bs/biosequenceIO

  (bs-reader [this]
    (->indexedOBOReader (:index this)
                        (bs/open-index-reader (:path this))
                        (bs/bs-path this)))

  bs/biosequenceFile

  (bs-path [this]
    (absolute-path (:path this)))

  (empty-instance [this path]
    (init-indexed-obo path)))

(defn init-indexed-obo [file]
  (->indexedOBOFile {} file))

(defmethod print-method clj_biosequence.go.indexedOBOFile
  [this w]
  (bs/print-tagged-index this w))

;; graph

(defrecord goGraph [hash])

(defn branch?
  [node]
  true)

(defn children
  [node]
  (seq (map ) (is-a-relations node)))

(defn make-node
  [node children]
  )

(defn get-node
  [g id]
  (get (:hash g) id))

(defn- namespace?
  [n]
  (condp = (term-namespace n)
    "biological_process" "GO:0008150"
    "molecular_function" "GO:0003674"
    "cellular_component" "GO:0005575"))

(defn get-path
  [g id]
  (let [b (fn [n] true)
        c (fn [n] (seq (map #(get-node g %) (is-a-relations n))))
        mn (fn [n c]
             (->goStanza (:id n) (:type n) (assoc (:info n)
                                                :is_a (vec c))))
        z (zip/zipper b c mn (get-node g id))
        ns (namespace? (get-node g id))]
    (loop [l z
           a []]
      (cond (= (:id (zip/node l)) ns)
            (recur (zip/next l)
                   (conj a (conj (vec (map :id (zip/path l))) ns)))
            (zip/end? l)
            (if (= (:id (zip/node l)) ns)
              (conj a (conj (vec (map :id (zip/path l))) ns))
              a)
            :else
            (recur (zip/next l) a)))))

(defn reduce-to-level
  [paths level]
  (into {}
        (->> (map #(take level (reverse %)) paths)
             (filter #(= (count %) level))
             frequencies
             (map (fn [[k v]] (vector (nth k (- level 1)) v))))))

(defn id-to-term
  [g t]
  (term-name (get-node g t)))

(defn in-mem-go-graph
  [obo-reader]
  (->goGraph (into {} (map #(vector (bs/accession %) %)
                           (bs/biosequence-seq obo-reader)))))

;; (defn add-node
;;   [g n]
;;   (if (get (:hash g) (bs/accession n))
;;     g
;;     (assoc g :hash
;;            (assoc (:hash g) (bs/accession n) n))))

;; (defn remove-node
;;   [g n]
;;   (assoc g :hash))

;; (defn remove-is-a
;;   [g n id]
;;   (assoc g :hash (-> g
;;                      (add-node n)
;;                      :hash
;;                      (update-in [(bs/accession n) :info :is_a]
;;                                 (fn [v n] (vec (remove #(= % n) v)))
;;                                 id))))

;; (defn add-is-a
;;   [g n id]
;;   (assoc g :hash (-> g
;;                      (add-node n)
;;                      :hash
;;                      (update-in [(bs/accession n) :info :is_a]
;;                                 (fn [v n] (vec (conj v n))) id))))

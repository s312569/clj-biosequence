(ns clj-biosequence.store
  (:require [clj-biosequence.core :as bs]
            [clj-biosequence.persistence :as ps]))

(declare prep-obj)

(defrecord storeReader [cursor]

  bs/biosequenceReader

  (biosequence-seq [this]
    (ps/record-seq (:cursor this)))

  java.io.Closeable

  (close [this]
    nil))

(defrecord biosequenceStore [name]

  bs/biosequenceIO

  (bs-reader [this]
    (->storeReader (ps/find-all (:name this) "sequence"))))

(defn save-biosequences
  "Takes a list of biosequences and saves them to the biosequenceStore
  `s`."
  [lst s]
  (ps/save-records (pmap prep-obj lst) (:name s) "sequence"))

(defn index-biosequence-file
  "Takes a biosequence file, a project designation and a name for the
  index and returns a biosequenceStore."
  [file name]
  (let [s (->biosequenceStore name)]
    (do
      (with-open [r (bs/bs-reader file)]
        (save-biosequences (bs/biosequence-seq r) s))
      s)))

(defn init-store
  "Returns a new biosequenceStore with the specified name."
  [name]
  (->biosequenceStore name))

(defn update-biosequence
  "Updates a biosequence in the store `s`."
  [bs s]
  (ps/update-record (prep-obj bs) (:name s) "sequence"))

(defn get-biosequence
  "Returns a biosequence from the store, `s`, with the accession,
  `a`."
  [a s]
  (ps/get-record {:acc a} (:name s) "sequence"))

;; utilities

(defn- prep-obj
  [o]
  (hash-map :acc (bs/accession o) :src (pr-str o) :_id (:_id o)))

(ns clj-biosequence.store
  (:require [clj-biosequence.core :as bs]
            [clj-biosequence.persistence :as ps]))

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
    (->storeReader (ps/find-all "clj-biosequence" (:name this)))))

(defn init-store
  [name]
  (->biosequenceStore name))

(defn update-biosequence
  [bs s]
  (ps/update-record (bs/bs-save bs) "clj-biosequence" (:name s)))

(defn save-biosequences
  [lst s]
  (ps/save-records (pmap bs/bs-save lst) "clj-biosequence" (:name s)))

(defn get-biosequence
  [a s]
  (ps/get-record a "clj-biosequence" (:name s)))

(defn index-biosequence-file
  [file name]
  (let [s (->biosequenceStore name)]
    (do
      (with-open [r (bs/bs-reader file)]
        (save-biosequences (bs/biosequence-seq r) s))
      s)))

(defn bs-collections
  []
  (ps/get-collections "clj-biosequence"))

(defn delete-store
  [st]
  (ps/drop-collection "clj-biosequence" (:name st)))

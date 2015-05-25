(ns clj-biosequence.indexing
  (:require [taoensso.nippy :refer [freeze thaw]]
            [clojure.edn :as edn]
            [fs.core :refer [exists? delete]])
  (:import [java.io RandomAccessFile]))

(defn- write-and-position
  [obj writer acc]
  (let [o (freeze obj)
        off (.getFilePointer writer)]
    (.write writer o)
    (vector acc (list off (count o)))))

(defn- read-one
  [off len reader]
  (let [bb (byte-array len)]
    (.seek reader off)
    (.read reader bb)
    (thaw bb)))

(defn close-index-reader
  [reader]
  (.close ^RandomAccessFile reader))

(defn load-indexed-file
  [path]
  (edn/read-string (slurp (str path ".idx"))))

(defn index-objects
  [writer objs id-func]
  (->> (map #(write-and-position % writer (id-func %)) objs)
       (into {})))

(defn index-reader
  [path]
  (RandomAccessFile. (str path ".bin") "r"))

(defn index-writer
  [path]
  (let [w (RandomAccessFile. (str path ".bin") "rw")]
    (.seek w (.length w))
    w))

(defn object-seq
  [reader indexes]
  (map (fn [[o l]] (read-one o l reader)) indexes))

(defn save-index
  [path index]
  (binding [*print-length* false]
    (spit (str path ".idx") (pr-str index))))

(defn delete-index
  [path]
  (if (exists? (str path ".bin"))
    (delete (str path ".bin")))
  (if (exists? (str path ".idx"))
    (delete (str path ".idx"))))

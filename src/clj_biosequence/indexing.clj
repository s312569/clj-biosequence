(ns clj-biosequence.indexing
  (:require [taoensso.nippy :refer [freeze thaw]]
            [clojure.edn :as edn]
            [miner.tagged :refer [pr-tagged-record-on tag->factory some-tag-reader-fn]]
            [clojure.data.xml :as xml]
            [fs.core :as fs])
  (:import [java.io RandomAccessFile]))

(defrecord indexWriter [strm]

  java.io.Closeable
  
  (close [this]
    (.close ^java.io.RandomAccessFile (:strm this))))

(defn- init-index-writer
  [path]
  (->indexWriter (RandomAccessFile. (str path ".bin") "rw")))

(defn- write-and-position
  [obj writer acc]
  (let [o (freeze obj)
        off (.getFilePointer (:strm writer))]
    (.write (:strm writer) o)
    (vector acc (list off (count o)))))

(defn- read-one
  [off len reader]
  (let [bb (byte-array len)]
    (.seek reader off)
    (.read reader bb)
    (thaw bb)))

(defn- my-record-tag-reader
  [tag val]
  (when-let [factory (and (map? val)
                          (tag->factory tag))]
    (factory val)))

(def my-tagged-default-reader 
  (some-tag-reader-fn my-record-tag-reader))

(defn close-index-reader
  [reader]
  (.close ^RandomAccessFile reader))

(defn print-tagged
  [o w]
  (pr-tagged-record-on o w))

(defn load-indexed-file
  [path]
  (edn/read-string {:default my-tagged-default-reader
                    :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
                   (slurp (str path ".idx"))))

(defn index-objects
  [writer objs id-func]
  (->> (map #(write-and-position % writer (id-func %)) objs)
       (into {})))

(defn index-reader
  [path]
  (RandomAccessFile. (str path ".bin") "r"))

(defn index-writer
  [path]
  (init-index-writer path))

(defn object-seq
  [reader indexes func]
  (map (fn [[o l]] (func (read-one o l reader)))
       indexes))

(defn save-index
  [path index]
  (spit (str path ".idx") (pr-str index)))

(defn delete-index
  [path]
  (if (fs/exists? (str path ".bin"))
    (fs/delete (str path ".bin")))
  (if (fs/exists? (str path ".idx"))
    (fs/delete (str path ".idx"))))

(ns clj-biosequence.index
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [taoensso.nippy :as nip])
  (:import [java.io RandomAccessFile]))

(declare init-file-index write-and-position read-one)

;; reader

(defrecord indexReader [strm]

  java.io.Closeable
  
  (close [this]
    (.close ^java.io.RandomAccessFile (:strm this))))

(defn init-index-reader
  [strm]
  (->indexReader strm))

;; indexed file

(defrecord indexFile [path]

  bs/biosequenceFile

  (bs-path [this]
    (:path this))

  bs/biosequenceIO

  (bs-reader [this]
    (init-index-reader (RandomAccessFile. (:path this) "rw"))))

(defn init-index-file
  [path]
  (->indexFile path))

;; index

(defrecord fileIndex [index]

  bs/biosequenceReader

  (biosequence-seq [this]
    (map (fn [[o l]]
           (read-one o l (bs/bs-path (:file this))))
         (vals (:index this))))

  (get-biosequence [this accession]
    (let [[o l] (get (:index this) accession)]
      (if o
        (read-one o l (bs/bs-path (:file this)))))))

(defn save-index
  [idx]
  (let [o (str (bs/bs-path (:file idx)) ".idx")]
    (with-open [out (bs/bs-reader (init-index-file o))]
      (.write (:strm out) (bs/bs-freeze idx)))))

(defn load-index
  [file]
  (let [l (.length (io/file file))
        bb (byte-array l)]
    (with-open [r (RandomAccessFile. file "rw")]
      (.read r bb)
      (bs/bs-thaw bb))))

;; creation

(defn index-biosequence-file
  [file]
  (let [ofile (init-index-file (str (bs/bs-path file) ".bin"))]
    (with-open [o (bs/bs-reader ofile)]
      (with-open [i (bs/bs-reader file)]
        (assoc (->> (bs/biosequence-seq i)
                    (map #(write-and-position % (:strm o)))
                    (into {})
                    init-file-index)
          :file ofile)))))

(defn index-biosequence-multi-file
  [files out]
  (let [ofile (init-index-file (str (fs/absolute-path out) ".bin"))]
    (with-open [o (bs/bs-reader ofile)]
      (-> (apply merge (doall (map
                                #(with-open [r (bs/bs-reader %)]
                                   (->> (bs/biosequence-seq r)
                                        (map (fn [x] (write-and-position x (:strm o))))
                                        (into {})))
                                files)))
          init-file-index
          (assoc :file ofile)))))

(defn index-biosequence-list
  [lst outfile]
  (let [ofile (init-index-file (str (fs/absolute-path outfile) ".bin"))]
    (with-open [o (bs/bs-reader ofile)]
      (assoc (->> lst
                  (map #(write-and-position % (:strm o)))
                  (into {})
                  init-file-index)
        :file ofile))))

;; utility

(defn- init-file-index
  [h]
  (->fileIndex h))

(defn- read-one
  [off len file]
  (let [bb (byte-array len)]
    (with-open [r (RandomAccessFile. file "rw")]
      (.seek r off)
      (.read r bb)
      (bs/bs-thaw bb))))

(defn- write-and-position
  [obj strm]
  (let [o (nip/freeze obj)
        off (.getFilePointer strm)]
    (.write strm o)
    (vector (bs/accession obj) (list off (count o)))))



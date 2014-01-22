(ns clj-biosequence.index
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [clojure.edn :as edn]
            [taoensso.nippy :as nip])
  (:import [java.io RandomAccessFile]))

(declare init-file-index write-and-position read-one)

;; reader

(defrecord indexWriter [strm]

  java.io.Closeable
  
  (close [this]
    (.close ^java.io.RandomAccessFile (:strm this))))

(defn init-index-writer
  [strm]
  (->indexWriter strm))

;; indexed file

(defrecord indexFile [path]

  bs/biosequenceFile

  (bs-path [this]
    (:path this)))

(defn bs-writer [this]
  (init-index-writer (RandomAccessFile. (:path this) "rw")))

(defn init-index-file
  [path]
  (->indexFile (fs/absolute-path path)))

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

;; creation

(defn index-biosequence-file
  [file]
  (let [ofile (init-index-file (str (bs/bs-path file) ".bin"))]
    (with-open [o (bs-writer ofile)]
      (let [i (with-open [i (bs/bs-reader file)]
                (assoc (->> (bs/biosequence-seq i)
                            (map #(write-and-position % (:strm o)))
                            (into {})
                            init-file-index)
                  :file ofile))]
        (spit (str (bs/bs-path file) ".idx") (:index i))
        i))))

(defn index-biosequence-multi-file
  [files out]
  (let [ofile (init-index-file (str (fs/absolute-path out) ".bin"))]
    (with-open [o (bs-writer ofile)]
      (let [i (-> (apply merge (doall (map
                                       #(with-open [r (bs/bs-reader %)]
                                          (->> (bs/biosequence-seq r)
                                               (map (fn [x] (write-and-position x (:strm o))))
                                               (into {})))
                                       files)))
                  init-file-index
                  (assoc :file ofile))]
        (spit (str (fs/absolute-path out) ".idx") (:index i))
        i))))

(defn index-biosequence-list
  [lst outfile]
  (let [ofile (init-index-file (str (fs/absolute-path outfile) ".bin"))
        i (with-open [o (bs-writer ofile)]
            (assoc (->> lst
                        (map #(write-and-position % (:strm o)))
                        (into {})
                        init-file-index)
              :file ofile))]
    (spit (str (fs/absolute-path outfile) ".idx") (:index i))
    i))

(defn load-biosequence-index
  [ind]
  (assoc (init-file-index (edn/read-string (slurp (str (fs/absolute-path ind) ".idx"))))
    :file (init-index-file (str (fs/absolute-path ind) ".bin"))))

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

(defn write-and-position
  [obj strm]
  (let [o (nip/freeze obj)
        off (.getFilePointer strm)]
    (.write strm o)
    (vector (bs/accession obj) (list off (count o)))))



(ns clj-biosequence.sam
  (:require [clojure.string :refer [split]]
            [clojure.java.io :refer [reader]]
            [fs.core :refer [file?]]
            [clojure.zip :as zip]
            [clojure.edn :as edn]
            [clj-biosequence.core :as bs]
            [instaparse.core :as insta]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; alignment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord samAlignment [qname flag rname pos mapq cigar
                         rnext pnext tlen seq qual])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; IO
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- hex-byte-array
  [s]
  (mapv #(format "0x%x" (Integer/parseInt (apply str %) 16))
        (partition 2 s)))

(defn- parse-alignment
  [line]
  (let [f (split line #"\t")]
    (-> (apply ->samAlignment (take 11 f))
        (assoc :flag (Integer/parseInt (second f))
               :pos (Integer/parseInt (nth f 3))
               :mapq (Integer/parseInt (nth f 4))
               :pnext (Integer/parseInt (nth f 7))
               :tlen (Integer/parseInt (nth f 8))
               :optional
               (into {}
                     (map #(let [f (split % #":")
                                 v (last f)]
                             (vector (first f)
                                     (case (second f)
                                       "i" (Integer/parseInt v)
                                       "f" (Float/parseFloat v)
                                       "H" (hex-byte-array v)
                                       v)))
                          (drop 11 f)))))))

(defrecord samReader [strm]
  bs/biosequenceReader
  (biosequence-seq [this]
    (->> (line-seq (:strm this))
         (drop-while #(= \@ (first %)))
         (filter #(not (= "" %)))
         (map parse-alignment)))
  java.io.Closeable
  (close [this]
    (.close (:strm this)))
  bs/biosequenceParameters
  (parameters [this]
    (:header this)))

(defn init-sam-reader
  "Initializes a samReader."
  [strm]
  (->samReader strm))

(defn tag-value=
  [line tag value]
  (->> (rest (rest line))
       (filter #(= tag (second %)))
       (filter #(= value (last %)))
       seq))

(defn tag-value
  [tag line]
  (->> (rest (rest line))
       (filter #(= tag (second %)))
       (map last)))

(defn sorted
  [reader]
  (->> (bs/parameters reader)
       (filter #(= "HD" (second %)))
       first
       (tag-value "SO")
       first))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- get-header
  [f]
  (with-open [r (reader f)]
    (let [p (insta/parser
             "line= <#'@'> #'(HD|SQ|RG|PG)' tag+ |
                    <#'@'> #'(CO|[a-z][a-z])' text
              tag= <#'\\t'> #'[A-Za-z][A-Za-z0-9]' <#':'> #'[ -~]+'
              text= <#'\\t'> #'.*'")]
      (doall
       (->> (line-seq r)
            (take-while #(= \@ (first %)))
            (filter #(not (= "" %)))
            (map p))))))

(defrecord samFile [file])

(extend samFile
  bs/biosequenceIO
  {:bs-reader
   (fn [this]
     (assoc (->samReader (reader (bs/bs-path this)))
       :header (get-header (bs/bs-path this))))}
  bs/biosequenceFile
  bs/default-biosequence-file)

(defn init-sam-file
  "Initializes a samFile record"
  [file]
  {:pre [(file? file)]}
  (->samFile file))

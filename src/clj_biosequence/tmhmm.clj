(ns clj-biosequence.tmhmm
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bs]
            [clojure.string :as string]
            [clojure.pprint :as pp]
            [clojure.data.xml :as xml]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; sequence
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord tmhmmProtein [accession length expaa
                         first60 predhel topology])

(extend tmhmmProtein
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession
    (fn [this]
      (:accession this))))

(defn numb-tm-domains
  [tp]
  (:predhel tp))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; reader
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- make-tmhmm-result
  [line]
  (let [[accession length expaa first60 predhel topology]
        (map #(let [spl (string/split % #"=")]
                (if (second spl)
                  (second spl)
                  (first spl)))
             (string/split line #"\s+"))]
    (->tmhmmProtein accession
                    (Integer/parseInt length)
                    (Float/parseFloat expaa)
                    (Float/parseFloat first60)
                    (Integer/parseInt predhel)
                    topology)))

(defrecord tmhmmReader [strm]
  bios/biosequenceReader
  (biosequence-seq [this]
    (map #(make-tmhmm-result %)
         (line-seq strm)))
  java.io.Closeable
  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord tmhmmFile [file])

(extend tmhmmFile
  bs/biosequenceIO
  {:bs-reader
   (fn [this]
     (->tmhmmReader (io/reader (bs/bs-path this))))}
  bs/biosequenceFile
  bs/default-biosequence-file)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; runnin tmhmm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn tmhmm
  [bs & {:keys [outfile] :or {outfile (fs/temp-file "sp-out")}}]
  {:pre [(not (some false? (map bs/protein? bs)))]}
  (let [in (bios/biosequence->file bs
                                   (fs/temp-file "sp-in")
                                   :append false)]
    (with-open [out (io/output-stream outfile)]
      (try (let [sp @(exec/sh ["tmhmm" "-short"
                               (fs/absolute-path in)] {:out out})]
             (if (= 0 (:exit sp))
               (->tmhmmFile (fs/absolute-path outfile))
               (if (:err sp)
                 (throw (Throwable. (str "TMHMM error: "
                                         (:err sp))))
                 (throw (Throwable. (str "Exception: "
                                         (:exception sp)))))))
           (finally (fs/delete in))))))

;; private

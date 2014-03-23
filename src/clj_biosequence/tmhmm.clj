(ns clj-biosequence.tmhmm
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bios]
            [clojure.string :as string]
            [clojure.pprint :as pp]
            [clojure.data.xml :as xml]))

(declare make-tmhmm-result)

(defrecord tmhmmProtein [accession length expaa first60 predhel topology]

  bios/Biosequence

  (accession [this]
    (:accession this)))

(defn numb-tm-domains
  [tp]
  (:predhel tp))

;; results

(defrecord tmhmmReader [strm]

  bios/biosequenceReader

  (biosequence-seq [this]
    (map #(make-tmhmm-result %)
         (line-seq strm)))

  (parameters [this]
    ())

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord tmhmmFile [file]

  bios/biosequenceIO

  (bs-reader [this]
    (->tmhmmReader (io/reader (:file this)))))

(defn tmhmm
  [bs & {:keys [outfile] :or {outfile (fs/temp-file "sp-out")}}]
  (let [in (bios/biosequence->file bs (fs/temp-file "sp-in")
                                   :append false
                                   :func (fn [x] (if (bios/protein? x)
                                                  (bios/fasta-string x)
                                                  (throw (Throwable. "TMHMM only analyses proteins.")))))]
    (with-open [out (io/output-stream outfile)]
      (try (let [sp @(exec/sh ["tmhmm" "-short" (fs/absolute-path in)] {:out out})]
             (if (= 0 (:exit sp))
               (->tmhmmFile (fs/absolute-path outfile))
               (if (:err sp)
                 (throw (Throwable. (str "TMHMM error: " (:err sp))))
                 (throw (Throwable. (str "Exception: " (:exception sp)))))))
           (finally (fs/delete in))))))

;; private

(defn- make-tmhmm-result
  [line]
  (let [[accession length expaa first60 predhel topology]
        (map #(let [spl (string/split % #"=")]
                (if (second spl)
                  (second spl)
                  (first spl))) (string/split line #"\s+"))]
    (->tmhmmProtein accession (Integer/parseInt length) (Float/parseFloat expaa)
                      (Float/parseFloat first60) (Integer/parseInt predhel) topology)))

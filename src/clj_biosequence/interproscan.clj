(ns clj-biosequence.interproscan
  (:require [clojure.data.xml :as xml]
            [clojure.data.zip.xml :as zf]
            [clojure.java.io :as io]
            [clj-commons-exec :as exec]
            [clojure.zip :as zip]
            [clojure.string :as string]
            [clj-biosequence.core :as bs]
            [fs.core :as fs]))

(declare run-ips)

;; matches

(defrecord interproscan [src])

;; ips protein

(defrecord interproscanProtein [src]

  bs/Biosequence

  (accession [this]
    (zf/xml1-> (zip/xml-zip (:src this))
               :xref (zf/attr :id))))

(defn match-seq [protein]
  (map #()))

;; ips searchx

(defrecord interproscanReader [strms]

  bs/biosequenceReader

  (biosequence-seq [this]
    (mapcat #(->> (:content (xml/parse %))
                  (filter (fn [x] (= :protein (:tag x))))
                  (map (fn [x] (->interproscanProtein x))))
            (:strms this)))

  (parameters [this]
    ())

  java.io.Closeable

  (close [this]
    (doseq [r (:strms this)]
      (.close ^java.io.BufferedReader r))))

(defrecord interproscanResult [files]

  bs/biosequenceIO

  (bs-reader [this]
    (->interproscanReader (map io/reader (:files this)))))

;; ips run

(defn ips
  "Runs interproscan on a list of biosequences."
  [bsl & {:keys [outfile appl lookup goterms precalc pathways]
          :or {outfile (fs/temp-file "ips")
               appl '("pfam")
               lookup true
               goterms true
               precalc false
               pathways false}}]
  (->interproscanResult
   (doall (map #(let [i (bs/biosequence->file % (fs/temp-file "seq-") :append false)]
                  (try
                    (run-ips (fs/absolute-path i) (fs/absolute-path outfile)
                             seqtype appl precalc pathways lookup goterms)
                    (finally (fs/delete i))))
               (partition-all 10000 bsl)))))

;; pfam tigrfam - hmmer3-match
;; prodom
;; 

;; utilities

(defn ips-command
  [i o s a p path l g]
  (vec
   (remove nil? (-> (list
                     "interproscan.sh" "-i" i "-o" o "-seqtype" s "-f" "XML"
                     "-appl" (apply str (interpose "," a))
                     (if (not p) "-dp")
                     (if (not path) "-pa")
                     (if (or l g) "-iprlookup")
                     (if g "-goterms"))))))

(defn- run-ips
  [i o s a p path l g]
  (try
    (let [ips @(exec/sh (ips-command i o s a p path l g))]
      (if (= 0 (:exit ips))
        o
        (throw (Throwable. (str "Interproscan error: " (:err ips))))))
    (catch Exception e
      (throw e)
      (fs/delete o))))

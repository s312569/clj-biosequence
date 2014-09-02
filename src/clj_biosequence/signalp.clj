(ns clj-biosequence.signalp
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bs]
            [clojure.string :as string]
            [clojure.pprint :as pp]
            [clojure.data.xml :as xml]))

(declare signal-command make-signal-result)

;; signalp analysis

(defrecord signalpProtein [name cmax cpos ymax ypos smax spos smean D result Dmaxcut network]

  bs/Biosequence

  (accession [this]
    (:name this)))

(defn signalp?
  [sp]
  (= "Y" (:result sp)))

;; results

(defrecord signalpReader [strm]

  bs/biosequenceReader

  (biosequence-seq [this]
    (map #(make-signal-result %)
         (drop 2 (line-seq strm))))

  (parameters [this]
    ())

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord signalpFile [file]

  bs/biosequenceIO

  (bs-reader [this]
    (->signalpReader (io/reader (:file this))))

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this))))

(defn signalp
  [bs outfile & {:keys [params] :or {params {}}}]
  (let [in (bs/biosequence->file
            bs (fs/temp-file "sp-in")
            :append false
            :func (fn [x]
                    (if (bs/protein? x)
                      (bs/fasta-string x)
                      (throw (Throwable. "SignalP only analyses proteins.")))))]
    (try
      (with-open [out (io/output-stream outfile)]
        (let [sp @(exec/sh (signal-command params in) {:out out})]
          (if (and (= 0 (:exit sp)) (nil? (:err sp)))
            (->signalpFile (fs/absolute-path outfile))
            (if (:err sp)
              (throw (Throwable. (str "SignalP error: " (:err sp))))
              (throw (Throwable. (str "Exception: " (:exception sp))))))))
      (catch Exception e
        (fs/delete outfile)
        (fs/delete in)
        (println e))
      (finally (fs/delete in)))))

(defn- register-outfile
  [a]
  (first (swap! a #(cons (fs/temp-file "sp-") %))))

(defn filter-signalp
  [bsl & {:keys [trim params] :or {trim false params {}}}]
  (let [flist (atom ())
        sps (map #(signalp % (register-outfile flist) :params params)
                 (partition-all 10000 bsl))]
    (try (let [h (->> (doall
                       (pmap #(with-open [r (bs/bs-reader %)]
                                (into {}
                                      (map (fn [x]
                                             (vector (bs/accession x)
                                                     (list (:result x)
                                                           (:cpos x))))
                                           (bs/biosequence-seq r)))) sps))
                      (apply merge))]
           (remove nil?
                   (map #(if (= "Y" (first (get h (bs/accession %))))
                           (if trim
                             (bs/init-fasta-sequence
                              (bs/accession %)
                              (bs/def-line %)
                              :iupacAminoAcids
                              (subvec
                               (bs/bs-seq %)
                               (+ 1 (second
                                     (get h (bs/accession %))))))
                             %))
                        bsl)))
         (finally
           (doall (map #(fs/delete %) @flist))))))

;; private

(defn- make-signal-result
  [line]
  (let [[name cmax cpos ymax ypos smax spos smean D result Dmaxcut network]
        (string/split line #"\s+")]
    (->signalpProtein name (read-string cmax) (read-string cpos) (read-string ymax) (read-string ypos) (read-string smax) (read-string spos) (read-string smean) (read-string D) result (read-string Dmaxcut) network)))

(defn- signal-command
  [params infile]
  (conj (->> (concat
              ["-f" "short"        
               "-t" "euk"]
              (flatten (vec params)))
             (cons "signalp")
             vec)
        (fs/absolute-path infile)))

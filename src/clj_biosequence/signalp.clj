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

(defrecord signalpProtein [name cmax cpos ymax ypos
                           smax spos smean D result
                           Dmaxcut network]

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

(defn- protein-only-error
  [x]
  (if (bs/protein? x)
    (bs/fasta-string x)
    (throw (Throwable.
            "SignalP only analyses proteins."))))

(defn signalp
  "Runs signalp on a collection of biosequences and returns a signalp
  result file. Only the first 10,000 biosequences will be analysed. For
  greater numbers see `filter-signalp` or partition collections before
  running."
  [bs outfile & {:keys [params] :or {params {}}}]
  (let [in (bs/biosequence->file (take 10000 bs) (fs/temp-file "sp-in")
                                 :append false
                                 :func protein-only-error)]
    (try
      (with-open [out (io/output-stream outfile)]
        (let [sp @(exec/sh (signal-command params in) {:out out}
                           :close-err? false)]
          (if (and (= 0 (:exit sp)) (nil? (:err sp)))
            (->signalpFile (fs/absolute-path outfile))
            (if (:err sp)
              (throw (Throwable. (str "SignalP error: " (:err sp))))
              (throw (Throwable. (str "!Exception: " (:exception sp))))))))
      (catch Exception e
        (fs/delete outfile)
        (fs/delete in)
        (println e))
      (finally (fs/delete in)))))

(defn- trimmed-fasta
  [s c]
  (bs/init-fasta-sequence (bs/accession s) (bs/def-line s)
                          :iupacAminoAcids
                          (subvec (bs/bs-seq s) (+ 1 c))))

(defn- result-hash
  [file]
  (with-open [r (bs/bs-reader file)]
    (->> (map #(vector (bs/accession %)
                       (vector (:result %)
                               (:cpos %)))
              (bs/biosequence-seq r))
         (into {}))))

(defn filter-signalp
  [bsl & {:keys [trim params] :or {trim false params {}}}]
  (let [sps (map #(signalp % (fs/temp-file "sp-") :params params)
                 (partition-all 10000 bsl))]
    (try (let [h (->> (pmap result-hash sps)
                      (apply merge))]
           (->> (map #(let [r (get h (bs/accession %))]
                        (if (= "Y" (first r))
                          (if trim
                            (trimmed-fasta % (second r))
                            %)))
                     bsl)
                (remove nil?)))
         (finally
           (doall (map #(fs/delete (bs/bs-path %)) sps))))))

;; private

(defn- make-signal-result
  [line]
  (let [[name cmax cpos ymax ypos smax spos smean D result Dmaxcut network]
        (string/split line #"\s+")]
    (->signalpProtein name (Float/parseFloat cmax)
                      (Integer/parseInt cpos)
                      (Float/parseFloat ymax)
                      (Integer/parseInt ypos)
                      (Float/parseFloat smax)
                      (Integer/parseInt spos)
                      (Float/parseFloat smean)
                      (Float/parseFloat D) result
                      (Float/parseFloat Dmaxcut)
                      network)))

(defn- signal-command
  [params infile]
  (let [ak ["-s" "-t" "-u" "-U" "-M" "-c"]
        up (flatten (vec (select-keys params ak)))]
    (doseq [w (keys (apply dissoc params ak))]
      (println (str "warning: " w " is not a signalp argument - ignored.")))
    (-> (cons "signalp" up)
        vec
        (conj (fs/absolute-path infile)))))

(ns clj-biosequence.signalp
  (:require [clojure.java.io :refer [writer reader output-stream]]
            [fs.core :as fs]
            [clj-commons-exec :refer [sh]]
            [clj-biosequence.core :as bs]
            [clojure.string :refer [split]]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; signalp analysis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord signalpProtein [name cmax cpos ymax ypos
                           smax spos smean D result
                           Dmaxcut network])

(extend signalpProtein
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession
    (fn [this]
      (:name this))))

(defn signalp?
  [sp]
  (= "Y" (:result sp)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; reader
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- make-signal-result
  [line]
  (let [[name cmax cpos ymax ypos smax spos
         smean D result Dmaxcut network]
        (split line #"\s+")]
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

(defrecord signalpReader [strm]
  bs/biosequenceReader
  (biosequence-seq [this]
    (map #(make-signal-result %)
         (drop 2 (line-seq strm))))
  java.io.Closeable
  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord signalpFile [file])

(extend signalpFile
  bs/biosequenceIO
  {:bs-reader
   (fn [this]
     (->signalpReader (reader (:file this))))}
  bs/biosequenceFile
  bs/default-biosequence-file)

(defn init-signalp-result
  [file]
  (->signalpFile file))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; run signalp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- signal-command
  [params infile]
  (let [ak ["-s" "-t" "-u" "-U" "-M" "-c"]
        up (flatten (vec (select-keys params ak)))]
    (doseq [w (keys (apply dissoc params ak))]
      (println (str "warning: " w
                    " is not a signalp argument - ignored.")))
    (-> (cons "signalp" up)
        vec
        (conj (fs/absolute-path infile)))))

(defn signalp
  "Runs signalp on a collection of biosequences and returns a signalp
  result file. Only the first 10,000 biosequences will be
  analysed. For greater numbers see `filter-signalp` or partition
  collections before running."
  [bs outfile & {:keys [params] :or {params {}}}]
  {:pre [(not (some false? (map bs/protein? bs)))]}
  (let [in (bs/biosequence->file (take 10000 bs)
                                 (fs/temp-file "sp-in")
                                 :append false)]
    (try
      (with-open [out (output-stream outfile)]
        (let [sp @(sh (signal-command params in) {:out out}
                           :close-err? false)]
          (if (and (= 0 (:exit sp)) (nil? (:err sp)))
            (->signalpFile (fs/absolute-path outfile))
            (if (:err sp)
              (throw (Throwable. (str "SignalP error: "
                                      (:err sp))))
              (throw (Throwable. (str "!Exception: "
                                      (:exception sp))))))))
      (catch Exception e
        (fs/delete outfile)
        (fs/delete in)
        (println e)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; filter
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- result-hash
  [file]
  (with-open [r (bs/bs-reader file)]
    (->> (map #(vector (bs/accession %)
                       (vector (:result %)
                               (:cpos %)))
              (bs/biosequence-seq r))
      (into {}))))

(defn- sp-outfile
  [p]
  (if p
    (let [c (atom 0)]
      (fn []
        (str p "-" (swap! c inc))))
    (fn [] (fs/temp-file "sp-"))))

(defn- sp-combine-results
  [outfile sps]
  (with-open [w (writer outfile :append true)]
    (with-open [r (reader (bs/bs-path (first sps)))]
      (doseq [l (take 2 (line-seq r))]
        (.write w (str l "\n"))))
    (doseq [f (map bs/bs-path sps)]
      (with-open [r (reader f)]
        (dorun (map #(.write w (str % "\n"))
                    (filter #(not (= \# (first %)))
                            (line-seq r))))))))

(defn filter-signalp
  [bsl & {:keys [trim params outfile] :or {trim false params {} outfile nil}}]
  (let [of (sp-outfile outfile)
        sps (pmap #(signalp % (of) :params params)
                  (partition-all 10000 bsl))]
    (try (let [h (->> (pmap result-hash sps)
                      (apply merge))]
           (->> (map #(let [r (get h (bs/accession %))]
                        (if (= "Y" (first r))
                          (if trim
                            (bs/sub-bioseq % (second r))
                            %)))
                     bsl)
                (remove nil?)))
         (finally
           (if outfile
             (sp-combine-results outfile sps))
           (dorun (map #(fs/delete (bs/bs-path %)) sps))))))

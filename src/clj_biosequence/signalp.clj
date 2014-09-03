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
  (let [flist (atom ())
        sps (map #(signalp % (register-outfile flist) :params params)
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
           (doall (map #(fs/delete %) @flist))))))

;; private

(defn- make-signal-result
  [line]
  (let [[name cmax cpos ymax ypos smax spos smean D result Dmaxcut network]
        (string/split line #"\s+")]
    (->signalpProtein name (Float/parseFloat cmax) (Integer/parseInt cpos) (Float/parseFloat ymax)
                      (Integer/parseInt ypos) (Float/parseFloat smax) (Integer/parseInt spos)
                      (Float/parseFloat smean) (Float/parseFloat D) result (Float/parseFloat Dmaxcut)
                      network)))

(defn- signal-command
  [params infile]
  (let [up (flatten (vec (dissoc params "-f" "-g" "-k" "-m" "-n" "-T" "-w" "-l" "-v" "-V" "-h")))]
    (-> (cons "signalp" up)
        vec
        (conj (fs/absolute-path infile)))))


;; allowed
  ;; -s   Signal peptide networks to use ('best' or 'notm'). Default: 'best'
  ;; -t   Organism type> (euk, gram+, gram-). Default: 'euk'
  ;; -u   user defined D-cutoff for noTM networks
  ;; -U   user defined D-cutoff for TM networks
  ;; -M   Minimal predicted signal peptide length. Default: [10]
  ;; -c   truncate to sequence length - 0 means no truncation. Default '70'


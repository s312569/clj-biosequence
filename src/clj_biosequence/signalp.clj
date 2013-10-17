(ns clj-biosequence.signalp
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bios]
            [clojure.string :as string]
            [clojure.pprint :as pp]
            [clojure.data.xml :as xml]))

(declare signal-command make-signal-result)

;; signalp analysis

(defrecord signalpProtein [name cmax cpos ymax ypos smax spos smean D result Dmaxcut network]

  bios/Biosequence

  (accession [this]
    (:name this)))

;; results

(defrecord signalpReader [strm]

  bios/biosequenceReader

  (biosequence-seq [this]
    (map #(make-signal-result %)
         (drop 2 (line-seq strm))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord signalpResult [file]

  bios/biosequenceIO

  (bs-reader [this]
    (->signalpReader (io/reader (:file this)))))

(defn signalp
  [bs & {:keys [params outfile] :or {params {} outfile (fs/temp-file "sp")}}]
  (let [in (bios/fasta->file bs (fs/temp-file "sp-in")
                             :append false :func (fn [x] (if (bios/protein? x)
                                                          (bios/fasta-string x)
                                                          (throw (Throwable. "SignalP only analyses proteins.")))))]
    (with-open [out (io/output-stream outfile)]
      (let [sp @(exec/sh (signal-command params in) {:out out})]
        (if (= 0 (:exit sp))
          (->signalpResult (fs/absolute-path outfile))
          (if (:err sp)
            (throw (Throwable. (str "SignalP error: " (:err sp))))
            (throw (Throwable. (str "Exception: " (:exception sp))))))))))

(defn signalp?
  "Returns true if sequence contains signal sequence, false otherwise. Really slow."
  [prot]
  (let [spr (signalp (list prot))]
    (try
      (with-open [r (bios/bs-reader spr)]
        (= (:result (signalp (list prot))) "Y"))
      (finally (fs/delete (:file spr))))))

;; private

(defn- make-signal-result
  [line]
  (let [[name cmax cpos ymax ypos smax spos smean D result Dmaxcut network]
        (string/split line #"\s+")]
    (->signalpProtein name (read-string cmax) (read-string cpos) (read-string ymax) (read-string ypos) (read-string smax) (read-string spos) (read-string smean) (read-string D) result (read-string Dmaxcut) network)))

(defn- signal-command
  [params infile]
  (-> (into []
            (merge
             {"-f" "short" ;Setting the output format ('short', 'long', 'summary' or 'all')
              "-g" "Off" ;Graphics 'png' or 'png+eps'. Default: 'Off'
              "-k" "Off" ;Keep temporary directory. Default: 'Off'
              "-s" "best" ;Signal peptide networks to use ('best' or 'notm'). Default: 'best'
              "-t" "euk" ;Organism type> (euk, gram+, gram-). Default: 'euk'
              "-m" "Off" ;Make fasta file with mature sequence. Default: 'Off'
              "-c" "70" ; truncate to sequence length - 0 means no truncation. Default '70'
              "-l" "STDERR" ;Logfile if -v is defined. Default: 'STDERR'
              "-v" "Off"    ;Verbose. Default: 'Off'
              } params))
      flatten
      (conj "signalp")
      vec
      (conj (fs/absolute-path infile))))

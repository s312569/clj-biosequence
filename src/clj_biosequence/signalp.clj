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

; to do batch signalp

(declare signal-default-parameters make-signal-result)

(defrecord signalpResult [name cmax cpos ymax ypos smax spos smean D result Dmaxcut network])

(defn signalp
  "Returns a hash-map of the results of analysing a biosequence with signalP. The
   return hash-map keys are: :name, :cmax, :cpos, :ymax, :ypos, :smax, :spos, :smean, 
   :D, :result, :Dmaxcut and :network. All are standard outputs from signalP. Params 
   is a hash-map mapping signalP parameters to values. They are the same as the command
    version of the program. The only one that should need changing is the organism type:
   {'-t' 'euk' (euk, gram+, gram-; default euk)}"
  ([bioseq] (signalp bioseq {}))
  ([bioseq params]
     (if (bios/protein? bioseq)
      (let [in (fs/temp-file "signalp-")
            defs (flatten (vec (signal-default-parameters params)))]
        (spit in (bios/fasta-string bioseq))
        (let [args (conj (vec (cons "signalp" defs))
                         (fs/absolute-path in))
              sp @(exec/sh args)]
          (if (= 0 (:exit sp))
            (make-signal-result (string/split
                                 (first (rest (rest (string/split-lines (:out sp)))))
                                 #"\s+"))
            (if (:err sp)
              (throw (Throwable. (str "SignalP error: " (:err sp))))
              (throw (Throwable. (str "Exception: " (:exception sp))))))))
      (throw (Throwable. "Only protein sequences can be analysed for signal sequence.")))))

(defn signalp?
  "Returns true if sequence contains signal sequence, false otherwise."
  [prot]
  (= (:result (signalp prot)) "Y"))

;; private

(defn- make-signal-result
  [[name cmax cpos ymax ypos smax spos smean D result Dmaxcut network]]
  (->signalpResult (second (re-find #"^[^|]+\|([^|]+)\|" name)) (read-string cmax) (read-string cpos) (read-string ymax) (read-string ypos) (read-string smax) (read-string spos) (read-string smean) (read-string D) result (read-string Dmaxcut) network))

(defn- signal-default-parameters
  [params]
  (merge {"-f" "short" ;Setting the output format ('short', 'long', 'summary' or 'all')
          "-g" "Off"   ;Graphics 'png' or 'png+eps'. Default: 'Off'
          "-k" "Off"   ;Keep temporary directory. Default: 'Off'
          "-s" "best" ;Signal peptide networks to use ('best' or 'notm'). Default: 'best'
          "-t" "euk"  ;Organism type> (euk, gram+, gram-). Default: 'euk'
          "-m" "Off" ;Make fasta file with mature sequence. Default: 'Off'
          "-n" "Off" ;Make gff file of processed sequences. Default: 'Off'
          "-c" "70" ; truncate to sequence length - 0 means no truncation. Default '70'
          "-l" "STDERR"       ;Logfile if -v is defined. Default: 'STDERR'
          "-v" "Off"          ;Verbose. Default: 'Off'
          }
         params))

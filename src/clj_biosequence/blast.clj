(ns clj-biosequence.blast
  (:require [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bs]
            [clojure.pprint :as pp]
            [clojure.string :refer [split]]
            [clojure.data.xml :as xml]
            [clj-biosequence.alphabet :as ala])
  (:import [clj_biosequence.core fastaSequence]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blast hsp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn get-hsp-value
  "Takes a blastHsp object and returns the value corresponding to key.
   Keys are the keyword version of the XML nodes in the BLAST xml
   output.  All values are returned as strings. Typical BLAST HSP
   keys are:
    :Hsp_bit-score
    :Hsp_score
    :Hsp_evalue
    :Hsp_query-from
    :Hsp_query-to
    :Hsp_hit-from
    :Hsp_hit-to
    :Hsp_positive
    :Hsp_identity
    :Hsp_gaps
    :Hsp_hitgaps
    :Hsp_querygaps
    :Hsp_qseq
    :Hsp_hseq
    :Hsp_midline
    :Hsp_align-len
    :Hsp_query-frame
    :Hsp_hit-frame
    :Hsp_num
    :Hsp_pattern-from
    :Hsp_pattern-to
    :Hsp_density"
  [this key]
  (if (= key :Hsp_midline)
    (->> (:content (:src this))
         (filter #(= (:tag %) :Hsp_midline))
         first :content first)
    (zf/xml1-> (zip/xml-zip (:src this)) key zf/text)))

(defrecord blastHsp [src]
  bs/biosequenceTranslation
  (frame [this]
    (Integer/parseInt (get-hsp-value this :Hsp_query-frame))))

(defn- residue-counter [start end]
  (let [c (atom start)]
    {:increment (if (< end start)
                  (fn [x] (reset! c (- @c x)))
                  (fn [x] (reset! c (+ @c x))))
     :value (fn [] (deref c))}))

(defn- formatter
  [start end]
  (let [c (residue-counter start end)]
    (fn [x]
      (let [n (count (remove #{\-} x))]
        (vector (str ((:value c)))
                (apply str x)
                (str (- ((:increment c) n) 1)))))))

(defn- get-lines
  [hsp keys]
  (let [f (partial get-hsp-value hsp)
        args (map f keys)
        form (apply formatter (map #(Integer/parseInt %)
                                   (drop 1 args)))]
    (map form (partition-all 52 (first args)))))

(defn hsp-alignment
  "Returns the alignment from a HSP as a string."
  [hsp]
  (let [l (interleave
           (get-lines hsp '(:Hsp_qseq :Hsp_query-from :Hsp_query-to))
           (map #(vector "" (apply str %) "")
                (partition-all 52 (get-hsp-value hsp :Hsp_midline)))
           (get-lines hsp '(:Hsp_hseq :Hsp_hit-from :Hsp_hit-to)))
        m (apply max (mapcat #(list (count (first %))
                                    (count (nth % 2)))
                             l))
        b (fn [s] (apply str
                         (repeat (+ 2 (- m (count s))) \space)))
        pf (fn [x] (str (first x) (b (first x)) (second x) "  "
                        (last x) "\n"))]
    (apply str (interpose "\n" (map #(apply str %)
                                    (partition-all 3 (map pf l)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blast hit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord blastHit [src])

(defn get-hit-value
  "Takes a blastHit object and returns the value corresponding to key.
  Keys are the keyword version of the XML nodes in the BLAST xml
  output.  All values are returned as strings. Typical BLAST Hit
  values are:
   :Hit_id
   :Hit_len
   :Hit_accession
   :Hit_def
   :Hit_num"
  [this key]
  (zf/xml1-> (zip/xml-zip (:src this)) key zf/text))

(defn hit-accession
  "Returns the accession of the Blast hit."
  [hit]
  (get-hit-value hit :Hit_accession))

(defn hit-def
  "Returns the definition line of a Blast hit."
  [hit]
  (get-hit-value hit :Hit_def))

(defn hsp-seq
  "Takes a blastHit object and returns a lazy list of the blastHsp 
   objects contained in the hit."
  [this]
  (map #(->blastHsp (zip/node %))
       (zf/xml-> (zip/xml-zip (:src this))
                 :Hit_hsps
                 :Hsp)))

(defn hit-bit-scores
  "Takes a blastHit object and returns a list of floats corresponding
  to the bit scores of the HSPs composing the hit."
  [hit]
  (map #(Float/parseFloat (get-hsp-value % :Hsp_bit-score))
       (hsp-seq hit)))

(defn hit-e-values
  "Takes a blastHit object and returns a list of floats corresponding
  to the e-values of the HSPs composing the hit."
  [hit]
  (map #(Float/parseFloat (get-hsp-value % :Hsp_evalue))
       (hsp-seq hit)))

(defn hit-frames
  "Takes a blastHit object and returns a list of frames from each of
  the HSPs."
  [hit]
  (map #(Integer/parseInt (get-hsp-value % :Hsp_query-frame))
       (hsp-seq hit)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blast iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord blastIteration [src])

(extend blastIteration
  bs/biosequenceID
  (assoc bs/default-biosequence-id
    :accession (fn [this]
                 (-> (zf/xml1-> (zip/xml-zip (:src this))
                                :Iteration_query-def zf/text)
                     (split #"\s")
                     (first)))
    :accessions (fn [this]
                  (list (bs/accession this)))))

(defn iteration-query-length
  "Takes a blastIteration object and returns the query length."
  [iteration]
  (Integer/parseInt
   (zf/xml1-> (zip/xml-zip (:src iteration))
              :Iteration_query-len zf/text)))

(defn hit-seq
  "Returns a lazy list of blastHit objects from a blastIteration
  object."
  [this]
  (map #(->blastHit (zip/node %))
       (zf/xml-> (zip/xml-zip (:src this))
                 :Iteration_hits :Hit)))

(defn significant-hit-seq
  "Returns a list of blastHit objects from a blastIteration object
  that have a bit score equal to or greater than that specified (or
  default of 50). Measure argument accepts :bits or :evalue."
  [iteration score measure]
  (if (not (#{:bits :evalue} measure))
    (throw (Throwable. "Only :bits or :evalue allowable arguments for :measure keyword.")))
  (filter #(if (= measure :bits)
             (some (partial <= score) (hit-bit-scores %))
             (some (partial >= score) (hit-e-values %)))
          (hit-seq iteration)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord blastParameters [src])

(defn blast-parameter-value
  "Returns the value of a blast parameter from a blastParameters. Key
   denotes parameter keys used in the blast xml. All values returned
   as strings. Typical keys include:
   :Parameters_matrix
   :Parameters_expect
   :Parameters_include
   :Parameters_sc-match
   :Parameters_sc-mismatch
   :Parameters_gap-open
   :Parameters_gap-extend
   :Parameters_filter"
  [p key]
  (zf/xml1-> (zip/xml-zip (:src p))
             key
             zf/text))

(defn blast-evalue
  "Returns the evalue used from blast parameter record."
  [param]
  (Integer/parseInt (zf/xml1-> (zip/xml-zip (:src param))
                               :Parameters_expect zf/text)))

(defn blast-matrix
  "Returns the matrix used from a blast parameter record."
  [param]
  (zf/xml1-> (zip/xml-zip (:src param))
             :Parameters_matrix zf/text))

(defn blast-filter
  "Returns the filter used from a blast parameter record."
  [param]
  (zf/xml1-> (zip/xml-zip (:src param)) :Parameters_filter zf/text))

(defn blast-database
  "Returns the database used from a blast parameter record."
  [param]
  (:database (:src param)))

(defn blast-version
  "Returns the blast version from a blast parameter record."
  [param]
  (:version (:src param)))

(defn blast-program
  "Returns the program used from a blast parameter record."
  [param]
  (:program (:src param)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blast reader
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord blastReader [strm parameters]
  bs/biosequenceReader
  (biosequence-seq [this]
    (->> (:content (xml/parse (:strm this)))
         (some #(if (= :BlastOutput_iterations (:tag %)) %))
         :content
         (filter #(= :Iteration (:tag %)))
         (map #(->blastIteration %))))
  (get-biosequence [this acc]
    (first (filter #(= (bs/accession %) acc)
                   (bs/biosequence-seq this))))
  bs/biosequenceParameters
  (parameters [this] (:parameters this))
  java.io.Closeable
  (close [this] (.close ^java.io.BufferedReader (:strm this))))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blast search
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defrecord blastSearch [file opts])

(extend blastSearch
  bs/biosequenceIO
  {:bs-reader
   (fn [this]
     (let [p (with-open [r (apply bs/bioreader (bs/bs-path this)
                                  (:opts this))]
               (let [x (xml/parse r)
                     pa (->> (:content x)
                             (filter #(= :BlastOutput_param (:tag %)))
                             first
                             :content
                             first)]
                 (->blastParameters
                  (assoc pa
                    :database
                    (->> (:content x)
                         (filter #(= :BlastOutput_db (:tag %)))
                         first
                         :content
                         first)
                    :version
                    (->> (:content x)
                         (filter #(= :BlastOutput_version (:tag %)))
                         first
                         :content
                         first)
                    :program
                    (->> (:content x)
                         (filter #(= :BlastOutput_program (:tag %)))
                         first
                         :content
                         first)))))
           r (apply bs/bioreader (bs/bs-path this) (:opts this))]
       (->blastReader r p)))}
  bs/biosequenceFile
  bs/default-biosequence-file)

(defn init-blast-search
  "Initialises a blast search with a blast ouput file in xml format."
  [file & opts]
  {:pre [(fs/file? file)]}
  (->blastSearch (fs/absolute-path file) opts))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blast db
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- get-sequence-from-blast-db [db id]
  (let [s @(exec/sh (list "blastdbcmd" "-entry"
                          id "-db" (bs/bs-path db)))]
    (if (= 0 (:exit s))
      (:out s)
      (if (:err s)
        (throw (Throwable. (str "Blast error: " (:err s))))
        (throw (Throwable. (str "Exception: " (:exception s))))))))

(defrecord blastDB [file alphabet])

(extend blastDB
  bs/biosequenceFile
  bs/default-biosequence-file
  bs/biosequenceReader
  {:biosequence-seq
   (fn [_]
     (throw (Throwable. "Can't open a stream on a blast database.")))
   :get-biosequence
   (fn [this acc]
     (let [fs (-> (get-sequence-from-blast-db this acc)
                  (bs/init-fasta-string (:alphabet this)))]
       (with-open [r (bs/bs-reader fs)]
         (first (bs/biosequence-seq r)))))})

(defn init-blast-db
  "Initialises a blastDB object with the path and name of a blast
  database (omitting the indexed file extensions)."
  [file alphabet]
  {:pre [(fs/file? file)]}
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet.")))
  (->blastDB file alphabet))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; blasting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- blast-default-params
  [params in-file out-file db]
  (doall
   (remove #(nil? %)
           (flatten (seq (merge {"-evalue" "10"
                                 "-outfmt" "5"
                                 "-max_target_seqs" "3"
                                 "-query"
                                 in-file
                                 "-out"
                                 out-file
                                 "-db" db}
                                params))))))

(defn- run-blast 
  [prog db in out params]
  (let [defs (blast-default-params params
                                   in
                                   out
                                   (bs/bs-path db))]
    (let [bl @(exec/sh (cons prog defs))]
      (if (= 0 (:exit bl))
        (->blastSearch out nil)
        (if (:err bl)
          (throw (Throwable.
                  (str "Blast error: " (:err bl))))
          (throw (Throwable.
                  (str "Exception: " (:exception bl)))))))))

(defn blast
  [bs program db outfile & {:keys [params] :or {params {}}}]
  {:pre [(not (fs/file? outfile))]}
  (let [i (bs/biosequence->file bs (fs/temp-file "seq-")
                                :append false)]
    (try
      (run-blast program db
                 (fs/absolute-path i)
                 (fs/absolute-path outfile)
                 params)
      (finally (fs/delete i)))))


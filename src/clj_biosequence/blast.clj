(ns clj-biosequence.blast
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bs]
            [clj-biosequence.store :as st]
            [clojure.pprint :as pp]
            [clojure.string :refer [split]]
            [clojure.data.xml :as xml]
            [clj-biosequence.alphabet :as ala])
  (:import [clj_biosequence.core fastaSequence]))

(import '(java.io BufferedReader StringReader))

(declare blastp-defaults run-blast get-sequence-from-blast-db blast-default-params split-hsp-align iteration-query-id init-blast-collection get-hit-value)

;; blast hsp

(defrecord blastHsp [src])

(defn get-hsp-value
  "Takes a blastHsp object and returns the value corresponding to key.
     Keys are the keyword version of the XML nodes in the BLAST xml output. 
    All values are returned as strings. Typical BLAST HSP values are:
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
         first :content)
    (zf/xml1-> (zip/xml-zip (:src this)) key zf/text)))

(defn frame
  "Returns the frame of query match (if there is one)."
  [hsp]
  (Integer/parseInt (get-hit-value hsp :Hsp_query-frame)))

(defn print-alignment [hsp]
  (let [ss (fn [x] (map (partial apply str) (partition-all 52 x)))
        l (apply interleave (map ss (list (get-hsp-value hsp :Hsp_qseq)
                                          (get-hsp-value hsp :Hsp_midline)
                                          (get-hsp-value hsp :Hsp_hseq))))]
    (doseq [s (partition 3 l)]
      (doseq [f s]
        (println f))
      (println))))

;; blast hit

(defrecord blastHit [src])

(defn get-hit-value
  "Takes a blastHit object and returns the value corresponding to key. 
   Keys are the keyword version of the XML nodes in the BLAST xml output. 
   All values are returned as strings. Typical BLAST Hit values are:
   :Hit_id
   :Hit_len
   :Hit_accession
   :Hit_def
   :Hit_num"
  [this key]
  (zf/xml1-> (zip/xml-zip (:src this)) key zf/text))

(defn hit-accession
  [hit]
  (get-hit-value hit :Hit_accession))

(defn hit-def
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
  (let [v (map #(Float/parseFloat (get-hsp-value % :Hsp_bit-score)) (hsp-seq hit))]
    (if (not (empty? v)) v (list 0))))

(defn hit-e-value
  "Takes a blastHit object and returns a list of floats corresponding
  to the e-values of the HSPs composing the hit."
  [hit]
  (let [v (map #(Float/parseFloat (get-hsp-value % :Hsp_evalue)) (hsp-seq hit))]
    (if (not (empty? v)) v (list 1000000))))

(defn remove-hit-duplicates
  "Needs TESTING"
  [l]
  (map (fn [[k v]]
         (apply max-key #(first (hit-bit-scores %)) v))
       (seq (group-by #(get-hit-value % :Hit_accession) l))))

;; blast iteration

(defrecord blastIteration [src]
  
  bs/Biosequence
  
  (accession [this]
    (iteration-query-id this))

  st/mongoBSRecordIO

  (mongo-bs-save [this pname cname]
    (let [s (hash-map :acc (bs/accession this) :element "sequence"
                      :pname pname :cname cname
                      :type "biosequence/blast-iteration"
                      :src (bs/bs-freeze this))]
      (if (:_id this)
        (assoc s :_id (:_id this))
        s))))

(defn init-blast-iteration
  [src]
  (->blastIteration src))

(defn iteration-query-id
  "Takes a blastIteration object and returns the query ID."
  [this]
  (-> (zf/xml1-> (zip/xml-zip (:src this)) :Iteration_query-def zf/text)
      (split #"\s")
      (first)))

(defn hit-seq
  "Returns a (lazy) list of blastHit objects from a blastIteration object."
  [this]
  (map #(->blastHit (zip/node %))
       (zf/xml-> (zip/xml-zip (:src this)) :Iteration_hits :Hit)))

(defn top-hit
  "Returns the highest scoring blastHit object from a blastIteration object."
  [this]
  (or (first (hit-seq this)) (->blastHit nil)))

(defn top-hsp
  "Returns the highest scoring hsp from the highest scoring hit in a blast iteration."
  [it]
  (or (->> it hit-seq first hsp-seq first)
      (->blastHsp nil)))

;; parameters

(defrecord blastParameters [src]

  st/mongoBSRecordIO

  (mongo-bs-save [this pname cname]
    (let [s (hash-map :acc (bs/accession this) :element "parameter"
                      :pname pname :cname cname
                      :type "biosequence/blast-parameters"
                      :src (bs/bs-freeze this))]
      (if (:_id this)
        (assoc s :_id (:_id this))
        s))))

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

(defn blast-database
  [param]
  (:database (:src param)))

(defn blast-version
  [param]
  (:version (:src param)))

(defn blast-program
  [param]
  (:program (:src param)))

(defn- init-blast-params
  [src]
  (->blastParameters src))

;; blastSearch

(defprotocol blastSearchAccess
  (result-by-accession [this accession]
    "Returns the blast search for the specified protein."))

(defrecord blastReader [strm parameters]

  bs/biosequenceReader

  (biosequence-seq [this]
    (->> (:content (xml/parse (:strm this)))
         (some #(if (= :BlastOutput_iterations (:tag %)) %))
         :content
         (filter #(= :Iteration (:tag %)))
         (map init-blast-iteration)))
  
  (parameters [this]
    (:parameters this))
  
  java.io.Closeable
  
  (close [this]
    (.close ^java.io.BufferedReader (:strm this)))
  
  blastSearchAccess
  
  (result-by-accession [this accession]
    (some #(if (= accession (iteration-query-id %))
             %)
          (bs/biosequence-seq this))))

(defn init-blast-reader
  [s p]
  (->blastReader s p))

(defrecord blastSearch [file]

  bs/biosequenceIO

  (bs-reader [this]
    (let [p (with-open [r (io/reader (:file this))]
              (let [p (->> (:content (xml/parse r))
                           (filter #(= :BlastOutput_param (:tag %)))
                           first
                           :content
                           first)]
                (list
                 (init-blast-params (assoc p
                                      :database
                                      (->> (:content (:xml this))
                                           (filter #(= :BlastOutput_db (:tag %)))
                                           first
                                           :content
                                           first)
                                      :version
                                      (->> (:content (:xml this))
                                           (filter #(= :BlastOutput_version (:tag %)))
                                           first
                                           :content
                                           first)
                                      :program
                                      (->> (:content (:xml this))
                                           (filter #(= :BlastOutput_program (:tag %)))
                                           first
                                           :content
                                           first))))))
          r (io/reader (:file this))]
      (init-blast-reader r p)))

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:file this))))

(defn init-blast-search
  [file]
  (->blastSearch (fs/absolute-path file)))


(defn delete-blast-search
  [search]
  (fs/delete (:file search)))

;; blast db

(defrecord blastDB [path alphabet]

  bs/biosequenceFile

  (bs-path [this]
    (fs/absolute-path (:path this))))

(defn blast-get-sequence
  "Returns the specified sequence from a blastDB object as a fastaSequence object."
  [db id]
  (if id
    (let [fs (-> (get-sequence-from-blast-db db id)
                 (bs/init-fasta-string (:alphabet db)))]
      (with-open [r (bs/bs-reader fs)]
        (first (bs/biosequence-seq r))))))

(defn init-blast-db
  "Initialises a blastDB object."
  [path alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet."))
    (if (and (not (nil? path)) (fs/file? path))
      (->blastDB path alphabet)
      (throw (Throwable. (str "File not found: " path))))))

;; blasting

(defn blast
  ([bs program db outfile] (blast bs program db outfile {}))
  ([bs program db outfile params]
     (let [i (bs/biosequence->file bs (fs/temp-file "seq-") :append false)]
       (try
         (run-blast program db
                    (fs/absolute-path i)
                    (fs/absolute-path outfile)
                    params)
         (finally (fs/delete i))))))

;; helpers

(defn- get-sequence-from-blast-db [db id]
  (let [s @(exec/sh (list "blastdbcmd" "-entry" id "-db" (:path db)))]
    (if (= 0 (:exit s))
      (:out s)
      (if (:err s)
        (throw (Throwable. (str "Blast error: " (:err s))))
        (throw (Throwable. (str "Exception: " (:exception s))))))))

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
  "Need timeout"
  (let [defs (blast-default-params params
                                   in
                                   out
                                   (bs/bs-path db))]
    (let [bl @(exec/sh (cons prog defs))]
      (if (= 0 (:exit bl))
        (->blastSearch out)
        (if (:err bl)
          (throw (Throwable. (str "Blast error: " (:err bl))))
          (throw (Throwable. (str "Exception: " (:exception bl)))))))))

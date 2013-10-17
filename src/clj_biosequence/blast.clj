(ns clj-biosequence.blast
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bios]
            [clojure.pprint :as pp]
            [clojure.string :as st]
            [clojure.data.xml :as xml]
            [clj-biosequence.alphabet :as ala])
  (:import [clj_biosequence.core fastaSequence]))

(import '(java.io BufferedReader StringReader))

(declare blastp-defaults run-blast get-sequence-from-blast-db blast-default-params split-hsp-align iteration-query-id)

;; blast hsp

(defrecord blastHSP [src])

(defmethod print-method clj_biosequence.blast.blastHSP
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

(defn get-hsp-value
  "Takes a blastHSP object and returns the value corresponding to key.
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
    (first (:content
            (first
             (filter #(= (:tag %) :Hsp_midline)
                     (:content (:src this))))))
    (zf/xml1-> (zip/xml-zip (:src this))
               key
               zf/text)))

(defn alignment-string
  "Takes a blastHSP object and returns a string of its alignment."
  [this]
  (pp/cl-format nil "~{~A~%~A~%~A~%~%~}"
                (interleave
                 (split-hsp-align
                  (get-hsp-value this :Hsp_qseq)
                  (read-string (get-hsp-value this :Hsp_query-from)))
                 (map #(pp/cl-format nil "~7@< ~>~A~7< ~>"
                                     (apply str %))
                      (partition-all 58 (get-hsp-value this :Hsp_midline)))
                 (split-hsp-align
                  (get-hsp-value this :Hsp_hseq)
                  (read-string (get-hsp-value this :Hsp_hit-from))))))

;; blast hit

(defrecord blastHit [src])

(defmethod print-method clj_biosequence.blast.blastHit
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

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
  (zf/xml1-> (zip/xml-zip (:src this))
             key
             zf/text))

(defn hsp-seq
  "Takes a blastHit object and returns a lazy list of the blastHSP 
   objects contained in the hit."
  [this]
  (map #(->blastHSP (zip/node %))
       (zf/xml-> (zip/xml-zip (:src this))
                 :Hit_hsps
                 :Hsp)))

(defn top-hsp
  "Takes a blastHit object and returns the top scoring blastHSP 
   object in the hit. Returns an empty blastHSP if blastHit was empty."
  [this]
  (or (first (hsp-seq this))
      (->blastHSP nil)))

(defn hit-bit-score
  "Takes a blastHit object and returns the bit score of the top scoring HSP
   in the hit. Returns 0 if the blastHit was empty."
  [this]
  (if (:src this)
    (Float/parseFloat (get-hsp-value (top-hsp this) :Hsp_bit-score))
    0))

(defn hit-e-value
  "Takes a blastHit object and returns the bit score of the top scoring HSP
   in the hit. Returns 0 if the blastHit was empty."
  [this]
  (if (:src this)
    (Float/parseFloat (get-hsp-value (top-hsp this) :Hsp_evalue))))

(defn hit-string
  "A convenience function that takes a blastHit object and returns a formatted
     summary of the top scoring HSP in the hit. Includes accession, bit score
     (of the top scoring hit), definition of the hit and the alignment. Returns 
    'No hits in search' if empty blastHit. Example output:

   Accession: sp|Q8HY10|CLC4M_NOMCO
   Bit score: 50.0617822382917
   Def: C-type lectin domain family 4 member M OS=Nomascus concolor GN=CLEC
   Alignment:
   
   72     QAQQRDIEKEIESQKTSLTESWKKIIAEDIENRTNR-----SELKMEGQLSDLQEALT    124
          +++Q++I +E+   K ++ E  +K   ++I     R      EL  + +   + + LT       
   198    KSKQQEIYQELTRLKAAVGELPEKSKQQEIYQELTRLKAAVGELPDQSKQQQIYQELT    255"
  [this]
  (if (:src this)
    (let [hsp (top-hsp this)]
      (pp/cl-format nil "Accession: ~A~%Bit score: ~A~%Def: ~A~%Alignment:~%~%~A"
                    (get-hit-value this :Hit_id)
                    (get-hsp-value hsp :Hsp_bit-score)
                    (first (map #(apply str %)
                                (partition-all 67 (get-hit-value this :Hit_def))))
                    (alignment-string hsp)))
    "No hits in search.\n"))

;; blast iteration

(defrecord blastIteration [src]

  bios/Biosequence

  (accession [this]
    (iteration-query-id this))

  (bs-save [this]
    (assoc this :src (pr-str (dissoc this :_id))
           :acc (bios/accession this))))

(defmethod print-method clj_biosequence.blast.blastIteration
  [this ^java.io.Writer w]
  (bios/print-tagged this w))

(defn iteration-query-id
  "Takes a blastIteration object and returns the query ID."
  [this]
  (-> (zf/xml1-> (zip/xml-zip (:src this))
                 :Iteration_query-def
                 zf/text)
      (st/split #"\s")
      (first)))

(defn hit-seq
  "Returns a (lazy) list of blastHit objects from a blastIteration object."
  [this]
  (map #(->blastHit (zip/node %))
       (zf/xml-> (zip/xml-zip (:src this))
                 :Iteration_hits
                 :Hit)))

(defn top-hit
  "Returns the highest scoring blastHit object from a blastIteration object."
  [this]
  (or (first (hit-seq this))
      (->blastHit nil)))

;; blastSearch

(defrecord blastReader [strm]

  bios/biosequenceReader

  (biosequence-seq [this]
    (->> (:content (xml/parse (:strm this)))
         (filter #(= :BlastOutput_iterations (:tag %)))
         first
         :content
         (filter #(= :Iteration (:tag %)))
         (map #(->blastIteration %))))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this))))

(defrecord blastSearch [file]

  bios/biosequenceIO

  (bs-reader [this]
    (->blastReader (io/reader (:file this)))))

(defn init-blast-search
  [file]
  (->blastSearch (fs/absolute-path file)))

(defn get-iteration-by-id
  "Returns the blastIteration object for the specified biosequence from
     a blastSearch object."
  [this accession]
  (with-open [^java.io.BufferedReader r (bios/bs-reader (:src this))]
    (some #(if (= accession (iteration-query-id %))
             %)
          (bios/biosequence-seq r))))

(defn get-parameter-value
  "Returns the value of a blast parameter from a blastSearch object. Key 
     denotes parameter keys used in the blast xml. All values returned as
     strings. Typical keys include:
   :Parameters_matrix
   :Parameters_expect
   :Parameters_include
   :Parameters_sc-match
   :Parameters_sc-mismatch
   :Parameters_gap-open
   :Parameters_gap-extend
   :Parameters_filter"
  [this key]
  (with-open [rdr (io/reader (:src this))]
    (zf/xml1-> (zip/xml-zip (xml/parse rdr))
               :BlastOutput_param
               :Parameters
               key
               zf/text)))

(defn database-searched
  "Returns the path (as a string) of the database used in a blast search
     from a blastSearch object."
  [this]
  (with-open [rdr (io/reader (:src this))]
    (zf/xml1-> (zip/xml-zip (xml/parse rdr))
               :BlastOutput_db
               zf/text)))

(defn program-used
  "Returns the program used in a blast search from a blastSearch object."
  [this]
  (with-open [rdr (io/reader (:src this))]
    (zf/xml1-> (zip/xml-zip (xml/parse rdr))
               :BlastOutput_program
               zf/text)))

;; blast db

(defrecord blastDB [path alphabet])

(defn get-sequence
  "Returns the specified sequence from a blastDB object as a fastaSequence object."
  [db id]
  (if id
    (let [fs (-> (get-sequence-from-blast-db db id)
                 (bios/init-fasta-string (:alphabet db)))]
      (with-open [r (bios/bs-reader fs)]
        (first (bios/biosequence-seq r))))))

(defn init-blast-db
  "Initialises a blastDB object."
  [path alphabet]
  (if-not (ala/alphabet? alphabet)
    (throw (Throwable. "Unrecognised alphabet."))
    (if (fs/file? path)
      (->blastDB path alphabet)
      (throw (Throwable. (str "File not found: " path))))))

;; blasting

(defn blast
  [bs program db & {:keys [outfile params] :or {outfile (fs/temp-file "blast")
                                                params {}}}]
  (let [i (bios/fasta->file bs (fs/temp-file "seq-") :append false)]
    (try
      (run-blast program db
                 (fs/absolute-path i)
                 (fs/absolute-path outfile)
                 params)
      (finally (fs/delete i)))))

;; helpers

(defn- split-hsp-align
  "Helper function that takes a list of blast query or hit alignment strings
   and produces a list of strings formatted for the 'alignment-string' function.
   Basically it calculates the sequence numbers and sticks them on the start and end."
  [string begin]
  (loop [s (partition-all 58 string)
         b begin
         st []]
    (if (empty? s)
      st
      (let [l (count (remove #(= % \-) (first s)))]
        (recur (rest s)
               (+ b l)
               (conj st
                     (pp/cl-format nil "~7@<~A~>~A~7<~A~>" 
                                   b (apply str (first s)) (+ b (- l 1)))))))))

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
                                 "-max_target_seqs" "1"
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
                                   (:path db))]
    (let [bl @(exec/sh (cons prog defs))]
      (if (= 0 (:exit bl))
        (->blastSearch out)
        (if (:err bl)
          (throw (Throwable. (str "Blast error: " (:err bl))))
          (throw (Throwable. (str "Exception: " (:exception bl)))))))))

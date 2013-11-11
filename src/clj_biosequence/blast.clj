(ns clj-biosequence.blast
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-commons-exec :as exec]
            [clojure.data.zip.xml :as zf]
            [clojure.zip :as zip]
            [clj-biosequence.core :as bios]
            [clj-biosequence.store :as sto]
            [clojure.pprint :as pp]
            [clojure.string :refer [split]]
            [clojure.data.xml :as xml]
            [clj-biosequence.alphabet :as ala])
  (:import [clj_biosequence.core fastaSequence]))

(import '(java.io BufferedReader StringReader))

(declare blastp-defaults run-blast get-sequence-from-blast-db blast-default-params split-hsp-align iteration-query-id init-blast-collection)

;; blast hsp

(defrecord blastHSP [src])

(defmethod print-method clj_biosequence.blast.blastHSP
  [this ^java.io.Writer w]
  (bios/print-biosequence this w))

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
  (zf/xml1-> (zip/xml-zip (:src this)) key zf/text))

;; blast hit

(defrecord blastHit [src])

(defmethod print-method clj_biosequence.blast.blastHit
  [this ^java.io.Writer w]
  (bios/print-biosequence this w))

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

(defn hit-bit-scores
  "Takes a blastHit object and returns a list of floats corresponding
  to the bit scores of the HSPs composing the hit."
  [hit]
  (map #(Float/parseFloat (get-hsp-value % :Hsp_bit-score)) (hsp-seq hit)))

(defn hit-e-value
  "Takes a blastHit object and returns a list of floats corresponding
  to the e-values of the HSPs composing the hit."
  [hit]
  (map #(Float/parseFloat (get-hsp-value % :Hsp_evalue)) (hsp-seq hit)))

;; blast iteration

(defrecord blastIteration [src]

  bios/Biosequence

  (accession [this]
    (iteration-query-id this))

  (save-rep [this]
    (hash-map :acc (bios/accession this)
              :src (pr-str this)
              :element) "sequence"))

(defmethod print-method clj_biosequence.blast.blastIteration
  [this ^java.io.Writer w]
  (bios/print-biosequence this w))

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

;; parameters

(defrecord blastParameters [src])

(defn parameter-value
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

(defn- init-blast-params
  [src]
  (->blastParameters src))

(defmethod print-method clj_biosequence.blast.blastParameters
  [this ^java.io.Writer w]
  (bios/print-biosequence this w))

;; blastSearch

(defprotocol blastSearchAccess
  (result-by-accession [this accession]
    "Returns the blast search for the specified protein.")
  (parameters [this]
    "returns a blastParameters from the search.")
  (database [this]
    "Returns the path (as a string) of the database used in a blast
     search from a blastSearch object.")
  (version [this]
    "returns the version of the blast used in the search.")
  (program [this]
    "Returns the program used in a blast search from a blastSearch object."))

(defrecord blastReader [strm xml]

  bios/biosequenceReader

  (biosequence-seq [this]
    (map #(->blastIteration (zip/node %))
          (zf/xml-> (zip/xml-zip (:xml this))
                    :BlastOutput_iterations
                    :Iteration)))

  java.io.Closeable

  (close [this]
    (.close ^java.io.BufferedReader (:strm this)))

  blastSearchAccess

  (result-by-accession [this accession]
    (some #(if (= accession (iteration-query-id %))
             %)
          (bios/biosequence-seq this)))

  (parameters [this]
    (->> (:content (:xml this))
         (filter #(= :BlastOutput_param (:tag %)))
         first
         :content
         first
         init-blast-params))

  (database [this]
    (->> (:content (:xml this))
         (filter #(= :BlastOutput_db (:tag %)))
         first
         :content
         first))

  (version [this]
    (->> (:content (:xml this))
         (filter #(= :BlastOutput_version (:tag %)))
         first
         :content
         first))

  (program [this]
    (->> (:content (:xml this))
         (filter #(= :BlastOutput_program (:tag %)))
         first
         :content
         first)))

(defrecord blastSearch [file]

  bios/biosequenceIO

  (bs-reader [this]
    (let [s (io/reader (:file this))]
      (->blastReader s (xml/parse s))))

  sto/storeCollectionIO

  (mongo-save-file [this project name]
    (let [i (init-blast-collection name (:name project) "biosequence/blast")]
      (with-open [r (bios/bs-reader this)]
        (sto/save-list (concat
                        (list (hash-map :acc (str (java.util.UUID/randomUUID))
                                        :src (pr-str (database r))
                                        :element "db")
                              (hash-map :acc (str (java.util.UUID/randomUUID))
                                        :src (pr-str (program r))
                                        :element "program")
                              (hash-map :acc (str (java.util.UUID/randomUUID))
                                        :src (pr-str (version r))
                                        :element "version")
                              (hash-map :acc (str (java.util.UUID/randomUUID))
                                        :src (pr-str (parameters r))
                                        :element "parameters"))
                        (map bios/save-rep (bios/biosequence-seq r)))
                       i)))))

(defn init-blast-search
  [file]
  (->blastSearch (fs/absolute-path file)))

;; store

(defrecord blastResultCollection [name pname type]

  blastSearchAccess

  (result-by-accession [this accession]
    (first (sto/get-record this :acc accession :element "sequence")))

  (parameters [this]
    (first (sto/get-record this :element "parameters")))

  (database [this]
    (first (sto/get-record this :element "db")))

  (version [this]
    (first (sto/get-record this :element "version")))

  (program [this]
    (first (sto/get-record this :element "program")))

  sto/storeCollectionAccess

  (collection-seq [this]
    (sto/get-record this :element "sequence")))

(defmethod print-method clj_biosequence.blast.blastResultCollection
  [this ^java.io.Writer w]
  (bios/print-biosequence this w))

(defn init-blast-collection
  [name pname type]
  (->blastResultCollection name pname type))

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

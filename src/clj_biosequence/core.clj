(ns clj-biosequence.core
  (:require [clojure.java.io :refer [writer reader input-stream]]
            [fs.core :refer [file? absolute-path extension exists?]]
            [clj-http.client :as client]
            [clojure.string :refer [trim split upper-case]]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.indexing :as ind]
            [iota :as iot]
            [clj-time.format :refer [formatter parse unparse]]
            [clojure.core.reducers :as r])
  (:import
   (org.apache.commons.compress.compressors.bzip2
    BZip2CompressorInputStream)))

(declare init-fasta-sequence)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro get-text
  "Low level macro for retrieving data from xml elements."
  [obj & keys]
  `(zf/xml1-> (clojure.zip/xml-zip (:src ~obj))
              ~@keys clojure.data.zip.xml/text))

(defmacro get-list
  "Low level macro for retrieving data from xml elements."
  [obj & keys]
  `(zf/xml-> (clojure.zip/xml-zip (:src ~obj))
            ~@keys))

(defmacro get-one
  "Low level macro for retrieving data from xml elements."
  [obj & keys]
  `(zf/xml1-> (clojure.zip/xml-zip (:src ~obj))
              ~@keys))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; date
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn make-date
  [str f]
  (parse f str))

(defn make-date-format
  [str]
  (formatter str))

(defn parse-date
  [d]
  (if d
    (unparse (formatter "yyyy-MM-dd") d)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; protocols 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(load "protocols")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; network
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def bioseq-proxy (atom {}))

(defn set-bioseq-proxy!
  [params]
  (reset! bioseq-proxy params))

(defn get-req
  ([a] (get-req a {}))
  ([a param]
   (client/get a (merge param @bioseq-proxy))))

(defn post-req
  [a param]
  (client/post a (merge param @bioseq-proxy)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn bioseq->string
  "Returns the sequence of a biosequence as a string with newlines
  every 80 chars."
  [bs]
  (->> (:sequence bs)
       (partition-all 80)
       (map #(apply str %))
       (interpose "\n")
       (apply str)))

(defn fasta-string
  [bioseq]
  "Returns the biosequence as a fasta formatted string."
  (str ">" (accession bioseq) " "
       (description bioseq) "\n"
       (bioseq->string bioseq) "\n"))

(defn sub-bioseq
  "Returns a new fasta sequence object with the sequence corresponding
  to 'beg' (inclusive) and 'end' (exclusive) of 'bs'. If no 'end'
  argument returns from 'start' to the end of the sequence. Unlike
  vectors assumes NOT a zero based index to correlate with common
  intervals in annotation files."
  ([bs beg] (sub-bioseq bs beg nil))
  ([bs beg end]
     (init-fasta-sequence (accession bs)
                          (str (description bs)
                               " [" beg " - "
                               (if end end "End") "]")
                          (alphabet bs)
                          (->> (if end
                                 (subvec (bs-seq bs) (- beg 1) end)
                                 (subvec (bs-seq bs) (- beg 1)))
                               (apply str)))))

(defn interval-complete?
  "Returns the interval if it has a start and end value or a point
  value, nil otherwise."
  [interval]
  (and (or (and (start interval) (end interval))
           (point interval))
       interval))

(defn get-interval-sequence
  "Returns a fasta sequence corresponding to the provided interval."
  [interval bs]
  (cond
    ;; circular DNA spanning origin not comp
    (and (> (start interval) (end interval))
         (not (comp? interval)))
    (let [o (sub-bioseq bs (start interval))
          t (sub-bioseq bs 0 (end interval))]
      (assoc o :sequence (apply str (concat (bs-seq o) (bs-seq t)))
             :description (str (second (re-find #"^(.+)\s[^\[]+\]$"))
                               " [" (start interval) " - " (end interval) "]")))
    ;; circular DNA spanning origin comp direction
    (and (comp? interval)
         (> (end interval) (start interval)))
    (let [o (sub-bioseq bs (end interval))
          t (sub-bioseq bs 0 (start interval))]
      (assoc o :sequence (apply str (concat (bs-seq o) (bs-seq t)))
             :description (str (second (re-find #"^(.+)\s[^\[]+\]$"))
                               " [" (end interval) " - " (start interval) "]")))
    ;; otherwise as usual
    :else
    (let [s (if (comp? interval)
              (end interval)
              (start interval))
          e (if (comp? interval)
              (start interval)
              (end interval))]
      (sub-bioseq bs s e))))

;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;

(defn biosequence->file
  "Takes a collection of biosequences and prints them to file. To
  append to an existing file use `:append true` and the `:func`
  argument can be used to pass a function that will be used to prepare
  the printed output, the default is `fasta-string` which will print
  the biosequences to the file in fasta format. Returns the file."
  [bs file & {:keys [append func] :or {append true func fasta-string}}]
  (with-open [w (writer file :append append)]
    (dorun (map #(let [n (func %)]
                   (if n
                     (.write w (str n "\n"))))
                bs)))
  file)

(defn- clean-sequence
  [s ala]
  (let [w #{\space \newline}]
    (loop [l s a []]
      (if-not (seq l)
        (apply str a)
        (let [c (first l)]
          (if (not (w c))
            (recur (rest l) (conj a (if (ala/allowed-character? ala c) c \X)))
            (recur (rest l) a)))))))

(defn clean-sequences
  "Removes spaces and newlines and checks that all characters are
  legal characters for the supplied alphabet. Replaces non-valid
  characters with \\X. Only performs cleaning if supplied alphabet is
  a checked alphabet."
  [alphabet coll]
  (let [a (ala/get-alphabet alphabet)]
    (if (ala/checked? a)
      (map #(assoc % :sequence
                   (clean-sequence (:sequence %) a)) coll)
      coll)))

(defn object->file
  "Spits an object to file after making sure *print-length* is
  temporarily set to false."
  [obj file]
  (binding [*print-length* false]
    (spit file (pr-str obj))))

;; helper files

(load "dna")
(load "protein")
(load "bioreader")
(load "serialising")
(load "mapping")
(load "fasta")

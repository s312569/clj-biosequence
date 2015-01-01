(ns clj-biosequence.core
  (:require [clojure.java.io :refer [writer reader input-stream]]
            [fs.core :refer [file? absolute-path extension]]
            [clj-http.client :as client]
            [clojure.string :refer [trim split upper-case]]
            [clj-biosequence.alphabet :as ala]
            [clj-biosequence.indexing :as ind]
            [clojure.edn :as edn]
            [clojure.data.xml :as xml]
            [miner.tagged :as tag]
            [iota :as iot]
            [clj-time.format :refer [formatter parse]]
            [clojure.core.reducers :as r])
  (:import
   (org.apache.commons.compress.compressors.bzip2
    BZip2CompressorInputStream)))

(declare init-fasta-store init-fasta-sequence translate init-index-file)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; protocols 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(load "protocols")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro get-text
  [obj & keys]
  `(zf/xml1-> (zip/xml-zip (:src ~obj))
              ~@keys zf/text))

(defmacro get-list
  [obj & keys]
  `(zf/xml-> (zip/xml-zip (:src ~obj))
            ~@keys))

(defmacro get-one
  [obj & keys]
  `(zf/xml1-> (zip/xml-zip (:src ~obj))
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

;; network

(def bioseq-proxy (atom {}))

(defn set-bioseq-proxy!
  [params]
  (reset! bioseq-proxy params))

(defn get-req
  [a param]
  (client/get a (merge param @bioseq-proxy)))

(defn post-req
  [a param]
  (client/post a (merge param @bioseq-proxy)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- bioseq->string
  "Returns the sequence of a biosequence as a string with newlines
  every 80 chars."
  [bs]
  (apply str (interpose "\n"
                        (map #(apply str %)
                             (partition-all 80 (bs-seq bs))))))

(defn fasta-string
  [bioseq]
  "Returns the biosequence as a fasta formatted string."
  (str ">" (accession bioseq) " "
       (description bioseq) "\n" (bioseq->string bioseq) "\n"))

(defn sub-bioseq
  "Returns a new fasta sequence object with the sequence corresponding
   to 'beg' (inclusive) and 'end' (exclusive) of 'bs'. If no 'end'
   argument returns from 'start' to the end of the sequence. Zero
   based index."
  ([bs beg] (sub-bioseq bs beg nil))
  ([bs beg end]
     (init-fasta-sequence (accession bs)
                          (str (description bs)
                               " [" beg " - "
                               (if end end "End") "]")
                          (alphabet bs)
                          (if end
                            (apply str (subvec (bs-seq bs) beg end))
                            (apply str (subvec (bs-seq bs) beg))))))

(defn get-interval-sequence
  "Returns a fasta sequence corresponding to the provided interval."
  [interval bs]
  (cond
   ;; circular DNA spanning origin not comp
   (and (> (start interval) (end interval))
        (not (comp? interval))) 
   (let [o (sub-bioseq bs (- (start interval) 1))
         t (sub-bioseq bs 0 (end interval))]
     (assoc o :sequence (vec (concat (bs-seq o) (bs-seq t)))
            :description (str (second (re-find #"^(.+)\s[^\[]+\]$"))
                              " [" (start interval) " - " (end interval) "]")))
   ;; circular DNA spanning origin comp direction
   (and (comp? interval)
        (> (end interval) (start interval)))
   (let [o (sub-bioseq bs (- (end interval) 1))
         t (sub-bioseq bs 0 (start interval))]
     (assoc o :sequence (vec (concat (bs-seq o) (bs-seq t)))
            :description (str (second (re-find #"^(.+)\s[^\[]+\]$"))
                              " [" (end interval) " - " (start interval) "]")))
   ;; otherwise as usual
   :else
   (let [s (if (comp? interval)
             (- (end interval) 1)
             (- (start interval) 1))
         e (if (comp? interval)
             (start interval)
             (end interval))]
     (sub-bioseq bs s e))))

;; (defn get-feature-sequence
;;   "Returns a fastaSequence object containing the sequence specified in
;;   a feature object from a biosequence."
;;   [feature bs]
;;   (let [intervals (intervals feature)]
;;     (init-fasta-sequence
;;      (accession bs)
;;      (str (description bs) " - Feature: " (feature-type feature)
;;           " - [" (start (first intervals)) "-" (end (last intervals)) "]")
;;      (alphabet bs)
;;      (vec (mapcat #(if (comp? %)
;;                      (apply str (subvec (ala/revcom (bs-seq bs))
;;                                         (- (end %) 1)
;;                                         (start %)))
;;                      (apply str (subvec (bs-seq bs)
;;                                         (- (start %) 1)
;;                                         (end %)))) intervals)))))

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

(defn clean-sequence
  "Removes spaces and newlines and checks that all characters are
   legal characters for the supplied alphabet. Replaces non-valid
   characters with \\X. If `a' is not a defined alphabet throws an
   exception."
  [s a]
  (let [k (ala/get-alphabet a)
        w #{\space \newline}
        t (if (string? s) (upper-case s) s)]
    (loop [l t a []]
      (if-not (seq l)
        a
        (let [c (first l)]
          (if (not (w c))
            (recur (rest l) (conj a (if (k c) c \X)))
            (recur (rest l) a)))))))

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

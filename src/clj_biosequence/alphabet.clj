(ns clj-biosequence.alphabet
  (:require [clojure.edn :as edn]
            [clojure.java.io :refer [resource]]))

;; alphabets

(defprotocol biosequenceAlphabet
  (alphabet-is-protein? [this]
    "Returns true if alphabet is a protein alphabet.")
  (allowed-character? [this char]
    "True if character argument is permitted in alphabet.")
  (checked? [this]
    "Specifies if sequence alphabet restrictions should be enforced.")
  (data [this char]
    "Returns meta data for corresponding characters."))

(def ^{:doc "Default alphabet parameters."}
  default-biosequence-alphabet
  {:allowed-character? (fn [this char]
                         (contains? (.alphabet this) char))
   :checked? (fn [_] true)
   :data (fn [this char] (get (.alphabet this) char))})

(deftype iupacAminoAcids [alphabet])
(deftype iupacNucleicAcids [alphabet])
(deftype signalpAminoAcids [alphabet])
(deftype uncheckedDNA [alphabet])

(extend iupacAminoAcids
  biosequenceAlphabet
  (assoc default-biosequence-alphabet
    :alphabet-is-protein? (fn [this] true)))

(extend iupacNucleicAcids
  biosequenceAlphabet
  (assoc default-biosequence-alphabet
         :alphabet-is-protein? (fn [_] false)))

(extend signalpAminoAcids
  biosequenceAlphabet
  (assoc default-biosequence-alphabet
         :alphabet-is-protein? (fn [_] true)))

(extend uncheckedDNA
  biosequenceAlphabet
  (assoc default-biosequence-alphabet
         :alphabet-is-protein? (fn [_] false)
         :checked? (fn [_] false)))

(def alphabets {:iupacAminoAcids
                (iupacAminoAcids.
                 (edn/read-string (slurp (resource "aa-iupacAA.clj"))))
                :iupacNucleicAcids
                (iupacNucleicAcids.
                 (edn/read-string (slurp (resource "dna-iupacdna.clj"))))
                :signalpAminoAcids
                (signalpAminoAcids.
                 (edn/read-string (slurp (resource "aa-signalp.clj"))))
                :uncheckedDNA
                (uncheckedDNA. {})})

;; codon tables

(def codon-tables
  (edn/read-string (slurp (resource "codon-table.clj"))))

;; functions

(defn get-alphabet
  "If a keyword corresponding to a defined alphabet or, if an object
  that satisifies biosequenceAlphabet protocol, the object. Throws an
  exception otherwise."
  [k]
  (or (alphabets k)
      (if (satisfies? biosequenceAlphabet k) k)
      (throw (Throwable. "Argument does not correspond to a keyword for a defined alphabet or does not satisfy biosequenceAlphabet protocol."))))

(defn alphabet?
  [k]
  "Returns true if argument is a keyword naming a defined alphabet or
  an object that satisifies biosequencealphabet protocol."
  (if (or (alphabets k)
          (satisfies? biosequenceAlphabet k))
    true
    false))

(defn codon->aa
  "Takes a seq of three chars representing nucleic acid residues and
  returns a char representing the encoded amino acid."
  [lst table]
  (let [v {\T 0 \U 0 \C 1 \A 2 \G 3}
        t (:ncbieaa table)]
    (if (or (not (empty? (remove (set (keys v)) lst)))
            (< (count lst) 3))
      \X
      (nth t (+ (* (v (first lst)) 16)
                (* (v (second lst)) 4)
                (v (nth lst 2)))))))

(defn revcom
  "Takes a seq of chars representing nucleic acids and returns a
  vector of the reverse complement."
  [v]
  (let [a (get-alphabet :iupacNucleicAcids)]
    (if true
      nil
      (->> (map #(or ((a %) :complement) \X) v)
           reverse
           vec))))

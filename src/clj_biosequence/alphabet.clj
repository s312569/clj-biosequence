(ns clj-biosequence.alphabet
  (:require [clojure.edn :as edn]
            [clojure.java.io :refer [resource]]))

(def iupacNucleicAcids
  (edn/read-string (slurp (resource "dna-iupacdna.clj"))))

(def iupacAminoAcids
  (edn/read-string (slurp (resource "aa-iupacAA.clj"))))

(def signalpAminoAcids
  (edn/read-string (slurp (resource "aa-signalp.clj"))))

(def codon-tables
  (edn/read-string (slurp (resource "codon-table.clj"))))

;; functions

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
  (vec (reverse (map #(or (get (iupacNucleicAcids %) :complement) \X)
                     v))))

(defn alphabet?
  [k]
  "Takes a keyword and returns the keyword if it is a recognised
   alphabet. Returns nil otherwise."
  (#{:iupacNucleicAcids :iupacAminoAcids :signalpAminoAcids} k))

(defn get-alphabet
  "Takes a keyword and returns the corresponding alphabet hash."
  [k]
  ({:iupacAminoAcids iupacAminoAcids
    :iupacNucleicAcids iupacNucleicAcids
    :signalpAminoAcids signalpAminoAcids} k))

(defn alphabet-is-protein
  "Returns true if alphabet is a protein alphabet."
  [alphabet]
  (if (#{:iupacAminoAcids :signalpAminoAcids} alphabet)
    true false))

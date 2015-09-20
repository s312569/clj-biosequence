(ns clj-biosequence.alphabet
  (:require [clojure.edn :as edn]
            [clojure.java.io :refer [resource]]))

;; alphabets

(def iupacAminoAcids
  {:protein? true
   :checked? true
   :alphabet (edn/read-string (slurp (resource "aa-iupacAA.clj")))})

(def iupacNucleicAcids
  {:protein? false
   :checked? true
   :alphabet (edn/read-string (slurp (resource "dna-iupacdna.clj")))})

(def signalpAminoAcids
  {:protein? true
   :checked? true
   :alphabet (edn/read-string (slurp (resource "aa-signalp.clj")))})

(def uncheckedDNA
  {:protein? false
   :checked? false})

(def uncheckedProtein
  {:protein? true
   :checked? false})

(def alphabets
  {:iupacAminoAcids iupacAminoAcids
   :iupacNucleicAcids iupacNucleicAcids
   :signalpAminoAcids signalpAminoAcids
   :uncheckedDNA uncheckedDNA
   :uncheckedProtein uncheckedProtein})

(defn alphabet-is-protein?
  [ala]
  (:protein ala))

(defn allowed-character?
  [ala char]
  (contains? (:alphabet ala) char))

(defn data
  [ala char]
  ((:alphabet ala) char))

(defn checked?
  [ala]
  (ala :checked?))

;; codon tables

(def codon-tables
  (edn/read-string (slurp (resource "codon-table.clj"))))

;; functions

(defn get-alphabet
  [k]
  (alphabets k))

(defn alphabet?
  [k]
  "Returns true if argument is a keyword naming a defined alphabet or
  an object that satisifies biosequencealphabet protocol."
  (contains? alphabets k))

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
  (let [a (:alphabet (get-alphabet :iupacNucleicAcids))]
    (->> (map #(or ((a %) :complement) \X) v)
         reverse)))

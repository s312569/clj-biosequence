(ns clj-biosequence.alphabet
  (:require [clojure.edn :as edn]
            [clojure.java.io :refer [resource]]))

(def iupacNucleicAcids
  {\A {:longname "Adenine" :complement \T}
   \C {:longname "Cytosine" :complement \G}
   \G {:longname "Guanine" :complement \C}
   \T {:longname "Thymine" :complement \A}
   \U {:longname "Uracil" :complement \A}
   \R {:longname "Purine (A or G)" :complement \Y}
   \Y {:longname "Pyrimidine (C, T, or U)" :complement \R}
   \M {:longname "C or A" :complement \K}
   \K {:longname "T, U, or G" :complement \M}
   \W {:longname "T, U, or A" :complement \W}
   \S {:longname "C or G" :complement \S}
   \B {:longname "C, T, U, or G (not A)" :complement \V}
   \D {:longname "A, T, U, or G (not C)" :complement \H}
   \H {:longname "A, T, U, or C (not G)" :complement \D}
   \V {:longname "A, C, or G (not T, not U)" :complement \B}
   \N {:longname "Any base (A, C, G, T, or U)" :complement \N}
   \X {:longname "Any base (A, C, G, T, or U)" :complement \X}})

(def iupacAminoAcids
  {\A {:three "Ala" :longname "Alanine" :monoisotopic 71.03711 :average 71.0788}
   \R {:three "Arg" :longname "Arginine" :monoisotopic 156.10111 :average 156.1875}
   \N {:three "Asn" :longname "Asparagine" :monoisotopic 114.04293 :average 114.1038}
   \D {:three "Asp" :longname "Aspartic acid" :monoisotopic 115.02694 :average 115.0886}
   \C {:three "Cys" :longname "Cysteine" :monoisotopic 103.00919 :average 103.1388}
   \Q {:three "Gln" :longname "Glutamine" :monoisotopic 128.05858 :average 128.1307}
   \E {:three "Glu" :longname "Glutamic acid" :monoisotopic 129.04259 :average 129.1155}
   \G {:three "Gly" :longname "Glycine" :monoisotopic 57.02146 :average 57.0519}
   \H {:three "His" :longname "Histidine" :monoisotopic 137.05891 :average 137.1411}
   \I {:three "Ile" :longname "Isoleucine" :monoisotopic 113.08406 :average 113.1594}
   \L {:three "Leu" :longname "Leucine" :monoisotopic 113.08406 :average 113.1594}
   \K {:three "Lys" :longname "Lysine" :monoisotopic 128.09496  :average 128.1741}
   \M {:three "Met" :longname "Methionine" :monoisotopic 131.04049 :average 131.1926}
   \F {:three "Phe" :longname "Phenylalanine" :monoisotopic 147.06841 :average 147.1766}
   \P {:three "Pro" :longname "Proline" :monoisotopic 97.05276 :average 97.1167}
   \S {:three "Ser" :longname "Serine" :monoisotopic 87.03203 :average 87.0782}
   \T {:three "Thr" :longname "Threonine" :monoisotopic 101.04768 :average 101.1051}
   \W {:three "Trp" :longname "Tryptophan" :monoisotopic 186.07931 :average 186.2132}
   \Y {:three "Tyr" :longname "Tyrosine" :monoisotopic 163.06333 :average 163.1760}
   \V {:three "Val" :longname "Valine" :monoisotopic 99.06841 :average 99.1326}
   \B {:three "Asx" :longname "Aspartic acid or Asparagine" :monoisotopic nil :average nil}
   \Z {:three "Glx" :longname "Glutamine or Glutamic acid" :monoisotopic nil :average nil}
   \J {:three "Xle" :longname "Leucine or Isoleucine" :monoisotopic 113.08406
       :average 113.1594}
   \X {:three "Xaa" :longname "Any amino acid" :monoisotopic 100.0 :average 100.0}
   \U {:three "Sec" :longname "Selenocysteine" :monoisotopic 150.953636 :average 150.0388}
   \O {:three "Pyl" :longname "Pyrrolysine" :monoisotopic 237.147727 :average 237.3018}
   \* {:three "Stop" :longname "Stop" :monoisotopic 0 :average 0}})

(def signalpAminoAcids
  {\A {:three "Ala" :longname "Alanine" :monoisotopic 71.03711 :average 71.0788 :pka 0}
   \R {:three "Arg" :longname "Arginine" :monoisotopic 156.10111 :average 156.1875 :pka 12.4}
   \N {:three "Asn" :longname "Asparagine" :monoisotopic 114.04293 :average 114.1038 :pka 0}
   \D {:three "Asp" :longname "Aspartic acid" :monoisotopic 115.02694 :average 115.0886 :pka 3.86}
   \C {:three "Cys" :longname "Cysteine" :monoisotopic 103.00919 :average 103.1388 :pka 8.33}
   \Q {:three "Gln" :longname "Glutamine" :monoisotopic 128.05858 :average 128.1307 :pka 0}
   \E {:three "Glu" :longname "Glutamic acid" :monoisotopic 129.04259 :average 129.1155 :pka 4.25}
   \G {:three "Gly" :longname "Glycine" :monoisotopic 57.02146 :average 57.0519 :pka 0}
   \H {:three "His" :longname "Histidine" :monoisotopic 137.05891 :average 137.1411 :pka 6.0}
   \I {:three "Ile" :longname "Isoleucine" :monoisotopic 113.08406 :average 113.1594 :pka 0}
   \L {:three "Leu" :longname "Leucine" :monoisotopic 113.08406 :average 113.1594 :pka 0}
   \K {:three "Lys" :longname "Lysine" :monoisotopic 128.09496  :average 128.1741 :pka 10.5}
   \M {:three "Met" :longname "Methionine" :monoisotopic 131.04049 :average 131.1926 :pka 0}
   \F {:three "Phe" :longname "Phenylalanine" :monoisotopic 147.06841 :average 147.1766 :pka 0}
   \P {:three "Pro" :longname "Proline" :monoisotopic 97.05276 :average 97.1167 :pka 0}
   \S {:three "Ser" :longname "Serine" :monoisotopic 87.03203 :average 87.0782 :pka 0}
   \T {:three "Thr" :longname "Threonine" :monoisotopic 101.04768 :average 101.1051 :pka 0}
   \W {:three "Trp" :longname "Tryptophan" :monoisotopic 186.07931 :average 186.2132 :pka 0}
   \Y {:three "Tyr" :longname "Tyrosine" :monoisotopic 163.06333 :average 163.1760 :pka 10.0}
   \V {:three "Val" :longname "Valine" :monoisotopic 99.06841 :average 99.1326 :pka 0}
   \X {:three "Xaa" :longname "Any amino acid" :monoisotopic nil :average nil :pka 0}})

;; ; average mod
;; 'I (+ cys 57.0513)
;;                'M (+ cys 45.988)
;; 'O (+ cys 1.007947)

;; ; mono mod
;; 'I (+ cys 57.021464)
;;                'M (+ cys 45.988)
;;                'O (+ cys 1.007947)

(def codon-tables
  (edn/read-string (slurp (resource "codon-table.clj"))))

;; functions

(defn codon->aa
  "Takes a seq of three chars representing nucleic acid residues and returns a
   char representing the encoded amino acid."
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
  "Takes a seq of chars representing nucleic acids and returns a vector of the
   reverse complement."
  [v]
  (vec (reverse (map #(or (get (iupacNucleicAcids %) :complement) \X) v))))

(defn alphabet?
  [k]
  "Takes a keyword and returns the keyword if it is a recognised alphabet. Returns
   nil otherwise."
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

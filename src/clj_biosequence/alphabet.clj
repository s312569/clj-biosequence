(ns clj-biosequence.alphabet)

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
  {\A {:three "Ala" :longname "Alanine"}
   \R {:three "Arg" :longname "Arginine"}
   \N {:three "Asn" :longname "Asparagine"}
   \D {:three "Asp" :longname "Aspartic acid"}
   \C {:three "Cys" :longname "Cysteine"}
   \Q {:three "Gln" :longname "Glutamine"}
   \E {:three "Glu" :longname "Glutamic acid"}
   \G {:three "Gly" :longname "Glycine"}
   \H {:three "His" :longname "Histidine"}
   \I {:three "Ile" :longname "Isoleucine"}
   \L {:three "Leu" :longname "Leucine"}
   \K {:three "Lys" :longname "Lysine"}
   \M {:three "Met" :longname "Methionine"}
   \F {:three "Phe" :longname "Phenylalanine"}
   \P {:three "Pro" :longname "Proline"}
   \S {:three "Ser" :longname "Serine"}
   \T {:three "Thr" :longname "Threonine"}
   \W {:three "Trp" :longname "Tryptophan"}
   \Y {:three "Tyr" :longname "Tyrosine"}
   \V {:three "Val" :longname "Valine"}
   \B {:three "Asx" :longname "Aspartic acid or Asparagine"}
   \Z {:three "Glx" :longname "Glutamine or Glutamic acid"}
   \J {:three "Xle" :longname "Leucine or Isoleucine"}
   \X {:three "Xaa" :longname "Any amino acid"}
   \U {:three "Sec" :longname "Selenocysteine"}
   \O {:three "Pyl" :longname "Pyrrolysine"}
   \* {:three "Stop" :longname "Stop"}})

(def codon-tables
  {1 {:name "Standard"
      :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \* \W \L \L \L \L \P
                 \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                 \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                 \G]
      :sncbieaa [\- \- \- \M \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \-
                 \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \-]}
   2 {:name "Vertebrate Mitochondrial"
      :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \L \L \L \L \P
                 \P \P \P \H \H \Q \Q \R \R \R \R \I \I \M \M \T \T \T \T \N \N
                 \K \K \S \S \* \* \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                 \G]
      :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \M \M \M \M \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                 \-]}
   3 {:name "Yeast Mitochondrial"
      :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \T \T \T \T \P
                 \P \P \P \H \H \Q \Q \R \R \R \R \I \I \M \M \T \T \T \T \N \N
                 \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                 \G]
      :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \- \- \M \M \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \-]}
   4 {:name "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"
      :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \L \L \L \L \P
                 \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                 \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                 \G]
      :sncbieaa [\- \- \M \M \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \-
                 \- \- \- \- \- \- \- \- \- \- \- \M \M \M \M \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                 \-]}
   5 {:name "Invertebrate Mitochondrial"
      :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \L \L \L \L \P
                 \P \P \P \H \H \Q \Q \R \R \R \R \I \I \M \M \T \T \T \T \N \N
                 \K \K \S \S \S \S \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                 \G]
      :sncbieaa [\- \- \- \M \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \M \M \M \M \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                 \-]}
   6 {:name "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"
      :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \Q \Q \C \C \* \W \L \L \L \L \P
                 \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                 \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                 \G]
      :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \-]}
   9 {:name "Echinoderm Mitochondrial; Flatworm Mitochondrial"
      :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \L \L \L \L \P
                 \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                 \N \K \S \S \S \S \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                 \G]
      :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                 \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                 \-]}
   10 {:name "Euplotid Nuclear" ,
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \C \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   11 {:name "Bacterial, Archaeal and Plant Plastid"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \* \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \M \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \-
                  \- \- \- \- \- \- \- \- \- \- \- \M \M \M \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   12 {:name "Alternative Yeast Nuclear"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \* \W \L \L \L \S \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   13 {:name "Ascidian Mitochondrial"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \M \M \T \T \T \T \N \N
                  \K \K \S \S \G \G \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \M \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \M \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   14 {:name "Alternative Flatworm Mitochondrial"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \Y \* \C \C \W \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \N \K \S \S \S \S \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   15 {:name "Blepharisma Macronuclear"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \Q \C \C \* \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   16 {:name "Chlorophycean Mitochondrial"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \L \C \C \* \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   21 {:name "Trematode Mitochondrial"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \M \M \T \T \T \T \N \N
                  \N \K \S \S \S \S \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   22 {:name "Scenedesmus obliquus Mitochondrial"
       :ncbieaa  [\F \F \L \L \S \S \* \S \Y \Y \* \L \C \C \* \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   23 {:name "Thraustochytrium Mitochondrial"
       :ncbieaa  [\F \F \* \L \S \S \S \S \Y \Y \* \* \C \C \* \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \R \R \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \- \- \M \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                  \-]}
   24 {:name "Pterobranchia Mitochondrial"
       :ncbieaa  [\F \F \L \L \S \S \S \S \Y \Y \* \* \C \C \W \W \L \L \L \L \P
                  \P \P \P \H \H \Q \Q \R \R \R \R \I \I \I \M \T \T \T \T \N \N
                  \K \K \S \S \S \K \V \V \V \V \A \A \A \A \D \D \E \E \G \G \G
                  \G]
       :sncbieaa [\- \- \- \M \- \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \-
                  \- \- \- \- \- \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \-
                  \- \- \- \- \- \- \- \- \- \M \- \- \- \- \- \- \- \- \- \- \-
                  \-]}})

;; functions

(defn codon->aa
  "Takes a seq of three chars representing nucleic acid residues and returns a
   char representing the encoded amino acid."
  [lst table]
  (let [v {\T 0 \U 0 \C 1 \A 2 \G 3}]
    (if (or (not (empty? (remove (set (keys v)) lst)))
            (< (count lst) 3))
      \X
      (nth (:ncbieaa table) (+ (* (v (first lst)) 16)
                               (* (v (second lst)) 4)
                               (v (nth lst 2)))))))

(defn revcom
  "Takes a seq of chars representing nucleic acids and returns a vector of the
   reverse complement."
  [v]
  (vec (reverse (map #(let [c (iupacNucleicAcids %)]
                        (if c
                          (c :complement)
                          \X))
                     v))))

(defn alphabet?
  [k]
  "Takes a keyword and returns the keyword if it is a recognised alphabet. Returns
   nil otherwise."
  (#{:iupacNucleicAcids :iupacAminoAcids} k))

(defn get-alphabet
  "Takes a keyword and returns tha corresponding alphabet hash."
  [k]
  ({:iupacAminoAcids iupacAminoAcids
    :iupacNucleicAcids iupacNucleicAcids} k))

(defn alphabet-chars
  [a]
  (let [k (keys (get-alphabet a))]
    (set (concat k (map #(Character/toLowerCase %) k)))))


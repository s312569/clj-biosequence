(in-ns 'clj-biosequence.core)

(defn protein-charge
  "Calculates the theoretical protein charge at the specified
  pH (default 7). Uses pKa values set out in the protein alphabets
  from cl-biosequence.alphabet. Considers Lys, His, Arg, Glu, Asp, Tyr
  and Cys residues only and ignores all other amino acids. The number
  of disulfides can be specified and 2 times this figure will be
  deducted from the number of Cys residues used in the calculation.
  Values used for the pKa of the N-term and C-term are 9.69 and 2.34
  respectively."
  [p & {:keys [ph disulfides] :or {ph 7 disulfides 0}}]
  {:pre [(protein? p)]}
  (let [a (ala/get-alphabet :signalpAminoAcids)
        ncalc (fn [x n]
                (* n
                   (/ (Math/pow 10 x)
                      (+ (Math/pow 10 ph)
                         (Math/pow 10 x)))))
        ccalc (fn [x n]
                (* n
                   (/ (Math/pow 10 ph)
                      (+ (Math/pow 10 ph)
                         (Math/pow 10 x)))))
        freq (frequencies (bs-seq p))
        cys (let [c (freq \C)
                  f (if c (- c (* 2 disulfides)) 0)]
              (if (>= f 0) f
                  (throw (Throwable.
                          "More disulfides specified than Cys residues."))))]
    (- (reduce +
               (cons (ncalc 9.69 1)
                     (map (fn [[x n]] (ncalc (:pka (a x)) n))
                          (select-keys freq [\K \H \R]))))
       (+ (ccalc (:pka (a \C)) cys)
          (reduce + (cons (ccalc 2.34 1)
                          (map (fn [[x n]] (ccalc (:pka (a x)) n))
                               (select-keys freq [\E \D \Y]))))))))

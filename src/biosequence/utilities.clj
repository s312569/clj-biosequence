(ns biosequence.utilities
  (:require [clojure.java.io :as io]
            [clojure.string :as string]
            [fs.core :as fs]
            [clojure.pprint :as pp]))


;; protein mass

(defn aa-average-masses
  "Returns a map hash of average amino acid masses. The optional argument c-mod specifies Cys modifications. Default is no modification and current permissable values are:
   I - Iodoacetamide - 160.1608
   M - MMTS - 149.1268"
  ([] (aa-average-masses nil))
  ([c-mod]
     (let [cys 103.1388
           c (condp = c-mod
               'I (+ cys 57.0513)
               'M (+ cys 45.988)
               'O (+ cys 1.007947)
               cys)]
       (merge {\C c}
              {\A 71.0788
               \R 156.1875
               \N 114.1038
               \D 115.0886
               \E 129.1155
               \Q 128.1307
               \G 57.0519
               \H 137.1411
               \I 113.1594
               \L 113.1594
               \K 128.1741
               \M 131.1926
               \F 147.1766
               \P 97.1167
               \S 87.0782
               \T 101.1051
               \W 186.2132
               \Y 163.1760
               \V 99.1326
               \* 0
               \X 110}))))

(defn aa-monoisotopic-masses
  "Returns a map hash of amino acid masses. The optional argument c-mod specifies Cys modifications. Default is no modification and current permissable values are:
   I - Iodoacetamide - 160.1608
   M - MMTS - 149.1268"
  ([] (aa-monoisotopic-masses nil))
  ([c-mod]
     (let [cys 103.00919
           c (condp = c-mod
               'I (+ cys 57.021464)
               'M (+ cys 45.988)
               'O (+ cys 1.007947)
               cys)]
       (merge {\C c}
              {\A 71.03711
               \R 156.10111
               \N 114.04293
               \D 115.02694
               \E 129.04259
               \Q 128.05858
               \G 57.02146
               \H 137.05891
               \I 113.08406
               \L 113.08406
               \K 128.09496
               \M 131.04049
               \F 147.06841
               \P 97.05276
               \S 87.03203
               \T 101.04768
               \W 186.07931
               \Y 163.06333
               \V 99.06841
               \* 0
               \X 110}))))

(defn string-mass
  "Returns the mass of a string of amino acids in single letter format. Monoisotopic or average can be specified using the 'average' argument, false or nil will give the monoisotopic values. Default is false. A Cys modification can be specified using the 'c-mod' argument. Default is no modification and current permissable values are:
   I - Iodoacetamide - +57.022
   M - MMTS - +45.988
   O - oxidised - +1.007947"
  ([st] (string-mass st nil nil))
  ([st average] (string-mass st average nil))
  ([st average c-mod]
     (let [aas (if average
                 (aa-average-masses c-mod)
                 (aa-monoisotopic-masses c-mod))]
       (+ 17.0074
          1.007947
          (reduce + (map #(get aas %) (string/upper-case st)))))))

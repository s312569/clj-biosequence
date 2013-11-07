(in-ns 'clj-biosequence.core)

(defn if-string-int
  "If a string an integer is parsed, if not returns e. Will throw an
   expception if no integer can be parsed if `error?' is true. Used
   only when parsing optional fields in files where value could be nil
   or a string representation of an integer."
  ([e] (if-string-int e true))
  ([e error?]
     (try
       (if (string? e) (Integer/parseInt e) e)
       (catch NumberFormatException e
         (if error?
           (throw e))))))

;; sequence

(defn clean-sequence
  "Removes spaces and newlines and checks that all characters are
   legal characters for the supplied alphabet. Replaces non-valid
   characters with \\X. If `a' is not a defined alphabet throws an
   exception."
  [s a]
  (let [k (complement (ala/alphabet-chars a))
        w #{\space \newline}]
    (vec (remove nil? (map #(cond (k %) \X
                                  (w %) nil
                                  :else %) (vec s))))))

;; store utilities

(defn biosequence-save
  [this project name type]
  (let [i (st/init-biosequence-collection name (:name project) type)]
    (with-open [r (bs-reader this)]
      (st/save-list (pmap #(hash-map :acc (accession %)
                                  :src (pr-str %))
                       (biosequence-seq r))
                 i))
    i))

;; printing

(defn print-biosequence
  [biosequence w]
  (print-tagged biosequence w))

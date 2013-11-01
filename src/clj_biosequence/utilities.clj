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

;; for print-method defs

(defn print-tagged
  "Used for printing objects tagged so that edn/read-string can read
  them in."
  [obj w]
  (tag/pr-tagged-record-on obj w))


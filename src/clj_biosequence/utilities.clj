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

;; printing objects

(defn print-tagged
  "Used for printing objects tagged so that edn/read-string can read
  them in."
  [obj w]
  (time (tag/pr-tagged-record-on obj w)))

(defn my-tag->factory
  "Returns the map-style record factory for the `tag` symbol.  Returns nil if `tag` does not
  refer to a record."
  [tag]
  (when (namespace tag)
    (resolve (symbol (str (namespace tag) "/map->" (name tag))))))

(defn default-reader
  [tag value]
  (if-let [factory (and (map? value)
                     (my-tag->factory tag))]
    (factory value)
    (throw (Throwable. (str "Record not supported: " tag)))))

(defn bs-read
  [h]
  (if-let [o (ed/read-string
              {:default #'default-reader
               :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
              (:src h))]
    (merge o (dissoc h :src))))

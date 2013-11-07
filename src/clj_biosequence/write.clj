(ns clj-biosequence.write
  (:require [miner.tagged :as tag]
            [clojure.data.xml :as xml]
            [clojure.edn :as ed]))

(defn print-tagged
  "Used for printing objects tagged so that edn/read-string can read
  them in."
  [obj w]
  (tag/pr-tagged-record-on obj w))

(defn my-tag->factory
  "Returns the map-style record factory for the `tag` symbol. Returns
  nil if `tag` does not refer to a record."
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
  [s]
  (ed/read-string
   {:default #'default-reader
    :readers {'clojure.data.xml.Element clojure.data.xml/map->Element}}
   s))

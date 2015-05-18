(defproject clj-biosequence "0.3.5"
  :description "Library for the manipulation of biological sequences."
  :url ""
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :repositories [["biojava" "http://www.biojava.org/download/maven/"]]
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [org.clojure/data.xml "0.0.8"]
                 [org.clojure/data.zip "0.1.1"]
                 [org.apache.commons/commons-compress "1.9"]
                 [com.taoensso/nippy "2.7.1"]
                 [com.velisco/tagged "0.3.4"]
                 [clj-http "1.0.1"]
                 [iota "1.1.2"]
                 [clj-time "0.9.0"]
                 [fs "1.3.3"]
                 [org.clojars.hozumi/clj-commons-exec "1.1.0"]]
  :resource-paths ["shared" "resources"]
  :plugins [[codox "0.8.10"]]
  :codox {:src-dir-uri "https://github.com/s312569/clj-biosequence/blob/master/"
          :src-linenum-anchor-prefix "L"}
  :repl-options {:init (set! *print-length* 100)}
  :jvm-opts ["-Xmx1000M"])

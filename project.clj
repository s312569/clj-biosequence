(defproject clj-biosequence "0.1.4-SNAPSHOT"
  :description "Library for the manipulation of biological sequences."
  :url ""
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :repositories [["biojava" "http://www.biojava.org/download/maven/"]]
  :dependencies [[org.clojure/clojure "1.6.0"]
                 [org.clojure/data.xml "0.0.7"]
                 [org.clojure/data.zip "0.1.1"]
                 [org.apache.commons/commons-compress "1.7"]
                 [com.novemberain/monger "1.7.0"]
                 [com.taoensso/nippy "2.5.2"]
                 [com.velisco/tagged "0.3.4"]
                 [clj-http "1.0.0"]
                 [iota "1.1.2"]
                 [clj-time "0.8.0"]
                 [fs "1.3.3"]
                 [org.clojars.hozumi/clj-commons-exec "1.0.7"]]
  :resource-paths ["shared" "resources"]
  :plugins [[codox "0.6.4"]]
  :codox {:src-dir-uri "https://github.com/s312569/clj-biosequence/blob/master/"
          :src-linenum-anchor-prefix "L"}
  :repl-options {:init (set! *print-length* 100)}
  :jvm-opts ["-Xmx1000M"])

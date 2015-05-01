(in-ns 'clj-biosequence.core)

(defrecord gzipped [file])
(defrecord zipped [file])
(defrecord bzipped [file])
(defrecord uncompressed [file])

(defn- compressed?
  [f]
  (condp = (extension f)
      ".gz" (->gzipped f)
      ".zip" (->zipped f)
      ".bz2" (->bzipped f)
      (->uncompressed f)))

(defprotocol bioReader
  (make-reader [x opts]))

(def default-reader
  {:make-reader (fn [x o] (apply reader x o))})

(extend String
  bioReader
  (assoc default-reader
         :make-reader (fn [x o] (make-reader (compressed? x) o))))

(extend java.io.File
  bioReader
  (assoc default-reader
         :make-reader (fn [x o] (make-reader (compressed? x) o))))

(extend Object
  bioReader
  default-reader)

(extend uncompressed
  bioReader
  (assoc default-reader
         :make-reader (fn [x o] (apply reader (:file x) o))))

(extend gzipped
  bioReader
  (assoc default-reader
         :make-reader (fn [x o] (apply reader
                                       (java.util.zip.GZIPInputStream.
                                        (input-stream (:file x)))
                                       o))))

(extend zipped
  bioReader
  (assoc default-reader
         :make-reader (fn [x o] (apply reader
                                       (java.util.zip.ZipInputStream.
                                        (input-stream (:file x)))
                                       o))))

(extend bzipped
  bioReader
  (assoc default-reader
         :make-reader (fn [x o] (apply reader
                                       (BZip2CompressorInputStream.
                                        (input-stream (:file x)))
                                       o))))

(defn bioreader
  [x & opts]
  (make-reader x opts))

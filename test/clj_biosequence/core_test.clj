(ns clj-biosequence.core-test
  (:require [clojure.java.io :as io]
            [fs.core :as fs]
            [clj-time.format :as ti]
            [clojure.string :as st])
  (:use clojure.test
        clj-biosequence.core
        clj-biosequence.uniprot
        clj-biosequence.signalp
        clj-biosequence.genbank
        clj-biosequence.blast
        clj-biosequence.fastq
        clj-biosequence.entrezgene))

(def fasta-lower-nuc (with-open [f (bs-reader
                              (init-fasta-file
                               (io/resource "test-files/lowercase-nuc.fasta")
                               :iupacNucleicAcids))]
                 (first (biosequence-seq f))))

(def fasta-lower-aa (with-open [f (bs-reader
                              (init-fasta-file
                               (io/resource "test-files/lowercase-AA.fasta")
                               :uncheckedProtein))]
                 (first (biosequence-seq f))))

(def fasta-nuc (with-open [f (bs-reader
                              (init-fasta-file
                               (io/resource "test-files/nuc-sequence.fasta")
                               :iupacNucleicAcids))]
                 (first (biosequence-seq f))))

(def un-seq (with-open [f (bs-reader
                           (init-uniprot-file
                            (io/resource
                             "test-files/uniprot-s-mansoni-20121217.xml")))]
              (first (biosequence-seq f))))

(deftest fasta
  (testing "Fasta"
    (is (= "gi|116025203|gb|EG339215.1|EG339215"
           (accession fasta-nuc)))
    (is (= (list "gi|116025203|gb|EG339215.1|EG339215")
           (accessions fasta-nuc)))
    (is (= (vec "KAAN-aaa29f08.b1")
           (subvec (vec (description fasta-nuc)) 0 16)))
    (is (= [\G \T \A \C \A \A \A]
           (subvec (bs-seq fasta-nuc) 0 7)))
    (is (= false (protein? fasta-nuc)))
    (is (= :iupacNucleicAcids (alphabet fasta-nuc)))
    (is (= [\G \T \A \C \A \A \A]
           (bs-seq (sub-bioseq fasta-nuc 1 7))))
    (is (= [\G \T \A]
           (first (partition-all 3 (bs-seq fasta-nuc)))))
    (is (= [\T \A \A \T]
           (bs-seq (sub-bioseq (reverse-comp fasta-nuc) 1 4))))
    (is (= [\A \T \T \A]
           (bs-seq (sub-bioseq (reverse-seq fasta-nuc) 1 4))))
    ;; lower-case tests
    (is (= [\G \T \A \C \A \A \A]
           (bs-seq (sub-bioseq fasta-lower-nuc 1 7))))
    (is (= "anact"
           (subs (:sequence fasta-lower-aa) 0 5)))
    (is (= "accession"
           (with-open [s
                       (bs-reader (init-fasta-string
                                   ">accession desc\n sequence"
                                   :iupacNucleicAcids))]
             (accession (first (biosequence-seq s)))))))
  (testing "Fasta indexing"
    (let [i (index-biosequence-file
             (init-fasta-file
              (io/resource "test-files/nuc-sequence.fasta")
              :iupacNucleicAcids))]
      (try
        (with-open [r (bs-reader i)]
          (is (= "gi|114311762|gb|EE738912.1|EE738912"
                 (accession (first (biosequence-seq r)))))
          (is (= "gi|114314166|gb|EE741316.1|EE741316"
                 (accession
                  (get-biosequence
                   r
                   "gi|114314166|gb|EE741316.1|EE741316")))))
        (let [ii (load-biosequence-index
                  (fs/absolute-path
                   (io/resource "test-files/nuc-sequence.fasta")))]
          (with-open [rr (bs-reader ii)]
            (is (= "gi|114311762|gb|EE738912.1|EE738912"
                   (accession (first (biosequence-seq rr)))))
            (is (= "gi|114314166|gb|EE741316.1|EE741316"
                   (accession
                    (get-biosequence
                     rr "gi|114314166|gb|EE741316.1|EE741316"))))))
        (finally
          (delete-indexed-biosequence i))))))

(deftest uniprot
  (testing "Uniprot"
    (is (= "C4PYP8" (accession un-seq)))
    (is (= (list "C4PYP8" "G4VQR6" "DRE2_SCHMA")
           (accessions un-seq)))
    (is (= "Anamorsin homolog [Schistosoma mansoni]"
           (description un-seq)))
    (is (= [\M \E \Q \C \V \A]
           (subvec (bs-seq un-seq) 0 6)))
    (is (= true (protein? un-seq)))
    (is (= :iupacAminoAcids (alphabet un-seq)))
    (is (= [\M \E \Q \C \V \A] (bs-seq (sub-bioseq un-seq 1 6))))
    (is (= [\M \E \Q] (first (partition-all 3 (bs-seq un-seq)))))
    (is (= "Schistosoma mansoni" (-> (tax-refs un-seq)
                                     first tax-name)))
    (is (= "Eukaryota" (-> (tax-refs un-seq)
                           first
                           lineage
                           (st/split #";")
                           first)))
    (is (= "DRE2_SCHMA" (-> (accessions un-seq) last)))
    (is (= "Anamorsin homolog" (-> (names un-seq) first)))
    (is (= 29662.0 (calc-mol-wt un-seq)))
    (is (= "Smp_207000" (-> (genes un-seq) first orf)))
    (testing "Citations"
      (let [c (first (citations un-seq))]
        (is (= "The genome of the blood fluke Schistosoma mansoni."
               (title c)))
        (is (= "Nature" (journal c)))
        (is (= 2009 (year c)))
        (is (= "460" (volume c)))
        (is (= "352" (pstart c)))
        (is (= "358" (pend c)))
        (is (= "Berriman M." (first (authors c))))))
    (is (= "function" (-> (comments un-seq)
                          first
                          obj-type)))
    (is (= "chain" (-> (feature-seq un-seq)
                       first
                       obj-type)))
    (is (= 1 (-> (feature-seq un-seq) first
                 intervals first start))))
  (testing "Uniprot indexing"
    (let [i (index-biosequence-file
             (init-uniprot-file
              (io/resource
               "test-files/uniprot-s-mansoni-20121217.xml")))]
      (try
        (with-open [r (bs-reader i)]
          (is (= "P35661"
                 (accession (first (biosequence-seq r)))))
          (is (= "P35661"
                 (accession (get-biosequence r "P35661")))))
        (let [ii (-> (io/resource
                      "test-files/uniprot-s-mansoni-20121217.xml")
                     fs/absolute-path
                     load-biosequence-index)]
          (with-open [rr (bs-reader ii)]
            (is (= "P35661"
                   (accession (first (biosequence-seq rr)))))
            (is (= "P35661"
                   (accession (get-biosequence rr "P35661"))))))
        (finally
          (delete-indexed-biosequence i))))))

(deftest entrez-gene
  (testing "Entrezgene"
    (let [egf (with-open [r (bs-reader
                             (init-entrezgene-file
                              (io/resource
                               "test-files/entrez-gene.xml")))]
                (first (biosequence-seq r)))]
      (is (= 3875 (accession egf)))
      (is (= (list 3875) (accessions egf)))
      (is (= '("K18" "CYK18") (synonyms egf)))
      (is (= "keratin 18" (description egf)))
      (is (= "HGNC" (-> (get-db-refs egf)
                        first
                        database-name)))
      (is (= "HGNC:6430" (-> (get-db-refs egf)
                             first
                             object-id)))
      (is (= "1" (trans-table egf)))
      (is (= "KRT18 encodes the type I intermediate filament chain keratin 18. Keratin 18, together with its filament partner keratin 8, are perhaps the most commonly found members of the intermediate filament gene family. They are expressed in single layer epithelial tissues of the body. Mutations in this gene have been linked to cryptogenic cirrhosis. Two transcript variants encoding the same protein have been found for this gene. [provided by RefSeq, Jul 2008]"
             (summary egf)))
      (is (= :iupacNucleicAcids (alphabet egf)))
      (is (= "6" (moltype egf)))
      (is (not (protein? egf)))
      (is (= :iupacNucleicAcids (alphabet egf)))
      (is (= "Homo sapiens" (-> (tax-refs egf)
                                first
                                tax-name)))
      (is (= "man" (-> (tax-refs egf)
                       first
                       synonyms
                       first)))
      (is (= "Eukaryota" (-> (tax-refs egf)
                             first
                             lineage
                             (st/split #"; ")
                             first)))
      (is (= "human" (-> (tax-refs egf) first common-name)))
      (is (= "9606" (-> (tax-refs egf)
                        first
                        get-db-refs
                        first
                        object-id)))
      (is (= "12" (-> (entrez-biosource egf)
                      subsource
                      first
                      subsource-name)))
      (is (= ["chromosome" "1"] (-> (entrez-biosource egf)
                                    subsource
                                    first
                                    subsource-type)))
      (is (= ["natural" 1] (-> (entrez-biosource egf)
                               origin)))
      (is (= "keratin, type I cytoskeletal 18"
             (-> (proteins egf)
                 first
                 synonyms
                 first)))
      (is (= "keratin, type I cytoskeletal 18"
             (-> (proteins egf)
                 first
                 description)))
      (is (= "KRT18" (locus egf)))
      (is (= "12q13" (map-location egf)))
      (is (= "PIG46" (locus-tag egf)))
      (let [c (first (comments egf))]
        (is (= 254 (obj-type c)))
        (is (= "comment" (obj-value c)))
        (is (= "RefSeq Status" (obj-heading c)))
        (is (= "REVIEWED" (obj-label c))))
      (is (= 23 (->> (comments egf)
                     (filter #(seq (comments %)))
                     first
                     obj-type)))
      (is (= "genomic" (-> (comments egf)
                           (nth 5)
                           product-comments
                           first
                           obj-value)))
      (is (= 254 (-> (property-comments egf)
                     first
                     obj-type)))
      (is (= "NC_000012" (->> (comments egf)
                              (some accession))))
      (is (= '("NC_000012") (->> (comments egf)
                                 (map accessions)
                                 (some seq))))
      (is (= "12" (->> (comments egf)
                       (some version))))
      (is (= '(["UNIGENE" "Hs.406013"])
             (-> (walk-comments egf
                                obj-heading "Additional Links"
                                obj-text "UniGene")
                 first
                 xtra-terms)))
      (is (= "24660549" (->> (walk-comments egf obj-type 254)
                             (map pmids)
                             (some seq)
                             first)))
      (let [i (->> (comments egf)
                   (mapcat seq-locations)
                   first
                   intervals
                   first)]
        (is (= 52948870 (start i)))
        (is (= 52952900 (end i)))
        (is (not (comp? i))))
      (is (->> (walk-comments egf obj-type 254
                              obj-type 1)
               (mapcat seq-locations)
               (mapcat intervals)
               (some comp?)))
      (is (= "Review record(s) in Gene"
             (->> (comments egf)
                  (mapcat other-sources)
                  first
                  anchor)))
      (is (= "70 found"
             (->> (comments egf)
                  (mapcat other-sources)
                  first
                  pre-text)))
      (is (= "(e-PCR)" (->> (walk-comments egf
                                           obj-type 254
                                           obj-type 254)
                            (mapcat other-sources)
                            (some post-text))))
      (is (= "http://www.ncbi.nlm.nih.gov/sutils/evv.cgi?taxid=9606&contig=NT_029419.13&gene=KRT18&lid=3875&from=15713618&to=15717648"
             (->> (walk-comments egf
                                 obj-type 254
                                 obj-type 254)
                  (mapcat other-sources)
                  (some url))))
      (is (= "HuGENet" (->> (comments egf)
                            (mapcat other-sources)
                            (mapcat get-db-refs)
                            first
                            database-name)))
      (is (= "3875.16143128" (->> (comments egf)
                                  (mapcat other-sources)
                                  (mapcat get-db-refs)
                                  first
                                  object-id)))
      (is (= "3875" (-> (entrez-gene-source egf)
                        source-id)))
      (is (= "LocusLink" (-> (entrez-gene-source egf)
                             source-name)))
      (is (= "3875" (-> (entrez-gene-source egf)
                        string2)))
      (is (= "12q13" (-> (entrez-location egf)
                         first
                         display-string)))
      (is (= "cyto" (-> (entrez-location egf)
                        first
                        map-method-type)))
      (is (= "1991-09-20" (parse-date (creation-date egf))))
      (is (= "2014-12-30" (parse-date (update-date egf))))
      (is (= ["live" 0] (status egf))))))

(deftest genbank-test
  (testing "Genbank"
    (let [gsn (with-open [r (bs-reader
                             (init-genbank-file
                              (io/resource "test-files/nucleotide-gb.xml")
                              :iupacNucleicAcids))]
                (second (biosequence-seq r)))]
      (is (= "GAAZ01003035" (accession gsn)))
      (is (= "TSA: Crotalus horridus Chorr_CTL-10 mRNA sequence")
          (description gsn))
      (is (= false (protein? gsn)))
      (is (= ">gb|GAAZ01003035|gb|GAAZ01003035.1||gnl|TSA:GAAZ01|Chorr_CTL-10|gi|521752463| TSA: Crotalus horridus Chorr_CTL-10 mRNA sequence\nctccattcccctcataggaggaaagcctggagttgcctctgagcagacttgctacctgtggaggccgaggaacagttctctctgcagggaaggaaagaacgccatggggcgattcatcttcgtgagcttcaacttgctggtcgtgttcctctccctaagtggaactctagctgatttggaatgtccctccggttggtcttcctatgatcggtattgctacaagcccttcaaacaagagatgacctgggccgatgcagagaggttctgctcggagcaggcgaagggcgggcatctcctctctgtcgaaaccgccctagaagcatcctttgtggacaatgtgctctatgcgaacaaagagtacctcacacgttatatctggattggactgagggttcaaaacaaaggacagccatgctccagcatcagttatgagaacctggttgacccatttgaatgttttatggtgagcagagacacaaggcttcgtgagtggtttaaagttgactgtgaacaacaacattctttcatatgcaagttcacgcgaccacgttaagatccggctgtgtgaagtctggagaagcaaggaagccccccacctctccccaccccccacctgccgcaatctctg\n") (fasta-string gsn))
      (is (= :iupacNucleicAcids (alphabet gsn)))
      (is (= (parse-date (creation-date gsn)) "2013-07-08"))
      (is (= (parse-date (update-date gsn)) "2013-07-08"))
      (is (= 1 (version gsn)))
      (is (= "35024" (-> (tax-refs gsn)
                         first
                         get-db-refs
                         first
                         object-id)))
      (is (= "Eukaryota" (-> (tax-refs gsn)
                             first
                             lineage
                             (st/split #";\s")
                             first)))
      (is (= "Crotalus horridus" (-> (tax-refs gsn)
                                     first
                                     tax-name)))
      (is (= "mRNA" (moltype gsn)))
      (is (= "GAAZ01003035" (locus gsn)))
      (is (= "GAAZ01003035"
             (with-open [r (bs-reader
                            (init-genbank-connection '("GAAZ01003035")
                                                     :nucest :xml
                                                     :iupacNucleicAcids))]
               (accession (first (biosequence-seq r))))))
      (is (= "Crotalus horridus" (-> (tax-refs gsn)
                                     first
                                     tax-name)))
      (is (= "source" (-> (feature-seq gsn) first obj-type)))
      (is (= "1..628" (-> (feature-seq gsn) first feature-location)))
      (is (= 1 (-> (feature-seq gsn) first intervals first start)))))
  (testing "Genbank indexing"
    (let [i (index-biosequence-file
             (init-genbank-file
              (io/resource "test-files/nucleotide-gb.xml")
              :iupacNucleicAcids))]
      (try
        (with-open [r (bs-reader i)]
          (is (= "GAAZ01003041"
                 (accession (first (biosequence-seq r)))))
          (is (= "GAAZ01003035"
                 (accession (get-biosequence r "GAAZ01003035")))))
        (let [ii (load-biosequence-index
                  (fs/absolute-path
                   (io/resource "test-files/nucleotide-gb.xml")))]
          (with-open [rr (bs-reader ii)]
            (is (= "GAAZ01003041"
                   (accession (first (biosequence-seq rr)))))
            (is (= "GAAZ01003035"
                   (accession (get-biosequence rr "GAAZ01003035"))))))
        (finally
          (delete-indexed-biosequence i))))))

(deftest blast-test
  (testing "Blast"
    (let [toxindb (init-blast-db
                   (io/resource "test-files/toxins.fasta")
                   :iupacAminoAcids)
          toxins (init-fasta-file
                  (io/resource "test-files/toxins.fasta")
                  :iupacAminoAcids)
          tox-bl (with-open [r (bs-reader toxins)]
                   (blast (take 20 (biosequence-seq r))
                          "blastp"
                          toxindb
                          "/tmp/blast.xml"))]
      (is (= "B3EWT5"
             (with-open [r (bs-reader tox-bl)]
               (first
                (->> (biosequence-seq r)
                     (filter #(>= (-> (hit-seq %)
                                      second hit-bit-scores first) 50))
                     (map #(-> (hit-seq %) second hit-accession)))))))
      (testing "Blast indexing"
        (let [i (index-biosequence-file tox-bl)]
          (try
            (with-open [r (bs-reader i)]
              (is (= "sp|P0CE81|A1HB1_LOXIN"
                     (accession (get-biosequence
                                 r
                                 "sp|P0CE81|A1HB1_LOXIN")))))
            (let [ii (load-biosequence-index (bs-path tox-bl))]
              (with-open [rr (bs-reader ii)]
                (is (= "sp|P0CE81|A1HB1_LOXIN"
                       (accession (get-biosequence
                                   rr
                                   "sp|P0CE81|A1HB1_LOXIN"))))))))))))

(deftest fastq-test
  (testing "Fastq"
    (let [fs (with-open [r (bs-reader (init-fastq-file
                                       (io/resource "test-files/fastq-test.fastq")))]
               (first (biosequence-seq r)))]
      (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:2860:1993 1:N:0:ACAGTG"
             (accession fs)))
      (is (= [\C \G \C] (first (partition-all 3 (bs-seq fs)))))
      (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:2860:1993 1:N:0:ACAGTG"
             (description fs)))
      (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:2860:1993 1:N:0:ACAGTG\nCGCTGTACTCGATTATATGTCGATGTAACTTTTTCTGTACGTTTTAACTGGCATGTTTTTCTATATTAGATCGTGCGAGAATCACAGTACCTTAGTGGGG\n+\nCCCFFFFFHHHFHIIJJJJJIJJJIHJJIJJJIJJIIIJJJHIJIJIJJJJJJIIGCHIJJJHEIIIJJDGGIHGEFEFCDEEEDDDCDDCCDDC@CBDB\n"
             (fastq->string fs)))))
  (testing "Fastq indexing"
    (let [i (index-biosequence-file
             (init-fastq-file
              (io/resource "test-files/fastq-test.fastq")))]
      (try
        (with-open [r (bs-reader i)]
          (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG"
                 (accession (get-biosequence r "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG")))))
        (let [ii (-> (io/resource "test-files/fastq-test.fastq")
                     fs/absolute-path
                     load-biosequence-index)]
          (with-open [rr (bs-reader ii)]
            (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG"
                   (accession
                    (get-biosequence
                     rr
                     "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG"))))))
        (finally
          (delete-indexed-biosequence i))))))

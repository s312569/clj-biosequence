(ns clj-biosequence.core-test
  (:require [clojure.java.io :as io]
            [fs.core :as fs])
  (:use clojure.test
        clj-biosequence.core
        clj-biosequence.uniprot
        clj-biosequence.signalp
        clj-biosequence.genbank
        clj-biosequence.blast
        clj-biosequence.fastq
        clj-biosequence.citation))

(def fasta-nuc (with-open [f (bs-reader
                              (init-fasta-file
                               (io/resource "test-files/nuc-sequence.fasta")
                               :iupacNucleicAcids))]
                 (first (biosequence-seq f))))

(def un-seq (with-open [f (bs-reader
                           (init-uniprotxml-file
                            (io/resource
                             "test-files/uniprot-s-mansoni-20121217.xml")))]
              (first (biosequence-seq f))))

(deftest fasta
  (testing "Fasta"
    (is (= "gi|116025203|gb|EG339215.1|EG339215" (accession fasta-nuc)))
    (is (= (list "gi|116025203|gb|EG339215.1|EG339215") (accessions fasta-nuc)))
    (is (= (vec "KAAN-aaa29f08.b1") (subvec (vec (def-line fasta-nuc)) 0 16)))
    (is (= [\G \T \A \C \A \A \A] (subvec (bs-seq fasta-nuc) 0 7)))
    (is (= false (protein? fasta-nuc)))
    (is (= :iupacNucleicAcids (alphabet fasta-nuc)))
    (is (= java.lang.String (class (bioseq->string fasta-nuc))))
    (is (= [\G \T \A \C \A \A \A] (bs-seq (sub-bioseq fasta-nuc 0 7))))
    (is (= [\G \T \A] (first (partition-all 3 (bs-seq fasta-nuc)))))
    (is (= [\T \A \A \T] (bs-seq (sub-bioseq (reverse-comp fasta-nuc) 0 4))))
    (is (= [\A \T \T \A] (bs-seq (sub-bioseq (reverse-seq fasta-nuc) 0 4))))
    (is (= "accession" (with-open [s
                                   (bs-reader (init-fasta-string
                                               ">accession desc\n sequence"
                                               :iupacNucleicAcids))]
                         (accession (first (biosequence-seq s)))))))
  (testing "Fasta indexing"
    (let [i (index-biosequence-file
             (init-fasta-file (io/resource "test-files/nuc-sequence.fasta")
                              :iupacNucleicAcids))]
      (try
        (with-open [r (bs-reader i)]
          (is (= "gi|116025203|gb|EG339215.1|EG339215"
                 (accession (first (biosequence-seq r)))))
          (is (= "gi|114314166|gb|EE741316.1|EE741316"
                 (accession (get-biosequence
                             r
                             "gi|114314166|gb|EE741316.1|EE741316")))))
        (let [ii (load-biosequence-index
                  (fs/absolute-path
                   (io/resource "test-files/nuc-sequence.fasta")))]
          (with-open [rr (bs-reader ii)]
            (is (= "gi|116025203|gb|EG339215.1|EG339215"
                   (accession (first (biosequence-seq rr)))))
            (is (= "gi|114314166|gb|EE741316.1|EE741316"
                   (accession (get-biosequence
                               rr "gi|114314166|gb|EE741316.1|EE741316"))))))
        (finally
          (delete-indexed-biosequence i))))))

(deftest uniprot
  (testing "Uniprot"
    (is (= "C4PYP8" (accession un-seq)))
    (is (= (list "C4PYP8" "G4VQR6") (accessions un-seq)))
    (is (= (vec "Anamorsin") (subvec (vec (def-line un-seq)) 0 9)))
    (is (= [\M \E \Q \C \V \A] (subvec (bs-seq un-seq) 0 6)))
    (is (= true (protein? un-seq)))
    (is (= :iupacAminoAcids (alphabet un-seq)))
    (is (= java.lang.String (class (bioseq->string un-seq))))
    (is (= [\M \E \Q \C \V \A] (bs-seq (sub-bioseq un-seq 0 6))))
    (is (= [\M \E \Q] (first (partition-all 3 (bs-seq un-seq)))))
    (is (= "C4PYP8"
           (with-open [f (-> (io/resource "test-files/uniprot-s-mansoni-20121217.xml")
                             slurp
                             init-uniprot-string
                             bs-reader)]
             (accession (first (biosequence-seq f))))))
    (is (= nil
           (with-open [f (-> (io/resource "test-files/uniprot-s-mansoni-20121217.xml")
                             init-uniprotxml-file
                             bs-reader)]
             (parameters f))))
    (is (= "Schistosoma mansoni" (get (organism un-seq) "scientific")))
    (is (= "Eukaryota" (first (lineage un-seq))))
    (is (= "DRE2_SCHMA" (prot-name un-seq)))
    (is (= "Anamorsin homolog" (recommended-name un-seq)))
    (is (= 29662.0 (:mass (sequence-info un-seq))))
    (is (= "Smp_207000" (:gene (first (gene un-seq)))))
    (testing "Citations"
      (let [c (first (references un-seq))]
        (is (= "The genome of the blood fluke Schistosoma mansoni." (title c)))
        (is (= "Nature" (journal c)))
        (is (= "2009" (year c)))
        (is (= "460" (volume c)))
        (is (= "352" (pstart c)))
        (is (= "358" (pend c)))
        (is (= "Berriman M." (first (authors c))))))
    (is (= "function" (:type (:attrs (first (comments un-seq))))))
    (is (= "chain" (feature-type (first (feature-seq un-seq)))))
    (is (= 1 (-> (feature-seq un-seq) first interval-seq first start))))
  (testing "Uniprot indexing"
    (let [i (index-biosequence-file
             (init-uniprotxml-file
              (io/resource "test-files/uniprot-s-mansoni-20121217.xml")))]
      (try
        (with-open [r (bs-reader i)]
          (is (= "C4PYP8"
                 (accession (first (biosequence-seq r)))))
          (is (= "P35661"
                 (accession (get-biosequence r "P35661")))))
        (let [ii (-> (io/resource
                      "test-files/uniprot-s-mansoni-20121217.xml")
                     fs/absolute-path
                     load-biosequence-index)]
          (with-open [rr (bs-reader ii)]
            (is (= "C4PYP8"
                   (accession (first (biosequence-seq rr)))))
            (is (= "P35661"
                   (accession (get-biosequence rr "P35661"))))))
        (finally
          (delete-indexed-biosequence i))))))

(deftest signalp-test
  (testing "Signalp"
    (let [bs (with-open [r (bs-reader
                            (init-fasta-file
                             (io/resource "test-files/toxins.fasta")
                             :iupacAminoAcids))]
               (doall (take 20 (biosequence-seq r))))]
      (is (= "sp|C1IC47|3FN3_WALAE"
             (accession (first (filter-signalp bs))))))))

(deftest genbank-test
  (testing "Genbank"
    (let [gsn (with-open [r (bs-reader
                             (init-genbank-file
                              (io/resource "test-files/nucleotide-gb.xml")))]
                (second (biosequence-seq r)))]
      (is (= "GAAZ01003035" (accession gsn)))
      (is (= "TSA: Crotalus horridus Chorr_CTL-10 mRNA sequence") (def-line gsn))
      (is (= false (protein? gsn)))
      (is (= ">gb|GAAZ01003035|gb|GAAZ01003035.1||gnl|TSA:GAAZ01|Chorr_CTL-10|gi|521752463| TSA: Crotalus horridus Chorr_CTL-10 mRNA sequence\nctccattcccctcataggaggaaagcctggagttgcctctgagcagacttgctacctgtggaggccgaggaacagttctctctgcagggaaggaaagaacgccatggggcgattcatcttcgtgagcttcaacttgctggtcgtgttcctctccctaagtggaactctagctgatttggaatgtccctccggttggtcttcctatgatcggtattgctacaagcccttcaaacaagagatgacctgggccgatgcagagaggttctgctcggagcaggcgaagggcgggcatctcctctctgtcgaaaccgccctagaagcatcctttgtggacaatgtgctctatgcgaacaaagagtacctcacacgttatatctggattggactgagggttcaaaacaaaggacagccatgctccagcatcagttatgagaacctggttgacccatttgaatgttttatggtgagcagagacacaaggcttcgtgagtggtttaaagttgactgtgaacaacaacattctttcatatgcaagttcacgcgaccacgttaagatccggctgtgtgaagtctggagaagcaaggaagccccccacctctccccaccccccacctgccgcaatctctg\n") (fasta-string gsn))
      (is (= :iupacNucleicAcids (alphabet gsn)))
      (is (= (created gsn) "08-JUL-2013"))
      (is (= (modified gsn) "08-JUL-2013"))
      (is (= 1 (version gsn)))
      (is (= 35024 (taxid gsn)))
      (is (= "Eukaryota" (first (taxonomy gsn))))
      (is (= "Crotalus horridus" (org-scientific-name gsn)))
      (is (= "mRNA" (moltype gsn)))
      (is (= "GAAZ01003035" (gb-locus gsn)))
      (is (= "GAAZ01003035"
             (with-open [r (bs-reader
                            (init-genbank-connection '("GAAZ01003035")
                                                     :nucest :xml))]
               (accession (first (biosequence-seq r))))))
      (is (= "source" (-> (feature-seq gsn) first feature-type)))
      (is (= "1..628" (-> (feature-seq gsn) first feature-location)))
      (is (= 1 (-> (feature-seq gsn) first interval-seq first start)))
      (is (= "organism" (-> (feature-seq gsn) first qualifier-seq first qualifier-name)))
      (is (= "Crotalus horridus" (-> (feature-seq gsn) first qualifier-seq first qualifier-value)))))
  (testing "Genbank indexing"
    (let [i (index-biosequence-file
             (init-genbank-file
              (io/resource "test-files/nucleotide-gb.xml")))]
      (try
        (with-open [r (bs-reader i)]
          (is (= "KE373594"
                 (accession (first (biosequence-seq r)))))
          (is (= "GAAZ01003035"
                 (accession (get-biosequence r "GAAZ01003035")))))
        (let [ii (load-biosequence-index
                  (fs/absolute-path
                   (io/resource "test-files/nucleotide-gb.xml")))]
          (with-open [rr (bs-reader ii)]
            (is (= "KE373594"
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
             (def-line fs)))
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

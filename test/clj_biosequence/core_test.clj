(ns clj-biosequence.core-test
  (:require [clojure.java.io :as io]
            [fs.core :as fs])
  (:use clojure.test
        clj-biosequence.core
        clj-biosequence.uniprot
        clj-biosequence.signalp
        clj-biosequence.genbank
        clj-biosequence.blast
        clj-biosequence.fastq))

(def fasta-nuc (with-open [f (bs-reader (init-fasta-file
                                         (io/resource "test-files/nuc-sequence.fasta")
                                         :iupacNucleicAcids))]
                 (first (biosequence-seq f))))

(def un-seq (with-open [f (bs-reader
                           (init-uniprotxml-file
                            (io/resource "test-files/uniprot-s-mansoni-20121217.xml")))]
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
    (is (= "accession" (with-open [s (bs-reader (init-fasta-string
                                                 ">accession\n desc\n sequence"
                                                 :iupacNucleicAcids))]
                         (accession (first (biosequence-seq s)))))))
  (testing "Fasta indexing"
    (try
      (let [i (index-file (init-fasta-file (io/resource "test-files/nuc-sequence.fasta")
                                           :iupacNucleicAcids))]
        (is (= "gi|116025203|gb|EG339215.1|EG339215"
               (accession (first (biosequence-seq i)))))
        (is (= "gi|114314166|gb|EE741316.1|EE741316"
               (accession (get-biosequence i "gi|114314166|gb|EE741316.1|EE741316")))))
      (let [i (load-indexed-file (io/resource "test-files/nuc-sequence.fasta"))]
        (is (= "gi|116025203|gb|EG339215.1|EG339215"
               (accession (first (biosequence-seq i)))))
        (is (= "gi|114314166|gb|EE741316.1|EE741316"
               (accession (get-biosequence i "gi|114314166|gb|EE741316.1|EE741316")))))
      (finally
        (fs/delete (str (fs/absolute-path (io/resource "test-files/nuc-sequence.fasta"))
                        ".idx"))
        (fs/delete (str (fs/absolute-path (io/resource "test-files/nuc-sequence.fasta"))
                        ".bin"))))))

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
      (let [c (first (citations un-seq))]
        (is (= "journal article" (ref-type c)))
        (is (= "The genome of the blood fluke Schistosoma mansoni." (title c)))
        (is (= "Nature" (journal c)))
        (is (= 2009 (year c)))
        (is (= 460 (volume c)))
        (is (= 352 (pstart c)))
        (is (= 358 (pend c)))
        (is (= "Berriman M." (first (authors c))))))
    (is (= "function" (:type (:attrs (first (comments un-seq))))))
    (is (= "chain" (feature-type (first (feature-seq un-seq)))))
    (is (= 1 (-> (feature-seq un-seq) first interval-seq first start))))
  (testing "Uniprot indexing"
    (try
      (let [i (index-file (init-uniprotxml-file
                           (io/resource "test-files/uniprot-s-mansoni-20121217.xml")))]
        (is (= "C4PYP8"
               (accession (first (biosequence-seq i)))))
        (is (= "P35661"
               (accession (get-biosequence i "P35661")))))
      (let [i (load-indexed-file (io/resource "test-files/uniprot-s-mansoni-20121217.xml"))]
        (is (= "C4PYP8"
               (accession (first (biosequence-seq i)))))
        (is (= "P35661"
               (accession (get-biosequence i "P35661")))))
      (finally
        (fs/delete (str (fs/absolute-path
                         (io/resource "test-files/uniprot-s-mansoni-20121217.xml"))
                        ".idx"))
        (fs/delete (str (fs/absolute-path
                         (io/resource "test-files/uniprot-s-mansoni-20121217.xml"))
                        ".bin"))))))

(deftest signalp-test
  (testing "Signalp"
    (let [bs (with-open [r (bs-reader (init-fasta-file
                                       (io/resource "test-files/toxins.fasta")
                                       :iupacAminoAcids))]
               (doall (take 20 (biosequence-seq r))))]
      (is (= "sp|C1IC47|3FN3_WALAE" (accession (first (filter-signalp bs))))))))

(deftest genbank-test
  (testing "Genbank"
    (let [gsn (with-open [r (bs-reader (init-genbank-file
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
      (is (= "GAAZ01003035" (with-open [r (bs-reader
                                           (init-genbank-connection '("GAAZ01003035")
                                                                    :nucest :xml))]
                              (accession (first (biosequence-seq r))))))
      (is (= "source" (-> (feature-seq gsn) first feature-type)))
      (is (= "1..628" (-> (feature-seq gsn) first feature-location)))
      (is (= 1 (-> (feature-seq gsn) first interval-seq first start)))
      (is (= "organism" (-> (feature-seq gsn) first qualifier-seq first qualifier-name)))
      (is (= "Crotalus horridus" (-> (feature-seq gsn) first qualifier-seq first qualifier-value)))))
  (testing "Genbank indexing"
    (try
      (let [i (index-file (init-genbank-file
                           (io/resource "test-files/nucleotide-gb.xml")))]
        (is (= "KE373594"
               (accession (first (biosequence-seq i)))))
        (is (= "GAAZ01003035"
               (accession (get-biosequence i "GAAZ01003035")))))
      (let [i (load-indexed-file (io/resource "test-files/nucleotide-gb.xml"))]
        (is (= "KE373594"
               (accession (first (biosequence-seq i)))))
        (is (= "GAAZ01003035"
               (accession (get-biosequence i "GAAZ01003035")))))
      (finally
        (fs/delete (str (fs/absolute-path
                         (io/resource "test-files/nucleotide-gb.xml"))
                        ".idx"))
        (fs/delete (str (fs/absolute-path
                         (io/resource "test-files/nucleotide-gb.xml"))
                        ".bin"))))))

(deftest blast-test
  (testing "Blast"
    (let [toxindb (init-blast-db (io/resource "test-files/toxins.fasta")
                                 :iupacAminoAcids)
          toxins (init-fasta-file (io/resource "test-files/toxins.fasta") :iupacAminoAcids)
          tox-bl (with-open [r (bs-reader toxins)]
                   (blast (take 20 (biosequence-seq r))
                          "blastp"
                          toxindb
                          "/tmp/blast.xml"))]
      (is (= "B3EWT5" (with-open [r (bs-reader tox-bl)]
                        (first
                         (->> (biosequence-seq r)
                              (filter #(>= (-> (hit-seq %) second hit-bit-scores first) 50))
                              (map #(-> (hit-seq %) second hit-accession)))))))
      (testing "Blast indexing"
        (try
          (let [i (index-file tox-bl)]
            (is (= "sp|P84001|29C0_ANCSP" (accession (first (biosequence-seq i)))))
            (is (= "sp|P0CE81|A1HB1_LOXIN"
                   (accession (get-biosequence i "sp|P0CE81|A1HB1_LOXIN"))))
            (is (= 10 (blast-evalue (parameters i)))))
          (let [i (load-indexed-file (bs-path tox-bl))]
            (is (= "sp|P84001|29C0_ANCSP" (accession (first (biosequence-seq i)))))
            (is (= "sp|P0CE81|A1HB1_LOXIN"
                   (accession (get-biosequence i "sp|P0CE81|A1HB1_LOXIN"))))
            (is (= 10 (blast-evalue (parameters i)))))
          (finally
            (fs/delete "/tmp/blast.xml")
            (fs/delete "/tmp/blast.xml.idx")
            (fs/delete "/tmp/blast.xml.bin")))))))

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
    (try
      (let [i (index-file (init-fastq-file
                           (io/resource "test-files/fastq-test.fastq")))]
        (println (accession (second (biosequence-seq i))))
        (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:6408:1991 1:N:0:ACAGTG"
               (accession (first (biosequence-seq i)))))
        (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG"
               (accession (get-biosequence i "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG")))))
      (let [i (load-indexed-file (io/resource "test-files/fastq-test.fastq"))]
        (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:6408:1991 1:N:0:ACAGTG"
               (accession (first (biosequence-seq i)))))
        (is (= "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG"
               (accession (get-biosequence i "@HWI-ST1213:141:C17PWACXX:6:1101:5381:1994 1:N:0:ACAGTG")))))
      (finally
        (fs/delete (str (fs/absolute-path
                         (io/resource "test-files/fastq-test.fastq"))
                        ".idx"))
        (fs/delete (str (fs/absolute-path
                         (io/resource "test-files/fastq-test.fastq"))
                        ".bin"))))))

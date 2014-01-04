(ns clj-biosequence.core-test
  (:require [clojure.java.io :as io]
            [fs.core :as fs])
  (:use clojure.test
        clj-biosequence.core
        clj-biosequence.store
        clj-biosequence.uniprot
        clj-biosequence.signalp))

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
    (is (= 143 (get (residue-frequencies fasta-nuc) \G)))
    (is (= [\G \T \A \C \A \A \A] (bs-seq (sub-bioseq fasta-nuc 0 7))))
    (is (= [\G \T \A] (first (partition-bioseq fasta-nuc))))
    (is (= [\T \A \A \T] (bs-seq (sub-bioseq (reverse-comp fasta-nuc) 0 4))))
    (is (= [\A \T \T \A] (bs-seq (sub-bioseq (reverse-seq fasta-nuc) 0 4))))
    (is (= "accession" (with-open [s (bs-reader (init-fasta-string
                                                 ">accession\n desc\n sequence"
                                                 :iupacNucleicAcids))]
                         (accession (first (biosequence-seq s))))))
    (is (= nil
           (with-open [f (bs-reader (init-fasta-file
                                     (io/resource "test-files/nuc-sequence.fasta")
                                     :iupacNucleicAcids))]
             (parameters f))))))

(deftest uniprot
  (testing "Uniprot"
    (is (= "C4PYP8" (accession un-seq)))
    (is (= (list "C4PYP8" "G4VQR6") (accessions un-seq)))
    (is (= (vec "Anamorsin") (subvec (vec (def-line un-seq)) 0 9)))
    (is (= [\M \E \Q \C \V \A] (subvec (bs-seq un-seq) 0 6)))
    (is (= true (protein? un-seq)))
    (is (= :iupacAminoAcids (alphabet un-seq)))
    (is (= java.lang.String (class (bioseq->string un-seq))))
    (is (= 16 (get (residue-frequencies un-seq) \A)))
    (is (= [\M \E \Q \C \V \A] (bs-seq (sub-bioseq un-seq 0 6))))
    (is (= [\M \E \Q] (first (partition-bioseq un-seq))))
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
    (is (= '("Schistosoma mansoni") (organism un-seq "scientific")))
    (is (= "Eukaryota" (first (lineage un-seq))))
    (is (= "DRE2_SCHMA" (prot-name un-seq)))
    (is (= '("Anamorsin homolog") (nomenclature un-seq "recommendedName")))
    (is (= '("29662") (sequence-info un-seq "mass")))
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
    (is (= "function" (first (comments un-seq))))
    (is (= (vec "May be required") (subvec
                                    (vec (:text (first (comment-value un-seq "function"))))
                                    0 15)))))

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
                (first (biosequence-seq r)))]
      (is (= (created gsn) "08-JUL-2013"))
      )))

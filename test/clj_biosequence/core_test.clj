(ns clj-biosequence.core-test
  (:require [clojure.java.io :as io])
  (:use clojure.test
        clj-biosequence.core
        clj-biosequence.store
        clj-biosequence.uniprot))

(deftest fasta
  (testing "Fasta:"
    (let [nuc (init-fasta-sequence "gi|196122522|gb|EU910552.1|" "Viola baoshanensis cyclotide precursor 7c mRNA, complete cds" :iupacNucleicAcids "ATGAAGATGTTTGTTGCCCTTGTGCTCTTTGCAGCCTTCGCTCTCCCCGCTGCCTTTGCAGCTCAGCAAGATGTCATCACTCTCCGAGCTTATGAGGAGCTTCTCAAGAATGGGGCTGCTAATGGAATGACCAAAACTGTCATCTCAAGCCCTGTTCTTGAAGAGGCTCTCGTCTCTTACTCCAAGAACAAGCTCGGCGGTCTTCCTGTCTGCGGAGAGACTTGCGTCGGTGGAACATGCAACACCCCTGGGTGCGGCTGCAGCTGGCCAGTCTGCACCAGAAACTCTCTTGAAAGCACCAAATCTGCAAACCCTCTCCTTGAAGAGGCACTCACCACGTTCGCCAAGAAAGGTCTTGGCGGTCTTCCCGTCTGCGGTGAGACTTGCGTCGGCGGAACATGCAACACCCCTGGGTGCACTTGCAGCTGGCCTGTTTGCACCAGAAATGCTCTTGAGACTCAGAAACCCAACCATTTGCTCGAAGAGGCACTCATTGCATTCGCCAAGAAAGGTAACCTCGGTGGTCTTCCCGTATGCGGCGAAACTTGCGTTGGTGGAACATGCAACACCCCTGGGTGCTCCTGCAGCTGGCCTGTTTGCACCAGAAACTCTCTTGCCATGTAAGACTGAGATTTAGGCATATGGCTGGACAATGCTCGATTTGTTGTTTGTCGTTCTCGTTCTTGTTGTTTGTCGGAGAGGAGGCTATACATGCTGGTGACATATGAAATTACCTAAAATTAAAATGTAAAGCATGCCTCCTTTCCCT")
          prot (init-fasta-sequence "tr|Q93X48|Q93X48_LENER" "Lectin OS=Lens ervoides GN=lectin PE=4 SV=2" :iupacAminoAcids "MASLQTQMISFYLIFLSILLTTIFFFKVNSTETTSFSITKFSPDQQNLIFQGDGYTTKEKLTLTKAVKNTVGRALYSTPIHIWDRDTGNVANFVTSFTFVINAPNSYNVADGFTFFIAPVDTKPQTGGGYLGVFNSKDYDKTSQTVAVEFDTFYNAAWDPSNKDRHIGIDVNSIKSVSTKSWNLQNGERANVVIAFNAATNVLTVTLTYPNSLEEENVTSYTLNEVVPMKDVLPEWVRIGFSATTGAEFAAHEVLSWSFHSELGGTSSSKQAADA")]
      (testing "Translation:"
        (is (= "MKMFVALVLFAAFALPAAFAAQQDVITLRAYEELLKNGAANGMTKTVISSPVLEEALVSYSKNKLGGLPVCGETCVGGTCNTPGCGCSWPVCTRNSLESTKSANPLLEEALTTFAKKGLGGLPVCGETCVGGTCNTPGCTCSWPVCTRNALETQKPNHLLEEALIAFAKKGNLGGLPVCGETCVGGTCNTPGCSCSWPVCTRNSLAM*D*DLGIWLDNARFVVCRSRSCCLSERRLYMLVTYEIT*N*NVKHASFPX"
               (bioseq->string (translate nuc 1))))
        (is (= "RERRHALHFNFR*FHMSPACIASSPTNNKNENDKQQIEHCPAICLNLSLTWQESFWCKQASCRSTQGCCMFHQRKFRRIREDHRGYLSWRMQ*VPLRANGWVSESQEHFWCKQASCKCTQGCCMFRRRKSHRRREDRQDLSWRTW*VPLQGEGLQIWCFQESFWCRLASCSRTQGCCMFHRRKSLRRQEDRRACSWSKRREPLQEQGLR*QFWSFH*QPHS*EAPHKLGE**HLAELQRQRGERRLQRAQGQQTSSX"
               (bioseq->string (nth (six-frame-translation nuc) 3)))))
      (testing "Fasta string"
        (is (= ">gi|196122522|gb|EU910552.1| Viola baoshanensis cyclotide precursor 7c mRNA, complete cds\nATGAAGATGTTTGTTGCCCTTGTGCTCTTTGCAGCCTTCGCTCTCCCCGCTGCCTTTGCAGCTCAGCAAGATGTCATCACTCTCCGAGCTTATGAGGAGCTTCTCAAGAATGGGGCTGCTAATGGAATGACCAAAACTGTCATCTCAAGCCCTGTTCTTGAAGAGGCTCTCGTCTCTTACTCCAAGAACAAGCTCGGCGGTCTTCCTGTCTGCGGAGAGACTTGCGTCGGTGGAACATGCAACACCCCTGGGTGCGGCTGCAGCTGGCCAGTCTGCACCAGAAACTCTCTTGAAAGCACCAAATCTGCAAACCCTCTCCTTGAAGAGGCACTCACCACGTTCGCCAAGAAAGGTCTTGGCGGTCTTCCCGTCTGCGGTGAGACTTGCGTCGGCGGAACATGCAACACCCCTGGGTGCACTTGCAGCTGGCCTGTTTGCACCAGAAATGCTCTTGAGACTCAGAAACCCAACCATTTGCTCGAAGAGGCACTCATTGCATTCGCCAAGAAAGGTAACCTCGGTGGTCTTCCCGTATGCGGCGAAACTTGCGTTGGTGGAACATGCAACACCCCTGGGTGCTCCTGCAGCTGGCCTGTTTGCACCAGAAACTCTCTTGCCATGTAAGACTGAGATTTAGGCATATGGCTGGACAATGCTCGATTTGTTGTTTGTCGTTCTCGTTCTTGTTGTTTGTCGGAGAGGAGGCTATACATGCTGGTGACATATGAAATTACCTAAAATTAAAATGTAAAGCATGCCTCCTTTCCCT\n"
               (fasta-string nuc))))
      (testing "Core biosequence functions"
        (is (= [\A \G \G \G \A \A \A \G \G \A \G]
               (bs-seq (sub-bioseq (reverse-comp nuc) 0 11))))
        (is (= [\M \K \M \F \V \A \L \V \L \F]
               (bs-seq (sub-bioseq (translate nuc 1) 0 10))))
        (is (= [\G \K \E \A \C \F \T \F]
               (bs-seq (sub-bioseq (nth (six-frame-translation nuc) 5) 0 8)))))
      (testing "File access"
        (let [f (io/resource "test-files/nuc-sequence.fasta")
              ff (init-fasta-file f :iupacNucleicAcids)
              fs (init-fasta-string (slurp f) :iupacNucleicAcids)]
          (with-open [r (bs-reader ff)]
            (is (= [\G \T \A \C \A \A \A]
                   (bs-seq (sub-bioseq (first (biosequence-seq r)) 0 7)))))
          (with-open [r (bs-reader fs)]
            (is (= [\G \T \A \C \A \A \A]
                   (bs-seq (sub-bioseq (first (biosequence-seq r)) 0 7)))))
          (testing "Store"
            (let [p (init-project "clj-test")
                  c (mongo-connect)
                  i (mongo-save-file ff p "testing")]
              (is (= 6 (count (collection-seq i))))
              (is (= "gi|116025203|gb|EG339215.1|EG339215"
                     (accession (first (collection-seq i)))))
              (drop-project p)
              (mongo-disconnect)))))
      (testing "Mapping"
        (is (= "B5B3Z7" ((id-convert '("196122523") "P_GI" "ACC"
                                     "jason.mulvenna@gmail.com")
                         "196122523")))))))

(deftest uniprot
  (testing "Uniprot"
    (let [uf (io/resource "test-files/uniprot-s-mansoni-20121217.xml")
          uff (init-uniprotxml-file uf)
          up (with-open [r (bs-reader uff)]
               (first (biosequence-seq r)))
          us (init-uniprot-string (slurp uf))
          uc (init-uniprot-connection '("P56871" "P56879" "P84641" "P84642")
                                      :xml "jason.mulvenna@gmail.com")]
      (is (= [\X \M \E \Q \C \V] (bs-seq (sub-bioseq up 0 6))))
      (is (= "C4PYP8" (accession up)))
      (is (= "P35661" (with-open [r (bs-reader us)]
                        (accession (second (biosequence-seq r))))))
      (is (= "P56871"
             (first (wget-uniprot-search "cyclotide" "jason.mulvenna@gmail.com"))))
      (is (= "P84641"
             (with-open [r (bs-reader uc)]
               (accession (nth (biosequence-seq r) 2))))))))

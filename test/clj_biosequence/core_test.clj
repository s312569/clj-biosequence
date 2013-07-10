(ns clj-biosequence.core-test
  (:use clojure.test
        clj-biosequence.core
        clj-biosequence.blast
        clj-biosequence.uniprot)
  (:require [fs.core :as fs]
            [clj-biosequence.genbank :as gb]))

(deftest core
  (testing "Core:"
    (testing "Sequence initialisation"
      (let [nuc (init-fasta-sequence "gi|196122522|gb|EU910552.1|" "Viola baoshanensis cyclotide precursor 7c mRNA, complete cds" :nucleotide "ATGAAGATGTTTGTTGCCCTTGTGCTCTTTGCAGCCTTCGCTCTCCCCGCTGCCTTTGCAGCTCAGCAAGATGTCATCACTCTCCGAGCTTATGAGGAGCTTCTCAAGAATGGGGCTGCTAATGGAATGACCAAAACTGTCATCTCAAGCCCTGTTCTTGAAGAGGCTCTCGTCTCTTACTCCAAGAACAAGCTCGGCGGTCTTCCTGTCTGCGGAGAGACTTGCGTCGGTGGAACATGCAACACCCCTGGGTGCGGCTGCAGCTGGCCAGTCTGCACCAGAAACTCTCTTGAAAGCACCAAATCTGCAAACCCTCTCCTTGAAGAGGCACTCACCACGTTCGCCAAGAAAGGTCTTGGCGGTCTTCCCGTCTGCGGTGAGACTTGCGTCGGCGGAACATGCAACACCCCTGGGTGCACTTGCAGCTGGCCTGTTTGCACCAGAAATGCTCTTGAGACTCAGAAACCCAACCATTTGCTCGAAGAGGCACTCATTGCATTCGCCAAGAAAGGTAACCTCGGTGGTCTTCCCGTATGCGGCGAAACTTGCGTTGGTGGAACATGCAACACCCCTGGGTGCTCCTGCAGCTGGCCTGTTTGCACCAGAAACTCTCTTGCCATGTAAGACTGAGATTTAGGCATATGGCTGGACAATGCTCGATTTGTTGTTTGTCGTTCTCGTTCTTGTTGTTTGTCGGAGAGGAGGCTATACATGCTGGTGACATATGAAATTACCTAAAATTAAAATGTAAAGCATGCCTCCTTTCCCT")
            prot (init-fasta-sequence "tr|Q93X48|Q93X48_LENER" "Lectin OS=Lens ervoides GN=lectin PE=4 SV=2" :protein "MASLQTQMISFYLIFLSILLTTIFFFKVNSTETTSFSITKFSPDQQNLIFQGDGYTTKEKLTLTKAVKNTVGRALYSTPIHIWDRDTGNVANFVTSFTFVINAPNSYNVADGFTFFIAPVDTKPQTGGGYLGVFNSKDYDKTSQTVAVEFDTFYNAAWDPSNKDRHIGIDVNSIKSVSTKSWNLQNGERANVVIAFNAATNVLTVTLTYPNSLEEENVTSYTLNEVVPMKDVLPEWVRIGFSATTGAEFAAHEVLSWSFHSELGGTSSSKQAADA")
            ffile (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa" :protein)
            store (index-fasta-file ffile)
            sseq (with-biosequences [l store]
                   (first l))]
       (is (= clj_biosequence.core.fastaSequence
              (class nuc)))
       (is (= 6 (count (six-frame-translation nuc (codon-tables :standard)))))
       (is (= clj_biosequence.core.fastaSequence
              (class prot)))
       (is (= clj_biosequence.core.fastaFile
              (class ffile)))
       (is (= "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa"
              (file-path ffile)))
       (is (= nil (with-biosequences-in-file [l ffile]
                    (println (first l))
                    (println (last l)))))
       (is (= nil
              (fasta-to-stream nuc *out*)))
       (is (= clj_biosequence.core.fastaStore
              (class store)))
       (is (= nil (with-biosequences [l store]
                    (println (first l))
                    (println (last l)))))
       (is (= clj_biosequence.core.fastaSequence (class sseq)))
       (is (= clj_biosequence.core.fastaSequence
              (class (with-connection-to-store [store]
                       (update-object (assoc sseq :sequence "jason"))))))
       (is (= "jason" (with-connection-to-store [store]
                        (sequence-string (get-biosequence (accession sseq))))))
       (is (= clj_biosequence.core.fastaStore
              (class (load-fasta-store "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa.ind" :protein))))
       (fs/delete-dir (fs/parent (:file store)))))))

(deftest blast
  (testing "Blast"
    (testing "Blastp"
      (let [seq (init-fasta-sequence "tr|Q93X48|Q93X48_LENER" "Lectin OS=Lens ervoides GN=lectin PE=4 SV=2" :protein "MASLQTQMISFYLIFLSILLTTIFFFKVNSTETTSFSITKFSPDQQNLIFQGDGYTTKEKLTLTKAVKNTVGRALYSTPIHIWDRDTGNVANFVTSFTFVINAPNSYNVADGFTFFIAPVDTKPQTGGGYLGVFNSKDYDKTSQTVAVEFDTFYNAAWDPSNKDRHIGIDVNSIKSVSTKSWNLQNGERANVVIAFNAATNVLTVTLTYPNSLEEENVTSYTLNEVVPMKDVLPEWVRIGFSATTGAEFAAHEVLSWSFHSELGGTSSSKQAADA")
            bdb (init-blast-db
                 "/Users/jason/Dropbox/clj-biosequence/test-files/toxins.fasta"
                 :protein)
            file (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa" :protein)
            fbl (blast-file file bdb "blastp")
            store (index-fasta-file file)
            sbl (blast-store store bdb "blastp")
            bl (blast-biosequence seq bdb "blastp")]
        (is (= clj_biosequence.blast.blastDB (class bdb)))
        (is (= clj_biosequence.core.fastaFile (class file)))
        (is (= clj_biosequence.core.fastaStore (class store)))
        (is (= clj_biosequence.blast.blastSearch (class fbl)))
        (is (= clj_biosequence.blast.blastIteration (class bl)))
        (is (= 24.6386102277464 (hit-bit-score (first (hit-seq bl)))))
        (fs/delete-dir (fs/parent (:file store)))
        (fs/delete (:src fbl))
        (fs/delete (:src sbl))))))

(deftest uniprot
  (testing "Uniprot"
    (let [ufile (init-uniprotxml-file "/Users/jason/Dropbox/clj-biosequence/test-files/uniprot-s-mansoni-20121217.xml")
          useq (with-biosequences-in-file [l ufile]
                 (first l))
          ustore (index-uniprotxml-file ufile)
          ustore2 (load-uniprot-store (fs/parent (:file ustore)))
          bdb (init-blast-db
                 "/Users/jason/Dropbox/clj-biosequence/test-files/toxins.fasta"
                 :protein)
          ubl (blast-biosequence useq bdb "blastp")]
      (is (= clj_biosequence.uniprot.uniprotXmlFile
             (class ufile)))
      (is (= clj_biosequence.uniprot.uniprotProtein
             (class useq)))
      (is (= '("C4PYP8" "G4VQR6")
             (accessions useq)))
      (is (= nil (println (fasta-string useq))))
      (is (= nil (fasta-to-stream useq *out*)))
      (is (= "Schistosoma mansoni" (org-scientific-name useq)))
      (is (= "2010-03-23" (created useq)))
      (is (= "2012-10-03" (modified useq)))
      (is (= 20 (version useq)))
      (is (= "Swiss-Prot" (database useq)))
      (is (= '("Eukaryota" "Metazoa" "Platyhelminthes" "Trematoda" "Digenea" "Strigeidida" "Schistosomatoidea" "Schistosomatidae" "Schistosoma") (taxonomy useq)))
      (is (= 6183 (taxid useq)))
      (is (= {:ncbi-taxid 6183, :abbreviation (), :synonym (), :full (), :common '("Blood fluke"), :scientific '("Schistosoma mansoni")} (organism-name useq)))
      (is (= "DRE2_SCHMA" (prot-name useq)))
      (is (= {:mass 29662, :checksum "11A38702F71B84D8", :modified "2010-03-23", :version 2, :amino-acids "MEQCVADCLNSDDCVMIVWSGEVQEDVMRGLQVAVSTYVKKLQFENLEKFVDSSAVDSQLHECSVILCGWPNSISVNILKLGLLSNLLSCLRPGGRFFGRDLITGDWDSLKKNLTLSGYINPYQLSCENHLIFSASVPSNYTQGSSVKLPWANSDVEAAWENVDNSSDANGNIINTNTLLQKSDLKTPLSVCGKEAATDSVGKKKRACKNCTCGLAEIEAAEEDKSDVPISSCGNCYLGDAFRCSTCPYRGLPPFKPGERILIPDDVLRADL"} (amino-acids useq)))
      (is (= {:alternate '({:fullname "Fe-S cluster assembly protein DRE2 homolog", :shortname nil, :ecnumber nil}), :recommended '({:fullname "Anamorsin homolog", :shortname nil, :ecnumber nil})} (nomenclature useq)))
      (is (= '({:name "Smp_207000", :type "ORF"}) (gene useq)))
      (is (= {:mitochondrion (), :chloroplast (), :apicoplast (), :nucleomorph (), :plasmid (), :hydrogenosome (), :cyanelle (), :plastid (), :organellar-chromatophore (), :non-photosynthetic-plastid ()} (gene-location useq)))
      (is (= '({:country nil, :last "358", :date "2009", :pubmed "19606141", :institute nil, :name "Nature", :first "352", :title "The genome of the blood fluke Schistosoma mansoni.", :city nil, :scope ("NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]"), :type "journal article", :consortium (), :number nil, :authors ("Berriman M." "Haas B.J." "LoVerde P.T." "Wilson R.A." "Dillon G.P." "Cerqueira G.C." "Mashiyama S.T." "Al-Lazikani B." "Andrade L.F." "Ashton P.D." "Aslett M.A." "Bartholomeu D.C." "Blandin G." "Caffrey C.R." "Coghlan A." "Coulson R." "Day T.A." "Delcher A." "DeMarco R." "Djikeng A." "Eyre T." "Gamble J.A." "Ghedin E." "Gu Y." "Hertz-Fowler C." "Hirai H." "Hirai Y." "Houston R." "Ivens A." "Johnston D.A." "Lacerda D." "Macedo C.D." "McVeigh P." "Ning Z." "Oliveira G." "Overington J.P." "Parkhill J." "Pertea M." "Pierce R.J." "Protasio A.V." "Quail M.A." "Rajandream M.A." "Rogers J." "Sajid M." "Salzberg S.L." "Stanke M." "Tivey A.R." "White O." "Williams D.L." "Wortman J." "Wu W." "Zamanian M." "Zerlotini A." "Fraser-Liggett C.M." "Barrell B.G." "El-Sayed N.M."), :source {:tissue (), :transposon (), :plasmid (), :strain ("Puerto Rican")}, :editors (), :publisher nil, :volume "460", :db nil}) (citation useq)))
      (is (= '({:comments nil, :orientation nil, :topology nil, :location "Cytoplasm", :evidence "by similarity"}) (subcellular-location useq)))
      (is (= () (alternative-products useq)))
      (is (= () (mass-spectroscopy useq)))
      (is (= {:similarity "Belongs to the anamorsin family.", :function "May be required for the maturation of extramitochondrial Fe/S proteins (By similarity). Has anti-apoptotic effects in the cell (By similarity)."} (comments useq)))
      (is (= '({:data {"molecule type" "Genomic_DNA", "protein sequence ID" "CCD81528.1"}, :db "EMBL", :id "HE601631"} {:data {"evidence" "IEA:UniProtKB-SubCell", "term" "C:cytoplasm"}, :db "GO", :id "GO:0005737"} {:data {"evidence" "IEA:UniProtKB-KW", "term" "P:apoptotic process"}, :db "GO", :id "GO:0006915"} {:data {"entry name" "CIAPIN1"}, :db "InterPro", :id "IPR007785"} {:data {"match status" "1", "entry name" "PTHR13273"}, :db "PANTHER", :id "PTHR13273"} {:data {"match status" "1", "entry name" "CIAPIN1"}, :db "Pfam", :id "PF05093"}) (db-references useq)))
      (is (= '("inferred from homology") (existence useq)))
      (is (= {:KW-1185 "Reference proteome", :KW-0963 "Cytoplasm", :KW-0181 "Complete proteome", :KW-0053 "Apoptosis"} (keywords useq)))
      (is (= '(({:id "PRO_0000392331", :description "Anamorsin homolog", :type "chain", :begin 1, :end 272, :position nil}) (features useq))))
      (is (= nil
             (println (with-biosequences-in-file [l ufile]
                        (println (fasta-string (first l)))
                        (println (fasta-string (second l)))))))
      (is (= "/Users/jason/Dropbox/clj-biosequence/test-files/uniprot-s-mansoni-20121217.xml" (file-path ufile)))
      (is (= clj_biosequence.uniprot.uniprotStore
             (class ustore)))
      (is (= nil (with-biosequences [l ustore]
                   (println (fasta-string (first l)))
                   (println (fasta-string (last l))))))
      (is (= clj_biosequence.uniprot.uniprotStore
             (class ustore2)))
      (fs/delete-dir (fs/parent (:file ustore))))))

(deftest genbank
  (testing "Genbank"
    (let [pgfile (gb/init-genbank-file "/Users/jason/Dropbox/clj-biosequence/test-files/protein-gb.xml")
          pgseq (with-biosequences-in-file [l pgfile]
                  (first l))
          pgstore (gb/index-genbank-file pgfile)
          pgstorep (with-connection-to-store [pgstore]
                     (get-biosequence "B3EWH7"))
          pgfeat  (nth (gb/feature-seq pgseq) 2)
          pgint (first (gb/interval-seq pgfeat))]
      (is (= clj_biosequence.genbank.genbankFile (class pgfile)))
      (is (= clj_biosequence.genbank.genbankSequence (class pgseq)))
      (is (= clj_biosequence.genbank.genbankStore (class pgstore)))
      (is (= '("sp|B3EWH7.1|CYCM_PETHY" "gi|408407585")
             (accessions pgstorep)))
      (is (= 8 (with-biosequences-in-file [l pgfile]
                 (count l))))
      (is (= "B3EWH7" (accession pgseq)))
      (is (= '("sp|B3EWH7.1|CYCM_PETHY" "gi|408407585") (accessions pgseq)))
      (is (= "03-OCT-2012" (created pgseq)))
      (is (= "06-MAR-2013" (modified pgseq)))
      (is (= 1 (version pgseq)))
      (is (= 4102 (taxid pgseq)))
      (is (= '("Eukaryota" " Viridiplantae" " Streptophyta" " Embryophyta" " Tracheophyta" " Spermatophyta" " Magnoliophyta" " eudicotyledons" " core eudicotyledons" " asterids" " lamiids" " Solanales" " Solanaceae" " Petunioideae" " Petunia")
             (taxonomy pgseq)))
      (is (= "Petunia x hybrida" (org-scientific-name pgseq)))
      (is (= true (protein? pgseq)))
      (is (= ">gb|B3EWH7|sp|B3EWH7.1|CYCM_PETHYgi|408407585| RecName: Full=Acyclotide phyb-M\nQSISCAESCVWIPCATSLIGCSCVNSRCIYSK"
             (fasta-string pgseq)))
      (is (= nil (fasta-to-stream pgseq *out*)))
      (is (= "CYCM_PETHY" (gb/gb-locus pgseq)))
      (is (= 12 (count (gb/feature-seq pgseq))))
      (is (= "Region" (gb/feature-type pgfeat)))
      (is (= #clj_biosequence.core.fastaSequence{:accession "B3EWH7", :description "RecName: Full=Acyclotide phyb-M - Feature: Region - [1-32]", :type :protein, :sequence "QSISCAESCVWIPCATSLIGCSCVNSRCIYSK"} (gb/get-feature-sequence pgfeat pgseq)))
      (is (= 1 (count (gb/interval-seq pgfeat))))
      (is (= clj_biosequence.genbank.genbankInterval (class pgint)))
      (is (= 1 (gb/start pgint)))
      (is (= 32 (gb/end pgint)))
      (is (= false (gb/comp? pgint)))
      (is (= #clj_biosequence.core.fastaSequence{:accession "B3EWH7", :description "RecName: Full=Acyclotide phyb-M [1-32]", :type :protein, :sequence "QSISCAESCVWIPCATSLIGCSCVNSRCIYSK"} (gb/get-interval-sequence pgint pgseq)))
      (is (= "152031586" (first (gb/wget-genbank-search "cyclotide" :protein))))
      (is (= #clj_biosequence.core.fastaSequence{:accession "gi|152031586|sp|P58440.2|CYO8_VIOOD", :description "RecName: Full=Cycloviolacin-O8; AltName: Full=Cyclotide c1; Flags: Precursor", :type :protein, :sequence "MEMKNMVVGLFLIAAFALPALATNFEKDFITHETVQAILKKVGPSSNGMLDEQTISALTGKTIISNPVLEEALLTHSNSINALGGTLPCGESCVWIPCISSVVGCSCKSKVCYKNSLA"}
             (first (gb/wget-genbank-sequence "152031586" :protein :fasta))))
      (is (= '("sp|P58440.2|CYO8_VIOOD" "gi|152031586")
             (accessions
              (first
               (gb/wget-genbank-sequence "152031586" :protein :xml)))))
      (fs/delete-dir (fs/parent (:file pgstore))))))

(deftest genbank-nuc
  (testing "Genbank nucleotide"
    (let [ngfile (gb/init-genbank-file "/Users/jason/Dropbox/clj-biosequence/test-files/nucleotide-gb.xml")
          ngseq (with-biosequences-in-file [l ngfile]
                  (second l))
          ngstore (gb/index-genbank-file ngfile)
          ngstoren (with-connection-to-store [ngstore]
                     (get-biosequence "GAAZ01003035"))
          ngfeat  (nth (gb/feature-seq ngseq) 1)
          ngint (first (gb/interval-seq ngfeat))]
      (is (= clj_biosequence.genbank.genbankFile (class ngfile)))
      (is (= clj_biosequence.genbank.genbankSequence (class ngseq)))
      (is (= clj_biosequence.genbank.genbankStore (class ngstore)))
      (is (= '("gb|GAAZ01003035.1|" "gnl|TSA:GAAZ01|Chorr_CTL-10" "gi|521752463")
             (accessions ngstoren)))
      (is (= 8 (with-biosequences-in-file [l ngfile]
                 (count l))))
      (is (= "GAAZ01003035" (accession ngseq)))
      (is (= '("gb|GAAZ01003035.1|" "gnl|TSA:GAAZ01|Chorr_CTL-10" "gi|521752463")
             (accessions ngseq)))
      (is (= "08-JUL-2013" (created ngseq)))
      (is (= "08-JUL-2013" (modified ngseq)))
      (is (= 1 (version ngseq)))
      (is (= 35024 (taxid ngseq)))
      (is (= '("Eukaryota" " Metazoa" " Chordata" " Craniata" " Vertebrata" " Euteleostomi" " Lepidosauria" " Squamata" " Scleroglossa" " Serpentes" " Colubroidea" " Viperidae" " Crotalinae" " Crotalus")
             (taxonomy ngseq)))
      (is (= "Crotalus horridus" (org-scientific-name ngseq)))
      (is (= false (protein? ngseq)))
      (is (= ">gb|GAAZ01003035|gb|GAAZ01003035.1|gnl|TSA:GAAZ01|Chorr_CTL-10gi|521752463| TSA: Crotalus horridus Chorr_CTL-10 mRNA sequence\nCTCCATTCCCCTCATAGGAGGAAAGCCTGGAGTTGCCTCTGAGCAGACTTGCTACCTGTGGAGGCCGAGGAACAGTTCTCTCTGCAGGGAAGGAAAGAACGCCATGGGGCGATTCATCTTCGTGAGCTTCAACTTGCTGGTCGTGTTCCTCTCCCTAAGTGGAACTCTAGCTGATTTGGAATGTCCCTCCGGTTGGTCTTCCTATGATCGGTATTGCTACAAGCCCTTCAAACAAGAGATGACCTGGGCCGATGCAGAGAGGTTCTGCTCGGAGCAGGCGAAGGGCGGGCATCTCCTCTCTGTCGAAACCGCCCTAGAAGCATCCTTTGTGGACAATGTGCTCTATGCGAACAAAGAGTACCTCACACGTTATATCTGGATTGGACTGAGGGTTCAAAACAAAGGACAGCCATGCTCCAGCATCAGTTATGAGAACCTGGTTGACCCATTTGAATGTTTTATGGTGAGCAGAGACACAAGGCTTCGTGAGTGGTTTAAAGTTGACTGTGAACAACAACATTCTTTCATATGCAAGTTCACGCGACCACGTTAAGATCCGGCTGTGTGAAGTCTGGAGAAGCAAGGAAGCCCCCCACCTCTCCCCACCCCCCACCTGCCGCAATCTCTG"
             (fasta-string ngseq)))
      (is (= nil (fasta-to-stream ngseq *out*)))
      (is (= "GAAZ01003035" (gb/gb-locus ngseq)))
      (is (= 2 (count (gb/feature-seq ngseq))))
      (is (= "CDS" (gb/feature-type ngfeat)))
      (is (= #clj_biosequence.core.fastaSequence{:accession "GAAZ01003035", :description "TSA: Crotalus horridus Chorr_CTL-10 mRNA sequence - Feature: CDS - [104-553]", :type :nucleotide, :sequence "ATGGGGCGATTCATCTTCGTGAGCTTCAACTTGCTGGTCGTGTTCCTCTCCCTAAGTGGAACTCTAGCTGATTTGGAATGTCCCTCCGGTTGGTCTTCCTATGATCGGTATTGCTACAAGCCCTTCAAACAAGAGATGACCTGGGCCGATGCAGAGAGGTTCTGCTCGGAGCAGGCGAAGGGCGGGCATCTCCTCTCTGTCGAAACCGCCCTAGAAGCATCCTTTGTGGACAATGTGCTCTATGCGAACAAAGAGTACCTCACACGTTATATCTGGATTGGACTGAGGGTTCAAAACAAAGGACAGCCATGCTCCAGCATCAGTTATGAGAACCTGGTTGACCCATTTGAATGTTTTATGGTGAGCAGAGACACAAGGCTTCGTGAGTGGTTTAAAGTTGACTGTGAACAACAACATTCTTTCATATGCAAGTTCACGCGACCACGTTAA"}
             (gb/get-feature-sequence ngfeat ngseq)))
      (is (= 1 (count (gb/interval-seq ngfeat))))
      (is (= clj_biosequence.genbank.genbankInterval (class ngint)))
      (is (= 104 (gb/start ngint)))
      (is (= 553 (gb/end ngint)))
      (is (= false (gb/comp? ngint)))
      (is (= #clj_biosequence.core.fastaSequence{:accession "GAAZ01003035", :description "TSA: Crotalus horridus Chorr_CTL-10 mRNA sequence [104-553]", :type :nucleotide, :sequence "ATGGGGCGATTCATCTTCGTGAGCTTCAACTTGCTGGTCGTGTTCCTCTCCCTAAGTGGAACTCTAGCTGATTTGGAATGTCCCTCCGGTTGGTCTTCCTATGATCGGTATTGCTACAAGCCCTTCAAACAAGAGATGACCTGGGCCGATGCAGAGAGGTTCTGCTCGGAGCAGGCGAAGGGCGGGCATCTCCTCTCTGTCGAAACCGCCCTAGAAGCATCCTTTGTGGACAATGTGCTCTATGCGAACAAAGAGTACCTCACACGTTATATCTGGATTGGACTGAGGGTTCAAAACAAAGGACAGCCATGCTCCAGCATCAGTTATGAGAACCTGGTTGACCCATTTGAATGTTTTATGGTGAGCAGAGACACAAGGCTTCGTGAGTGGTTTAAAGTTGACTGTGAACAACAACATTCTTTCATATGCAAGTTCACGCGACCACGTTAA"}
             (gb/get-interval-sequence ngint ngseq)))
      (fs/delete-dir (fs/parent (:file ngstore))))))

(deftest genbank-genome
  (testing "Genbank genome"
    (let [gfile (gb/init-genbank-file "/Users/jason/Dropbox/clj-biosequence/test-files/akata-sequence.xml")
          gseq (with-biosequences-in-file [l gfile]
                 (first l))
          gfeat (first (filter #(= (gb/feature-type %) "CDS")
                               (gb/feature-seq gseq)))]
      (is (= 1 (with-biosequences-in-file [l gfile]
                 (count l))))
      (is (= 402 (count (gb/feature-seq gseq))))
      (let [f (filter #(= (gb/feature-type %)
                          "CDS") (gb/feature-seq gseq))]
        (doseq [feat f]
          (is (= (gb/qualifier-extract feat "translation")
                 (apply str (remove #(= % \*)
                                    (sequence-string
                                     (translate-biosequence
                                      (gb/get-feature-sequence feat gseq)
                                      1))))))))
      (is (= '("MGSLEMVPMGAGPPSPGGDPDGDDGGNNSQYPSASGSSGNTPTPPNDEERESNEEPPPPYEDPYWGNGDRHSDYQPLGTQDQSLYLGLQHDGNDGLPPPPYSPRDDSSQHIYEEAGRGX"
"MNPVCLPVIVAPYLFWLAAIAASCFTASVSTVVTATGLALSLLLLAAVASSYAAAQRKLLTPVTVLTAVVTX"
"FAICLTWRIEDPPFNSLLFALLAAAGGLQGIYX"
"LVMLVLLILAYRRRWRRLTVCGGIMFLACVLVLIVDAVLQLSPLLGAVTVVSMTLLLLAFVLWLSSPGGLGTLGAALLTLAAX"
"LALLASLILGTLNLTTMFLLMLLWTLX"
"VLLICSSCSSCPLTKILLARLFLYALALLLLASALIAGGSILQTNFKSLSSTEFIPX"
"LFCMLLLIVAGILFILAILTEWGSGNRTYGPVFMCLGGLLTMVAGAVWLTVMTNTLLSAWILTAGFLIFLIX"
"FALFGVIRCCRYCCYYCLTLESEERPPTPYRNTV*")
             (map #(sequence-string
                    (translate-biosequence
                     (gb/get-interval-sequence % gseq)
                     (:frame %)))
                  (gb/interval-seq gfeat)))))))
# `clj-biosequence`

`clj-biosequence` is a library designed to make working with biological sequence data easier. Currently it supports working with Genbank, Uniprot and generic fasta sequences. Basic functionality includes:

- Parses and accessors for Genbank and Uniprot XML and fasta format.
- A wrapper for BLAST.
- A wrapper for signalP.
- A very simple persistence mechanism.
- Mechanisms for lazy processing of sequences from very large sequence files.
- Interfaces for search and retrieval of sequences from online databases.
- Translation functions for DNA and RNA sequences.
- ID mapping functionality using the Uniprot's ID mapping tool.

Written by biologists and not computer scientists so improvements and
suggestions welcome.

[API docs](http://s312569.github.io/clj-biosequence/)

## Installation

`clj-biosequence` is available from [Clojars](https://clojars.org/clj-biosequence). For the current version add the following to your project.clj file:

```clojure
[clj-biosequence "0.1.4-SNAPSHOT"]
```

To use in your namespace:

```clojure
(ns my-app.core
  (:require [clj-biosequence.core :as cbs] ;; for base functionality and fasta
  	    [clj-biosequence.uniprot :as up] ;; for Uniprot functionality
	    [clj-biosequence.genbank :as gb] ;; for Genbank functionality
	    [clj-biosequence.blast :as bl] ;; for BLAST functionality
	    [clj-biosequence.signalp :as sp] ;; for a wrapper to signalp
	    [clj-biosequence.utilities :as ut] ;; some utilities))
```

## Usage

### Core functionality

`clj-biosequence` provides access to sequences from a variety of sources, including files, persistent stores and websites. The `clj-biosequence.core` library provides the core functionality as well as an implementation for handling fasta sequences. So if working with a fasta formatted file in the REPL, a workflow could go something like this:

```clojure
user> (use 'clj-biosequence.core)

;; initialise the file containing the sequences

user> (def f-file (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa" :protein))

;; The macro `with-biosequences` provides a handle to a lazy list of biosequence objects, from
;; a variety of sequence sources including biosequence files, stores and collections of biosequences.
;; In this case a list of `fastaSequence` objects is returned from a `fastaFile`.

user> (with-biosequences [list f-file]
        (first list))
-> \#clj_biosequence.core.fastaSequence{:accession "sp|P84001|29C0_ANCSP", :description "U3-ctenitoxin-Asp1a (Fragment) OS=Ancylometes sp. PE=1 SV=1", :type :protein, :sequence "ANACTKQADCAEDECCLDNLFFKRPYCEMRYGAGKRCAAASVYKEDKDLY"}
```

This workflow can be used for all sequence types defined in `clj-biosequence` (see below for working with Uniprot and Genbank sequences). The core library also implements a protocol that all sequences objects satisfy. The protocol provides the following functions that work on all sequence types (if only to return nil for some of them in the case of fasta sequences):

- `accession` - the primary accession.
- `accessions` - other accessions associated with the sequence.
- `sequence-string` - the sequence of the string.
- `def-line` - the definition line of a sequence.
- `protein?` - boolean
- `fasta-string` - a string of the sequence in fasta format.
- `org-scientific-name` - scientific name of the organism from which the sequence is derived
- `created` - date sequence was created.
- `modified` - date sequence was modified.
- `version` - version of the sequence.
- `database` - database the sequence is from.
- `taxonomy` - returns a list containing all the taxa of the organism from which the sequence is derived.
- `taxid` - the NCBI identification of the organism.

All sequence objects are records so extra information can be associated with a sequence, for example BLAST results, using `assoc`. This can be handy when using the presistence capabilities described below.

### Uniprot

`clj-biosequence.uniprot` provides an interface to Uniprot sequences in the UniprotXML format.

```clojure
user> (use 'clj-biosequence.uniprot)

;; initialise a uniprotxml file
user> (def ufile (init-uniprotxml-file "/Users/jason/Dropbox/clj-biosequence/test-files/uniprot-s-mansoni-20121217.xml"))

;; use `with-biosequences` as above to access a lazy list of uniprotProtein
;; objects.

user> (with-biosequences [l ufile]
          (accession (first l)))
"C4PYP8"
```

`uniprotProtein` objects implement the protocol described above as well as a few specific functions for extracting information from a uniprot biosequence. All of the following return either a `clojure.data.xml.Element` or a list of these elements. Information can then be accessed using the `clojure.data.xml`, `clojure.data.zip.xml` and `clojure.zip` libraries.

- `organism` - extensive organism information
- `prot-name` - the name of the protein
- `amino-acids` - extensive sequence information
- `nomenclature` - extensive naming information
- `gene` - gene information
- `gene-location` - gene location information
- `citation` - citations
- `subcellular-location`
- `alternative-products`
- `interactions`
- `mass-spectroscopy`
- `comments`
- `db-references`
- `existence` - experimental evidence for the existence of the protein
- `keywords`
- `features` - records describing the sequence's features

`clj-biosequence.uniprot` also provides a web interface to uniprot comprised of two functions. A search function `wget-uniprot-search`:

```clojure
;; `wget-uniprot-search` returns a non-lazy list of uniprot accession numbers
;; corresponding to a search term.
;; The search term has the same format as that used at the Uniprot website
;; (http://www.uniprot.org/).
;; So to retrieve accessions from Schistosoma mansoni reference proteome set 
;; that are intrinsic to the membrane:

user> (take 10 (wget-uniprot-search "taxonomy:6183 AND keyword:1185 AND go:0031224"))

("Q26597" "Q8MZK8" "C4Q533" "C4PYI6" "C4PYZ0" "Q27779" "Q26586" "Q86D97" "C4PY08" "Q26579")

;; or all reviewed human sequences

user> (count (wget-uniprot-search "reviewed:yes AND organism:9606"))
20264
```
A lazy list of sequences can be obtained from uniprot using the macro `with-wget-uniprot-sequence`. This macro provides a handle to a lazy list of biosequences from Uniprot corresponding to those specified in the list of accessions provided as an argument. Sequences can be returned as `uniprotProtein` or `fastaSequence` objects:

```clojure
;; `with-wget-uniprot-sequence` takes a collection of Uniprot accession numbers
;;  and returns a handle to a lazy list of sequence objects from Uniprot.

user> (with-wget-uniprot-sequence [l (take 10 (wget-uniprot-search "taxonomy:6183 AND keyword:1185 AND go:0031224")) :xml "jason.mulvenna@gmail.com"]
        (doseq [seq l]
          (println (accession seq))))
Q26597
Q8MZK8
C4Q533
C4PYI6
C4PYZ0
Q27779
Q26586
Q86D97
C4PY08
Q26579

user> (with-wget-uniprot-sequence [l '("Q26597" "Q8MZK8") :fasta "jason.mulvenna@gmail.com"]
        (doseq [seq l]
          (println (class seq))))
clj_biosequence.core.fastaSequence
clj_biosequence.core.fastaSequence
```

Downloaded sequences can be streamed to a file or a persistent store.

### GenBank

`clj-biosequence.uniprot` provides an interface to Genbank sequences in the GenbankXML format.

The same core functionality is available for Genbank as is available for Uniprot. For example:

```clojure

;; initialise a genbank XML file

user> (def gb (init-genbank-file "/Users/jason/Dropbox/clj-biosequence/test-files/protein-gb.xml"))

;; access sequences using `with-biosequences`

user> (with-biosequences [l gb]
        (println (version (first l))))
1
```

In addition to these, the functions `moltype` and `gb-locus` are defined and return the molecule type and locus of the biosequence.

To access features of a Genbank sequence `feature-seq` is defined and returns a lazy list of feature objects. To retrieve information from feature objects `feature-type` and `qualifier-extract` can be used. `feature-type` returns the type of featuer --- protein, Region, Site, CDS etc. `qualifier-extract` can be used to retrieve values from the 'GBFeature_quals' elements of a feature by supplying a qualifier name. Otherwise, further information can be obtained from the contents of `:src`, which is an xml element containing the feature and which can be accessed as described above.

The sequence corresponding to a feature can be returned using `get-feature-sequence` which takes a feature and its parent sequence and returns a fastaSequence object containing the seqeunce of the feature.

Intervals of a feature object can be accessed using `interval-seq` which provides a non-lazy list of intervals. `get-interval-sequence` can be used to obtain interval sequences in the same way as described above for `get-feature-sequence`.

```clojure
;; get a genbank sequence from the file defined above

user> (def gbs (with-biosequences [l gb]
                 (first l)))

;; feature types in the sequence

user> (doseq [f (feature-seq gbs)]
        (println (feature-type f)))
source
Protein
Region
Site
Bond
Bond
Bond

;; disulfide bonds

user> (doseq [b (filter #(= (feature-type %) "Bond") (feature-seq gbs))]
        (if (= "disulfide" (qualifier-extract b "bond_type"))
          (println (feature-location b))))
bond(5,21)
bond(9,23)
bond(14,28)

;; intervals

user> (doseq [f (feature-seq gbs)]
        (print (str (feature-type f) ": "))
          (doseq [i (interval-seq f)]
		    (println (str (start i) " - " (end i)))))
source: 1 - 32
Protein: 1 - 32
Region: 1 - 32
Site: 1 - 1
unsure: 3 - 3
unsure: 12 - 12
unsure: 18 - 18
unsure: 19 - 19
unsure: 29 - 29
```

The function and macro `wget-genbank-search` and `with-wget-genbank-sequences` provide access to Genbank sequences over the web in the same way as described above for Uniprot.

### BLAST

The `clj-bioseqeunce.blast` package provides one function and two macros for running BLAST and providing access to the results. The core function `blast` takes a biosequence file, store (see below) or collection of biosequences and will BLAST every biosequence against the specified database using the specified parameters. Results of the search are written to the specified outfile and a blast Search object is returned initialised with the file path.

The parameters argument is a hash-map with two mandatory keys, :db and :program, which specifiy the database to be searched and the program (blastp etc) to search it with. An optional key :options contains another hash-map of BLAST parameters for changing the default parameters. They keys being BLAST command line parameters and the values the value of the argument.

```clojure

;; initialise a blast database

user> (def bdb (init-blast-db "/Users/jason/Dropbox/clj-biosequence/test-files/toxins.fasta" :protein))

;; initialise a fasta file

user> (def ffile (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa" :protein))

;; run blast

user> (blast ffile "/Users/jason/Dropbox/clj-biosequence/test-files/blast.txt" {:db bdb :program "blastp" :options {"-evalue" "10"}})
#clj_biosequence.blast.blastSearch{:src "/Users/jason/Dropbox/clj-biosequence/test-files/blast.txt"}
```

The macro `with-iterations-in-search` provides a handle to a lazy list of blast results from a blastSearch object. One iteration corresponding to one of the biosequences from the file, store or collection of biosequences. Which biosequence a search belongs to can be determined using `iteration-query-id`. BLAST hits and high-scoring segment pairs can be accessed using `hit-seq` on iteration objects and `hsp-seq` on blastHit objects respectively.

```clojure

;; produce a blastResult object

user> (def bs (blast ffile "/Users/jason/Dropbox/clj-biosequence/test-files/blast.txt" {:db bdb :program "blastp" :options {"-evalue" "10"}})

;; access iterations using `with-iterations-in-search`

user> (with-iterations-in-search [l bs]
        (doseq [i l]
          (println (iteration-query-id i))))
sp|P84001|29C0_ANCSP
sp|P84027|32C7_ANCSP
sp|C1IC47|3FN3_WALAE
sp|C1IC48|3FN4_WALAE
sp|C1IC49|3FN5_WALAE
sp|P84028|45C1_ANCSP
sp|C0JAT4|A1H1_LOXHI
sp|P0CE78|A1H1_LOXRE
sp|C0JAT9|A1H1_LOXSP
sp|C0JAT5|A1H2_LOXHI

;; access blast hits using `hit-seq`. The function `get-hit-value` provides access
;; to various hit values corresponding to values in blast xml. Common keys are
;; :Hit_id, :Hit_len, :Hit_accession, :Hit_def, :Hit_num

user> (with-iterations-in-search [l bs]
        (doseq [i l]
          (println (get-hit-value (second (hit-seq i)) :Hit_id))))
sp|B3EWT5|TXC20_CUPSA
sp|P33034|TOG3A_AGEAP
sp|Q5UFR8|NXS1_ACAPE
sp|C1IC47|3FN3_WALAE
sp|Q53B61|TXW72_OPHHA
sp|O76199|TX25_PHONI
sp|C0JAT6|A1H3_LOXHI
sp|P0CE79|A1H2_LOXRE
sp|C0JAU1|A1H2_LOXSP
sp|C0JAT6|A1H3_LOXHI
```

The macro `with-blast-results` provides a handle to a lazy list of search results without the worry of specifying an out-file. Once the macro exits the results are deleted so they should be stored elsewhere, for example in a biosequence store, if they will be needed again.

## SignalP

If you have the program `signalp`
([SignalP](http://www.cbs.dtu.dk/services/SignalP/)) in your path,
`clj-biosequence.core` provides a wrapper consisting of two
functions. The first, `signalp`, takes any type of biosequence and an
optional parameter hash-map and returns a `signalpResult` object
containing the results of the analysis. The function `signalp?`
function returns true or false depending on the value of the `:result`
parameter from a signalP search. Throws an exception if the
biosequence is not a protein (ie. the biosequence returns false from
the function `protein?`). Example:

```clojure

;; initialise a fasta file

user> (def ffile (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/toxins.fasta" :protein))

#'user/ffile

;; define a protein sequence

user> (def prot (with-biosequences [l ffile]
                  (first l)))
#'user/prot

;; run signalp on the protein

user> (signalp prot)

#clj_biosequence.signalp.signalpResult{:name "P84001", :cmax 0.11, :cpos 38, :ymax 0.116, :ypos 3, :smax 0.134, :spos 2, :smean 0.132, :D 0.125, :result "N", :Dmaxcut 0.45, :network "SignalP-noTM"}

;; does the sequence have a signal seqeunce

user> (signalp? prot)
false
```

## Translation

Any nucleotide biosequence can be translated using the
`translate-biosequence` function in `clj-biosequence.core`. The
function takes a biosequence and an integer frame argument and returns
a new `fastaSequence` object containing the translated sequence (using
the standard translation table) and annotated definition line etc. The
function also takes an optional argument which is a nested list
representing a codon table.

The function `codon-table` returns pre-prepared codon tables according
to a keyword argument. At the moment the only key it accepts is
`:standard` which returns the standard codon table.

```clojure

;; initialise a fasta file

user> (def ffile (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/nuc-sequences.fasta" :protein))

#'user/ffile

;; define a nucleotide

user> (def nucl (with-biosequences [l ffile]
                  (first l)))
#'user/nucl

;; show the nucleotide biosequence

user> nucl

#clj_biosequence.core.fastaSequence{:accession "gi|116025203|gb|EG339215.1|EG339215", :description "KAAN-aaa29f08.b1 Platypus_EST_Cell_line_1.0-4.0kb Ornithorhynchus anatinus cDNA similar to ref|NP_005715.1| tetraspan 3; tetraspanin TM4-A; tetraspan TM4SF; transmembrane 4 superfamily, member 8; tetraspanin 3 [Homo sapiens] sp|O60637|T4S8_HUMAN Transmembrane 4 superfamily, member 8 (Tetraspanin 3) (Tspan-3) (Tetraspanin TM4-A) pir|A592, mRNA sequence", :type :nucleotide, :sequence "GTACAAAAAAGTTGGCCCAGGCAGGACCGGCAGCAACAGGAGGAGGAGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCGCCATTTCACCCCGCGCCGCCGCCGCCGCTGCCGCCGCCGCCATGGGCCAGTGCGGCATCACCTCGTCCAAGACGGTGCTGGTTTTCCTCAACCTCATCTTCTGGGCGGCAGCAGGCATCCTGTGCTACGTGGGAGCCTATGTCTTCATTACTTACGACGACTACGATCATTTCTTTGAAGACGTGTACACGCTCATTCCCGCCGTTGTGATCATAGCTGTAGGCACGCTCCTTTTCATCATCGGACTTATCGGATGCTGTGCCACAATCCGAGAAAGCCGTTGTGGACTAGCAACGTTTGTGATCATCCTCCTCTTGGTTTTTGTCACAGAAGTTGTCGTCGTGGTTTTGGGTTACATTTACCGAGCAAAGGTGGAAAACGAGGTCGACCGCAGCATCGAGAAGGTGTACAGGGCCTATAATGAAACCTCCTCTGACGCTGCCCGTCTAGCTATTGATTA"}

;; translate it in third frame

user> (translate-biosequence nucl 3)

#clj_biosequence.core.fastaSequence{:accession "gi|116025203|gb|EG339215.1|EG339215-3", :description "KAAN-aaa29f08.b1 Platypus_EST_Cell_line_1.0-4.0kb Ornithorhynchus anatinus cDNA similar to ref|NP_005715.1| tetraspan 3; tetraspanin TM4-A; tetraspan TM4SF; transmembrane 4 superfamily, member 8; tetraspanin 3 [Homo sapiens] sp|O60637|T4S8_HUMAN Transmembrane 4 superfamily, member 8 (Tetraspanin 3) (Tspan-3) (Tetraspanin TM4-A) pir|A592, mRNA sequence - Translated frame: 3", :type :protein, :sequence "TKKLAQAGPAATGGGAAAAAAAAAAAISPRAAAAAAAAAMGQCGITSSKTVLVFLNLIFWAAAGILCYVGAYVFITYDDYDHFFEDVYTLIPAVVIIAVGTLLFIIGLIGCCATIRESRCGLATFVIILLLVFVTEVVVVVLGYIYRAKVENEVDRSIEKVYRAYNETSSDAARLAIDX"}
```


## License

Copyright © 2013

Distributed under the Eclipse Public License, the same as Clojure.

# `clj-biosequence`

`clj-biosequence` is a library designed to make working with
biological sequence data easier. Basic functions include:

- Parses and accessors for Genbank, Uniprot XML, fasta and fastq formats.
- A wrapper for BLAST.
- A wrapper for signalP.
- A wrapper for TMHMM.
- Indexing of files for random access.
- A simple MongoDB interface.
- Mechanisms for lazy processing of sequences from very large files.
- Interfaces for search and retrieval of sequences from online databases.
- Translation functions for DNA and RNA sequences.
- ID mapping functionality using the Uniprot's ID mapping tool.

## Installation

Available from [Clojars](https://clojars.org/clj-biosequence). For the
current version add the following to your project.clj file:

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
            [clj-biosequence.fastq :as fq] ;; for fastq functionality
            [clj-biosequence.index :as ind] ;; for indexing functionality
            [clj-biosequence.interproscan :as ips] ;; for interproscan functionality
	        [clj-biosequence.signalp :as sp] ;; for a wrapper for signalp
            [clj-biosequence.store :as st] ;; for a mongoDB interface
            [clj-biosequence.tmhmm :as tm] ;; for a wrapper for TMHMM))
```

## Basic usage

`clj-biosequence` provides a reader and sequence mechanism for the
lazy access of biosequences in a variety of formats. For example, if
working with fasta sequences a typical session in the REPL could go
something like:

```clojure

;; import core and fasta functions

user> (use 'clj-biosequence.core)

;; to use test files included in library use clojure.java.io.
;; Otherwise string or java file object can be used.
;; For fasta an alphabet is also required to initialise a file.

user> (def fa-file (init-fasta-file (resource "test-files/nuc-sequence.fasta") :iupacNucleicAcids))
#'user/fa-file

;; then `bs-reader` can be used with `with-open` and `biosequence-seq`
;; to get access to a lazy sequence of fasta sequences in the file.

user> (with-open [r (bs-reader fa-file)]
                 (realized? (biosequence-seq r)))
false
user> (with-open [r (bs-reader fa-file)]
                 (count (biosequence-seq r)))
6
```

And thats just about it. The same pattern is used for all sequence
formats supported (at the moment this includes Genbank xml, Uniprot
xml, fasta, fastq, bed and arf formats) and each format has a number
of accessor functions providing access to information contained in the
format. `clj-biosequence.core` also defines a set of functions
supported by all formats as outlined here.

Some examples:

```clojure

;; a lazy sequence of translations in six reading frames

user> (with-open [r (bs-reader fa-file)]
                 (->> (biosequence-seq r)
                      (mapcat #(six-frame-translation %))
                      realized?))
false
user> (with-open [r (bs-reader fa-file)]
                 (->> (biosequence-seq r)
                      (mapcat #(six-frame-translation %))
                      count))
36

;; `fasta-string` can be used to convert biosequences to fasta strings

user> (use 'clj-biosequence.uniprot)
nil
user> (def uniprot-f (init-uniprotxml-file
                       (resource "test-files/uniprot-s-mansoni-20121217.xml")))
#'user/uniprot-f
user> (with-open [r (bs-reader uniprot-f)]
                 (println (->> (biosequence-seq r) first fasta-string)))

>sp|C4PYP8|DRE2_SCHMA Anamorsin homolog | Fe-S cluster assembly protein DRE2 homolog [Schistosoma mansoni]
MEQCVADCLNSDDCVMIVWSGEVQEDVMRGLQVAVSTYVKKLQFENLEKFVDSSAVDSQLXHECSVILCGWPNSISVNILK
LGLLSNLLSCLRPGGRFFGRDLITGDWDSLKKNLTLSGYIXNPYQLSCENHLIFSASVPSNYTQGSSVKLPWANSDVEAAW
ENVDNSSDANGNIINTNTLLXQKSDLKTPLSVCGKEAATDSVGKKKRACKNCTCGLAEIEAAEEDKSDVPISSCGNCYLGD
XAFRCSTCPYRGLPPFKPGERILIPDDVLRADL

;; filters can be implemented pretty easily

user> (with-open [r (bs-reader fa-file)]
                 (->> (biosequence-seq r)
                      (filter #(second (re-find #"(Mus musculus)" (def-line %))))
                      first
                      accession))
"gi|114311762|gb|EE738912.1|EE738912"

;; The function `biosequence->file` sends biosequences to a file and
;; also accepts a function argument to transform the biosequence
;; before writing (the default is `fasta-string`).

;; a Uniprot to fasta converter is thus:

user> (with-open [r (bs-reader uniprot-f)]
                 (biosequence->file (biosequence-seq r) "/tmp/fasta.fa"))
"/tmp/fasta.fa"
user> (with-open [r (bs-reader (init-fasta-file "/tmp/fasta.fa" :iupacAminoAcids))]
                 (count (biosequence-seq r)))
2
user> (with-open [r (bs-reader (init-fasta-file "/tmp/fasta.fa" :iupacAminoAcids))]
                 (class (first (biosequence-seq r))))
clj_biosequence.core.fastaSequence

;; sequences can be filtered to file using this function
;; for eg. filter Cytoplasmic proteins to file in fasta format

user> (with-open [r (bs-reader uniprot-f)]
                 (biosequence->file
                  (->> (biosequence-seq r)
                       (filter #(some (partial = "Cytoplasm")
                           (map :text (subcellular-location %)))))
                  "/tmp/fasta.fa"))
"/tmp/fasta.fa"
user> (with-open [r (bs-reader (init-fasta-file "/tmp/fasta.fa" :iupacAminoAcids))]
                 (count (biosequence-seq r)))
1
user> (with-open [r (bs-reader (init-fasta-file "/tmp/fasta.fa" :iupacAminoAcids))]
                 (println (fasta-string (first (biosequence-seq r)))))
>sp|C4PYP8|DRE2_SCHMA Anamorsin homolog | Fe-S cluster assembly protein DRE2 homolog [Schistosoma mansoni]
MEQCVADCLNSDDCVMIVWSGEVQEDVMRGLQVAVSTYVKKLQFENLEKFVDSSAVDSQLXHECSVILCGWPNSISVNILK
LGLLSNLLSCLRPGGRFFGRDLITGDWDSLKKNLTLSGYIXNPYQLSCENHLIFSASVPSNYTQGSSVKLPWANSDVEAAW
ENVDNSSDANGNIINTNTLLXQKSDLKTPLSVCGKEAATDSVGKKKRACKNCTCGLAEIEAAEEDKSDVPISSCGNCYLGD
XAFRCSTCPYRGLPPFKPGERILIPDDVLRADL
```
For strings containing fasta, Uniprot XML or Genbank XML formatted sequences the functions
`init-fasta-string`, `init-uniprot-string` and `init-genbank-string` allow the use of
strings with the `with-open` idiom. For Uniprot and Genbank connection initialisation
functions provide the same capability with remotely stored sequences from the relevant
servers (see below).

## Indexing

It can get tedious using `with-open` when using sequence files and random access to
particular biosequences is also slow as it relies on filtering the lazy sequences for
accession numbers. So `clj-biosequence.index` provides a simple mechanism for indexing
biosequence files.

Typical usage as follows:

```clojure
;; calling `index-biosequence-file` on any biosequence file returns a
;; biosequence index. This index can be used in calls to `biosequence-seq`
;; and `get-biosequence` (amongst others) without the `with-open` construct

user> (use 'clj-biosequence.index)
nil
user> (def fasta-in (index-biosequence-file fa-file))
#'user/fasta-in
user> (count (biosequence-seq fasta-in))
6
user> (first (biosequence-seq fasta-in))
#clj_biosequence.core.fastaSequence{:acc "gi|116025203|gb|EG339215.1|EG339215", :description "KAAN-aaa29f08.b1 ... etc"

;; Indexed files also offer random access to the biosequences,
;; this is many times faster than using `get-biosequence` with
;; a biosequence file opened with `bs-reader`.

user> (accession (get-biosequence fasta-in "gi|114311762|gb|EE738912.1|EE738912"))
"gi|114311762|gb|EE738912.1|EE738912"

;; when a file is indexed two additional files are created with the same
;; base-name as the biosequence file but with the extensions `.bin` and `.idx`.
;; The `.bin` file is compressed sequences and the `.idx` file is a 
;; text file containing the index. The `.idx` file is readable with
;; edn/read-string. To load an index use `load-biosequence-index` with the
;; path and basename of the index files.

user> (def fa-ind-2 (load-biosequence-index "/Users/jason/Dropbox/clj-biosequence/resources/test-files/nuc-sequence.fasta"))
#'user/fa-ind-2
user> (accession (get-biosequence fa-ind-2 "gi|114311762|gb|EE738912.1|EE738912"))
"gi|114311762|gb|EE738912.1|EE738912"

;; biosequence collections can be indexed using `index-biosequence-list` but the
;; base name of the index needs to be supplied

user> (def fa-ind-3 (with-open [r (bs-reader fa-file)]
                               (index-biosequence-list (biosequence-seq r)
                                                       "/tmp/fasta-ind")))
#'user/fa-ind-3
user> (accession (get-biosequence fa-ind-3 "gi|114311762|gb|EE738912.1|EE738912"))
"gi|114311762|gb|EE738912.1|EE738912"

;; this can be handy when filtering biosequences. For example secreted proteins
;; can be filtered into their own index

user> (def secreted (with-open [r (bs-reader toxins)]
                               (index-biosequence-list (-> (take 20 (biosequence-seq r))
                                                           (filter-signalp :trim true))
                                                       "/tmp/secreted-ind")))
#'user/secreted
user> (count (biosequence-seq secreted))
6

;; Finally, multi biosequence files can be merged into a single index using.
;; `index-biosequence-multi-file`. Once again the path and basename of the index
;; files needs to be supplied.
```
## BLAST

`clj-biosequence` supports most forms of BLAST, with the exception of PSI-BLAST. As
with other parts of clj-biosequence the BLAST functions seek to be as lazy and composable
as possible.

Typical usage as follows:

```clojure
;; initialise a BLAST db by passing the basename of the indexes to `init-blast-db`

user> (use 'clj-biosequence.blast)
nil
user> (def toxindb (init-blast-db (resource "test-files/toxins.fasta") :iupacAminoAcids))
#'user/toxindb

;; The function`blast` takes a list of biosequence objects and blasts them against
;; a blast database. It returns a blast search result which is a pointer to the
;; blast result file. This can then be opened using `bs-reader` and results
;; accessed using `biosequence-seq`

user> (use 'clj-biosequence.blast)
nil
user> (def toxindb (init-blast-db (resource "test-files/toxins.fasta") :iupacAminoAcids))
#'user/toxindb
user> (def tox-bl (with-open [r (bs-reader toxins)]
                             (blast (take 20 (biosequence-seq r))
                                    "blastp"
                                    toxindb
                                    "/tmp/blast.xml")))
#'user/tox-bl
user> (with-open [r (bs-reader tox-bl)]
                 (count (biosequence-seq r)))
20

;; BLAST results can be accessed using the accessors defined in the package
;; and the functions `hit-seq` and `hsp-seq`. For example to filter all
;; proteins in `tox-bl` that had a hit with a bit-score greater than
;; 50 and report their accession (note the use of second to avoid hits
;; to themselves):

user> (with-open [r (bs-reader tox-bl)]
                 (count (biosequence-seq r)))
20
user> (with-open [r (bs-reader tox-bl)]
                 (doall (->> (biosequence-seq r)
                             (filter #(>= (-> (hit-seq %) second hit-bit-scores first) 50))
                             (map #(-> (hit-seq %) second hit-accession))))
("B3EWT5" "Q5UFR8" "C1IC47" "Q53B61" "O76199" "C0JAT6" "P0CE79" "C0JAU1" "C0JAT6" "P0CE78"
"C0JAT9" "C0JAT5" "C0JAT9" "C0JAT6" "P0CE81" "P0CE80" "P0CE81" "P0CE82" "P0CE81")

;; Or a hash-map of the query id and  hit id of hits with a bit score greater than 50
;; (note that calling `accession` on a BLAST iteraton returns the query accession):

user> (with-open [r (bs-reader tox-bl)]
                 (->> (biosequence-seq r)
                      (filter #(>= (-> (hit-seq %) second hit-bit-scores first) 50))
                      (map #(vector (accession %)
                                    (-> (hit-seq %) second hit-accession)))
                      (into {}))))
{"sp|P84001|29C0_ANCSP" "B3EWT5", "sp|P0CE81|A1HB1_LOXIN" "P0CE80", "sp|C0JAT9|A1H1_LOXSP"
"C0JAU1", "sp|P0CE82|A1HB2_LOXIN" "P0CE81", "sp|P0CE80|A1HA_LOXIN" "P0CE81",
"sp|C0JAT8|A1H4_LOXHI" "C0JAT6", "sp|C0JAT5|A1H2_LOXHI" "C0JAT6", "sp|C0JAT6|A1H3_LOXHI"
"C0JAT5", "sp|C0JAT4|A1H1_LOXHI" "C0JAT6", "sp|C0JAU1|A1H2_LOXSP" "C0JAT9",
"sp|C0JAU2|A1H3_LOXSP" "C0JAT9", "sp|Q4VDB5|A1H_LOXGA" "P0CE82", "sp|C1IC47|3FN3_WALAE"
"Q5UFR8", "sp|C1IC48|3FN4_WALAE" "C1IC47", "sp|C1IC49|3FN5_WALAE" "Q53B61",
"sp|P84028|45C1_ANCSP" "O76199", "sp|Q56JA9|A1H_LOXSM" "P0CE81", "sp|P0CE78|A1H1_LOXRE"
"P0CE79", "sp|P0CE79|A1H2_LOXRE" "P0CE78"}

;; This can be combined with indexes or biosequence files to obtain the original
;; query biosequences.

user> (with-open [r (bs-reader tox-bl)]
                 (->> (biosequence-seq r)
                      (filter #(>= (-> (hit-seq %) second hit-bit-scores first) 50))
                      (map #(get-biosequence toxin-index (accession %)))
                      first))
#clj_biosequence.core.fastaSequence{:acc "sp|P84001|29C0_ANCSP", :description
"U3-ctenitoxin-Asp1a (Fragment) OS=Ancylometes sp. PE=1 SV=1", :alphabet :iupacAminoAcids,
:sequence [\A \N \A \C \T \K \Q \A \D \C \A \E \D \E \C \C \L \D \N \L \F \F \K \R \P \Y
\C \E \M \R \Y \G \A \G \K \R \C \A \A \A \S \V \Y \K \E \D \K \D \L \Y]}

;; or sent off to a file.

user> (with-open [r (bs-reader tox-bl)]
                 (biosequence->file
                  (->> (biosequence-seq r)
                       (filter #(>= (-> (hit-seq %) second hit-bit-scores first) 50))
                       (map #(get-biosequence toxin-index (accession %))))
                  "/tmp/blast.fa"))
"/tmp/blast.fa"

;; BLAST readers also provide access to the parameters used and these
;; can be accessed by calling `parameters` on the reader. This will return a
;; blast parameters object with accessors defined in the package.


;; As the entire chain is lazy these methods will work with as big a file as
;; can be thrown at them (hopefully). So one could annotate a large fasta file
;; with something like:

user> (with-open [r (bs-reader tox-bl)]
                 (biosequence->file
                  (->> (biosequence-seq r)
                       (filter #(>= (-> (hit-seq %) second hit-bit-scores first) 50))
                       (map #(let [s (get-biosequence toxin-index (accession %))
                                   h (-> (hit-seq %) first)]
                               (assoc s :description
                                      (str (def-line s) " - "
                                           "Similar to " (hit-def h) " - "
                                           (first (hit-bit-scores h)))))))
                  "/tmp/annotated-sequeunces.fa"))
"/tmp/annotated-sequeunces.fa"
user> (with-open [r (bs-reader (init-fasta-file "/tmp/annotated-sequeunces.fa"
                                                :iupacAminoAcids))]
                 (println (def-line (first (biosequence-seq r)))))
U3-ctenitoxin-Asp1a (Fragment) OS=Ancylometes sp. PE=1 SV=1 - Similar to Toxin CSTX-20 \
OS=Cupiennius salei PE=1 SV=1 - 89.737335

;; Although this is getting a bit complicated for the REPL and should probably
;; be a function of itself.

;; BLAST searches can be indexed like any other biosequence file. In which case
;; the index are keyed to the query accession.

user> (def blast-ind (index-biosequence-file tox-bl))
#'user/blast-ind
user> (-> (get-biosequence blast-ind "sp|Q56JA9|A1H_LOXSM") hit-seq first hit-accession)
"P0CE82"


```








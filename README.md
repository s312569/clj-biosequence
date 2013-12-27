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

```

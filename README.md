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
;; before writing (the default is `fasta-string`):



```

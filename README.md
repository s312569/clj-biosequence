# `clj-biosequence`

`clj-biosequence` is a library designed to make working with biological sequence data easier. Currently it supports working with Genbank, Uniprot and generic fasta sequences. Basic functionality includes:

- Parses and accessors for Genbank and Uniprot XML and fasta format.
- A wrapper for BLAST.
- A wrapper for signalP.
- A very simple persistence mechanism.
- Mechanisms for lazy processing of sequences from very large sequence files.
- Interfaces for search and retrieval of sequences from online databases.

Written by biologists and not computer scientists so improvements and
suggestions welcome.

## Installation

`clj-biosequence` is available from [Clojars](https://clojars.org/clj-biosequence). For the current version add the following to your project.clj file:

```clojure
[clj-biosequence "0.1.2-SNAPSHOT"]
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

`clj-biosequence` provides access to sequences from a variety of sources, including files, persistent stores and web-based websites. The clj-biosequence.core library provides the core functionality as well as an implementation for handling fasta sequences. So, for example, working with a fasta formatted file in the REPL could go something like this:

```clojure
user> (use 'clj-biosequence.core)

;; initialise the file containing the sequences

user> (def f-file (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa" :protein))

;; to access the sequences use the macro 'with-biosequences-in-file'. This provides
;; a handle to a lazy list of sequence objects, the type depending on the file
;; type. In this case a list of fastaSequence objects is returned.

user> (with-biosequences-in-file [list f-file]
        (first list))
\#clj_biosequence.core.fastaSequence{:accession "sp|P84001|29C0_ANCSP", :description "U3-ctenitoxin-Asp1a (Fragment) OS=Ancylometes sp. PE=1 SV=1", :type :protein, :sequence "ANACTKQADCAEDECCLDNLFFKRPYCEMRYGAGKRCAAASVYKEDKDLY"}
```

This workflow is applicable for all sequence types defined in `clj-biosequence` (see below for working with Uniprot and Genbank sequences). The core library also implements a protocol that all sequences objects satisfy. The protocol provides the following functions that work on all sequence types (if only to return nil in the case of fasta sequences):

- accession - the primary accession.
- accessions - other accessions associated with the sequence.
- sequence-string - the sequence of the string.
- def-line - the definition line of a sequence.
- protein? - boolean
- fasta-string - a string of the sequence in fasta format.
- fasta-to-stream - outputs a sequence in fasta format to a stream.
- org-scientific-name - scientific name of the organism from which the sequence is derived
- created - date sequence was created.
- modified - date sequence was modified.
- version - version of the sequence.
- database - database the sequence is from.
- taxonomy - returns a list containing all the taxa of the organism from which the sequence is derived.
- taxid - the NCBI identification of the organism.


## License

Copyright © 2013 FIXME

Distributed under the Eclipse Public License, the same as Clojure.

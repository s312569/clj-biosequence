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

`clj-biosequence` provides access to sequences from a variety of sources, including files, persistent stores and web-based websites. The `clj-biosequence.core` library provides the core functionality as well as an implementation for handling fasta sequences. So, for example, working with a fasta formatted file in the REPL could go something like this:

```clojure
(use 'clj-biosequence.core)

;; initialise the file containing the sequences

(def f-file (init-fasta-file "/Users/jason/Dropbox/clj-biosequence/test-files/bl-test.fa" :protein))

;; to access the sequences use the macro `with-biosequences-in-file`. This provides
;; a handle to a lazy list of sequence objects, the type depending on the file
;; type. In this case a list of `fastaSequence` objects is returned.

(with-biosequences-in-file [list f-file]
  (first list))
-> \#clj_biosequence.core.fastaSequence{:accession "sp|P84001|29C0_ANCSP", :description "U3-ctenitoxin-Asp1a (Fragment) OS=Ancylometes sp. PE=1 SV=1", :type :protein, :sequence "ANACTKQADCAEDECCLDNLFFKRPYCEMRYGAGKRCAAASVYKEDKDLY"}
```

This workflow can be used for all sequence types defined in `clj-biosequence` (see below for working with Uniprot and Genbank sequences). The core library also implements a protocol that all sequences objects satisfy. The protocol provides the following functions that work on all sequence types (if only to return nil for some of them in the case of fasta sequences):

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

All sequence objects are records so extra information can be associated with a sequence, for example BLAST results, using 'assoc'. This can be handy when using the presistence capabilities described below.

### Uniprot

`clj-biosequence.uniprot` provides an interface to Uniprot sequences in the UniprotXML format.

```clojure
(use 'clj-biosequence.uniprot)

;; initialise a uniprotxml file
(def ufile (init-uniprotxml-file "/Users/jason/Dropbox/clj-biosequence/test-files/uniprot-s-mansoni-20121217.xml"))

;; use `with-biosequences-in-file` as above to access a lazy list of uniprotProtein
;; objects.

(with-biosequences-in-file [l ufile]
        (accession (first l)))
-> "C4PYP8"
```

`uniprotProtein` objects implement the protocol described above as well as a few specific functions (see the docs for detailed information):

- `organism-name` - extensive organism information
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
- `features` - objects describing the sequences features

`clj-biosequence.uniprot` also provides a web interface to uniprot comprised of two functions. A search function `wget-uniprot-search`:

```clojure
;; `wget-uniprot-search` returns a lazy'ish (results are obtained from the server
;; 1000 at a time) list of uniprot accession numbers corresponding to a search term.
;; The search term has the same format as that used at the Uniprot website
;; (http://www.uniprot.org/).
;; So to retrieve accessions from Schistosoma mansoni reference proteome set 
;; that are intrinsic to the membrane:

(take 10 (wget-uniprot-search "taxonomy:6183 AND keyword:1185 AND go:0031224"))

->("Q26597" "Q8MZK8" "C4Q533" "C4PYI6" "C4PYZ0" "Q27779" "Q26586" "Q86D97" "C4PY08" "Q26579")

;; or all reviewed human sequences

(count (wget-uniprot-search "reviewed:yes AND organism:9606"))
->20264
```
The retrieval function, `wget-uniprot-sequence`, retrieves sequences from Uniprot. Sequences can be returned as `uniprotProtein` or `fastaSequence` objects:

```clojure
;; `wget-uniprot-sequence` takes a collection of Uniprot accession numbers
;; and returns a lazy'ish (1000 retrieved at a time) list of sequence
;; objects

(doseq [s (wget-uniprot-sequence '("Q26597" "Q8MZK8") :fasta)]
        (println (class s)))
->clj_biosequence.core.fastaSequence
  clj_biosequence.core.fastaSequence

(doseq [s (wget-uniprot-sequence '("Q26597" "Q8MZK8") :xml)]
        (println (:mass (amino-acids s))))
->58602
 84159
```

Downloaded sequences can be streamed to a file or a persistent store.

## License

Copyright © 2013 FIXME

Distributed under the Eclipse Public License, the same as Clojure.

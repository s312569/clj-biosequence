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
	    [clj-biosequence.signalp :as sp] ;; for a wrapper to signalp))
```

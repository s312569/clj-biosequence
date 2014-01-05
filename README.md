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

It can get tedious using `with-open` and random access to particular biosequences
is also slow as it relies on filtering the lazy sequences for accession numbers.
So `clj-biosequence.index` provides a simple mechanism for indexing biosequence files.

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

`clj-biosequence` supports most forms of BLAST, with the exception of
PSI-BLAST. As with other parts of clj-biosequence the BLAST functions
seek to be as lazy and composable as possible. To work the various
BLAST+ programs from the NCBI need to be in your path.

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

user> (def tox-bl (with-open [r (bs-reader toxins)]
                             (blast (take 20 (biosequence-seq r))
                                    "blastp"
                                    toxindb
                                    "/tmp/blast.xml")))
#'user/tox-bl
user> (with-open [r (bs-reader tox-bl)]
                 (count (biosequence-seq r)))
20

;; Addiitonal parameters can be passed to `blast` using the `:param` keyword
;; argument. Format is a hash-map with keys strings of the command line switches
;; with the desired value as a string. For example:

user> (def tox-bl (with-open [r (bs-reader toxins)]
                             (blast (take 20 (biosequence-seq r))
                                    "blastp"
                                    toxindb
                                    "/tmp/blast.xml"
                                    :params {"-max_target_seqs" "3"
                                             "-evalue" "1"})))
#'user/tox-bl

;; BLAST results can be accessed using the accessors defined in the package
;; and the functions `hit-seq` and `hsp-seq`. For example to filter all
;; proteins in `tox-bl` that had a hit with a bit-score greater than
;; 50 and report their accession (note the use of second to avoid hits
;; to themselves):

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

;; As the entire chain is lazy these methods will work with as big a file as
;; can be thrown at them (hopefully). So one could annotate a large fasta file
;; starting with a fasta index and a blast DB by:

user> (def toxin-index (index-biosequence-file toxins))
#'user/toxin-index
user> (with-open [r (bs-reader (blast (biosequence-seq toxin-index) "blastp" toxindb
                                      "/tmp/outfile.xml"))]
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
;; be a function(s) of itself (and the blast outfile might need to be deleted).

;; BLAST readers also provide access to the parameters used and these
;; can be accessed by calling `parameters` on the reader. This will return a
;; blast parameters object with accessors defined in the package.

user> (with-open [r (bs-reader tox-bl)]
                 (blast-database (parameters r)))
"/Users/jason/Dropbox/clj-biosequence/resources/test-files/toxins.fasta"
user> (with-open [r (bs-reader tox-bl)]
                 (blast-version (parameters r)))
"BLASTP 2.2.24+"
user> (with-open [r (bs-reader tox-bl)]
                 (blast-evalue (parameters r)))
"10"
user> (with-open [r (bs-reader tox-bl)]
                 (blast-filter (parameters r)))
"F"

;; BLAST searches can be indexed like any other biosequence file. In which case
;; the index is keyed to the query accession. Although, the parameter information
;; is lost.

user> (def blast-ind (index-biosequence-file tox-bl))
#'user/blast-ind
user> (-> (get-biosequence blast-ind "sp|Q56JA9|A1H_LOXSM") hit-seq first hit-accession)
"P0CE82"

```

## SignalP

SignalP works in a similar way as BLAST. If you have signalp in your
path it can be applied to collections of bioseqeunces using the
function `signalp` (which returns a signalp result object) or a
SignalP output file in short form format can be initialised as a
result object.

Basic usage as follows:

```clojure
;; running signalp

user> (use 'clj-biosequence.signalp)
nil
user> (def sr (signalp (take 20 (biosequence-seq toxin-index)) "/tmp/signalp.txt"))
#'user/sr
user> (with-open [r (bs-reader sr)]
                 (first (biosequence-seq r)))
#clj_biosequence.signalp.signalpProtein{:name "sp|P58809|CTX_CONMR", :cmax 0.105,
:cpos 7, :ymax 0.147, :ypos 1, :smax 0.208, :spos 1, :smean 0.0, :D 0.068, :result "N",
:Dmaxcut 0.45, :network "SignalP-noTM"}
user> (with-open [r (bs-reader sr)]
                 (accession (first (biosequence-seq r))))
"sp|P58809|CTX_CONMR"

;; `signalp?` can be used to determine if a result is positive or not

user> (with-open [r (bs-reader sr)]
                 (signalp? (first (biosequence-seq r))))
false
user> (with-open [r (bs-reader sr)]
                 (-> (filter signalp? (biosequence-seq r))
                     first))
#clj_biosequence.signalp.signalpProtein{:name "sp|Q9BP63|O3611_CONPE", :cmax 0.51,
:cpos 21, :ymax 0.696, :ypos 21, :smax 0.982, :spos 12, :smean 0.952, :D 0.834,
:result "Y", :Dmaxcut 0.45, :network "SignalP-noTM"}

;; a convenience function `filter-signalp` filters a collection of biosequence
;; proteins and returns only proteins containing a signal sequence. If the
;; keyword argument `:trim` is true the returned biosequences will have the
;; signal sequence trimmed from the sequence

user> (->> (filter-signalp (take 20 (biosequence-seq toxin-index)))
           first
           bioseq->string)
"MSRLGIMVLTLLLLVFIVTSHQDAGEKQATQRDAINFRWRRSLIRRTATEECEEYCEDEEKTCCGLEDGEPVCATTCLG"
user> (->> (filter-signalp (take 20 (biosequence-seq toxin-index)) :trim true)
           first
           bioseq->string)
"DAGEKQATQRDAINFRWRRSLIRRTATEECEEYCEDEEKTCCGLEDGEPVCATTCLG"

;; SignalP result objects can be indexed in the same manner as BLAST ie.
;; the query sequence accession becomes the index keys.

user> (def si (index-biosequence-file sr))
#'user/si
user> (accession (get-biosequence si "sp|P58809|CTX_CONMR"))
"sp|P58809|CTX_CONMR"
user> (signalp? (get-biosequence si "sp|P58809|CTX_CONMR"))
false

;; Search parameters can be passed to `signalp` and `filter-signalp` as hash-maps
;; using the `:param` keyword argument.

user> (def sr (signalp (take 20 (biosequence-seq toxin-index)) "/tmp/signalp.txt"
                       :params {"-s" "best" "-t" "gram+"}))
#'user/sr
user> (with-open [r (bs-reader sr)]
                 (first (biosequence-seq r)))
#clj_biosequence.signalp.signalpProtein{:name "sp|P58809|CTX_CONMR", :cmax 0.101,
:cpos 2, :ymax 0.119, :ypos 2, :smax 0.139, :spos 1, :smean 0.139, :D 0.127,
:result "N", :Dmaxcut 0.45, :network "SignalP-TM"}
```

## Accession mapping

`clj-biosequence` provides a facility for mapping accessions from one
database to another. It is provided in the core package and uses the
Uniprot mapping service so needs an active internet connection.

Basic usage:

```clojure
;; `id-convert` converts accessions. Its arguments are a list of accessions
;; to be converted, a 'from' database, a 'to' database and an email (required
;; by Uniprot). The 'from' and 'to' arguments are strings corresponding to
;; to the database codes used by the Uniprot mapping tool (full list at
;; http://www.uniprot.org/faq/28#id_mapping_examples and a partial list in
;; the doc string of `id-convert`.

;; `id-convert` returns a hash-map of query accessions and search results. If
;; an ID returned no result it is not in the result hash-map. There is a
;; 100,000 limit on individual queries imposed by Uniprot.

;; For example, to convert a list of Uniprot accessions to NCBI Genbank ids, 
;; using the previously defined toxin protein index which has accessions in
;; the format "sp|xxx|xxxx":

user> (map accession (take 5 (biosequence-seq toxin-index)))
("sp|P58809|CTX_CONMR" "sp|P61792|TXU2_HETVE" "sp|P86259|CT2X_CONTE"
"sp|Q9BP63|O3611_CONPE" "sp|A0SE59|CA13_CONMR")

user> (require '[clojure.string :as st])
nil
user> (-> (map #(second (st/split (accession %) #"\|"))
                (take 5 (biosequence-seq toxin-index)))
          (id-convert "ACC" "P_GI" "jason.mulvenna@gmail.com"))
{"P58809" "20454877", "P61792" "48428590", "P86259" "229485330", "Q9BP63" "74848505",
"A0SE59" "83657225"}
```

## Sequence retrieval

Sequences can be retrieved from both Genabnk and Uniprot using
`init-uniprot-connection` and `init-genbank-connection`. Both
functions take a list of accession numbers and a return type argument.
Uniprot also needs and email argument and Genbank a database argument.
Both functions can be used in conjunction with the search functions,
`genbank-search` and `uniprot-search`.

Basic usage:

```clojure
;; To generate a list of accessions search for Uniprot accessions (note this
;; generates a non-lazy list). Search syntax is exactly the same as Uniprot
;; search syntax (described at http://www.uniprot.org/help/text-search and
;; summarised in the doc string of `uniprot-search`).

;; For example, to get accessions of all proteins in the Schistosoma mansoni
;; reference proteome set:

user> (use 'clj-biosequence.uniprot)
nil
user> (def sm-prot (uniprot-search "organism:6183 AND keyword:1185" "jason.mulvenna@gmail.com"))
#'user/sm-prot
user> (count sm-prot)
11711
user> (first sm-prot)
"C4PYP8"

;; A lazy sequence of biosequences can be retrieved from Uniprot using
;; `init-uniprot-connection` and `bs-reader`. Sequences can be retrieved as
;; fasta or full Uniprot entries.

user> (def up-conn (init-uniprot-connection (take 10 sm-prot) :fasta "jason.mulvenna@gmail.com"))
#'user/up-conn
user> (with-open [r (bs-reader up-conn)]
                 (first (biosequence-seq r)))
#clj_biosequence.core.fastaSequence{:acc "sp|C4PYP8|DRE2_SCHMA", :description\
"Anamorsin homolog OS=Schistosoma mansoni GN=Smp_207000 PE=3 SV=2", :alphabet\
:iupacAminoAcids, :sequence [\M \E \Q \C \V \A \D \C \L \N \S \D \D \C \V \M ... etc

;; Uniprot

user> (def up-conn (init-uniprot-connection (take 10 sm-prot) :xml "jason.mulvenna@gmail.com"))
#'user/up-conn
user> (with-open [r (bs-reader up-conn)]
                 (class (first (biosequence-seq r))))
clj_biosequence.uniprot.uniprotProtein

;; Although sequences are downloaded as a compressed stream large sequence
;; downloads can take a long time ...

;; Genbank works exactly the same way. Search syntax is the same as Genbank query
;; format (see http://www.ncbi.nlm.nih.gov/books/NBK3837/ and a summary in the doc
;; doc string of `genbank-search`). A database also neds to specified and may be one
;; of :protein, :nucest, :nuccore, :nucgss or :popset.

;; So to get all Schistosoma mansoni proteins from Genbank

user> (use 'clj-biosequence.genbank)
nil
user> (def sm-prots (genbank-search "txid6183[Organism:noexp]" :protein))
#'user/sm-prots
user> (first sm-prots)
"566601372"
user> (with-open [r (bs-reader (init-genbank-connection (take 10 sm-prots) :protein :fasta))]
                 (second (biosequence-seq r)))
#clj_biosequence.core.fastaSequence{:acc "gi|566601352|gb|AHC70335.1|", :description
"nicotinic acetylcholine receptor [Schistosoma mansoni]", :alphabet :iupacAminoAcids,
:sequence [\M ... etc

user> (with-open [r (bs-reader (init-genbank-connection (take 10 sm-prots) :protein :xml))]
                 (class (second (biosequence-seq r))))
clj_biosequence.genbank.genbankSequence
```

## Supported formats

### Uniprot

```clojure

;; initialise uniprot files using `init-uniprotxml-file` and access
;; sequences as described above

user> (def up (init-uniprotxml-file (resource "test-files/uniprot-s-mansoni-20121217.xml")))
#'user/up
user> (with-open [r (bs-reader up)]
                 (count (biosequence-seq r)))
2

;; In addition to the accessors defined in the core module there are a
;; number of convenience functions (see documentation for full
;; details).

user> (with-open [r (bs-reader up)]
        (-> (biosequence-seq r) first sequence-info))
{:length 272, :mass 29662.0, :checksum "11A38702F71B84D8", :modified "2010-03-23", :version 2}

user> (with-open [r (bs-reader up)]
        (-> (biosequence-seq r) first lineage))
("Eukaryota" "Metazoa" "Platyhelminthes" "Trematoda" "Digenea"
"Strigeidida" "Schistosomatoidea" "Schistosomatidae" "Schistosoma")

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) first alternative-name))
"Fe-S cluster assembly protein DRE2 homolog"

;; citations can be extracted from the sequence in the form of
;; uniprotCitation objects that extend the citation protocol defined
;; in the core module (see documentation for details)

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) second citations first authors))
("Wright M.D." "Harrison R.A." "Melder A.M." "Newport G.R." "Mitchell G.F.")

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) second citations first ref-type))
"journal article"

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) second citations first journal))
"Mol. Biochem. Parasitol."

;; All comments can be accessed using the `comments` function which
;; returns the XML trees describing each comment.

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) first comments first))
#clojure.data.xml.Element{:tag :comment, :attrs {:type "function"},
 :content (#clojure.data.xml.Element{:tag :text, :attrs {}, :content
 ("May be required for the maturation of extramitochondrial Fe/S
 proteins (By similarity). Has anti-apoptotic effects in the cell (By
 similarity).")})}

;; A similar function, `db-references`, returns all database references.

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) first db-references first))
#clojure.data.xml.Element{:tag :dbReference, :attrs {:type "EMBL", :id
 "HE601631"}, :content (#clojure.data.xml.Element{:tag :property,
 :attrs {:type "protein sequence ID", :value "CCD81528.1"}, :content
 ()} #clojure.data.xml.Element{:tag :property, :attrs {:type "molecule
 type", :value "Genomic_DNA"}, :content ()})}

;; A convenience function, `go-terms`, uses `db-references` to extract
;; GO terms.

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) first go-terms))
("C:cytoplasm" "P:apoptotic process")

;; Features and intervals can be accessed through the protocols
;; defined in the core module.

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) second feature-seq first))
#clj_biosequence.uniprot.uniprotFeature{:src
 #clojure.data.xml.Element{:tag :feature, :attrs {:type "chain",
 :description "Glutathione S-transferase class-mu 26 kDa isozyme", :id
 "PRO_0000185811"}, :content (#clojure.data.xml.Element{:tag
 :location, :attrs {}, :content (#clojure.data.xml.Element{:tag
 :begin, :attrs {:position "1"}, :content ()}
 #clojure.data.xml.Element{:tag :end, :attrs {:position "195"},
 :content ()})})}}

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) second feature-seq first feature-type))
"chain"

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) second feature-seq first interval-seq first))
#clj_biosequence.uniprot.uniprotInterval{:src
 #clojure.data.xml.Element{:tag :location, :attrs {}, :content
 (#clojure.data.xml.Element{:tag :begin, :attrs {:position "1"},
 :content ()} #clojure.data.xml.Element{:tag :end, :attrs {:position
 "195"}, :content ()})}}

user> (with-open [r (bs-reader up)]
                 (-> (biosequence-seq r) second feature-seq first interval-seq first start))
1
```

### Genbank XML

`clj-biosequence` also supports Genbank XML.

```clojure
;; initialise a genbank file in the usual way

user> (use 'clj-biosequence.genbank)
nil
user> (def gbf (init-genbank-file (resource "test-files/nucleotide-gb.xml")))
#'user/gbf

;; Access to sequences is managed the same way as described above. As
;; with Uniprot sequences a number of convenience accessors are
;; defined (see the documentation for details)

user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first created))
"08-JUL-2013"
user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first taxid))
1268274
user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first gb-locus))
"KE373594"

;; The functions `feature-seq` and `interval-seq provide access to
;; features and intervals as described above for Uniprot sequences.

user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first feature-seq first feature-type))
"source"
user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first feature-seq first interval-seq first start))
1

;; The function `qualifier-seq` is also provided for genbank features
;; and it returns a lazy list of qualifiers.

user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first feature-seq first qualifier-seq first))
#clj_biosequence.genbank.genbankQualifier{:src
 #clojure.data.xml.Element{:tag :GBQualifier, :attrs {}, :content
 (#clojure.data.xml.Element{:tag :GBQualifier_name, :attrs {},
 :content ("organism")} #clojure.data.xml.Element{:tag
 :GBQualifier_value, :attrs {}, :content ("Blumeria graminis f. sp.
 tritici 96224")})}}
user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first feature-seq first
                     qualifier-seq first qualifier-name))
"organism"
user> (with-open [r (bs-reader gbf)]
                 (-> (biosequence-seq r) first feature-seq first
                     qualifier-seq first qualifier-value))
"Blumeria graminis f. sp. tritici 96224"

;; For Genbank files with large sequences, for example genome files,
;; `feature-seq` can also be called on a genbank reader to provide a
;; lazy sequence of features without loading the entire sequence into
;; memory. Note this will only access the first sequence in a file.

user> (def big-gb (init-genbank-file "/Users/jason/Dropbox/bl-ro1/hhv-4-annotation/akata/akata-sequence.xml"))
#'user/big-gb
user> (with-open [r (bs-reader big-gb)]
                 (-> (feature-seq r) first feature-type))
"source"
```






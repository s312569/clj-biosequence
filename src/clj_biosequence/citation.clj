(ns clj-biosequence.citation)

(defprotocol biosequenceCitation
  (title [this]
    "Returns the title of a citation object.")
  (journal [this]
    "Returns the journal of a citation object.")
  (year [this]
    "Returns the year of a citation object.")
  (volume [this]
    "Returns the volume of a citation object.")
  (pstart [this]
    "Returns the start page of a citation object.")
  (pend [this]
    "Returns the end page of a citation object.")
  (authors [this]
    "Returns the authors from a citation object.")
  (pubmed [this]
    "Returns the pubmed id of a reference if there is one.")
  (crossrefs [this]
    "Returns crossrefs - DOI, ISBN etc")
  (notes [this]
    "Returns any remarks in the citation.")
  (abstract [this]
    "Returns the abstract"))

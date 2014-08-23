(ns clj-biosequence.citation)

(defrecord biosequenceJournalCite
    [title journal year volume pstart
     pend authors citation pubmed crossrefs
     notes])



;; (defprotocol biosequenceCitation
;;   (ref-type [this]
;;     "Returns the citation type from a citation object.")
;;   (title [this]
;;     "Returns the title of a citation object.")
;;   (journal [this]
;;     "Returns the journal of a citation object.")
;;   (year [this]
;;     "Returns the year of a citation object.")
;;   (volume [this]
;;     "Returns the volume of a citation object.")
;;   (pstart [this]
;;     "Returns the start page of a citation object.")
;;   (pend [this]
;;     "Returns the end page of a citation object.")
;;   (authors [this]
;;     "Returns the authors from a citation object."))

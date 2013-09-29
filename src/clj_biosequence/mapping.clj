(in-ns 'clj-biosequence.core)

(defn id-convert
  "Takes either a single accession or a list of accessions and returns a hash-map
   mapping the accession numbers to the corresponding identification number in the 
   specified 'to' database. 'From' database also needs to be specified. If not 
   found returns an empty hash-map. Uses the Uniprot id mapping utility and a list 
   of supported databases is supplied at http://www.uniprot.org/faq/28#id_mapping_examples.
   Some common mappings include:
   DB Name                  Abbreviation     Direction
   UniProtKB AC/ID	    ACC+ID	     from
   UniProtKB AC	            ACC              to
   EMBL/GenBank/DDBJ	    EMBL_ID	     both
   EMBL/GenBank/DDBJ CDS    EMBL	     both
   Entrez Gene (GeneID)     P_ENTREZGENEID   both
   GI number	            P_GI	     both
   RefSeq Protein	    P_REFSEQ_AC	     both
   RefSeq Nucleotide	    REFSEQ_NT_ID     both
   WormBase	            WORMBASE_ID	     both

   There is a 100,000 limit on accessions in a single query imposed by Uniprot."
  [ids from to email]
  (let [i (if (seq? ids)
            ids (list ids))]
    (if (<= (count ids) 100000)
      (let [param {:from from :to to :format "tab" 
                   :query (apply str (doall (interpose "," (if (list? ids)
                                                             ids
                                                             (list ids)))))}
            address "http://www.uniprot.org/mapping/"
            r (client/post address 
                           {:client-params
                            {"http.useragent" (str "clj-http " email)}
                            :follow-redirects true
                            :force-redirects true
                            :form-params param})]
        (if-not (= 200 (:status r))
          (throw (Throwable. (str "Error in mapping request: " (:body r))))
          (into {} (map #(string/split % #"\t") (rest (string/split (:body r) #"\n"))))))
      (throw (Throwable. "No more than 100,000 mappings per query allowed.")))))

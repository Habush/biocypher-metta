# Author Abdulrahman S. Omar <xabush@singularitynet.io>
from biocypher_metta.adapters import Adapter
import pickle
import csv
import gzip

# Transcription factor - target gene relationships from TFLink

# UniprotID.TF	UniprotID.Target	NCBI.GeneID.TF	NCBI.GeneID.Target	Name.TF	Name.Target	Detection.method
# PubmedID	Organism	Source.database	Small-scale.evidence	TF.TFLink.ortho	TF.nonTFLink.ortho	Target.TFLink.ortho	Target.nonTFLink.ortho
# Q9H9S0	O94907	79923	22943	NANOG	DKK1	chromatin immunoprecipitation assay;inferred by curator	19148141;29087512;29126285;27924024	Homo sapiens	GTRD;ReMap;TRRUST	Yes	-	-	Dr:Q9PWH3;Dr:F1RBK0;Mm:O54908	Rn:D3Z9J1
# P37231	P10826	5468	5915	PPARG	RARB	chromatin immunoprecipitation assay;inferred by curator	17202159;12839938;29087512;27924024	Homo sapiens	GTRD;TRED;TRRUST	Yes	-	-	Mm:P22605;Rn:D3ZFD9	-
# P10242	P08047	4602	6667	MYB	SP1	chromatin immunoprecipitation assay;inferred by curator	29126285;27924024;17202159	Homo sapiens	GTRD;ReMap;TRED	Yes	Dr:F1QP24;Rn:A0A0G2K2A4	Mm:A0A087WPA7	Dr:F1QW97;Rn:Q01714	Mm:G3X8Q0

class TFLinkAdapter(Adapter):

    def __init__(self, filepath, entrez_to_ensemble_map):
        """
        Constructs TFLink adapter that returns edges between TFs and their target gene
        :param filepath: Path to the TSV file downloaded from tflink
        :param entrez_to_ensemble_map: file containing pickled dictionary mapping NCBI Entrez IDs to Ensemble IDs -
        this b/c we use Ensemble IDs to identify genes where TFLink uses Entrez Ids
        """
        self.filepath = filepath

        with open(entrez_to_ensemble_map, "rb") as f:
            self.entrez2ensemble = pickle.load(f)

        self.label = "regulates"
        self.source = "TFLink"
        self.source_url = "tflink.net"

    def get_edges(self):
        with gzip.open(self.filepath, 'rt') as fp:
            table = csv.reader(fp, delimiter="\t", quotechar='"')
            for row in table:
                tf_entrez_id = row[2]
                target_entrez_id = row[3]
                if tf_entrez_id in self.entrez2ensemble and target_entrez_id in self.entrez2ensemble:
                    tf_ensemble_id = self.entrez2ensemble[tf_entrez_id]
                    target_ensemble_id = self.entrez2ensemble[target_entrez_id]
                    _source = tf_ensemble_id
                    _target = target_ensemble_id
                    pubmed_ids = [f"pubmed:{i}" for i in row[7].split(";")]
                    sources = row[9].split(";")
                    small_scale_evidence = row[10]
                    if small_scale_evidence == "Yes":
                        evidence_type = "small_scale_evidence"
                    else:
                        evidence_type = "large_scale_evidence"
                    _props = {
                        "evidence": pubmed_ids,
                        "databases": sources,
                        "evidence_type": evidence_type,
                        "detection_method": row[6]
                    }
                    _id = f"{self.label}-{_source}-{_target}"
                    yield _id, _source, _target, self.label, _props

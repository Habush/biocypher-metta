
from biocypher_metta.adapters import Adapter
import pickle
import os

# https://coxpresdb.jp/download/Hsa-r.c6-0/coex/Hsa-r.v22-05.G16651-S235187.combat_pca.subagging.z.d.zip
# There is 16651 files. The file name is entrez gene id. The total genes annotated are 16651, one gene per file, each file contain logit score of other 16650 genes.
# There are two fields in each row: entrez gene id and logit score


class CoxpresdbAdapter(Adapter):

    def __init__(self, filepath, ensemble_to_entrez_path,
                 write_properties, add_provenance):  

        self.file_path = filepath
        self.ensemble_to_entrez_path = ensemble_to_entrez_path
        self.dataset = 'coxpresdb'
        self.label = 'coexpressed_with'
        self.source = 'CoXPresdb'
        self.source_url = 'https://coxpresdb.jp/'
        self.version = 'v8'

        assert os.path.isdir(self.file_path), "coxpresdb file path is not a directory"
        super(CoxpresdbAdapter, self).__init__(write_properties, add_provenance)

    def get_edges(self):
        # entrez_to_ensembl.pkl is generated using those two files:
        # gencode file: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz
        # Homo_sapiens.gene_info.gz file: https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
        # every gene has ensembl id in gencode file, every gene has hgnc id if available.
        # every gene has entrez gene id in gene_info file, every gene has ensembl id or hgcn id if available

        gene_ids = [f for f in os.listdir(self.file_path) if os.path.isfile(os.path.join(self.file_path, f))]

        with open(self.ensemble_to_entrez_path, 'rb') as f:
            entrez_ensembl_dict = pickle.load(f)
        for gene_id in gene_ids:
            gene_file_path = os.path.join(self.file_path, gene_id)
            entrez_id = gene_id
            ensembl_id = entrez_ensembl_dict.get(entrez_id)
            if ensembl_id:
                with open(gene_file_path, 'r') as input:
                    for line in input:
                        (co_entrez_id, score) = line.strip().split()
                        co_ensembl_id = entrez_ensembl_dict.get(co_entrez_id)
                        if co_ensembl_id:
                            _id = entrez_id + '_' + co_entrez_id + '_' + self.label
                            source = ensembl_id
                            target = co_ensembl_id
                            _props = {}
                            if self.write_properties:
                                _props['score'] = float(score)
                                if self.add_provenance:
                                    _props['source'] = self.source
                                    _props['source_url'] = self.source_url
                            yield source, target, self.label, _props

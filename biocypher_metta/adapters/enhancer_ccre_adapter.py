from biocypher_metta.adapters import Adapter

class EnhancerCCREAdapter(Adapter):
    def __init__(self, filepath, label):
        self.filepath = filepath
        self.label = label
        self.source = "EnhancerCCRE"
        self.source_url = "URL to dataset source"
        super(EnhancerCCREAdapter, self).__init__()

    def get_nodes(self):
        with open(self.filepath, 'r') as file:
            for line in file:
                if line.startswith("chr"):
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    score = float(fields[4])
                    gene_data = fields[3]

                    properties = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'score': score,
                        'gene_data': gene_data
                    }

                    node_id = f"{chrom}_{start}_{end}"

                    yield node_id, self.label, properties

    
    def get_edges(self):
        with open(self.filepath, 'r') as file:
            for line in file:
                if line.startswith("chr"):
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    gene_data = fields[3]

                    if gene_data != ".":
                        genes = gene_data.split(",")
                        for gene_info in genes:
                            if "_" in gene_info:
                                gene_id, _ = gene_info.split("_")
                            else:
                                gene_id = gene_info

                            enhancer_id = f"{chrom}_{start}_{end}"
                            yield enhancer_id, gene_id, "regulates", {}

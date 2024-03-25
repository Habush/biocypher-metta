from biocypher_metta.adapters import Adapter

# Example data from encodeCcreCombined.bed file:
# chr1	3119617	3119911	EM10E0431203	303	.	3119617	3119911	255,205,0	dELS	dELS	3.03931270232	enhD	E0431203	EM10E0431203 distal enhancer-like signature
# chr1	3119914	3120120	EM10E0431204	264	.	3119914	3120120	255,205,0	dELS	dELS	2.64144130291	enhD	E0431204	EM10E0431204 distal enhancer-like signature
# chr1	3120346	3120662	EM10E0431205	263	.	3120346	3120662	255,205,0	dELS	dELS	2.63858128228	enhD	E0431205	EM10E0431205 distal enhancer-like signature
# chr1	3292622	3292971	EM10E0431207	249	.	3292622	3292971	255,205,0	dELS	dELS	2.495384694	enhD	E0431207	EM10E0431207 distal enhancer-like signature


class EnhancerCCREAdapter(Adapter):
    def __init__(self, filepath, label):
        self.filepath = filepath
        self.label = label
        self.source = "EnhancerCCRE"
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

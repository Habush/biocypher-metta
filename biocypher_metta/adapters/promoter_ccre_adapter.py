from biocypher_metta.adapters import Adapter
import gzip

# Example data from the dataset encodeCcreCombined.bed:
# chr1	3119617	3119911	EM10E0431203	303	.	3119617	3119911	255,205,0	dELS	dELS	3.03931270232	enhD	E0431203	EM10E0431203 distal enhancer-like signature
# chr1	3119914	3120120	EM10E0431204	264	.	3119914	3120120	255,205,0	dELS	dELS	2.64144130291	enhD	E0431204	EM10E0431204 distal enhancer-like signature
# chr1	3120346	3120662	EM10E0431205	263	.	3120346	3120662	255,205,0	dELS	dELS	2.63858128228	enhD	E0431205	EM10E0431205 distal enhancer-like signature
# chr1	3292622	3292971	EM10E0431207	249	.	3292622	3292971	255,205,0	dELS	dELS	2.495384694	enhD	E0431207	EM10E0431207 distal enhancer-like signature

class PromoterCCREAdapter(Adapter):
    def __init__(self, filepath, write_properties, add_provenance, label="promoter", type="promoter"):
        self.filepath = filepath
        self.label = label
        self.type = type  
        self.source = "PromoterCCRE"
        self.source_url = "https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=mm10&g=encodeCcreCombined"

        super(PromoterCCREAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as file:
            for line in file:
                if line.startswith("chr"):
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    start = int(fields[1]) + 1
                    end = int(fields[2]) + 1

                    
                    try:
                        score = float(fields[4])
                        gene_data = fields[3]
                    except (IndexError, ValueError):
                        continue

                    region_type = fields[10]  
                

                    if region_type == "PLS":
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
        with gzip.open(self.filepath, 'rt') as file:
            for line in file:
                if line.startswith("chr"):
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    start = int(fields[1]) + 1
                    end = int(fields[2]) + 1 
                    gene_data = fields[3]
                    region_type = fields[10]  


                    if gene_data == ".":
                        continue

                    if region_type == "PLS":
                        promoter_id = f"{chrom}_{start}_{end}"
                        gene_id = gene_data.split("_")[0] if "_" in gene_data else gene_data
                        yield promoter_id, gene_id, self.label, {}
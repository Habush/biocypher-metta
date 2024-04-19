from biocypher_metta.adapters import Adapter
import gzip

# Example data from the dataset encodeCcreCombined.bed:
# chrom start end gene score strand chromStart chromEnd itemRgb signatures category signalValue type accessionId description
# chr1	3120346	3120662	EM10E0431205	263	.	3120346	3120662	255,205,0	dELS	dELS	2.63858128228	enhD	E0431205	EM10E0431205 distal enhancer-like signature
# chr1	3772420	3772769	EH38E1313009	595	.	3772420	3772769	255,0,0	PLS	PLS	5.957252591	prom	E1313009	EH38E1313009 promoter-like signature
# chr1	3772821	3773112	EH38E1313010	368	.	3772821	3773112	255,0,0	PLS,CTCF-bound	PLS	3.681354874	prom	E1313010	EH38E1313010 promoter-like signature
# chr1	3796383	3796733	EH38E1313041	536	.	3796383	3796733	255,0,0	PLS,CTCF-bound	PLS	5.367484968	prom	E1313041	EH38E1313041 promoter-like signature


class PromoterCCREAdapter(Adapter):
    def __init__(self, filepath, write_properties, add_provenance, label="promoter"):
        self.filepath = filepath
        self.label = label
        self.source = "PromoterCCRE"
        self.source_url = "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=regulation&hgta_track=encodeCcreCombined&hgta_table=encodeCcreCombined&hgta_doSchema=describe+table+schema"
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
                    except (IndexError, ValueError):
                        continue

                    region_type = fields[10]  
                

                    if region_type == "PLS":
                        # This block of code is executed if the region_type is "PLS" (promoter-Like Signature)
                        props = {}

                        if self.write_properties:
                            props['chr'] = chrom
                            props['start'] = start
                            props['end'] = end
                            props['score'] = score

                            if self.add_provenance:
                                props['source'] = self.source
                                props['source_url'] = self.source_url

                        node_id = f"{chrom}_{start}_{end}"
                        yield node_id, self.label, props

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

                    props = {}                                
                    if region_type == "PLS":
                    # This block of code is executed if the region_type is "PLS" (promoter-Like Signature)
                        if self.write_properties:
                            props['chr'] = chrom
                            props['start'] = start
                            props['end'] = end
                            props['gene'] = gene_data

                            if self.add_provenance:
                                props['source'] = self.source
                                props['source_url'] = self.source_url
                        promoter_id = f"{chrom}_{start}_{end}"
                        gene_id = gene_data.split("_")[0] if "_" in gene_data else gene_data
                        yield promoter_id, gene_id, self.label, {}
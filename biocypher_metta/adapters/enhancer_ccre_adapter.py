import gzip
from biocypher_metta.adapters import Adapter

# Example data from the dataset encodeCcreCombined.bed:
# chrom start end gene score strand chromStart chromEnd itemRgb signatures category signalValue type accessionId description
# chr1	181251	181601	EH38E1310153	488	.	181251	181601	255,167,0	pELS,CTCF-bound	pELS	4.88403259983	enhP	E1310153	EH38E1310153 proximal enhancer-like signature
# chr1	190865	191071	EH38E1310154	179	.	190865	191071	255,205,0	dELS,CTCF-bound	dELS	1.79282201562	enhD	E1310154	EH38E1310154 distal enhancer-like signature
# chr1	778562	778912	EH38E1310158	759	.	778562	778912	255,0,0	PLS,CTCF-bound	PLS	7.59852335807	prom	E1310158	EH38E1310158 promoter-like signature
# chr1	779086	779355	EH38E1310159	304	.	779086	779355	255,0,0	PLS,CTCF-bound	PLS	3.04663905401	prom	E1310159	EH38E1310159 promoter-like signature
# chr1	779727	780060	EH38E1310160	281	.	779727	780060	255,167,0	pELS,CTCF-bound	pELS	2.81743397307	enhP	E1310160	EH38E1310160 proximal enhancer-like signature

class EnhancerCCREAdapter(Adapter):
    def __init__(self, filepath,write_properties, add_provenance, label='enhancer'):
        self.filepath = filepath
        self.label = label
        self.source = "EnhancerCCRE"
        self.source_url = "https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=mm10&g=encodeCcreCombined"

        super(EnhancerCCREAdapter, self).__init__(write_properties, add_provenance)
    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as file:
            for line in file:
                if line.startswith("chr"):
                    fields = line.strip().split("\t")
                    chrom = fields[0]
                    start = int(fields[1]) + 1
                    end = int(fields[2])  + 1

                    try:
                        score = float(fields[4])
                    except (IndexError, ValueError):
                        continue

                    region_type = fields[10]  

                    if region_type in ["dELS", "pELS"]:
                        # This block of code is executed if the region_type is either "dELS" (distal Enhancer-Like Signature) or "pELS" (proximal Enhancer-Like Signature)
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
                    end = int(fields[2])  + 1
                    gene_data = fields[3]
                    region_type = fields[10]  

                    if gene_data == ".":
                        continue

                    props = {}
                    if region_type in ["dELS", "pELS"]:
                        if self.write_properties:
                            props['chr'] = chrom
                            props['start'] = start
                            props['end'] = end
                            props['gene'] = gene_data

                            if self.add_provenance:
                                props['source'] = self.source
                                props['source_url'] = self.source_url
                        
                        enhancer_id = f"{chrom}_{start}_{end}"
                        gene_id = gene_data.split("_")[0] if "_" in gene_data else gene_data
                        yield enhancer_id, gene_id, self.label, {}


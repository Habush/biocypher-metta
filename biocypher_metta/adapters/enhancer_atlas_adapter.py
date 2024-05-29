import gzip
import os
import pickle
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import build_regulatory_region_id, check_genomic_location

# Example enhancer atlas input file:

# enhancer_file
# chr1	875310	876520	11.4349429027718
# chr1	876690	877510	10.6950350518485
# chr1	902000	902080	10.429771917918
# chr1	902180	902820	10.3902750575485

# enhancer_gene_file
# chr1:874840-876520_ENSG00000225880$LINC00115$chr1$762902$-	1.104330
# chr1:876690-877510_ENSG00000225880$LINC00115$chr1$762902$-	1.264920
# chr1:902000-902080_ENSG00000225880$LINC00115$chr1$762902$-	1.663935
# chr1:955750-957580_ENSG00000225880$LINC00115$chr1$762902$-	1.512501
# chr1:874840-876520_ENSG00000237330$RNF223$chr1$1009687$-	1.065467


class EnhancerAtlasAdapter(Adapter):
    INDEX = {'chr': 0, 'coord_start': 1, 'coord_end': 2, 'snp': 7}

    def __init__(self, enhancer_filepath, enhancer_gene_filepath, tissue_to_ontology_filepath, 
                 write_properties, add_provenance, 
                 type='enhancer', label='enhancer',
                 chr=None, start=None, end=None):
        self.enhancer_filepath = enhancer_filepath
        self.enhancer_gene_filepath = enhancer_gene_filepath
        self.tissue_to_ontology_filepath = tissue_to_ontology_filepath
        self.chr = chr
        self.start = start
        self.end = end
        self.label = label
        self.type = type

        self.source = 'Enancer Atlas'
        self.version = '2.0'
        self.source_url = 'http://enhanceratlas.org/downloadv2.php'

        super(EnhancerAtlasAdapter, self).__init__(write_properties, add_provenance)
    
    def parse_enhancer_gene(self, info):
        enhancer_info = info.split('_')[0]
        chr = enhancer_info.split(':')[0]
        start = int(enhancer_info.split(':')[1].split('-')[0]) + 1 #+1 since it is 0-based genomic coordinate
        end = int(enhancer_info.split(':')[1].split('-')[1]) + 1
        gene = info.split('_')[1].split('$')[0]
        return chr, start, end, gene
    
    def get_nodes(self):
        with gzip.open(self.enhancer_filepath, 'rt') as f:
            for line in f:
                info = line.strip().split('\t')
                chr = info[EnhancerAtlasAdapter.INDEX['chr']]
                start = int(info[EnhancerAtlasAdapter.INDEX['coord_start']]) + 1 #+1 since it is 0-based genomic coordinate
                end = int(info[EnhancerAtlasAdapter.INDEX['coord_end']]) + 1
                enhancer_region_id = build_regulatory_region_id(chr, start, end)
                
                if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    props = {}
                    if self.write_properties:
                        props['chr'] = chr
                        props['start'] = start
                        props['end'] = end
                    
                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url

                    yield enhancer_region_id, self.label, props 

    def get_edges(self):
        tissues = [f for f in os.listdir(self.enhancer_gene_filepath) if os.path.isfile(os.path.join(self.enhancer_gene_filepath, f))]
        with open(self.tissue_to_ontology_filepath, 'rb') as f:
            tissues_ontology_map = pickle.load(f)
        
        for tissue in tissues:
            tissue_file_path = os.path.join(self.enhancer_gene_filepath, tissue)
            biological_context = tissues_ontology_map.get(tissue.replace('_EP.txt', ''))
            if biological_context:
                with open(tissue_file_path, 'r') as f:
                    for line in f:
                        info = line.strip().split('\t')
                        chr, start, end, gene = self.parse_enhancer_gene(line)
                        if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                            enhancer_region_id = build_regulatory_region_id(chr, start, end)
                            score = float(info[1])
                            props = {}
                            if self.write_properties:
                                props['biological_context'] = biological_context
                                props['score'] = score
                                if self.add_provenance:
                                    props['source'] = self.source
                                    props['source_url'] = self.source_url

                            yield enhancer_region_id, gene, self.label, props

from collections import defaultdict
import gzip
from biocypher_metta.adapters import Adapter

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

# enhancer_snps_file
# chr1	2165300	2167010	11.464638265294	chr1	2166404	2166408	rs61776614
# chr1	2165300	2167010	11.464638265294	chr1	2166807	2166811	rs76875252
# chr1	7530490	7533550	16.9883353758835	chr1	7532288	7532292	rs1149336
# chr1	8021850	8022240	13.3539912947246	chr1	8021971	8021975	rs35675666

class EnhancerAtlasAdapter(Adapter):
    INDEX = {'chr': 0, 'coord_start': 1, 'coord_end': 2, 'snp': 7}

    def __init__(self, filepath, enhancer_gene_filepath, enhancer_snps_filepath, write_properties, add_provenance, type='enhancer', label='enhancer'):
        self.filepath = filepath
        self.enhancer_gene_filepath = enhancer_gene_filepath
        self.enhancer_snps_filepath = enhancer_snps_filepath
        self.label = label
        self.type = type

        self.source = 'Enancer Atlas'
        self.version = '2.0'
        self.source_url = 'http://enhanceratlas.org/downloadv2.php'

        super(EnhancerAtlasAdapter, self).__init__(write_properties, add_provenance)
    
    def get_enhancer_id(self, chr, start, end):
        return chr+':'+start+'-'+end

    def get_nodes(self):
        enhancer_gene = defaultdict(list)
        with gzip.open(self.enhancer_gene_filepath, 'rt') as f:
            for line in f:
                info = line.strip().split('_')
                enhancer = info[0]
                gene = info[1].split('$')[0]
                enhancer_gene[enhancer].append(gene)
        
        enhancer_snps = defaultdict(list)
        with gzip.open(self.enhancer_snps_filepath, 'rt') as f:
            for line in f:
                info = line.strip().split('\t')
                chr = info[EnhancerAtlasAdapter.INDEX['chr']]
                start = info[EnhancerAtlasAdapter.INDEX['coord_start']]
                end = info[EnhancerAtlasAdapter.INDEX['coord_end']]
                enhancer = self.get_enhancer_id(chr, start, end)
                snp = info[EnhancerAtlasAdapter.INDEX['snp']]
                enhancer_snps[enhancer].append(snp)
        
        with gzip.open(self.filepath, 'rt') as f:
            for line in f:
                info = line.strip().split('\t')
                chr = info[EnhancerAtlasAdapter.INDEX['chr']]
                start = info[EnhancerAtlasAdapter.INDEX['coord_start']]
                end = info[EnhancerAtlasAdapter.INDEX['coord_end']]
                enhancer = self.get_enhancer_id(chr, start, end)
                genes = enhancer_gene.get(enhancer, [])
                snps = enhancer_snps.get(enhancer, [])
                
                props = {}
                if self.write_properties:
                    props['chr'] = chr
                    props['start'] = start
                    props['end'] = end
                    if genes:
                        props['genes'] = genes
                    if snps:
                        props['snps'] = snps
                
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield enhancer, self.label, props 

    def get_edges(self):
        with gzip.open(self.enhancer_gene_filepath, 'rt') as f:
            for line in f:
                info = line.strip().split('\t')
                enhancer = info[0].split('_')[0]
                gene = info[0].split('_')[1].split('$')[0]
                score = info[1]
                props = {}
                if self.write_properties:
                    props['score'] = score
                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url

                yield enhancer, gene, self.label, props

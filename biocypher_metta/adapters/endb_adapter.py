import gzip
import pickle

from biocypher_metta.adapters import Adapter
from liftover import get_lifter
# Example ENdb csv input file:
# Enhancer_id	Pubmed	Enhancer_symbol	Reference_genome	Chromosome	Start_position	End_position	Species	Biosample_name	Experiment_class	Enhancer_experiment	Enhancer_tar_ex_De	Enhancer_type	Target_gene	target_gene_strong_experiment	target_gene_weak_experiment	target_gene_experiment_description	target_gene_other_name	Disease	DO	mesh	Enhancer_function	Enhancer_function_experiment	En_function_ex_de	TF_name	TF_other_name	Experiment	TF_experiment_de	SNP_id	SNP_position	SNP_experiment	Distance_from_TSS	Distance_from_TSS_AB	enhancer_MID	TSS	n_Distance_from_TSS
# E_01_001	28652246	--	hg19	chr6	86160197	86176777	Human	Melanoma	Low+High throughput	"ChIP-seq,ChIP-qPCR"	A ChIP-qPCR approach revealed prominent binding of c-Jun at the AP-1 sites #1-2 and #5-6 (Fig. 3E) coinciding with enhancer regions (H3K27Ac) identified by the ENCODE ChIP-seq project (Fig. 3D).	Enhancer	NT5E	CRISPR/Cas9	"qPCR,ChIP-qPCR"	"Next, we devised a CRISPR/Cas9 approach to clarify the functional relevance of the individual c-Jun/AP-1 binding sites.c-Jun/AP-1 induces CD73 expression via binding to an intronic enhancer.Strong enrichment of AP-1 site mutagenesis frequencies in the CD73 low subfractions (52% versus 4% in CD73 lowest 10% versus highest 10%) was observed only in the cell population transfected with the sgRNA against the c-Jun/AP-1 site #5 (Fig. 3G-H) being consistent with reduced CD73 surface expression (Fig. S4A)."	"CALJA,CD73,E5NT,NT,NT5,NTE,eN,eNT"	Melanoma	DOID:1909	D008545	The intronic Enhancer of the CD73 gene controls its c-Jun/AP-1-dependent transcriptional activation downstream of mitogenic MAPK and inflammatory cytokine signaling. 	"ChIP-qPCR,CRISPR/Cas9"	"Taken together, we identified a critical c-Jun/AP-1 binding site in an intronic enhancer of the CD73 gene that controls its c-Jun/AP-1-dependent transcriptional activation downstream of mitogenic MAPK and inflammatory cytokine signaling."	JUN	"AP-1,AP1,c-Jun,cJUN,p39"	ChIP-qPCR	We looked for potential c-Jun/AP-1�Cbinding sites in the CD73 gene using ENCODE ChIP-Seq data.A ChIP-qPCR approach revealed prominent binding of c-Jun at the AP-1 sites #1�C2 and #5�C6 coinciding with Enhancer regions (H3K27Ac)identified by the ENCODE ChIP-seq project.	--	--	--	Intron				9186
# E_01_002	26751173	ER��enh588	hg19	chr11	69329684	69330954	Human	"MCF-7,T-47D"	Low+High throughput	"ChIP-seq,CRISPR/Cas9,Luciferase Reporter Assay,Western blot,Competitive Proliferation Assay"	"To demonstrate the generalizability of our screening approach, we designed a dropout screen to identify novel ER��-bound enhancers. 
# We started by cloning ER��enh588 WT region in a reporter vector and verified that it has strong ER���Cdependent transcription-enhancing activity, since mutations in the ER�� binding site completely abolish the enhancer activity and response to 17��-estradiol (Supplementary Fig. 4b). 

class ENdbAdapter(Adapter):
    ALLOWED_TYPES = ['enhancer']
    ALLOWED_LABELS = ['enhancer_regulates']
    INDEX = {'enhancer_id': 0, 'pubmed': 1,  'chr': 4, 'coord_start': 5, 'coord_end': 6, 'species' : 7, 'enhancer_type': 12, 'gene_id': 13, 'strong_experiment': 14, 'disease': 18, 'tf_name': 24, 'snp_id': 28}
    
    def __init__(self, filepath, hgnc_to_ensembl_map, type='enhancer', label='enhancer', delimiter='\t'):
        self.filepath = filepath
        self.hgnc_to_ensembl_map = pickle.load(open(hgnc_to_ensembl_map, 'rb'))
        self.type = type
        self.label = label 
        self.delimiter = delimiter
        self.genome_reference_converter = get_lifter('hg19', 'hg38')

        self.source = 'ENdb'
        self.version = ''
        self.source_url = 'https://bio.liclab.net/ENdb/Download.php'

    def parse_info_metadata(self, info):
        parsed_info = info.replace('"', '').split(',')
        return parsed_info
    
    def handle_target_genes(self, genes):
        parsed_genes = self.parse_info_metadata(genes)
        result = []
        for gene in parsed_genes:
            if gene in self.hgnc_to_ensembl_map:
                result.append(self.hgnc_to_ensembl_map[gene])
            elif gene.replace('-', '') in self.hgnc_to_ensembl_map:
                result.append(self.hgnc_to_ensembl_map[gene.replace('-', '')])
            elif gene[:-1] in self.hgnc_to_ensembl_map:
                result.append(self.hgnc_to_ensembl_map[gene[:-1]])
        if not result and parsed_genes:
            result.append(parsed_genes[0])
            
        return result
    
    def convert_to_hg38(self, chr,  pos):
        try:
            chr_no = chr.replace('chr', '').replace('ch', '')
            converted = self.genome_reference_converter.query(chr_no, pos)[0][1]
            return converted
        except:
            return False
    
    def get_nodes(self):
        self.temp = []
        with gzip.open(self.filepath, 'rt') as f:
            next(f) # Skip the header row
            for line in f:
                cur_data = line.strip().split(self.delimiter)
                first_val = cur_data[0]
                
                if first_val.startswith('E_'):
                    if self.temp:
                        data = self.temp.copy()
                    else:
                        self.temp = cur_data.copy()
                        continue
                else:
                    curr = cur_data.pop(0)
                    self.temp[-1] = self.temp[-1] + ' ' + curr
                    self.temp = self.temp + cur_data
                    continue

                species = data[ENdbAdapter.INDEX['species']]
                self.temp = cur_data
                
                if species != 'Human':
                    continue

                enhancer_id = data[ENdbAdapter.INDEX['enhancer_id']]
                chr = data[ENdbAdapter.INDEX['chr']]
                start_position = self.convert_to_hg38(chr, int(data[ENdbAdapter.INDEX['coord_start']]))
                end_position = self.convert_to_hg38(chr, int(data[ENdbAdapter.INDEX['coord_end']]))
                if not start_position or not end_position:
                    continue
                target_genes = data[ENdbAdapter.INDEX['gene_id']]
                enhancer_type = data[ENdbAdapter.INDEX['enhancer_type']]
                tf_name = data[ENdbAdapter.INDEX['tf_name']]
                snp_id = data[ENdbAdapter.INDEX['snp_id']]
                # disease = data[ENdbAdapter.INDEX['disease']]
                

                props = {
                    'chr': chr,
                    'start': start_position,
                    'end': end_position,
                    'enhancer_type': enhancer_type,
                    'genes': self.handle_target_genes(target_genes.upper()),
                    # Add other properties as needed based on your specific data interests
                }

                if tf_name != '--':
                    props['tf_name'] = self.parse_info_metadata(tf_name)
                if snp_id != '--':
                    props['snp_id'] = self.parse_info_metadata(snp_id)

                yield enhancer_id, self.label, props

    def get_edges(self):
        self.temp = []
        with gzip.open(self.filepath, 'rt') as f:
            next(f)
            for line in f:
                cur_data = line.strip().split(self.delimiter)
                first_val = cur_data[0]
                
                if first_val.startswith('E_'):
                    if self.temp:
                        data = self.temp.copy()
                    else:
                        self.temp = cur_data.copy()
                        continue
                else:
                    curr = cur_data.pop(0)
                    self.temp[-1] = self.temp[-1] + ' ' + curr
                    self.temp = self.temp + cur_data
                    continue

                species = data[ENdbAdapter.INDEX['species']]
                self.temp = cur_data

                if species != 'Human':
                    continue
                chr = data[ENdbAdapter.INDEX['chr']]
                start_position = self.convert_to_hg38(chr, int(data[ENdbAdapter.INDEX['coord_start']]))
                end_position = self.convert_to_hg38(chr, int(data[ENdbAdapter.INDEX['coord_end']]))
                if not start_position or not end_position:
                    continue

                enhancer_id = data[ENdbAdapter.INDEX['enhancer_id']]
                pubmed = 'pubmed:' + data[ENdbAdapter.INDEX['pubmed']]
                target_genes = data[ENdbAdapter.INDEX['gene_id']]
                target_genes_list = self.handle_target_genes(target_genes.upper())
                experiment = data[ENdbAdapter.INDEX['strong_experiment']]
                

                props = {'evidence': pubmed}
                if experiment != '--':
                    props['experiment'] = self.parse_info_metadata(experiment)

                for gene in target_genes_list:
                    yield enhancer_id, gene, self.label, props


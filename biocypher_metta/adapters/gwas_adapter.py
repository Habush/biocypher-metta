import csv
import gzip
from biocypher_metta.adapters import Adapter
# Example GWAS tsv input file:
# DATE ADDED TO CATALOG	PUBMEDID	FIRST AUTHOR	DATE	JOURNAL	LINK	STUDY	DISEASE/TRAIT	INITIAL SAMPLE SIZE	REPLICATION SAMPLE SIZE	REGION	CHR_ID	CHR_POS	REPORTED GENE(S)	MAPPED_GENE	UPSTREAM_GENE_ID	DOWNSTREAM_GENE_ID	SNP_GENE_IDS	UPSTREAM_GENE_DISTANCE	DOWNSTREAM_GENE_DISTANCE	STRONGEST SNP-RISK ALLELE	SNPS	MERGED	SNP_ID_CURRENT	CONTEXT	INTERGENIC	RISK ALLELE FREQUENCY	P-VALUE	PVALUE_MLOG	P-VALUE (TEXT)	OR or BETA	95% CI (TEXT)	PLATFORM [SNPS PASSING QC]	CNV
# 2019-11-04	30239722	Pulit SL	2018-09-14	Hum Mol Genet	www.ncbi.nlm.nih.gov/pubmed/30239722	Meta-analysis of genome-wide association studies for body fat distribution in 694,649 individuals of European ancestry.	Waist-to-hip ratio adjusted for BMI	315,284 European ancestry male individuals	NA	6p25.1	6	6745933	NR	LY86 - BTF3P7	ENSG00000112799	ENSG00000219986		90950	48982	rs1294436-C	rs1294436	0	1294436	intron_variant	1	0.5973	9E-48	47.045757490560675		0.0293	[0.025-0.033] unit increase	NR [~ 27400000] (imputed)	N
# 2019-11-04	30239722	Pulit SL	2018-09-14	Hum Mol Genet	www.ncbi.nlm.nih.gov/pubmed/30239722	Meta-analysis of genome-wide association studies for body fat distribution in 694,649 individuals of European ancestry.	Waist-to-hip ratio adjusted for BMI	315,284 European ancestry male individuals	NA	2p23.3	2	25263901	NR	DNMT3A			ENSG00000119772			rs12991495-T	rs12991495	0	12991495	intron_variant	0	0.6906	5E-13	12.301029995663981		0.0135	[0.0098-0.0172] unit increase	NR [~ 27400000] (imputed)	N
# 2019-11-04	30239722	Pulit SL	2018-09-14	Hum Mol Genet	www.ncbi.nlm.nih.gov/pubmed/30239722	Meta-analysis of genome-wide association studies for body fat distribution in 694,649 individuals of European ancestry.	Waist-to-hip ratio adjusted for BMI	315,284 European ancestry male individuals	NA	4q24	4	102267552	NR	SLC39A8			ENSG00000138821			rs13107325-T	rs13107325	0	13107325	missense_variant	0	0.0808	4E-19	18.397940008672037		0.0306	[0.024-0.037] unit decrease	NR [~ 27400000] (imputed)	N

class GWASAdapter(Adapter):
    INDEX = {'pubmed': 1, 'chr': 11, 'pos': 12, 'upstream_gene': 16, 'downstream_gene': 17, 'snp_gene': 18, 'snp': 21}

    def __init__(self, filepath, write_properties, add_provenance, label='snp'):
        self.filepath = filepath
        self.label = label

        self.source = 'GWAS'
        self.version = 'v1.0'
        self.source_url = 'https://www.ebi.ac.uk/gwas/docs/file-downloads'

        super(GWASAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as f:
            next(f)
            reader = csv.reader(f, delimiter='\t')
            
            for line in reader:
                snp_id = line[GWASAdapter.INDEX['snp']]
                if ',' in snp_id:
                    snp_id_list = snp_id.split(',')
                else:
                    snp_id_list = snp_id.split(';')
                chr = line[GWASAdapter.INDEX['chr']].split(';')
                pos = line[GWASAdapter.INDEX['pos']].split(';')
                if not pos:
                    continue
                for snp in snp_id_list:
                    for i in range(len(pos)):
                        try:    
                            if not pos[i]:
                                continue
                            props = {}
                            if self.write_properties:
                                props['chr'] = 'chr'+chr[i]
                                props['start'] = int(pos[i])
                                props['end'] = int(pos[i]) + 1
                                
                                if self.add_provenance:
                                    props['source'] = self.source
                                    props['source_url'] = self.source_url

                            yield snp.strip(), self.label, props
                        except:
                            pass

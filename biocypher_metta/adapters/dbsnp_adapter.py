import gzip
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import check_genomic_location
# Exaple dbSNP vcf input file:
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
# 1	10177	rs367896724	A	AC	.	.	RS=367896724;RSPOS=10177;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5747,0.4253;COMMON=1;TOPMED=0.76728147298674821,0.23271852701325178
# 1	10352	rs555500075	T	TA	.	.	RS=555500075;RSPOS=10352;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5625,0.4375;COMMON=1;TOPMED=0.86356396534148827,0.13643603465851172
# 1	10616	rs376342519	CCGCCGTTGCAAAGGCGCGCCG	C	.	.	RS=376342519;RSPOS=10617;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005040026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;KGPhase3;CAF=0.006989,0.993;COMMON=1

class DBSNPAdapter(Adapter):
    INDEX = {'chr': 0, 'pos': 1, 'id': 2, 'ref': 3, 'alt': 4, 'info': 7}
    def __init__(self, filepath, write_properties, add_provenance,
                 chr=None, start=None, end=None):
        self.filepath = filepath
        self.chr = chr
        self.start = start
        self.end = end
        self.label = 'snp'

        self.source = 'dbSNP'
        self.version = '2.0'
        self.source_url = 'https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/'
        super(DBSNPAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                data = line.strip().split('\t')
                rsid = data[DBSNPAdapter.INDEX['id']]
                chr = data[DBSNPAdapter.INDEX['chr']]
                pos = int(data[DBSNPAdapter.INDEX['pos']])
                ref = data[DBSNPAdapter.INDEX['ref']]
                alt = data[DBSNPAdapter.INDEX['alt']]

                if check_genomic_location(self.chr, self.start, self.end, chr, pos, pos):
                    props = {}
                    if self.write_properties:
                        props['chr'] = 'chr'+chr
                        props['pos'] = pos
                        props['ref'] = ref
                        props['alt'] = alt

                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url
                    
                    yield rsid, self.label, props
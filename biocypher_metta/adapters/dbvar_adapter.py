import gzip
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import check_genomic_location
# Example dbVar input file:
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
# 1	10000	nssv16889290	N	<DUP>	.	.	DBVARID=nssv16889290;SVTYPE=DUP;END=52000;SVLEN=42001;EXPERIMENT=1;SAMPLESET=1;REGIONID=nsv6138160;AC=1453;AF=0.241208;AN=6026
# 1	10001	nssv14768	T	<DUP>	.	.	DBVARID=nssv14768;SVTYPE=DUP;IMPRECISE;END=88143;CIPOS=0,0;CIEND=0,0;SVLEN=78143;EXPERIMENT=1;SAMPLE=NA12155;REGIONID=nsv7879
# 1	10001	nssv14781	T	<DUP>	.	.	DBVARID=nssv14781;SVTYPE=DUP;IMPRECISE;END=82189;CIPOS=0,0;CIEND=0,0;SVLEN=72189;EXPERIMENT=1;SAMPLE=NA18860;REGIONID=nsv7879

class DBVarVariantAdapter(Adapter):
    INDEX = {'chr': 0, 'coord_start': 1, 'id': 2, 'type': 4, 'info': 7}
    VARIANT_TYPES = {'<CNV>': 'copy number variation', '<DEL>': 'deletion', '<DUP>': 'duplication', '<INS>': 'insertion', '<INV>': 'inversion'}

    def __init__(self, filepath, write_properties, add_provenance, 
                 label='structural_variant', delimiter='\t',
                 chr=None, start=None, end=None):
        self.filepath = filepath
        self.delimiter = delimiter
        self.label = label
        self.chr = chr
        self.start = start
        self.end = end

        self.source = 'dbVar'
        self.version = ''
        self.source_url = 'https://www.ncbi.nlm.nih.gov/dbvar/content/ftp_manifest/'

        super(DBVarVariantAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                data = line.strip().split(self.delimiter)
                variant_id = data[DBVarVariantAdapter.INDEX['id']]
                variant_type_key = data[DBVarVariantAdapter.INDEX['type']]
                if variant_type_key not in DBVarVariantAdapter.VARIANT_TYPES:
                    continue
                variant_type = DBVarVariantAdapter.VARIANT_TYPES[variant_type_key]
                chr = 'chr' + data[DBVarVariantAdapter.INDEX['chr']]
                start = int(data[DBVarVariantAdapter.INDEX['coord_start']])
                info = data[DBVarVariantAdapter.INDEX['info']].split(';')
                end = start
                for i in range(len(info)):
                    if info[i].startswith('END='):
                        end = int(info[i].split('=')[1])
                        break
                
                if check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    props = {}

                    if self.write_properties:
                        props['chr'] = chr
                        props['start'] = start
                        props['end'] = end
                        props['variant_type'] = variant_type

                        if self.add_provenance:
                            props['source'] = self.source
                            props['source_url'] = self.source_url


                    yield variant_id, self.label, props

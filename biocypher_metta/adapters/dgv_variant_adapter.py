import gzip
from biocypher_metta.adapters import Adapter
from biocypher_metta.adapters.helpers import check_genomic_location
# Example dgv input file:
# variantaccession	chr	start	end	varianttype	variantsubtype	reference	pubmedid	method	platform	mergedvariants	supportingvariants	mergedorsample	frequency	samplesize	observedgains	observedlosses	cohortdescription	genes	samples
# dgv1n82	1	10001	22118	CNV	duplication	Sudmant_et_al_2013	23825009	Oligo aCGH,Sequencing			nsv945697,nsv945698	M		97	10	0		""	HGDP00456,HGDP00521,HGDP00542,HGDP00665,HGDP00778,HGDP00927,HGDP00998,HGDP01029,HGDP01284,HGDP01307
# nsv7879	1	10001	127330	CNV	gain+loss	Perry_et_al_2008	18304495	Oligo aCGH			nssv14786,nssv14785,nssv14773,nssv14772,nssv14781,nssv14771,nssv14775,nssv14762,nssv14764,nssv18103,nssv14766,nssv14770,nssv14777,nssv14789,nssv14782,nssv14788,nssv18117,nssv14790,nssv14791,nssv14784,nssv14776,nssv14787,nssv21423,nssv14783,nssv14763,nssv14780,nssv14774,nssv14768,nssv18113,nssv18093	M		31	25	1		""	NA07029,NA07048,NA10839,NA10863,NA12155,NA12802,NA12872,NA18502,NA18504,NA18517,NA18537,NA18552,NA18563,NA18853,NA18860,NA18942,NA18972,NA18975,NA18980,NA19007,NA19132,NA19144,NA19173,NA19221,NA19240
# nsv482937	1	10001	2368561	CNV	loss	Iafrate_et_al_2004	15286789	BAC aCGH,FISH			nssv2995976	M		39	0	1		""	

class DGVVariantAdapter(Adapter):
    INDEX = {'variant_accession': 0, 'chr': 1, 'coord_start': 2, 'coord_end': 3, 'type': 5, 'pubmedid': 7, 'genes': 17}

    def __init__(self, filepath, write_properties, add_provenance, 
                 label='structural_variant', delimiter='\t',
                 chr=None, start=None, end=None):
        self.filepath = filepath
        self.delimiter = delimiter
        self.label = label
        self.chr = chr
        self.start = start
        self.end = end

        self.source = 'dgv'
        self.version = ''
        self.source_url = 'http://dgv.tcag.ca/dgv/app/downloads'

        super(DGVVariantAdapter, self).__init__(write_properties, add_provenance)

    def get_nodes(self):
        with gzip.open(self.filepath, 'rt') as f:
            next(f)
            for line in f:
                data = line.strip().split(self.delimiter)
                variant_accession = data[DGVVariantAdapter.INDEX['variant_accession']]
                chr = data[DGVVariantAdapter.INDEX['chr']]
                start = int(data[DGVVariantAdapter.INDEX['coord_start']]) + 1
                end = int(data[DGVVariantAdapter.INDEX['coord_end']]) + 1
                variant_type = data[DGVVariantAdapter.INDEX['type']]
                pubmedid = data[DGVVariantAdapter.INDEX['pubmedid']]
                if not check_genomic_location(self.chr, self.start, self.end, chr, start, end):
                    continue
                props = {}

                if self.write_properties:
                    props['chr'] = 'chr' + chr
                    props['start'] = start
                    props['end'] = end
                    props['variant_type'] = variant_type
                    props['evidence'] = 'pubmed:'+pubmedid

                    if self.add_provenance:
                        props['source'] = self.source
                        props['source_url'] = self.source_url


                yield variant_accession, self.label, props

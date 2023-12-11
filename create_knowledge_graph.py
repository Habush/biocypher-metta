"""
Knowledge graph generation through BioCypher script
"""
from biocypher_metta.metta_writer import *
from biocypher._logger import logger
from biocypher_metta.adapters.uniprot_protein_adapter import UniprotProteinAdapter
from biocypher_metta.adapters.gencode_adapter import GencodeAdapter
from biocypher_metta.adapters.gencode_gene_adapter import GencodeGeneAdapter
from biocypher_metta.adapters.uniprot_adapter import UniprotAdapter
from biocypher_metta.adapters.reactome_pathway_adapter import ReactomePathwayAdapter
from biocypher_metta.adapters.reactome_adapter import ReactomeAdapter
from biocypher_metta.adapters.ontologies_adapter import OntologyAdapter
from biocypher_metta.adapters.gaf_adapter import GAFAdapter
from biocypher_metta.adapters.coxpresdb_adapter import CoxpresdbAdapter
from biocypher_metta.adapters.ccre_adapter import CCREAdapter
from biocypher_metta.adapters.encode_enhancer_gene_adapter import EncodeEnhancerGeneLinkAdapter
from biocypher_metta.adapters.tflink_adapter import TFLinkAdapter
from biocypher_metta.adapters.string_ppi_adapter import StringPPIAdapter
from biocypher_metta.adapters.tad_adapter import TADAdapter
from biocypher_metta.adapters.chromatin_state_adapter import ChromatinStateAdapter
from biocypher_metta.adapters.gtex_eqtl_adapter import GTExEQTLAdapter

ADAPTERS = {
    'gencode_gene': {'adapter': GencodeGeneAdapter(filepath='./samples/gencode_sample.gtf',
                                             gene_alias_file_path='./samples/Homo_sapiens.gene_info.gz'),
                     'outdir': 'gencode', 'nodes': True, 'edges': False},
    # 'gencode_transcripts': {
    #     'adapter': GencodeAdapter(filepath='./samples/gencode_sample.gtf', type='transcript',
    #                               label='transcript'),
    #     'outdir': 'gencode',
    #     'nodes': True,
    #     'edges': False
    # },
    #
    # 'transcribed_to': {
    #     'adapter': GencodeAdapter(filepath='./samples/gencode_sample.gtf', type='transcribed to',
    #                                  label='transcribed_to'),
    #     'outdir': 'gencode',
    #     'nodes': False,
    #     'edges': True
    # },
    # 'transcribed_from': {
    #     'adapter': GencodeAdapter(filepath='./samples/gencode_sample.gtf', type='transcribed from',
    #                               label='transcribed_from'),
    #     'outdir': 'gencode',
    #     'nodes': False,
    #     'edges': True
    # },
    # 'uniprotkb_sprot': {
    #     'adapter': UniprotProteinAdapter(filepath='./samples/uniprot_sprot_human_sample.dat.gz'),
    #     'outdir': 'uniprot',
    #     'nodes': True,
    #     'edges': False
    # },
    #
    # 'uniprotkb_translates_to': {
    #     'adapter': UniprotAdapter(filepath='./samples/uniprot_sprot_human_sample.dat.gz', type='translates to',
    #                               label='translates_to'),
    #     'outdir': 'uniprot',
    #     'nodes': False,
    #     'edges': True
    # },
    # 'uniprotkb_translation_of': {
    #     'adapter': UniprotAdapter(filepath='./samples/uniprot_sprot_human_sample.dat.gz', type='translation of',
    #                               label='translation_of'),
    #     'outdir': 'uniprot',
    #     'nodes': False,
    #     'edges': True
    # },
    #
    # 'pathway': {
    #     'adapter': ReactomePathwayAdapter('./samples/reactome/ReactomePathways.txt'),
    #     'outdir': 'pathway',
    #     'nodes': True,
    #     'edges': False
    # },
    #
    # 'genes_pathways': {
    #     'adapter': ReactomeAdapter('./samples/reactome/Ensembl2Reactome_All_Levels_sample.txt', 'genes_pathways'),
    #     'outdir': 'pathway',
    #     'nodes': False,
    #     'edges': True
    # },
    #
    # 'parent_pathway_of': {
    #     'adapter': ReactomeAdapter('./samples/reactome/ReactomePathwaysRelation.txt', 'parent_pathway_of'),
    #     'outdir': 'pathway',
    #     'nodes': False,
    #     'edges': True
    # },
    #
    # 'child_pathway_of': {
    #     'adapter': ReactomeAdapter('./samples/reactome/ReactomePathwaysRelation.txt', 'child_pathway_of'),
    #     'outdir': 'pathway',
    #     'nodes': False,
    #     'edges': True
    # },
    #
    # 'onotology_terms': {
    #     'adapter': OntologyAdapter(type='node', dry_run=True),
    #     'outdir': 'ontology',
    #     'nodes': True,
    #     'edges': False
    # },
    # 'onotology_relationships': {
    #     'adapter': OntologyAdapter(type='edge', dry_run=True),
    #     'outdir': 'ontology',
    #     'nodes': False,
    #     'edges': True
    # },

    # 'gaf': {
    #     'adapter': GAFAdapter(filepath='./samples/goa_human_sample.gaf.gz'),
    #     'outdir': 'gaf',
    #     'nodes': False,
    #     'edges': True
    # },
    #
    # 'coexpression': {
    #     'adapter': CoxpresdbAdapter(filepath='./samples/coxpresdb',
    #                                 ensemble_to_entrez_path='./aux_files/entrez_to_ensembl.pkl'),
    #     'outdir': 'coexpression',
    #     'nodes': False,
    #     'edges': True
    # },
    #
    # 'ccre': {
    #     'adapter': CCREAdapter(filepath='./samples/ccre_example.bed.gz'),
    #     'outdir': 'ccre',
    #     'nodes': True,
    #     'edges': False
    # },

    # 'encode_epiraction': {
    #     'adapter': EncodeEnhancerGeneLinkAdapter(filepath='./samples/epiraction_ENCFF712SUP.bed.gz',
    #                                              source='ENCODE_EpiRaction',
    #                                              source_url='https://www.encodeproject.org/annotations/ENCSR831INH/',
    #                                              biological_context='CL_0000765'),
    #
    #     'outdir': 'tf',
    #     'nodes': False,
    #     'edges': True
    # },

    'tflink': {
        'adapter': TFLinkAdapter(filepath='./samples/tflink_homo_sapiens_interactions.tsv',
                                entrez_to_ensemble_map='./aux_files/entrez_to_ensembl.pkl'),
        'outdir': 'tflink',
        'nodes': False,
        'edges': True
    },
    #
    # 'string': {
    #     'adapter': StringPPIAdapter(filepath="./samples/string_human_ppi_v12.0.txt",
    #                                 ensembl_to_uniprot_map="./aux_files/string_ensembl_uniprot_map.pkl"),
    #     "outdir": "ppi",
    #     "nodes": False,
    #     "edges": True
    # },

    # 'tad': {
    #     'adapter': TADAdapter(filepath="./samples/tad_sample.csv"),
    #     "outdir": "tad",
    #     "nodes": True,
    #     "edges": False
    # },

    # 'chromatin_state': {
    #     'adapter': ChromatinStateAdapter(filepath="./samples/roadmap",
    #                                      tissue_id_map="./aux_files/roadmap_epigenomics_tissue_id_map.pkl"),
    #     'outdir': "roadmap",
    #     "nodes": True,
    #     "edges": False
    # },

    'GTEx_eQTL': {
        'adapter': GTExEQTLAdapter(filepath="./samples/gtex"),
        'outdir': "gtex/eqtl",
        'nodes': False,
        'edges': True
    }



}
# Run build
def main():
    """
    Main function. Call individual adapters to download and process data. Build
    via BioCypher from node and edge data.
    """

    # Start biocypher

    bc = MeTTaWriter(schema_config="config/schema_config.yaml",
                     biocypher_config="config/biocypher_config.yaml",
                     output_dir="metta_out")

    # bc.show_ontology_structure()

    # Run adapters

    for k, v in ADAPTERS.items():
        logger.info("Running adapter: " + k)
        adapter = v['adapter']
        write_nodes = v['nodes']
        write_edges = v['edges']
        outdir = v['outdir']

        if write_nodes:
            nodes = adapter.get_nodes()
            bc.write_nodes(nodes, path_prefix=outdir)

        if write_edges:
            edges = adapter.get_edges()
            bc.write_edges(edges, path_prefix=outdir)



    logger.info("Done")

if __name__ == "__main__":
    main()

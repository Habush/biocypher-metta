"""
Knowledge graph generation through BioCypher script
"""

from biocypher_metta.metta_writer import *
from biocypher._logger import logger
from biocypher_metta.adapters.uniprot_adapter import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
    UniprotEdgeType,
    UniprotEdgeField,
)

from biocypher_metta.adapters.ppi_adapter import (
    PPI,
    IntactEdgeField,
    BiogridEdgeField,
    StringEdgeField,
)

from biocypher_metta.adapters.interpro_adapter import (
    InterPro,
    InterProNodeField,
    InterProEdgeField,
)

from biocypher_metta.adapters.go_adapter import (
    GO,
    GONodeType,
    GOEdgeType,
    GONodeField,
    GOEdgeField,
)

from biocypher_metta.adapters.drugbank_adapter import (
    DrugBank,
    DrugBankNodeField,
    DrugbankDTIEdgeField,
    DrugbankEdgeType,
    PrimaryNodeIdentifier,
)

from biocypher import BioCypher

# uniprot configuration
uniprot_node_types = [
    UniprotNodeType.PROTEIN,
    UniprotNodeType.GENE,
    UniprotNodeType.ORGANISM,
]

uniprot_node_fields = [
    # UniprotNodeField.PROTEIN_SECONDARY_IDS,
    UniprotNodeField.PROTEIN_LENGTH,
    UniprotNodeField.PROTEIN_MASS,
    UniprotNodeField.PROTEIN_ORGANISM,
    UniprotNodeField.PROTEIN_ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_PROTEOME,
    UniprotNodeField.PROTEIN_SUBCELLULAR_LOCATION,
    UniprotNodeField.PROTEIN_EC,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS,
    UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS,
    UniprotNodeField.PROTEIN_VIRUS_HOSTS,
    UniprotNodeField.PROTEIN_KEGG_IDS,
]

uniprot_edge_types = [
    UniprotEdgeType.PROTEIN_TO_ORGANISM,
    UniprotEdgeType.GENE_TO_PROTEIN,
]

uniprot_edge_fields = [
    UniprotEdgeField.GENE_ENSEMBL_GENE_ID,
]

# ppi configuration
intact_fields = [field for field in IntactEdgeField]
biogrid_fields = [field for field in BiogridEdgeField]
string_fields = [field for field in StringEdgeField]

# interpro (protein-domain edges and domain nodes) configuration
interpro_node_fields = [field for field in InterProNodeField]
interpro_edge_fields = [field for field in InterProEdgeField]

# GO configuration (Protein-GO, Domain-GO, GO-GO)
go_node_types = [GONodeType.PROTEIN, GONodeType.DOMAIN, GONodeType.BIOLOGICAL_PROCESS,
                 GONodeType.CELLULAR_COMPONENT, GONodeType.MOLECULAR_FUNCTION]
go_edge_types = [GOEdgeType.PROTEIN_TO_BIOLOGICAL_PROCESS, GOEdgeType.DOMAIN_TO_BIOLOGICAL_PROCESS,
                 GOEdgeType.DOMAIN_TO_CELLULAR_COMPONENT, GOEdgeType.DOMAIN_TO_MOLECULAR_FUNCTION] # There are more associations available, however current version of BioCypher probably only supports these
go_node_fields = [field for field in GONodeField]
go_edge_fields = [field for field in GOEdgeField]


# drugbank configuration
drugbank_node_fields = [field for field in DrugBankNodeField]
drugbank_dti_edge_fields = [field for field in DrugbankDTIEdgeField]
drugbank_edge_types = [DrugbankEdgeType.DRUG_TARGET_INTERACTION]
primary_drug_id = PrimaryNodeIdentifier.DRUGBANK


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

    logger.info("Downloading uniprot...")
    # Start uniprot adapter and load data
    uniprot_adapter = Uniprot(
        organism="9606",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        edge_types=uniprot_edge_types,
        edge_fields=uniprot_edge_fields,
        test_mode=True,
    )

    uniprot_adapter.download_uniprot_data(
        cache=True,
        retries=5,
    )

    logger.info("Downloaded uniprot data")

    # logger.info("Downloading ppi datasets...")
    # ppi_adapter = PPI(cache=True,
    #                   organism=9606,
    #                   intact_fields=intact_fields,
    #                   biogrid_fields=biogrid_fields,
    #                   string_fields=string_fields,
    #                   test_mode=True)
    #
    # # download and process intact data
    # ppi_adapter.download_intact_data()
    # ppi_adapter.intact_process()
    #
    # # download and process biogrid data
    # ppi_adapter.download_biogrid_data()
    # ppi_adapter.biogrid_process()
    #
    # # download and process string data
    # ppi_adapter.download_string_data()
    # ppi_adapter.string_process()
    #
    # # Merge all ppi data
    # ppi_adapter.merge_all()
    #
    # logger.info("Downloadeded ppi data")

    # interpro_adapter = InterPro(cache=False,
    #                             page_size=100,
    #                             organism="9606",
    #                             node_fields=interpro_node_fields,
    #                             edge_fields=interpro_edge_fields,
    #                             test_mode=True)
    #
    # # download domain data
    # interpro_adapter.download_domain_node_data()
    # interpro_adapter.download_domain_edge_data()
    #
    # # get interpro nodes and edge
    # interpro_adapter.get_interpro_nodes()
    # interpro_adapter.get_interpro_edges()
    #
    # go_adapter = GO(organism=9606, node_types=go_node_types, go_node_fields=go_node_fields,
    #                 edge_types=go_edge_types, go_edge_fields=go_edge_fields, test_mode=True)
    #
    # # download go data
    # go_adapter.download_go_data(cache=False)
    #
    # # get go nodes and go-protein, domain-go edges
    # go_adapter.get_go_nodes()
    # go_adapter.get_go_edges()
    #
    # drugbank_adapter = DrugBank(drugbank_user="drugbank_username", drugbank_passwd="drugbank_password", node_fields=drugbank_node_fields,
    #                             dti_edge_fields=drugbank_dti_edge_fields, edge_types=drugbank_edge_types, primary_node_id=primary_drug_id,
    #                             test_mode=True)
    #
    # # download drugbank data
    # drugbank_adapter.download_drugbank_data(cache=False)
    #
    # # get drug nodes and drug-target edges
    # drugbank_adapter.get_drug_nodes()
    # drugbank_adapter.get_dti_edges()

    # Write uniprot nodes and edges
    logger.info("Writing uniprot nodes & edges")
    bc.write_nodes(uniprot_adapter.get_nodes())
    bc.write_edges(uniprot_adapter.get_edges())

    # write ppi edges
    # logger.info("Writing uniprot edges")
    # bc.write_edges(ppi_adapter.get_ppi_edges())

    # # write interpro (domain) nodes
    # bc.write_nodes(interpro_adapter.node_list)
    #
    # # write interpro edges (protein-domain) edges
    # bc.write_edges(interpro_adapter.edge_list)
    #
    # # write GO nodes
    # bc.write_nodes(go_adapter.node_list)
    #
    # # write GO edges
    # bc.write_edges(go_adapter.edge_list)
    #
    # # write drug nodes
    # bc.write_nodes(drugbank_adapter.node_list)
    #
    # # write dti edges
    # bc.write_edges(drugbank_adapter.dti_edge_list)
    #
    # # Write import call and other post-processing
    # bc.write_import_call()
    # bc.summary()
    logger.info("Done")

if __name__ == "__main__":
    main()

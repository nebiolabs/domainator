"""NCBI taxonomy: the taxonomy database (:mod:`.ncbi`), and the taxonomy of records (:mod:`.records`)."""
from domainator.taxonomy.ncbi import NCBITaxonomy, default_ncbi_taxonomy_path

__all__ = ["NCBITaxonomy", "default_ncbi_taxonomy_path"]

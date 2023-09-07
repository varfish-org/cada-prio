import csv
import enum
import json
import typing
import warnings

import attrs
import cattrs
from logzero import logger
import pronto


@attrs.frozen(auto_attribs=True)
class GeneIds:
    """Mapping between gene IDs"""

    symbol: str
    hgnc_id: str
    ncbi_gene_id: str
    ensembl_gene_id: str


@enum.unique
class ClinSig(enum.Enum):
    BENIGN = "benign"
    LIKELY_BENIGN = "likely benign"
    UNCERTAIN = "uncertain significance"
    LIKELY_PATHOGENIC = "likely pathogenic"
    PATHOGENIC = "pathogenic"


@attrs.frozen(auto_attribs=True)
class GenePhenotypeRecord:
    """Full record from JSONL file"""

    #: RCV accession
    rcv: str
    #: SCV accession
    scv: str
    #: Clinical significance
    clinsig: typing.Optional[ClinSig]
    #: Submitter
    submitter: typing.Optional[str]
    #: Gene HGNC ID
    hgnc_ids: typing.List[str]
    #: Linked OMIM terms
    omim_terms: typing.List[str]
    #: Linked MONDO terms
    mondo_terms: typing.List[str]
    #: Linked HPO terms
    hpo_terms: typing.List[str]


@attrs.frozen(auto_attribs=True, order=True)
class LinkRecord:
    """Minimal link record"""

    #: Submitter
    submitter: str
    #: Gene HGNC ID
    hgnc_id: str
    #: Linked HPO terms
    hpo_terms: typing.Tuple[str, ...]


def load_gene_phenotype_records(path: str) -> typing.Iterable[GenePhenotypeRecord]:
    with open(path) as f:
        for line in f:
            yield cattrs.structure(json.loads(line), GenePhenotypeRecord)


def run(
    path_out: str,
    path_hgnc_json: str,
    path_gene_hpo_links: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
):
    logger.info("Loading HGNC info from %s...", path_hgnc_json)
    gene_ids = []
    ncbi_to_hgnc = {}
    with open(path_hgnc_json) as f:
        full_hgnc = json.load(f)
        for doc in full_hgnc["response"]["docs"]:
            if all(key in doc for key in ("symbol", "entrez_id", "ensembl_gene_id")):
                gene_id = GeneIds(
                    symbol=doc["symbol"],
                    hgnc_id=doc["hgnc_id"],
                    ncbi_gene_id=doc["entrez_id"],
                    ensembl_gene_id=doc["ensembl_gene_id"],
                )
                ncbi_to_hgnc[gene_id.ncbi_gene_id] = gene_id.hgnc_id
                gene_ids.append(gene_id)
    logger.info("... done loading %s records", len(gene_ids))

    logger.info("Loading ClinVar gene-phenotype links from %s...", path_gene_hpo_links)
    logger.info("  - parsing records ...")
    clinvar_gene_to_phenotype_set = {
        LinkRecord(
            submitter=record.submitter,
            hgnc_id=record.hgnc_ids[0],
            hpo_terms=tuple(record.hpo_terms),
        )
        for record in load_gene_phenotype_records(path_gene_hpo_links)
        if record.clinsig in (ClinSig.LIKELY_PATHOGENIC, ClinSig.PATHOGENIC)
        and len(record.hgnc_ids) == 1
        and record.submitter
    }
    logger.info("  - sorting ...")
    clinvar_gene_to_phenotype = list(sorted(clinvar_gene_to_phenotype_set))
    logger.info("... done loading %s records", len(clinvar_gene_to_phenotype_set))

    logger.info("Loading HPO gene-phenotype links from %s...", path_hpo_genes_to_phenotype)
    logger.info("  - parsing records ...")
    with open(path_hpo_genes_to_phenotype, "rt") as inputf:
        reader = csv.DictReader(inputf, delimiter="\t")
        hpo_by_gene = {}
        for record in reader:
            hpo_by_gene.setdefault(record["ncbi_gene_id"], set()).add(record["hpo_id"])

        hpo_gene_to_phenotype = {
            LinkRecord(
                submitter="HPO",
                hgnc_id=hgnc_id,
                hpo_terms=tuple(hpo_terms),
            )
            for hgnc_id, hpo_terms in hpo_by_gene.items()
        }
    logger.info("  - sorting ...")
    hpo_gene_to_phenotype = list(sorted(hpo_gene_to_phenotype))
    logger.info("... done loading %s records", len(hpo_gene_to_phenotype))

    logger.info("Loading HPO OBO from %s...", path_hpo_obo)
    logger.info("  - parsing ...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UnicodeWarning)
        hpo_ontology = pronto.Ontology(path_hpo_obo)
    logger.info("  - extracting term ids ...")
    hpo_id_from_alt: typing.Dict[str, str] = {}
    hpo_id_to_name: typing.Dict[str, str] = {}
    for term in hpo_ontology.terms():
        hpo_id_to_name[term.id] = term.name
        for alternate_id in term.alternate_ids:
            hpo_id_from_alt[alternate_id] = term.id
            hpo_id_to_name[alternate_id] = term.name
    logger.info("... done loading %s terms", len(hpo_ontology))

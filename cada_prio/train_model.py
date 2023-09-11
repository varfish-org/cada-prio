import csv
import enum
import gzip
import itertools
import json
import os
import pickle
import typing
import warnings

import attrs
import cattrs
from gensim.models import Word2Vec
from logzero import logger
import networkx as nx
import node2vec
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
class GeneToPhenotypeRecord:
    """Minimal link record"""

    #: Submitter
    submitter: str
    #: Gene HGNC ID
    hgnc_id: str
    #: Linked HPO terms
    hpo_terms: typing.Tuple[str, ...]


def load_gene_phenotype_records(path: str) -> typing.Iterable[GenePhenotypeRecord]:
    if path.endswith(".gz"):
        inputf = gzip.open(path, "rt")
    else:
        inputf = open(path, "rt")

    with inputf:
        for line in inputf:
            yield cattrs.structure(json.loads(line), GenePhenotypeRecord)


def load_hpo_ontology(
    path_hpo_obo,
) -> typing.Tuple[pronto.Ontology, typing.Dict[str, str], typing.Dict[str, str]]:
    logger.info("Loading HPO OBO from %s...", path_hpo_obo)
    logger.info("- parsing ...")
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UnicodeWarning)
        hpo_ontology = pronto.Ontology(path_hpo_obo)
    logger.info("- extracting term ids ...")
    hpo_id_from_alt: typing.Dict[str, str] = {}
    hpo_id_to_name: typing.Dict[str, str] = {}
    for term in hpo_ontology.terms():
        if not term.name or not term.id:
            continue
        hpo_id_to_name[term.id] = term.name
        for alternate_id in term.alternate_ids:
            hpo_id_from_alt[alternate_id] = term.id
            hpo_id_to_name[alternate_id] = term.name
    logger.info("... done loading %s terms", len(hpo_ontology))

    return hpo_ontology, hpo_id_from_alt, hpo_id_to_name


def load_hpo_gen2phen(
    path_hpo_genes_to_phenotype: str, ncbi_to_hgnc: typing.Dict[str, str]
) -> typing.List[GeneToPhenotypeRecord]:
    logger.info("Loading HPO gene-phenotype links from %s...", path_hpo_genes_to_phenotype)
    logger.info("- parsing records ...")
    warned_about = set()
    with open(path_hpo_genes_to_phenotype, "rt") as inputf:
        reader = csv.DictReader(inputf, delimiter="\t")
        hpo_by_gene: typing.Dict[str, typing.Set[str]] = {}
        for record in reader:
            ncbi_gene_id = record["ncbi_gene_id"]
            if ncbi_gene_id not in ncbi_to_hgnc:
                if ncbi_gene_id not in warned_about:
                    logger.warning("No HGNC entry for Entrez %s found", ncbi_gene_id)
                    warned_about.add(ncbi_gene_id)
                continue
            hgnc_id = ncbi_to_hgnc[ncbi_gene_id]
            hpo_by_gene.setdefault(hgnc_id, set()).add(record["hpo_id"])

        hpo_gene_to_phenotype_set = {
            GeneToPhenotypeRecord(
                submitter="HPO",
                hgnc_id=hgnc_id,
                hpo_terms=tuple(hpo_terms),
            )
            for hgnc_id, hpo_terms in hpo_by_gene.items()
        }
    logger.info("- sorting ...")
    hpo_gene_to_phenotype = list(sorted(hpo_gene_to_phenotype_set))
    logger.info("... done loading %s records", len(hpo_gene_to_phenotype))

    return hpo_gene_to_phenotype


def load_clinvar_gen2phen(path_gene_hpo_links: str) -> typing.List[GeneToPhenotypeRecord]:
    logger.info("Loading ClinVar gene-phenotype links from %s...", path_gene_hpo_links)
    logger.info("- parsing records ...")
    clinvar_gene_to_phenotype_set = {
        GeneToPhenotypeRecord(
            submitter=record.submitter,
            hgnc_id=record.hgnc_ids[0],
            hpo_terms=tuple(record.hpo_terms),
        )
        for record in load_gene_phenotype_records(path_gene_hpo_links)
        if record.clinsig in (ClinSig.LIKELY_PATHOGENIC, ClinSig.PATHOGENIC)
        and len(record.hgnc_ids) == 1
        and record.submitter
    }
    logger.info("- sorting ...")
    clinvar_gene_to_phenotype = list(sorted(clinvar_gene_to_phenotype_set))
    logger.info("... done loading %s records", len(clinvar_gene_to_phenotype_set))

    return clinvar_gene_to_phenotype


def load_hgnc_info(path_hgnc_json) -> typing.Tuple[typing.Dict[str, str], typing.List[GeneIds]]:
    logger.info("Loading HGNC info from %s...", path_hgnc_json)
    gene_ids = []
    ncbi_to_hgnc = {}
    with open(path_hgnc_json) as f:
        full_hgnc = json.load(f)
        for doc in full_hgnc["response"]["docs"]:
            if all(key in doc for key in ("symbol", "entrez_id")):
                gene_id = GeneIds(
                    symbol=doc["symbol"],
                    hgnc_id=doc["hgnc_id"],
                    ncbi_gene_id=doc["entrez_id"],
                    ensembl_gene_id=doc.get("ensembl_gene_id"),
                )
                ncbi_to_hgnc[gene_id.ncbi_gene_id] = gene_id.hgnc_id
                gene_ids.append(gene_id)
    logger.info("... done loading %s records", len(gene_ids))

    return ncbi_to_hgnc, gene_ids


#: type alias for an edge.
Edge = typing.Tuple[str, str]


def yield_hpo_edges(hpo_ontology: pronto.Ontology) -> typing.Iterator[Edge]:
    """Yield hierarchical edges from the HPO ontology"""
    for term in hpo_ontology.terms():
        parents = list(term.superclasses(distance=1))
        for parent in parents[1:]:
            yield (term.id, parent.id)


def yield_gene2phen_edges(links: typing.List[GeneToPhenotypeRecord]) -> typing.Iterator[Edge]:
    """Yield edges from link record list"""
    for link in links:
        for hpo_term in link.hpo_terms:
            yield (link.hgnc_id, hpo_term)


@attrs.frozen(auto_attribs=True)
class EmbeddingParams:
    """Parameters for the embedding"""

    #: The number of dimensions of feature vectors
    dimensions: int = 300
    #: The number of nodes in each random walk
    walk_length: int = 60
    #: Controls the probability for a walk to visit immediately back to the previous node
    p: float = 1.7987535798694703
    #: Controls the probability for a walk to visit previously unexplored neighborhoods in the
    #: graph
    q: float = 3.875406134463754
    #: Number of random walks to be generated from each node in the graph
    num_walks: int = 10
    #: Limit on the number of words in each context
    window: int = 4
    #: Set the min_count in the fitting
    min_count: int = 1
    #: Set the batch_words in the fitting
    batch_words: int = 4


def build_and_fit_model(*, clinvar_gen2phen, hpo_gen2phen, hpo_ontology, cpus: int = 1):
    # create graph edges combining HPO hierarchy and training edges from ClinVar
    logger.info("Constructing training graph ...")
    logger.info("- building edges ...")
    training_edges = list(
        itertools.chain(
            yield_hpo_edges(hpo_ontology),
            yield_gene2phen_edges(hpo_gen2phen),
            yield_gene2phen_edges(clinvar_gen2phen),
        )
    )
    logger.info("- graph construction")
    training_graph = nx.Graph()
    training_graph.add_edges_from(training_edges)
    logger.info("... done constructing training graph with %s edges", len(training_edges))

    logger.info("Computing the embedding / model fit...")
    logger.info("- embedding graph")
    embedding_params = EmbeddingParams()
    embedding = node2vec.Node2Vec(
        training_graph,
        dimensions=embedding_params.dimensions,
        walk_length=embedding_params.walk_length,
        num_walks=embedding_params.num_walks,
        p=embedding_params.p,
        q=embedding_params.q,
        workers=cpus,
    )
    logger.info("- fitting model")
    model = embedding.fit(
        window=embedding_params.window,
        min_count=embedding_params.min_count,
        batch_words=embedding_params.batch_words,
    )
    logger.info("... done computing the embedding")
    return training_graph, model


def write_graph_and_model(
    path_out, hgnc_info: typing.List[GeneIds], training_graph: nx.Graph, model: Word2Vec
):
    os.makedirs(path_out, exist_ok=True)

    path_hgnc_info = os.path.join(path_out, "hgnc_info.jsonl")
    logger.info("Saving HGNC info to %s...", path_hgnc_info)
    with open(path_hgnc_info, "w") as f:
        for record in hgnc_info:
            json.dump(cattrs.unstructure(record), f)
    logger.info("... done saving HGNC info")

    path_graph = os.path.join(path_out, "graph.gpickle")
    logger.info("Saving graph to %s...", path_graph)
    with open(path_graph, "wb") as outputf:
        pickle.dump(training_graph, outputf)
    logger.info("... done saving graph")

    logger.info("Saving embedding to %s...", path_out)
    path_embeddings = os.path.join(path_out, "embedding")
    logger.info("- %s", path_embeddings)
    model.wv.save_word2vec_format(path_embeddings)
    path_model = os.path.join(path_out, "model")
    logger.info("- %s", path_model)
    model.save(path_model)
    logger.info("... done saving embedding to")


def run(
    path_out: str,
    path_hgnc_json: str,
    path_gene_hpo_links: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
    cpus: int = 1,
):
    # load all data
    ncbi_to_hgnc, hgnc_info = load_hgnc_info(path_hgnc_json)
    clinvar_gen2phen = load_clinvar_gen2phen(path_gene_hpo_links)
    hpo_gen2phen = load_hpo_gen2phen(path_hpo_genes_to_phenotype, ncbi_to_hgnc)
    hpo_ontology, hpo_id_from_alt, hpo_id_to_name = load_hpo_ontology(path_hpo_obo)
    _, _ = hpo_id_from_alt, hpo_id_to_name

    # build and fit model
    training_graph, model = build_and_fit_model(
        clinvar_gen2phen=clinvar_gen2phen,
        hpo_gen2phen=hpo_gen2phen,
        hpo_ontology=hpo_ontology,
        cpus=cpus,
    )
    # write out graph and model
    write_graph_and_model(path_out, hgnc_info, training_graph, model)

"""Console script for CADA"""

import json
import logging
import sys
import tempfile
import typing

import cattr
import click
import logzero
import yaml

try:
    import optuna

    _ = optuna
    HAVE_OPTUNA = True
except ImportError:
    HAVE_OPTUNA = False

from cada_prio import _version, inspection, param_opt, predict, train_model


@click.group()
@click.version_option(_version.__version__)
@click.option("--verbose/--no-verbose", default=False)
@click.pass_context
def cli(ctx: click.Context, verbose: bool):
    """Main entry point for CLI via click."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    if verbose:
        logzero.loglevel(logging.DEBUG)
    else:
        logzero.loglevel(logging.INFO)


@cli.command("train")
@click.argument("path_out", type=str)
@click.option("--path-hgnc-json", type=str, help="path to HGNC JSON", required=True)
@click.option(
    "--path-gene-hpo-links", type=str, help="path to gene-HPO term links JSONL", required=True
)
@click.option(
    "--path-hpo-genes-to-phenotype",
    type=str,
    help="path to genes_to_phenotype.txt file",
    required=True,
)
@click.option("--path-hpo-obo", type=str, help="path HPO OBO file", required=True)
@click.option(
    "--path-embedding-params", type=str, help="optional path to JSON file with embedding parameters"
)
@click.option("--cpus", type=int, help="number of CPUs to use", default=1)
@click.pass_context
def cli_train(
    ctx: click.Context,
    path_out: str,
    path_hgnc_json: str,
    path_gene_hpo_links: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
    path_embedding_params: typing.Optional[str],
    cpus: int,
):
    """train model"""
    ctx.ensure_object(dict)
    train_model.run(
        path_out,
        path_hgnc_json,
        path_gene_hpo_links,
        path_hpo_genes_to_phenotype,
        path_hpo_obo,
        path_embedding_params,
        cpus,
    )


@cli.command("predict")
@click.argument("path_model", type=str)
@click.option(
    "--hpo-terms",
    type=str,
    help="comma-separate HPO terms or @file with space-separated ones",
    required=True,
)
@click.option(
    "--genes",
    type=str,
    help="comma-separated genes to restrict prediction to or @file with space-separated ones",
)
@click.pass_context
def cli_predict(
    ctx: click.Context,
    path_model: str,
    hpo_terms: str,
    genes: typing.Optional[str],
):
    """perform prediction/prioritication"""
    ctx.ensure_object(dict)

    if hpo_terms.startswith("@"):
        with open(hpo_terms[1:]) as f:
            hpo_term_list = [x.strip() for x in f.read().split()]
    else:
        hpo_term_list = [x.strip() for x in hpo_terms.split(",")]

    gene_list = None
    if genes:
        if genes.startswith("@"):
            with open(genes[1:]) as f:
                gene_list = [x.strip() for x in f.read().split()]
        else:
            gene_list = [x.strip() for x in genes.split(",")]

    if predict.run(path_model, hpo_term_list, gene_list) != 0:
        ctx.exit(1)


@cli.group("utils")
def cli_utils():
    """utilities"""


@cli_utils.command("dump-graph")
@click.option(
    "--hgnc-to-entrez/--no-hgnc-to-entrez", help="enable HGNC to Entrez mapping", default=True
)
@click.argument("path_graph", type=str)
@click.argument("path_hgnc_info", type=str)
@click.pass_context
def cli_dump_graph(
    ctx: click.Context,
    path_graph: str,
    path_hgnc_info: str,
    hgnc_to_entrez: bool,
):
    """dump graph edges for debugging"""
    ctx.ensure_object(dict)
    inspection.dump_graph(path_graph, path_hgnc_info, hgnc_to_entrez)


@cli_utils.command("dump-openapi-yaml")
@click.argument("path_out", type=str)
def cli_dump_openapi_yaml(path_out: str):
    """Dump OpenAPI YAML file"""
    from cada_prio import rest_server

    with open(path_out, "wt") as f:
        yaml.dump(rest_server.app.openapi(), f)


@cli.group("tune")
def cli_tune():
    """hyperparameter tuning"""


@cli_tune.command("train-eval")
@click.argument("path_out", type=str)
@click.option("--path-hgnc-json", type=str, help="path to HGNC JSON", required=True)
@click.option(
    "--path-hpo-genes-to-phenotype",
    type=str,
    help="path to genes_to_phenotype.txt file",
    required=True,
)
@click.option("--path-hpo-obo", type=str, help="path HPO OBO file", required=True)
@click.option(
    "--path-clinvar-phenotype-links",
    type=str,
    help="path to ClinVar phenotype links JSONL",
    required=True,
)
@click.option(
    "--fraction-links",
    type=float,
    help="fraction of links to add to the graph (conflicts with --path-validation-links)",
)
@click.option(
    "--path-validation-links",
    type=str,
    help="path to validation links JSONL (conflicts with --fraction-links)",
)
@click.option(
    "--path-embedding-params", type=str, help="path to JSON file with embedding params; optional"
)
@click.option(
    "--seed",
    type=int,
    help="seed for random number generator",
)
@click.option("--cpus", type=int, help="number of CPUs to use", default=1)
@click.pass_context
def cli_train_eval(
    ctx: click.Context,
    path_out: str,
    path_hgnc_json: str,
    path_hpo_genes_to_phenotype: str,
    path_hpo_obo: str,
    path_clinvar_phenotype_links: str,
    fraction_links: typing.Optional[float],
    path_validation_links: typing.Optional[str],
    path_embedding_params: typing.Optional[str],
    seed: typing.Optional[int],
    cpus: int,
):
    """train and evaluate model for one set of parameters"""
    if bool(fraction_links) == bool(path_validation_links):
        raise click.UsageError(
            "exactly one of --fraction-links and --path-validation-links must be given"
        )

    ctx.ensure_object(dict)
    param_opt.train_and_validate(
        path_out=path_out,
        path_hgnc_json=path_hgnc_json,
        path_hpo_genes_to_phenotype=path_hpo_genes_to_phenotype,
        path_hpo_obo=path_hpo_obo,
        path_clinvar_phenotype_links=path_clinvar_phenotype_links,
        fraction_links=fraction_links,
        path_validation_links=path_validation_links,
        path_embedding_params=path_embedding_params,
        seed=seed,
        cpus=cpus,
    )


@cli_tune.command("eval-only")
@click.argument("path_model", type=str)
@click.argument(
    "path_validation_links",
    type=str,
)
@click.pass_context
def cli_eval_only(
    ctx: click.Context,
    path_model: str,
    path_validation_links: str,
):
    """run a trained model on validation links"""
    ctx.ensure_object(dict)
    param_opt.run_validation(
        path_out=path_model,
        path_links_validation=path_validation_links,
    )


if HAVE_OPTUNA:

    @cli_tune.command("run-optuna")
    @click.argument("storage", type=str)
    @click.option("--n-trials", type=int, help="number of trials to run; default: 100", default=100)
    @click.option(
        "--study-name",
        type=str,
        help="name of Optuna study; default: cada-tune",
        default="cada-tune",
    )
    @click.option("--path-hgnc-json", type=str, help="path to HGNC JSON", required=True)
    @click.option(
        "--path-hpo-genes-to-phenotype",
        type=str,
        help="path to genes_to_phenotype.txt file",
        required=True,
    )
    @click.option("--path-hpo-obo", type=str, help="path HPO OBO file", required=True)
    @click.option(
        "--path-clinvar-phenotype-links",
        type=str,
        help="path to ClinVar phenotype links JSONL",
        required=True,
    )
    @click.option(
        "--fraction-links",
        type=float,
        help="fraction of links to add to the graph (conflicts with --path-validation-links)",
    )
    @click.option(
        "--path-validation-links",
        type=str,
        help="path to validation links JSONL (conflicts with --fraction-links)",
    )
    @click.option(
        "--seed",
        type=int,
        help="seed for random number generator",
    )
    @click.option("--cpus", type=int, help="number of CPUs to use", default=1)
    def cli_run_optuna(
        storage: str,
        study_name: str,
        n_trials: int,
        path_hgnc_json: str,
        path_hpo_genes_to_phenotype: str,
        path_hpo_obo: str,
        path_clinvar_phenotype_links: str,
        fraction_links: typing.Optional[float],
        path_validation_links: typing.Optional[str],
        seed: typing.Optional[int],
        cpus: int,
    ):
        """run hyperparameter tuning"""

        def objective(trial: optuna.trial.BaseTrial) -> float:
            with tempfile.TemporaryDirectory() as tmpdir, open(
                f"{tmpdir}/params.json", "wt"
            ) as tmpf:
                json.dump(
                    cattr.unstructure(
                        train_model.EmbeddingParams(
                            dimensions=trial.suggest_int("dimensions", low=100, high=500, step=10),
                            walk_length=trial.suggest_int("walk_length", low=1, high=100),
                            p=trial.suggest_float("p", low=0.1, high=2.5),
                            q=trial.suggest_float("q", low=0.0, high=1.0),
                            num_walks=trial.suggest_int("num_walks", low=10, high=50),
                            window=trial.suggest_int("window", low=4, high=8),
                            epochs=trial.suggest_int("epochs", low=3, high=6),
                            use_skipgram=trial.suggest_categorical("use_skipgram", [True, False]),
                            min_count=trial.suggest_int("min_count", low=1, high=5),
                            batch_words=trial.suggest_int("batch_words", low=3, high=6),
                        )
                    ),
                    tmpf,
                    indent=2,
                )
                tmpf.flush()
                result = param_opt.train_and_validate(
                    path_out=f"{tmpdir}/out",
                    path_hgnc_json=path_hgnc_json,
                    path_hpo_genes_to_phenotype=path_hpo_genes_to_phenotype,
                    path_hpo_obo=path_hpo_obo,
                    path_clinvar_phenotype_links=path_clinvar_phenotype_links,
                    fraction_links=fraction_links,
                    path_validation_links=path_validation_links,
                    path_embedding_params=f"{tmpdir}/params.json",
                    seed=seed,
                    cpus=cpus,
                )
                return result[100]

        optuna.logging.get_logger("optuna").addHandler(logging.StreamHandler(sys.stdout))
        study = optuna.create_study(study_name=study_name, storage=storage, load_if_exists=True)
        study.optimize(objective, n_trials=n_trials)

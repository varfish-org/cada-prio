"""Model quality tests"""

from cada_prio import param_opt


def test_quality(tmpdir):
    """Test quality of model built from 'classic' data from the original paper."""
    result = param_opt.train_and_validate(
        path_out=f"{tmpdir}/quality-model",
        path_hgnc_json="data/classic/hgnc_complete_set.json",
        path_hpo_genes_to_phenotype="data/classic/genes_to_phenotype.all_source_all_freqs_etc.txt",
        path_hpo_obo="data/classic/hp.obo",
        path_embedding_params=None,
        path_clinvar_phenotype_links="/dev/null",
        path_validation_links="data/classic/cases_validate.jsonl",
        fraction_links=None,
        seed=1,
        cpus=1,
    )

    assert result[100] >= 56.0, f"result={result}"

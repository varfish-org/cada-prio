from cada_prio import train_model


def test_train_smoke_test(tmpdir):
    train_model.run(
        f"{tmpdir}/train_out",
        "tests/data/train_smoke/hgnc_complete_set_2023-09-01.json",
        "tests/data/train_smoke/links.jsonl",
        "tests/data/train_smoke/genes_to_phenotype.txt",
        "tests/data/train_smoke/hp.obo",
    )

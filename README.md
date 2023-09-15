[![CI](https://github.com/bihealth/cada-prio/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/bihealth/cada-prio/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/bihealth/cada-prio/graph/badge.svg?token=HIBwaG4eYM)](https://codecov.io/gh/bihealth/cada-prio)
[![Documentation Status](https://readthedocs.org/projects/cada-prio/badge/?version=latest)](https://cada-prio.readthedocs.io/en/latest/?badge=latest)
[![Pypi](https://img.shields.io/pypi/pyversions/cada-prio.svg)](https://pypi.org/project/cada-prio)

# CADA: The Next Generation

This is a re-implementation of the [CADA](https://github.com/Chengyao-Peng/CADA) method for phenotype-similarity prioritization.

- Free software: MIT license
- Documentation: https://cada-prio.readthedocs.io/en/latest/
- Discussion Forum: https://github.com/bihealth/cada-prio/discussions
- Bug Reports: https://github.com/bihealth/cada-prio/issues

## Running Hyperparameter Tuning

Install with `tune` feature enabled:

```
pip install cada-prio[tune]
```

Create study with `optuna`:

```
optuna create-study --storage sqlite:///local_data/cada-tune.sqlite --study-name=cada-tune
```

Run tuning, e.g., on the "classic" model.
Thanks to [optuna](https://optuna.org/), you can run this in parallel as long as the database is shared.
Each run will use 4 CPUs in the example below and perform 1 trial.

```
cada-prio tune run-optuna \
    sqlite:///local_data/cada-tune.sqlite \
    --path-hgnc-json data/classic/hgnc_complete_set.json \
    --path-hpo-genes-to-phenotype data/classic/genes_to_phenotype.all_source_all_freqs_etc.txt \
    --path-hpo-obo data/classic/hp.obo \
    --path-clinvar-phenotype-links data/classic/cases_train.jsonl \
    --path-validation-links data/classic/cases_validate.jsonl \
    --n-trials 1 \
    --cpus=4
```

## Managing GitHub Project with Terraform

```
# export GITHUB_OWNER=bihealth
# export GITHUB_TOKEN=ghp_<thetoken>

# cd utils/terraform

# terraform init
# terraform import github_repository.cada-prio cada-prio
# terraform validate
# terraform fmt
# terraform plan
# terraform apply
```

[![CI](https://github.com/bihealth/cada-ng/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/bihealth/cada-ng/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/bihealth/cada-ng/graph/badge.svg?token=HIBwaG4eYM)](https://codecov.io/gh/bihealth/cada-ng)
[![Documentation Status](https://readthedocs.org/projects/cada-ng/badge/?version=latest)](https://cada-ng.readthedocs.io/en/latest/?badge=latest)
[![Pypi](https://img.shields.io/pypi/pyversions/cada-ng.svg)](https://pypi.org/project/cada-ng)

# CADA: The Next Generation

This is a re-implementation of the [CADA](https://github.com/Chengyao-Peng/CADA) method for phenotype-similarity prioritization.

- Free software: MIT license
- Documentation: https://cada-ng.readthedocs.io/en/latest/
- Discussion Forum: https://github.com/bihealth/cada-ng/discussions
- Bug Reports: https://github.com/bihealth/cada-ng/issues


## Managing GitHub Project with Terraform

```
# export GITHUB_OWNER=bihealth
# export GITHUB_TOKEN=ghp_<thetoken>

# cd utils/terraform

# terraform init
# terraform import github_repository.cada-ng cada-ng
# terraform validate
# terraform fmt
# terraform plan
# terraform apply
```

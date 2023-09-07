#!/bin/bash

# Utility script to start the Docker build process.

set -x
set -euo pipefail

ORG=bihealth
REPO=cada-prio
DOCKER_VERSION=${DOCKER_VERSION-adhoc}

sudo docker build . \
    --file utils/docker/Dockerfile \
    --pull \
    -t ghcr.io/$ORG/$REPO:$DOCKER_VERSION \
    "$@"

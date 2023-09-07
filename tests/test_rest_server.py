from unittest import mock

import pytest
from starlette.testclient import TestClient

from cada_prio import rest_server


@mock.patch("cada_prio.rest_server.PATH_DATA", "tests/data/model_smoke")
@pytest.mark.asyncio
async def test_version():
    with TestClient(rest_server.app) as client:
        response = client.get("/version")
        assert response.status_code == 200


@mock.patch("cada_prio.rest_server.PATH_DATA", "tests/data/model_smoke")
@pytest.mark.asyncio
async def test_predict_with_gene():
    with TestClient(rest_server.app) as client:
        response = client.get("/predict/?hpo_terms=HP:0008551&hpo=HP:0000007&gene=MKS1")
        assert response.status_code == 200


@mock.patch("cada_prio.rest_server.PATH_DATA", "tests/data/model_smoke")
@pytest.mark.asyncio
async def test_predict_without_gene():
    with TestClient(rest_server.app) as client:
        response = client.get("/predict/?hpo_terms=HP:0008551&hpo=HP:0000007")
        assert response.status_code == 200

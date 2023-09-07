"""REST API for CADA using FastAPI."""
from contextlib import asynccontextmanager
import os
import typing

from dotenv import load_dotenv
from fastapi import FastAPI, Query
from pydantic import BaseModel
from starlette.responses import Response

import cada_prio
from cada_prio import predict

# Load environment
env = os.environ
load_dotenv()

#: Debug mode
DEBUG = env.get("CADA_DEBUG", "false").lower() in ("true", "1")
#: Path to data / model
PATH_DATA = env.get("CADA_PATH_DATA", "/data/cada")

#: The CADA models, to be loaded on startup.
GLOBAL_STATIC = {}


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Load the models
    all_to_hgnc, hgnc_info_by_id = predict.load_hgnc_info(PATH_DATA)
    GLOBAL_STATIC["all_to_hgnc"] = all_to_hgnc
    GLOBAL_STATIC["hgnc_info_by_id"] = hgnc_info_by_id
    graph, model = predict.load_graph_model(PATH_DATA)
    GLOBAL_STATIC["graph"] = graph
    GLOBAL_STATIC["model"] = model

    yield

    GLOBAL_STATIC.clear()


app = FastAPI(lifespan=lifespan)


# Register endpoint for returning CADA version.
@app.get("/version")
async def version():
    return Response(content=cada_prio.__version__)


class PredictionResult(BaseModel):
    rank: int
    score: float
    gene_symbol: str
    ncbi_gene_id: str
    hgnc_id: str


# Register endpoint for the prediction
@app.get("/predict")
async def handle_predict(
    hpo_terms: typing.Annotated[typing.List[str], Query()],
    genes: typing.Annotated[typing.Optional[typing.List[str]], Query()] = [],
):
    _, sorted_scores = predict.run_prediction(
        hpo_terms,
        genes,
        GLOBAL_STATIC["all_to_hgnc"],
        GLOBAL_STATIC["graph"],
        GLOBAL_STATIC["model"],
    )
    hgnc_info_by_id = GLOBAL_STATIC["hgnc_info_by_id"]
    return [
        PredictionResult(
            rank=rank,
            score=score,
            gene_symbol=hgnc_info_by_id[hgnc_id].symbol,
            ncbi_gene_id=hgnc_info_by_id[hgnc_id].ncbi_gene_id,
            hgnc_id=hgnc_info_by_id[hgnc_id].hgnc_id,
        )
        for rank, (hgnc_id, score) in enumerate(sorted_scores.items(), start=1)
    ]

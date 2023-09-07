import typing

import gensim.models
import networkx as nx

class Node2Vec:
    def __init__(
        self,
        graph: nx.Graph,
        dimensions: int = 128,
        walk_length: int = 80,
        num_walks: int = 10,
        p: float = 1,
        q: float = 1,
        weight_key: str = "weight",
        workers: int = 1,
        sampling_strategy: typing.Optional[dict] = None,
        quiet: bool = False,
        temp_folder: typing.Optional[str] = None,
        seed: typing.Optional[int] = None,
    ): ...
    def fit(
        self, window: int, min_count: int, batch_words: int, **kwargs
    ) -> gensim.models.Word2Vec: ...

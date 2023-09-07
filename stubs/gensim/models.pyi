class KeyedVectors:
    def save_word2vec_format(self, fname: str): ...

class Word2Vec:
    wv: KeyedVectors

    def save(self, fname: str): ...

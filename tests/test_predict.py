from cada_prio import predict


def test_predict_smoke_test():
    predict.run("tests/data/model_smoke", ["HP:0008551"])

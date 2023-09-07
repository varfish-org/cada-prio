from cada_prio import predict


def test_predict_smoke_test(tmpdir):
    predict.run("tests/data/model_smoke", "HP:0008551")

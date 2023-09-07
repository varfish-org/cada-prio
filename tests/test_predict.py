from cada_prio import train_model


def test_predict_smoke_test(tmpdir):
    train_model.run("tests/data/model_smoke", "HP:0008551")

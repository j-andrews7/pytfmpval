from pytfmpval import tfmp


def test_score():
    m = tfmp.create_matrix("tests/MA0045.pfm")

    assert round(tfmp.pval2score(m, 0.00001), 2) == 8.77


def test_pval():
    m = tfmp.create_matrix("tests/MA0045.pfm")

    assert round(tfmp.score2pval(m, 8.7708), 5) == 0.00001

from spower.utils import getLogger
from playground import *
if __name__ == "__main__":
    d = PlayGround(getLogger(1))
    assert d.summary() == 0
    # set genotype by MAF
    assert True == L.GenotypeGenerator().apply(d.data)
    res = d.get(["haplotype1", "haplotype2"])
    assert res["haplotype1"] == (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
    assert res["haplotype2"] == (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)

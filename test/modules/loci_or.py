from spower.utils import getLogger
from playground import *
if __name__ == "__main__":
    d = PlayGround(getLogger(1))
    # Generate GF in cases, via OR model
    #- set OR model
    assert True == L.ORModel([0,1.78,0,0.9,1.1]).apply(d.data)
    res = d.get(["effect", 'par'])
    assert res['effect'] == (1.78, 1.78, 1.78, 1.78, 1.0, 1.1, 1.0, 0.9, 1.78, 1.78, 1.0, 1.78)
    almost_equal(res['par'], (0.0007527551788200869, 2.5975946043922594e-05, 0.00015583543875803612, 2.5975946043922594e-05, 0.0, 0.038674581525736176, 0.0, 3.700359640847335e-06, 0.0003116223926300634, 2.5975946043922594e-05, 0.0, 0.000726816919595037))
    #- update GF by OR method
    res0 = d.get(["gf0", "gf2"])
    assert True == L.ORGFUpdater(2.0, 0.01).apply(d.data)
    res = d.get(["gf0", "gf2"])
    #
    for x,y,z in zip(res0['gf0'], res['gf0'], d.get("direction")):
        if z == 'd':
            assert x > y
        else:
            assert round(x,5) <= round(y,5)

    for x,y,z in zip(res0['gf2'], res['gf2'], d.get("direction")):
        if z == 'd':
            assert x < y
        else:
            assert round(x,5) >= round(y,5)
    #
    almost_equal(res['gf0'], (0.9982956803933704, 0.9999411806025909, 0.9996471369540085, 0.9999411806025909, 0.9999000924155059, 0.616389730164512, 0.9341482876914409, 0.9999699971512321, 0.9992944017432116, 0.9999411806025909, 0.9982689749609894, 0.9983544003888292))
    almost_equal(res['gf2'], (5.901075301235597e-07, 7.021750290206489e-10, 2.5275076303345994e-08, 7.021750290206489e-10, 2.4955060214016003e-09, 0.04588573802725315, 0.0011213476914409, 2.223230920032367e-10, 1.0108487126230469e-07, 7.021750290206489e-10, 7.497609894787601e-07, 5.501261115497627e-07))
    #
    res = d.get(['wt_penetrance','heterozygotes_penetrance','homozygotes_penetrance','all_prevalence'])
    almost_equal(res["wt_penetrance"],(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01))
    almost_equal(res["heterozygotes_penetrance"],(0.01766223457035126, 0.01766223457035126, 0.01766223457035126, 0.01766223457035126, 0.01, 0.010989010989010992, 0.01, 0.00900900900900901, 0.01766223457035126, 0.01766223457035126, 0.01, 0.01766223457035126))
    almost_equal(res["homozygotes_penetrance"],(0.025324469140702525, 0.025324469140702525, 0.025324469140702525, 0.025324469140702525, 0.01, 0.011978021978021985, 0.01, 0.008018018018018021, 0.025324469140702525, 0.025324469140702525, 0.01, 0.025324469140702525))
    almost_equal(res["all_prevalence"], (0.010007400166226234, 0.010000255178156302, 0.010001531068937804, 0.010000255178156302, 0.01, 0.010394715076923078, 0.01, 0.00999996699667027, 0.010003062138488589, 0.010000255178156302, 0.01, 0.010007144987763445))
    # Generate disease status
    assert True == L.DiseaseEffectGenerator(0.01).apply(d.data) 
    assert round(d.get('loci_penetrance'),6) == round(0.8784655439599104,6)
    assert True == L.DiseaseStatusGenerator().apply(d.data) 
    assert d.get('phenotype') == 2.0

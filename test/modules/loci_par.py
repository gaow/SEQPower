from spower.utils import getLogger
from playground import *
if __name__ == "__main__":
    d = PlayGround(getLogger(1))
    # Generate GF in cases, via PAR model
    #- set PAR model
    # huge OR resulted
    assert True == L.PARModel([.25, .25, 0, 0],True).apply(d.data)
    res = d.get(["effect", 'par'])
    almost_equal(.5, sum(res["par"]))
    almost_equal(res['effect'], (3.69538420811764, 2446.0352249655757, 64.6018420615153, 2446.0352249655757, 1.0, 1.0, 1.0, 9.990065341916748e-05, 16.80005324201965, 2446.0352249655757, 1.0, 3.8916187359828944))
    almost_equal(res['par'], (0.002596439286207404, 0.07529673618148339, 0.01254945603024723, 0.07529673618148339, 0.0, 0.0, 0.0, 0.25, 0.006274726759048802, 0.07529673618148339, 0.0, 0.0026891693800463775))
    #- update GF by PAR method
    res0 = d.get(["gf0", "gf2"])
    assert True == L.PARGFUpdater(2).apply(d.data)
    res = d.get(["gf0", "gf2"])
    strictly_le(res["gf0"], res0["gf0"])
    strictly_ge(res["gf2"], res0["gf2"])
    almost_equal(res['gf0'],(0.9964586302149315, 0.9345310735676225, 0.9875652213478052, 0.9345310735676225, 0.9999000924155059, 0.64071956214016, 0.9341482876914409, 0.9999666969172785, 0.9934107737743848, 0.9345310735676225, 0.9982689749609894, 0.99640017414124))
    almost_equal(res['gf2'],(1.7189998721209188e-06, 1.2614814543107313e-06, 1.2737311635539153e-06, 1.2614814543107313e-06, 2.4955060214016003e-09, 0.03982036214015999, 0.0011213476914409, 2.772784468224e-10, 1.3332247368258924e-06, 1.2614814543107313e-06, 7.4976098947876e-07, 1.6874324919721214e-06))

from spower.utils import getLogger
from spower.simulator.sampler import Sample
if __name__ == "__main__":
    d = Sample(getLogger(1))
    d.seed(0,1)
    c = d.clone()
    c.seed(0,1)


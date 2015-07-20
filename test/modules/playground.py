from spower.simulator.sampler import * 
from spower.utils import almost_equal, strictly_le,strictly_ge
class PlayGround(Sample):
    def __init__(self, logger):
        Sample.__init__(self, logger)
        self.seed(911)
        self.maf = "0.0004828987	1.665168e-05	9.991008e-05	1.665168e-05	4.995504e-05	1.995504e-01	0.03348653	1.665168e-05	0.0001998202	1.665168e-05	0.0008658874	0.000466247".split()
        self.function_score = "0.002798582	0.002588026	1.8711e-08	3.713846e-10	0.02184875	0.05	0.0	-0.06091417	0.05147886	0.00377004	0.0	0.01586215".split()
        self.maf = list(map(float, self.maf))
        self.function_score = list(map(float, self.function_score))
        self.function_class = "1	1	1	1	0	1	0	1	1	1	0	1".split()
        self.position = "201	240	265	324	370	479	578	817	912	939	1073	1480".split()
        self.length = len(self.maf)
        self.rare_cutoff = 0.01
        self.reset()

    def reset(self):
        self.set(haplotype1 = [1] * self.length, haplotype2 = [1] * self.length, maf = self.maf)
        self.set(moi = 'A', gf0 = [(1 - x) ** 2 for x in self.maf], gf2 = [x ** 2 for x in self.maf])
        self.set(function_score = self.function_score, position = self.position)
        self.set(function_class=['ns' if x == '1' else 's' for x in self.function_class],
                 variant_class = ['r' if x < self.rare_cutoff else 'c' for x in self.maf],
                 direction = [x if y =='1' else 'n' for x, y in zip(['p' if x<0 else 'd' for x in self.function_score], self.function_class)],
                 missingness = ['NA'] * self.length)
        
    def summary(self):
        print(self.get(["maf", 'function_class', 'variant_class',
                         'direction', 'position', 'missingness', 'moi',
                         'haplotype1', 'haplotype2', 'function_score', 'gf0',
                         'gf2']))
        return 0

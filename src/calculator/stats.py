#!/usr/bin/env python
# -*- coding: utf-8 -*-
# $File: stats.py $
# $LastChangedDate:  $
# $Rev:  $
# This file is part of the SPower program
# Copyright (c) 2012, Gao Wang <ewanggao@gmail.com>
# GNU General Public License (http://www.gnu.org/licenses/gpl.html)
import sys
try:
    import numpy as np
    if sys.version_info.major == 2:
        from cstatgen.gsl_py2 import *
        from cstatgen.boostmath_py2 import *
    else:
        from cstatgen.gsl_py3 import *
        from cstatgen.boostmath_py3 import *
except ImportError as e:
    sys.exit(e)

from spower.utils import env

class RootFinder:
    ''' halex@stackoverflow.com:
    ... It's based on chapter 4.3 from the book Numerical Methods in Engineering with Python by Jaan Kiusalaas ...'''
    def __init__(self):
        '''does nothing here. I (Gao Wang) am just wrapping halex's functions'''
    def rootsearch(self,f,a,b,dx):
        x1 = a; f1 = f(a)
        x2 = a + dx; f2 = f(x2)
        while f1*f2 > 0.0:
            if x1 >= b:
                return None,None
            x1 = x2; f1 = f2
            x2 = x1 + dx; f2 = f(x2)
        return x1,x2

    def bisect(self,f,x1,x2,switch=0,epsilon=1.0e-9):
        f1 = f(x1)
        if f1 == 0.0:
            return x1
        f2 = f(x2)
        if f2 == 0.0:
            return x2
        if f1*f2 > 0.0:
            env.logger.info('Root is not bracketed')
            return None
        n = np.ceil(np.log(np.abs(x2 - x1)/epsilon)/np.log(2.0))
        for i in range(n):
            x3 = 0.5*(x1 + x2); f3 = f(x3)
            if (switch == 1) and (np.abs(f3) >np.abs(f1)) and (np.abs(f3) > np.abs(f2)):
                return None
            if f3 == 0.0:
                return x3
            if f2*f3 < 0.0:
                x1 = x3
                f1 = f3
            else:
                x2 =x3
                f2 = f3
        return (x1 + x2)/2.0

    def roots(self,f, a, b, eps=1e-6):
        '''Example:
        >>> f=lambda x:x*math.cos(x-4)
        >>> roots(f, -3, 3)
        '''
        env.logger.info('The roots on the interval [%f, %f] are:' % (a,b))
        while 1:
            x1,x2 = rootsearch(f,a,b,eps)
            if x1 != None:
                a = x2
                root = bisect(f,x1,x2,1)
                if root != None:
                    pass
                    env.logger.info('%f' % np.around(root,-int(np.log(eps, 10))))
            else:
                env.logger.info('Done')
                break


class CATTP:
    '''Power of Cochran-Armitage Trend Test'''
    def __init__(self, x, p, q, n1, n2):
        # data coding
        self.X = np.array(x)
        # sample size
        self.N1 = float(n1)
        self.N2 = float(n2)
        self.N = float(n1 + n2)
        # genotype frequency in cases
        self.p = np.array(p) 
        # genotype frequency in ctrls
        self.q = np.array(q)

    def getSigma0(self, sigma1, mu1):
       '''see proof in Appendix A of Freidlin et al 2002'''
       return np.sqrt(sigma1 ** 2 + mu1 ** 2)
       
    def getMu1(self, X, N1, N2, N, p, q):
        '''Equation 4 of Freidlin et al 2002'''
        return np.sum(X * (p - q) * N1 * N2 / N ** 2)

    def getSigma1(self, X, N1, N2, N, p, q):
        '''Equation 4 of Freidlin et al 2002.
        Matrix notation, Pfeiffer and Gail 2003'''
        Sp = self.correlationMatrix(p)
        Sq = self.correlationMatrix(q)
        sigma1 = (N1 * N2 / N ** 3) * \
            (N1 * np.dot(np.dot(X.T, Sp), X) + N2 * np.dot(np.dot(X.T, Sq), X)) 
        # equivalent non-matrix notation, easier to compute but less intuitive
        # sigma1 = (N1 * N2 / N ** 3) * \
        #     ( (N1 * (np.sum(X ** 2 * p) - np.sum(X * p) ** 2)) + \
        #           (N2 * (np.sum(X ** 2 * q) - np.sum(X * q) ** 2)) ) 
        return np.sqrt(sigma1) 

    def power(self, alpha, alternative = 'two-sided'):
        '''Equation 5 of Freidlin et al 2002'''
        sigma1 = self.getSigma1(self.X, self.N1, self.N2, self.N, self.p, self.q)
        mu1 = self.getMu1(self.X, self.N1, self.N2, self.N, self.p, self.q) 
        sigma0 = self.getSigma0(sigma1, mu1)
        z = gsl_cdf_ugaussian_Qinv(alpha / 2.0)
        power = 1 - gsl_cdf_ugaussian_P((z * sigma0 - np.sqrt(self.N) * mu1) / sigma1) + \
            gsl_cdf_ugaussian_P((-1.0 * z * sigma0 - np.sqrt(self.N) * mu1) / sigma1)
        return power

    ### ###
    
    def correlationMatrix(self, p):
        '''correlation matrix where m_ii = pi * (1-pi) and m_ik = -pi*pk'''
        m = np.empty((len(p),len(p),))
        for i in range(len(p)):
            for k in range(len(p)):
                if k == i:
                    m[i,k] = p[i] * (1 - p[i])
                else:
                    m[i,k] = p[i] * p[k] * -1.0
        return m

class HomogeneityTP:
    '''Power of Chisquare test for homogeneity, via non-centrality parameter
    H_0 is all input p's are equal to p0'''
    def __init__(self, p0, p, n):
        '''input is an array of proportions and sample sizes
        for binomial case (i = 2, j = 2) for example:
        binomial parameters under the null
        p0 (1 x j) = [0.1, 0.9]
        binomial parameters under the alternative
        for each population group
        p (i x j) = [[0.15,0.85],
                     [0.05,0.95]]
        number of observations
        n (1 x i) = [200, 200]
        where:
        assert len(p0) == len(p[0])
        assert len(p) == len(n)
        '''
        assert len(p0) == len(p[0])
        assert len(p) == len(n)
        self.p0 = p0
        self.p = p
        self.n = [float(x) for x in n]

    def power(self, alpha):
        '''non-centrality parameter based power calculation for homogeneity test
        - Power and Sample Size for Approximate Chi-Square Tests, William C. Guenther,
          The American Statistician, Vol. 31, No. 2 (May, 1977), pp. 83-85 
        - Zar, Jerrold H. Biostatistical Analysis. Upper Saddle River, NJ: Prentice Hall, 1999.
        - http://www.boost.org/doc/libs/1_52_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/stat_tut/weg/nccs_eg/nccs_power_eg.html
        '''
        # code up section 4 of William's notation
        n = np.sum(self.n)
        ncs = []
        for j, pj in enumerate(self.p0):
            cij = []
            pi = []
            for i, ni in enumerate(self.n):
                cij.append(np.sqrt(n) * (self.p[i][j] - pj))
                pi.append(ni/n)
            ncs.append(
                (np.sum([x**2 * y for x, y in zip(cij, pi)]) - np.sum([x*y for x, y in zip(cij, pi)]) ** 2) / pj
                )
        nc = np.sum(ncs)
        if nc == 0:
            # the boost nc chisq will result in negative prob. if nc = 0
            nc += 1.0 / np.finfo(np.float).max
        # power for chi-squared test, from my boostmath library 
        df = (len(self.n) - 1 ) * (len(self.p) - 1)
        cs = gsl_cdf_chisq_Qinv(alpha, df)
        return boost_cdf_non_central_chisq_Pinv(cs, nc, df)

class SimpleQtTP:
    '''sigma is the variance explained by additive effects at the marker of interest under an additive model.
    Test: linear model such as E(Yi) = μ + βG * Gi, with
    E(Yi) denoting the expected phenotypic value for each individual,
    μ denotine an overall mean, βG denoting the estimated per genotype effect,
    Gi denoting the observed genotype for each individual (coded as 0, 1 or 2)
    http://genome.sph.umich.edu/wiki/Power_Calculations:_Quantitative_Traits'''
    def __init__(self, n, sigma):
        self.nc = n * sigma
    def power(self, alpha):
        cs = gsl_cdf_chisq_Qinv(alpha, 1)
        return boost_cdf_non_central_chisq_Pinv(cs, nc, 1)

class TwoSampleTTP:
    '''two sample t-test for mean differences between populations 
       "p" is proportion of the 2nd group of the total populations'''
    def __init__(self, p, shift, sd1 = 1.0, sd2 = 1.0):
        self.p = p
        self.sd1 = sd1
        self.sd2 = sd2
        self.delta = shift * sd2

    def power(self, n, alpha, alternative = 'two-sided'):
        # use pooled variances here
        if alternative == 'two-sided':
            self.delta = np.absolute(self.delta)
            self.t = np.absolute(self.delta) / np.sqrt(self.sd1**2 / (n * self.p) + self.sd2**2 / (n * (1 - self.p)))
            z = gsl_cdf_ugaussian_Pinv(alpha/2) + self.t
            return gsl_cdf_ugaussian_P(z)
        else:
            raise NotImplementedError('not implemented for {0}'.format(alternative))
        # z = gsl_cdf_ugaussian_Qinv(alpha/2)
        # return 1.0 - gsl_cdf_ugaussian_P(z - self.t) + gsl_cdf_ugaussian_P(-1.0*z - self.t)

    def sample_size(self, power, alpha, alternative = 'two-sided'):
        if self.delta == 0.0:
            return float('nan')
        if alternative == 'two-sided':
            self.delta = np.absolute(self.delta)
            self.t = gsl_cdf_ugaussian_Pinv(power) - gsl_cdf_ugaussian_Pinv(alpha/2)
            return (self.sd1**2 / self.p + self.sd2**2 / (1 - self.p)) / (np.absolute(self.delta) / self.t) ** 2 
        else:
            raise NotImplementedError('not implemented for {0}'.format(alternative))

class OneSampleTTP:
    '''power and sample size for two-tailed test of the population mean, the one-sample t-test
    http://www.stata-journal.com/article.html?article=st0062
    T_df(.|t) is cumulative dist. of noncentral t-distribution with df degrees of freedom
    One sided test power: \( 1-\beta=T_{n-1}(t_{\alpha, n-1}|\frac{\delta\sqrt{n}}{\sigma}) \)
    One sided test sample size: \( n=\sigma^2 \frac{(t_\alpha + t_\beta)^2}{\delta^2} \)
    Two sided test power: \( 1-\beta=T_{n-1}(t_{\alpha/2, n-1}|\frac{|\delta|\sqrt{n}}{\sigma}) \)
    Two sided test sample size: \( n=\sigma^2 \frac{(t_{\alpha/2 }+ t_\beta)^2}{\delta^2} \)

    * Note: for two-sided power there is a strict interpretation instead of the formula above,
    which is also the implementation of this program
    
    http://en.wikipedia.org/wiki/Noncentral_t-distribution#Use_in_power_analysis

    '''
    def __init__(self, delta, sd = 1.0):
        self.delta = float(delta)
        self.sd = float(sd)

    def power(self, n, alpha, alternative = 'two-sided'):
        if alternative == 'two-sided':
            qu = gsl_cdf_tdist_Qinv(alpha/2.0, n - 1)
            return 1 + boost_cdf_non_central_tdist_Pinv(qu, np.sqrt(n) * self.delta / self.sd, n - 1) - \
                boost_cdf_non_central_tdist_Pinv(-1.0*qu, np.sqrt(n) * self.delta / self.sd, n - 1)
        else:
            raise NotImplementedError('not implemented for {0}'.format(alternative))

    def sample_size(self, power, alpha, alternative = 'two-sided'):
        if alternative == 'two-sided':
            # use normal approx. for now (otherwise have to use non-central tdist which my boostmath does not have for now)
            return (self.sd / self.delta) ** 2 * \
                (gsl_cdf_ugaussian_Pinv(alpha/2.0) + gsl_cdf_ugaussian_Pinv(1 - power)) ** 2
        else:
            raise NotImplementedError('not implemented for {0}'.format(alternative))

class OneSampleZTP:
    '''power and sample size for one-sided Z test, see TDTMendianPowerZtest.R
    '''
    def __init__(self):
        self.tside = 1
        self.sigma = 1

    def power(self, effect, alpha):
        '''
        power <- pnorm(qnorm(1-(alpha/tside))*sigma, mean = effect, sd = sigma, lower.tail = FALSE)
        '''
        z = gsl_cdf_ugaussian_Qinv(alpha/self.tside)
        # print(effect, alpha, z, gsl_cdf_ugaussian_Q(z*self.sigma - effect))
        return gsl_cdf_ugaussian_Q(z*self.sigma - effect)

    def effect(self, power, alpha):
        '''
        effect <- (qnorm(1-(alpha/tside))+qnorm(power))*sigma
        '''
        z = gsl_cdf_ugaussian_Qinv(alpha/self.tside)
        # print(power, alpha, z, gsl_cdf_ugaussian_Pinv(power))
        return (z + gsl_cdf_ugaussian_Pinv(power))*self.sigma


class ProportionNormalTP:
    '''Test for equal proportions between groups, via normal approx.
    '''
    def __init__(self, p1, p2, p):
        self.p1 = p1
        self.p2 = p2
        self.v1 = p1 * (1.0 - p1)
        self.v2 = p2 * (1.0 - p2)
        self.po = p
        # ratio of arm 2 to arm 1
        self.r = (1-p)/p
        self.p = (self.p1 + self.r * self.p2) / (self.r + 1)

    def power(self, n, alpha):
        '''
        Reference:
        modified from the "power.prop.test" function in R (which assumes self.po = 0.5)
        '''
        qu = gsl_cdf_ugaussian_Qinv(alpha/2)
        d = np.absolute(self.p1 - self.p2)
        pbar = self.p1 * self.po + self.p2 * (1 - self.po)
        wvbar = pbar * (1 - pbar) * ( 1 / (n * self.po) + 1 / (n * (1 - self.po)) )
        np.absolute(self.p1 - self.p2) / np.sqrt(self.v1 / (n * self.po) + self.v2 / (n * (1 - self.po)))
        s = np.sqrt(self.v1 / (n * self.po) + self.v2 / (n * (1 - self.po)))
        return gsl_cdf_ugaussian_P((d - qu * np.sqrt(wvbar)) / s) + gsl_cdf_ugaussian_Q((d + qu * np.sqrt(wvbar)) / s)
        # This test is from the pwr R package "pwr.2p2n.test"
        # h = np.absolute(2 * (np.arcsin(np.sqrt(self.p1)) - np.arcsin(np.sqrt(self.p2))))
        # return gsl_cdf_ugaussian_Q(gsl_cdf_ugaussian_Qinv(alpha/2) - h * np.sqrt(self.po * (1 - self.po) / n)) + \
        #     gsl_cdf_ugaussian_P(gsl_cdf_ugaussian_Pinv(alpha/2) - h * np.sqrt(self.po * (1 - self.po) / n))
 
    def sample_size(self, power, alpha):
        '''
        Reference:     

        Casagrande, Pike and Smith (1978) Biometrics 34: 483-486
        Fleiss, Tytun and Ury (1980) Biometrics 36: 343-346
        '''
        delta = np.absolute(self.p1-self.p2)
        za = gsl_cdf_ugaussian_Pinv(alpha/2)
        zb = gsl_cdf_ugaussian_Pinv(1 - power)
        mp = ( za * np.sqrt((self.r+1)*self.p*(1-self.p)) + zb * np.sqrt(self.r * self.v1 + self.v2) ) ** 2 / (self.r* delta**2)
        m = mp / 4.0 * (1 + np.sqrt(1 + 2 * (self.r + 1)/(self.r * mp * delta))) ** 2
        return (self.r+1) * m

      

if __name__ == '__main__':
    # t = CATTP([0,1,2], [0.1,0.4,0.5], [0.1,0.35,0.55], 200, 200)
    # print t.power(0.05)
    # t = HomogeneityTP([0.1, 0.9], [[0.15,0.85],[0.05,0.95]], [200,200])
    # print t.power(0.05)
    # t = HomogeneityTP([0.35, 0.65], [[0.3,0.7],[0.4,0.6], [0.35, 0.65]],[200,200, 200])
    # print t.power(0.05)
    # p = 0.34213693120684896
    # t = HomogeneityTP([p, 1-p], [[p,1-p],[p,1-p]],[500,500])
    # print t.power(0.05)
    # t = ProportionNormalTP(0.5,0.3,0.5)
    # print t.sample_size(0.8, 0.05)
    # print t.power(205, 0.05)
    # t = OneSampleTTP(0.3, 2.5)
    # print t.power(35, 0.05)
    # print t.sample_size(0.005, 0.05)
    #
    t = TwoSampleTTP(0.2, 0.3)
    print t.power(200, 0.05)
    # this should be the same as above: 0.1 = 0.2 * 0.5
    t = OneSampleTTP(0.06)
    print t.power(200, 0.05)

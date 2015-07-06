import logging
import numpy as np

logger=logging.getLogger(__file__)



class ExponentialDistribution:
    def __init__(self, lam, te):
        self.lam=lam
        self.te=te

    def sample(self, now, rng):
        return now+rng.exponential(self.lam)

    def hazard_integral(self, t0, t1):
        return self.lam*(t1-t0)

    def implicit_hazard_integral(self, xa, t0):
        return t0+xa/self.lam

    def enabling_time(self):
        return self.te

class WeibullDistribution:
    def __init__(self, lam, k, te, shift):
        self.lam=lam
        self.k=k
        self.te=te
        self.delta=shift

    def sample(self, now, rng):
        U=rng.uniform(0, 1)
        l=self.lam
        k=self.k
        d=now - (self.te - self.delta)
        value=0
        if d>0:
            value=l*np.pow(-np.log(1-U) + np.pow(d/l, k), 1/k)-d
        else:
            value=-d+l*np.pow(-np.log(1-U), 1/k)
        return now + value

    def hazard_integral(self, t0, t1):
        return ( np.pow((t1-self.te)/self.lam, self.k)
            -np.pow((t0-self.te)/self.lam, self.k) )

    def implicit_hazard_integral(self, xa, t0):
        return self.te+ self.lam*np.pow(xa+np.pow((t0-self.te)/self.lam,
                self.k), 1/self.k)

    def enabling_time(self):
        return self.te

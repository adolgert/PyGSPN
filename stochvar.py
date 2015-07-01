import logging
import math

logger=logging.getLogger(__file__)

class StochasticVariable:
    def __init__(self):
        pass

class Exponential:
    def __init__(self, lam):
        self.lam=lam
        self.te=None
    def build(self, te):
        x.te=te
        return x
    def sample(self, t0, rng):
        return te-math.log(uniform(rng))/self.lam
    def hazard_integral(self, t0, t1):
        return self.lam*(t1-t0)
    def implicit_hazard_integral(self, xa, t0):
        return t0+xa/self.lam
    def enabling_time(self):
        return self.te

class Weibull:
    def __init__(self, lam, k, shift):
        """
        lam is base hazard
        te is enabling
        shift moves weibull dist later in time.
        """
        self.lam=lam
        self.k=k
        self.shift=shift
        self.te=None
    def build(self, te):
        x.te=te
        return x
    def sample(self, t0, rng):
        d=t0-(te+self.shift)
        value=0
        U=uniform(rng)
        if d>0:
            value=self.lam*math.pow(-math.log(1-U)
                +math.pow(d/self.lam, self.k), 1/self.k)-d
        else:
            value=-d+self.lam*math.pow(-math.log(1-U), 1/self.k)
        return t0+value
    def implicit_hazard_integral(self, xa, t0):
        return te+self.lam*math.pow(xa
            +math.pow((t0-te)/self.lam, self.k), 1/self.k)
    def enabling_time(self):
        return self.te
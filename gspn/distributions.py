import logging
import numpy as np
import scipy.stats

logger=logging.getLogger(__file__)


def anderson_sample_tester(distribution, now, cnt, rng):
    """
    This checks whether hazard_integral and implicit_hazard_integral
    work by using them to sample a distribution.
    """
    samples=np.zeros(cnt)
    maxdiff=0.0
    for i in range(cnt):
        interval=-np.log(rng.uniform(0, 1))
        last_mod=now
        firing_time=distribution.implicit_hazard_integral(interval, last_mod)
        next_fire=firing_time
        stops=np.linspace(now, firing_time, num=rng.randint(0,5))
        for stop in stops[:-1]:
            time_penalty=distribution.hazard_integral(last_mod, stop)
            interval-=time_penalty
            next_fire=distribution.implicit_hazard_integral(interval, stop)
            last_mod=stop
        maxdiff=max(abs(next_fire-firing_time), maxdiff)
        samples[i]=next_fire
    return samples


class ExponentialDistribution(object):
    """
    This represents an exponential distribution.
    .. math::

        F(t) = 1-e^{-\int_0^t \lambda(s) ds}
    """
    def __init__(self, lam, te):
        self.lam=lam
        self.te=te

    def sample(self, now, rng):
        return now+rng.exponential(scale=1.0/self.lam)

    def hazard_integral(self, t0, t1):
        return self.lam*(t1-t0)

    def implicit_hazard_integral(self, xa, t0):
        return t0+xa/self.lam

    def loglikelihood(self, t0, tf):
        return -self.lam*(tf-t0)

    def enabling_time(self):
        return self.te


class WeibullDistribution(object):
    """
    This is a Weibull distribution.
    """
    def __init__(self, lam, k, te, shift):
        self.lam=lam
        self.k=k
        self.te=te
        self.delta=shift

    def sample(self, now, rng):
        logger.debug("WeibullDistribution.sample l={0}, k={1}, te={2}".format(
            self.lam, self.k, self.te))
        U=rng.uniform(0, 1)
        l=self.lam
        k=self.k
        d=now - (self.te - self.delta)
        value=0
        if d>0:
            value=l*np.power(-np.log(1-U) + np.power(d/l, k), 1/k)+self.te
        else:
            value=l*np.power(-np.log(1-U), 1/k)+self.te
        return now + value

    def hazard_integral(self, t0, t1):
        logger.debug("WeibullDistribution.hazard l={0}, k={1}, te={2}".format(
            self.lam, self.k, self.te))
        return ( np.power((t1-self.te)/self.lam, self.k)
            -np.power((t0-self.te)/self.lam, self.k) )

    def implicit_hazard_integral(self, xa, t0):
        t1=t0 + self.lam * xa**(1/self.k)
        logger.debug(("WeibullDistribution.implicit l={0}, k={1}, te={2} "+
            "xa={3} t0={4} t1={5}").format(
            self.lam, self.k, self.te, xa, t0, t1))
        return t1

    def enabling_time(self):
        return self.te


class GammaDistribution(object):
    """
    This is a gamma distribution with a shape and a rate,
    not a shape and a scale.
    Given a Gamma function,
    .. math::

        \Gamma(t)=\int_0^\infty x^{t-1}e^{-x}dx

    and the (lower) incomplete gamma function,
    .. math::

        \gamma(x;\alpha)=\int_0^x t^{\alpha-1}e^{-t}dt

    the CDF of the gamma distribution is
    .. math::

        F(x;\alpha, \beta)=\frac{\gamma(\alpha, \beta x)}{\Gamma(\alpha)}

    The PDF is
    .. math::

        f(x)=\frac{\beta^{\alpha}}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x}

    This is sampled with possible left censoring. 
    """
    def __init__(self, alpha, beta, te):
        self.alpha=alpha
        self.beta=beta
        self.te=te

    def sample(self, now, rng):
        """
        Sampling accounts for time shift and uses given random
        number generator.
        """
        U=rng.uniform(low=0, high=1, size=1)
        d=now-self.te
        if d>0:
            cumulative=scipy.stats.gamma.cdf(x=d, a=self.alpha,
                scale=1.0/self.beta, loc=0)
            return scipy.stats.gamma.isf(1-U*(1-cumulative)-cumulative,
                self.alpha, scale=1.0/self.beta, loc=0) + self.te
        else:
            return scipy.stats.gamma.isf(1-U,
                self.alpha, scale=1.0/self.beta, loc=0) + self.te

    def hazard_integral(self, t0, t1):
        """
        Our tools include

         - gammainc(a, x), normalized lower incomplete
         - gammaincinv(a, y), gammainc(a, x)=y
         - gammaincc(a, x), 1-gammainc(a, x)
         - gammainccinv(a, y), gammaincc(a, x)=y
        """
        return np.log(
            (1-scipy.special.gammainc(self.alpha, self.beta*(t0-self.te)))/
            (1-scipy.special.gammainc(self.alpha, self.beta*(t1-self.te)))
            )

    def implicit_hazard_integral(self, xa, t0):
        quad=1-np.exp(-xa)*(1-scipy.special.gammainc(self.alpha,
                self.beta*(t0-self.te)))
        return self.te+scipy.special.gammaincinv(self.alpha, quad)/self.beta

    def loglikelihood(self, t0, tf):
        t0e=t0-self.te
        logf=(np.log(np.power(self.beta, self.alpha)*
            np.power(t0e, self.alpha-1))
            -self.beta*t0e-scipy.special.gammaln(self.alpha))
        integrated=self.hazard_integral(self.te, t0)
        return logf + integrated

    def enabling_time(self):
        return self.te



class UniformDistribution(object):
    """
    Uniform distribution between a and b, offset by an enabling time te.
    """
    def __init__(self, a, b, te):
        """
        te is an absolute time.
        The earliest firing time is >te+a.
        The latest firing time is <=te+b.
        """
        self.a=a
        self.b=b
        self.te=te

    def sample(self, now, rng):
        """
        Sampling accounts for time shift and uses given random
        number generator.
        """
        if now<=self.te+self.a:
            return rng.uniform(low=self.te+self.a, high=self.te+self.b)
        elif now<=self.te+self.b:
            return rng.uniform(low=now, high=self.te+self.b)
        else:
            return np.float("nan")

    def hazard_integral(self, t0, t1):
        """
        Integrate the hazard, taking into account when the uniform
        interval starts and stops.
        """
        t0e=t0-self.te
        t1e=t1-self.te
        if t0e>self.b:
            return 0
        if t1e<self.a:
            return 0
        low=max(self.a, t0e)
        high=min(self.b, t1e)
        numerator=np.log(1-(t0e-self.a)/(self.b-self.a))
        denominator=np.log((self.b-t1e)/(self.b-self.a))
        return numerator-denominator

    def implicit_hazard_integral(self, xa, t0):
        Ft=1-np.exp(-xa)
        t0e=t0-self.te
        low=max(t0e, self.a)
        r=self.te+low*(1-Ft) + self.b*Ft
        return r

    def loglikelihood(self, t0, tf):
        t0e=t0-self.te
        tfe=tf-self.te
        if tfe<self.a or te>=self.b:
            return np.double("nan")
        t0e=max(t0e, self.a)
        ln_pdf=np.log((tf-self.a)/(self.b-self.a))
        int_hazard=self.hazard_integral(self.te, t0)
        return ln_pdf + int_hazard

    def enabling_time(self):
        return self.te

# LogLogistic
# Gaussian
# Histogram

class PiecewiseLinearDistribution(object):
    """
    This is a piecewise linear hazard, not a piecewise linear probability.
    """
    def __init__(self, times, hazards, enabling_time):
        assert(times[0]<1e-6)
        self.b=times
        self.w=hazards
        self.te=enabling_time
        self.partial_sum=np.zeros((len(times),), dtype=np.double)

        total=0
        for idx in range(len(times)-1):
            self.partial_sum[idx]=total
            total+=0.5*(hazards[idx]+hazards[idx+1])*(times[idx+1]-times[idx])
        self.partial_sum[len(times)-1]=total

    def sample(self, now, rng):
        """
        Sampling accounts for time shift and uses given random
        number generator.
        """
        raise RuntimeError("Cannot sample this yet.")

    def hazard_integral(self, t0, t1):
        """
        Integrate the hazard, taking into account when the uniform
        interval starts and stops.
        """
        t0e=t0-self.te
        t1e=t1-self.te
        if t1e<self.b[0]:
            return 0
        # Take the array of hazards at times and chop it off
        # on the left at t0e.
        b, w=(np.copy(self.b), np.copy(self.w))
        left=np.where(b<=t0e)[0][-1]
        b=b[left:]
        w=w[left:]
        b[0]=t0e
        # Chop off the right after t1e.
        right=np.where(b<t1e)[0][-1]
        b=b[:(right+1)]
        w=w[:(right+1)]
        # Add a point back for the end point.
        b=np.hstack([ b, np.array([t1e])])
        w=np.hstack([ w, np.array([w[-1]])])

        total=0.0
        for idx in range(b.shape[0]-1):
            total+=w[idx]*(b[idx+1]-b[idx])
        return total


    def implicit_hazard_integral(self, xa, t0):
        t0e=t0-self.te
        # Take the array of hazards at times and chop it off
        # on the left at t0e.
        b, w=(np.copy(self.b), np.copy(self.w))
        left=np.where(b<=t0e)[0][-1]
        b=b[left:]
        w=w[left:]
        b[0]=t0e
        logger.debug("Piecewise bleft {0} {1}".format(b, w))
        b-=t0e   # Now it starts at 0.
        logger.debug("Piecewise bleft {0} {1}".format(b, w))

        idx=0
        t1e=0
        while xa>0:
            x0=b[idx]
            if idx<b.shape[0]-1:
                x1=b[idx+1]
            else:
                x1=float("inf")
            if w[idx]>0:
                t1e=x0+xa/w[idx]
            else:
                t1e=float("inf")
            logger.debug("Piecewise midcalc x0 {0} x1 {1} t1e {2}".format(
                x0, x1, t1e))
            if t1e<x1:
                xa=0
            else:
                xa-=(x1-x0)*w[idx]
            idx+=1

        return t1e+t0

    def loglikelihood(self, t0, tf):
        return None
 
    def enabling_time(self):
        return self.te



class EmpiricalDistribution(object):
    """
    This distribution is used to collect samples and then
    compare the estimate of the sampled distribution with
    some known distribution.

    Wikipedia explains this. Kolmogorov-Smirnov test.
    https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
    """
    def __init__(self, samples):
        self.samples=samples

    def c_alpha(self):
        """
        The hypothesis is that the two distributions are different.
        The null hypothesis is rejected at level alpha if
        Dnn*sqrt(nn/(n+n))>c_alpha. This returns a table of c_alpha.
        The upshot is that you should look for values from the
        two comparison functions that are less than two.
        """
        return [[0.1, 1.22],
                [0.05, 1.36],
                [0.025, 1.48],
                [0.01, 1.63],
                [0.005, 1.73],
                [0.001, 1.95]]

    def compare_empirical(self, other):
        """
        Compare two empirical distributions using Kolmogorov-Smirnov.
        This returns the Kolmogorov-smirnov goodness-of-fit,
        sqrt(n)*D_n, where D_n is the Kolmogorov-smirnov statistic.
        """
        self.samples.sort()
        other.samples.sort()
        acnt=len(self.samples)
        aidx=-1
        bcnt=len(other.samples)
        bidx=-1
        maxdiff=0.0
        while aidx+1<acnt and bidx+1<bcnt:
            if bidx+1==bcnt:
                v=self.samples[aidx+1]
                while aidx+1!=acnt and self.samples[aidx+1]==v:
                    aidx+=1
                maxdiff=max(maxdiff, abs((aidx+1)/acnt - (bidx+1)/bcnt))
            elif aidx+1==acnt:
                v=other.samples[bidx+1]
                while bidx+1!=bcnt and other.samples[bidx+1]==v:
                    bidx+=1
                maxdiff=max(maxdiff, abs((aidx+1)/acnt - (bidx+1)/bcnt))
            elif self.samples[aidx+1]<other.samples[bidx+1]:
                v=self.samples[aidx+1]
                while aidx+1!=acnt and self.samples[aidx+1]==v:
                    aidx+=1
                maxdiff=max(maxdiff, abs((aidx+1)/acnt - (bidx+1)/bcnt))
            else:
                v=other.samples[bidx+1]
                while aidx+1!=acnt and self.samples[aidx+1]==v:
                    aidx+=1
                while bidx+1!=bcnt and other.samples[bidx+1]==v:
                    bidx+=1
                maxdiff=max(maxdiff, abs((aidx+1)/acnt - (bidx+1)/bcnt))
        return maxdiff*np.sqrt(acnt*bcnt/(acnt+bcnt))

    def compare_theoretical(self, cdf):
        """
        Compare this estimated CDF with a function which returns
        the CDF at a point. Use Kolmogorov-Smirnov to return a
        statistic which should be stable with increasing
        numbers of samples.
        This returns the Kolmogorov-smirnov goodness-of-fit,
        sqrt(n)*D_n, where D_n is the Kolmogorov-smirnov statistic.
        """
        self.samples.sort()
        acnt=len(self.samples)
        aidx=-1
        maxdiff=0.0
        while aidx!=acnt:
            v=self.samples[aidx+1]
            while self.samples[aidx+1]==v:
                aidx+=1
            maxdiff=max(maxdiff,
                abs(cdf(self.samples[aidx])-(aidx+1)/acnt))
        return maxdiff*np.sqrt(acnt)


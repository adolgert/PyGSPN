import logging
import itertools
#from enum import Enum
import numpy as np
import distributions
import llcp
import sample
import runner

logger=logging.getLogger(__file__)

class FarmState(object):
    susceptible=1
    latent=2
    subclinical=3
    clinical=4
    recovered=5

class FarmPlace(object):
    def __init__(self, farm):
        self.farm=farm
        self.state=FarmState.susceptible


class FarmABTransition:
    def __init__(self, farm_place, a, b, distribution):
        self.place=farm_place
        self.a=a
        self.b=b
        self.dist=distribution

    def depends(self):
        return [self.place]

    def affected(self):
        return [self.place]

    def enabled(self, now):
        if self.place.state==self.a:
            return (True, self.dist(now))
        else:
            return (False, None)

    def fire(self):
        self.place.state=self.b


class InfectiousIntensity:
    """
    This is the rate part of a transition.
    """
    def __init__(self, farm_place):
        self.place=place
    def depends(self):
        return [self.place]
    def affected(self):
        return []
    def intensity(self, now):
        if self.place.state in (FarmState.subclinical, FarmState.clinical):
            return 1
        return None

class InfectPartial:
    """
    This is the action-part of a transition.
    """
    def __init__(self, farm):
        self.farm=farm
    def depends(sefl):
        return [self.farm.place]
    def affected(self):
        return [self.farm.place]
    def enabled(self):
        if self.place.state in (FarmState.susceptible,):
            return True
        else:
            return False
    def fire(self, now):
        self.place.state=FarmState.latent

class Farm(object):
    """
    This is a scenario model for disease state within a farm.
    """
    def __init__(self, name):
        self.name=name
        self.place=FarmPlace(self)

    def write_places(self, writer):
        writer.add_place(self.place)

    def write_transitions(self, writer):
        t=FarmABTransition(self.place, FarmState.latent, FarmState.subclinical,
            lambda now: distributions.ExponentialDistribution(0.5, now))
        writer.add_transition(t)
        t=FarmABTransition(self.place, FarmState.subclinical, FarmState.clinical,
            lambda now: distributions.ExponentialDistribution(0.5, now))
        writer.add_transition(t)
        t=FarmABTransition(self.place, FarmState.clinical, FarmState.recovered,
            lambda now: distributions.ExponentialDistribution(0.5, now))
        writer.add_transition(t)

    def infectious_intensity(self):
        return InfectiousIntensity(self)

    def infection_partial(self):
        return InfectionPartial(self)

    def infectious(self):
        return self.place.state in (FarmState.subclinical,
            FarmState.clinical)

    def infectious_depends(self):
        return [self.place]

    def susceptible(self):
        return self.place.state in (FarmState.susceptible,)

    def susceptible_depends(self):
        return [self.place]

    def infect(self):
        self.place.state=FarmState.latent

    def infect_affected(self):
        return [self.place]


class InfectTransition(object):
    def __init__(self, intensity, action):
        self.intensity=intensity
        self.action=action

    def depends(self):
        deps=self.intensity.depends()
        deps.extend(self.action.depends())
        return deps

    def affected(self):
        return self.action.affected()

    def enabled(self, now):
        intensity=self.intensity.enabled()
        if intensity is not None and self.action.enabled():
            rate=0.5*intensity
            return (True, distributions.ExponentialDistribution(rate, now))
        else:
            return (False, None)

    def fire(self):
        self.action.infect()


class InfectNeighbor(object):
    """
    This is a scenario model for infection of one farm by another.
    """
    def __init__(self, farma, farmb):
        self.farma=farma
        self.farmb=farmb

    def write_places(self, writer):
        pass

    def write_transitions(self, writer):
        writer.add_transition(InfectTransition(
            self.farma.infectious_intensity(), self.farmb.infection_partial()))
        writer.add_transition(InfectTransition(
            self.farmb.infectious_intensity(), self.farma.infection_partial()))


def Build(individual_cnt):
    net=llcp.LLCP()
    farms=list()
    for add_idx in range(individual_cnt):
        farms.append(Farm(add_idx))

    for f in farms:
        f.write_places(net)
        f.write_transitions(net)

    for a, b in itertools.combinations(farms, 2):
        infect=InfectNeighbor(a, b)
        infect.write_places(net)
        infect.write_transitions(net)

    initial_idx=0
    farms[initial_idx].place.state=FarmState.latent
    return net


########################################
# This is the part that runs the SIR
def observer(transition, when):
    if isinstance(transition, FarmABTransition):
        print("AB {0} {1} {2} {3}".format(transition.place.farm.name,
            transition.a, transition.b, when))
    elif isinstance(transition, InfectTransition):
        print("Infect {0} {1} {2}".format(transition.intensity.place.farm.name,
            transition.action.place.farm.name,
            when))
    else:
        print("Unknown transition {0}".format(transition))
    return True

def test_farm():
    rng=np.random.RandomState()
    rng.seed(33333)
    net=Build(10)
    sampler=sample.NextReaction(net, rng)
    run=runner.RunnerFSM(sampler, observer)
    run.init()
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_farm()
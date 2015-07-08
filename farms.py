import logging
import itertools
#from enum import Enum
import numpy as np
import scipy.spatial.distance as distance
import distributions
import llcp
import sample
import runner
import point_process

logger=logging.getLogger(__file__)


##############################################################
# Within-farm disease model
##############################################################
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
    This is part of a transition which contributes to the hazard
    rate through a factor called an intensity.
    """
    def __init__(self, farm):
        self.farm=farm
        self.place=farm.place
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
    This is the action-part of a transition. It isn't a whole transition.
    """
    def __init__(self, farm):
        self.farm=farm
        self.place=farm.place
    def depends(self):
        return [self.farm.place]
    def affected(self):
        return [self.farm.place]
    def enabled(self):
        if self.place.state in (FarmState.susceptible,):
            return True
        else:
            return False
    def fire(self):
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
        return InfectPartial(self)


##############################################################
# Kernel-based neighbor infection
##############################################################
class InfectTransition(object):
    """
    This transition brings together pieces from different models
    into a full transition.
    """
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
        intensity=self.intensity.intensity(now)
        if intensity is not None and self.action.enabled():
            rate=0.5*intensity
            return (True, distributions.ExponentialDistribution(rate, now))
        else:
            return (False, None)

    def fire(self):
        self.action.fire()


class InfectNeighborModel(object):
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


##############################################################
# Movement restrictions
##############################################################
class RestrictionPlace(object):
    def __init__(self):
        self.restricted_date=None

class RestrictionTransition(object):
    def __init__(self, farm, restriction_place):
        self.farm=farm
        self.place=restriction_place
        self.detectable=farm.detectable_intensity()
        self.te=None

    def depends(self):
        dep=self.detectable.depends()
        dep.extend(self.place)
        return dep

    def affected(self):
        return [self.place]

    def enabled(self, now):
        if (self.detectable.intensity() is not None) and (self.restricted_date
                is not None):
            if self.te is None:
                self.te=now
            return (True, distributions.ExponentialDistribution(0.05, self.te))
        else:
            return (False, None)

    def fire(self, now):
        self.place.restricted_date=now

class RestrictionIntensity(object):
    def __init__(self, place):
        self.place=place
    def depends(self):
        return [self.place]
    def affected(self):
        return []
    def intensity(self, now):
        if self.restriction_place is None:
            return None
        else:
            # some kind of increasing restriction
            return (now-self.place.restricted_date)


class MovementRestrictionsModel(object):
    def __init__(self, landscape):
        self.landscape=landscape
        self.restriction_place=RestrictionPlace()

    def write_places(self, writer):
        writer.add_place(self.restriction_place)

    def write_transitions(self, writer):
        for farm in self.landscape.farms:
            t=RestrictionTransition(farm, self.restriction_place)
            writer.add_transition(t)
    def restriction_intensity(self):
        return RestrictionIntensity(self.restriction_place)



##############################################################
# Direct Movement Model
##############################################################
class DirectTransition(object):
    """
    Represents direct contact where a farm makes a certain
    number of shipments a day.
    """
    def __init__(self, farm, model):
        self.farm=farm
        self.model=model
    def depends(self):
        pass

class DirectModel(object):
    def __init__(self, landscape):
        self.landscape=landscape




class Landscape(object):
    def __init__(self):
        self.farm_locations=point_process.thomas_point_process_2D(
            5, 0.1, 5, (0, 1, 0, 1))
        individual_cnt=self.farm_locations.shape[0]
        self.distances=distance.squareform(distance.pdist(self.farm_locations,
            "euclidean"))
        self.farms=list()
        for add_idx in range(individual_cnt):
            self.farms.append(Farm(add_idx))


def Build():
    net=llcp.LLCP()
    landscape=Landscape()
    farms=landscape.farms

    for f in farms:
        f.write_places(net)
        f.write_transitions(net)

    for a, b in itertools.combinations(farms, 2):
        infect=InfectNeighborModel(a, b)
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
    net=Build()
    sampler=sample.NextReaction(net, rng)
    run=runner.RunnerFSM(sampler, observer)
    run.init()
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_farm()
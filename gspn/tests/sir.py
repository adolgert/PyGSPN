import logging
import numpy as np
import gspn

logger=logging.getLogger(__file__)


class CountPlace:
    def __init__(self, id):
        self.id=id
        self.count=0

class RecoverTransition:
    def __init__(self, i, r):
        self.i=i
        self.r=r

    def depends(self):
        return [self.i]

    def affected(self):
        return [self.i, self.r]

    def enabled(self, now):
        if self.i.count>0:
            return (True, gspn.ExponentialDistribution(1.0, now))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.i.count=0
        self.r.count=1


class InfectTransition:
    def __init__(self, source, susceptible, infected):
        self.i0=source
        self.s1=susceptible
        self.i1=infected

    def depends(self):
        return [self.i0, self.s1]

    def affected(self):
        return [self.s1, self.i1]

    def enabled(self, now):
        if self.i0.count>0 and self.s1.count>0:
            return (True, gspn.ExponentialDistribution(0.5, now))
        else:
            return (False, None)

    def fire(self, now, rng):
        self.s1.count=0
        self.i1.count=1


def BuildSIR(individual_cnt):
    net=gspn.LLCP()
    places=dict()
    for add_idx in range(individual_cnt):
        for disease_state in ['s', 'i', 'r']:
            p=CountPlace((add_idx, disease_state))
            places[p.id]=p
            net.add_place(p)
            p.count=0

    for internal_idx in range(individual_cnt):
        for internal_trans in ['r']:
            t=RecoverTransition(places[(internal_idx, 'i')],
                places[(internal_idx, 'r')])
            net.add_transition(t)

    for source_idx in range(individual_cnt):
        for target_idx in range(individual_cnt):
            if source_idx!=target_idx:
                t=InfectTransition(places[(source_idx, 'i')],
                    places[(target_idx, 's')], places[(target_idx, 'i')])
                net.add_transition(t)

    initial_idx=0
    for s_idx in range(individual_cnt):
        if s_idx==initial_idx:
            places[(s_idx, 'i')].count=1
        else:
            places[(s_idx, 's')].count=1
    net.init()
    return net


########################################
# This is the part that runs the SIR
def observer(transition, when):
    if isinstance(transition, RecoverTransition):
        print("Recover {0} {1}".format(transition.i.id, when))
    else:
        print("Infect {0} {1}".format(transition.s1.id, when))
    return True

def test_sir():
    rng=np.random.RandomState()
    rng.seed(33333)
    net=BuildSIR(10)
    sampler=gspn.NextReaction(net, rng)
    run=gspn.RunnerFSM(sampler, observer)
    run.init()
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_sir()

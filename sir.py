import logging
import numpy as np
import distributions
import llcp
import sample
import runner

logger=logging.getLogger(__file__)


class CountPlace:
    def __init__(self, id):
        self.id=id
        self.count=0

class RecoverTransition:
    def __init__(self):
        self.place=dict()
        self.dep=dict()

    def enable(self, now):
        if self.dep['i'].count>0:
            return (True, distributions.ExponentialDistribution(1.0, now))
        else:
            return (False, None)

    def fire(self):
        self.place['i'].count=0
        self.place['r'].count=1


class InfectTransition:
    def __init__(self):
        self.place=dict()
        self.dep=dict()

    def enable(self, now):
        if self.dep['i'].count>0 and self.dep['s'].count>0:
            return (True, distributions.ExponentialDistribution(0.5, now))
        else:
            return (False, None)

    def fire(self):
        self.place['s'].count=0
        self.place['n'].count=1


def BuildSIR(individual_cnt):
    net=llcp.LLCP()
    places=dict()
    for add_idx in range(individual_cnt):
        for disease_state in ['s', 'i', 'r']:
            p=CountPlace((add_idx, disease_state))
            places[p.id]=p
            net.add_place(p)
            p.count=0

    for internal_idx in range(individual_cnt):
        for internal_trans in ['r']:
            t=RecoverTransition()
            t.dep['i']=places[(internal_idx, 'i')]
            t.place['i']=places[(internal_idx, 'i')]
            t.place['r']=places[(internal_idx, 'r')]
            net.add_transition(t)

    for source_idx in range(individual_cnt):
        for target_idx in range(individual_cnt):
            if source_idx!=target_idx:
                t=InfectTransition()
                t.dep['s']=places[(target_idx, 's')]
                t.dep['i']=places[(source_idx, 'i')]
                t.place['s']=places[(target_idx, 's')]
                t.place['n']=places[(target_idx, 'i')]
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
    if 'i' in transition.place.keys():
        print("Recover {0} {1}".format(transition.place['i'].id, when))
    else:
        print("Infect {0} {1}".format(transition.place['s'].id, when))
    return True

def test_sir():
    rng=np.random.RandomState()
    rng.seed(33333)
    net=BuildSIR(10)
    sampler=sample.NextReaction(net, rng)
    run=runner.RunnerFSM(sampler, observer)
    run.init()
    run.run()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_sir()

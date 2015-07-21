from unittest import TestCase
import logging
import numpy as np
import sample
import sir
import runner

logger = logging.getLogger(__file__)

def observer(transition, when):
    if 'i' in transition.place.keys():
        print("Recover {0} {1}".format(transition.place['i'].id, when))
    else:
        print("Infect {0} {1}".format(transition.place['s'].id, when))
    return True

def test_sir():
    rng=np.random.RandomState()
    rng.seed(33333)
    net=sir.BuildSIR(10)
    sampler=sample.NextReaction(net, rng)
    run=runner.RunnerFSM(sampler, observer)
    run.init()
    run.run()


class TestSIR(TestCase):
    def test_sir_runs(self):
        logging.basicConfig(level=logging.DEBUG)
        test_sir()

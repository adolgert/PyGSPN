import logging
import numpy as np
import llcp
import sample

logger = logging.getLogger(__file__)

def test_sir():
    rng=np.random.RandomState()
    rng.seed(33333)
    net=llcp.BuildSIR(10)
    sampler=sample.NextReaction(net, rng)
    sampler.init()
    transition, when=sampler.next()
    while transition is not None:
        logger.debug("{0} {1}".format(transition, when))
        sampler.fire(transition, when)
        transition, when=sampler.next()


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    test_sir()

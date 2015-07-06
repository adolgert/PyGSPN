import numpy as np
import llcp
import sample



def test_sir():
    rng=np.random.RandomState()
    rng.seed(33333)
    net=llcp.BuildSIR(10)
    sampler=sample.NextReaction(net, rng)
    sampler.init()
    next=sampler.next()
    while next[0] is not None:
        print(next)
        net.fire(next[0], next[1])
        next=sampler.next()


if __name__ == "__main__":
    test_sir()

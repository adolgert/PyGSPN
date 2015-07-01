import logging
import heapq

logger=logging.getLogger(__file__)


class FirstReaction:
    """
    Gillespie's first reaction method.
    """
    def __init__(self):
        pass

    def next(self, system, rng):
        self.least=[None, float("inf")]
        self.rng=rng
        self.now=system.current_time()
        system.enabled_transitions(self._sample_trans)
        return self.least

    def _sample_trans(self, transition, distribution):
        trial_time=distribution(self.now, self.rng)
        if trial_time < self.least[1]:
            self.least[0]=transition
            self.least[1]=trial_time
        # else forget it


class NextReactionRecord:
    def __init__(self):
        self.remaining_exponential_interval=0.0
        self.last_modification_time=0.0
        self.heap_entry=[0, self]

class NextReaction:
    def __init__(self):
        self.heap=list()

    def next(self, system, rng):
        if len(self.heap)>0:
            return heap[0]
        return None

    def fire(self, system, transition):
        pass

    def enable(self, transition, distribution, now, rng):
        heapq.heappush(self.heap, entry)

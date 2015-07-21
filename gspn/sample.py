import logging
import math
import gspn.pairing_heap

logger=logging.getLogger(__file__)


class FirstReaction:
    """
    Gillespie's first reaction method.
    """
    def __init__(self, system, rng):
        self.system=system
        self.rng=rng

    def init(self):
        self.system.init()

    def next(self):
        self.least=[None, float("inf")]
        self.system.enabled_transitions(self._sample_trans)
        return self.least

    def fire(self, transition, when):
        self.system.fire(transition, when)

    def _sample_trans(self, transition, distribution, now):
        trial_time=distribution.sample(now, self.rng)
        if trial_time < self.least[1]:
            self.least[0]=transition
            self.least[1]=trial_time
        # else forget it


class NextReactionRecord:
    def __init__(self):
        self.remaining_exponential_interval=None
        self.last_modification_time=0.0
        self.heap_entry=None


class NextReaction:
    def __init__(self, system, rng):
        self.priority=gspn.pairing_heap.pairing_heap()
        self.system=system
        self.rng=rng

    def init(self):
        self.system.init(self._observe)

    def next(self):
        if not self.priority.empty():
            v=self.priority.peek()
            return (v[1], v[0])
        return (None, None)


    def fire(self, transition, when):
        self.system.fire(transition, when, self.rng, self._observe)


    def _observe(self, transition, olddist, newdist, now):
        if newdist is not None:
            if "_nr" in transition.__dict__:
                record=transition._nr
                if record.heap_entry is not None:
                    time_penalty=olddist.hazard_integral(
                        record.last_modification_time, now)
                    record.remaining_exponential_interval-=time_penalty
                # else: transition was disabled previously
                when_fire=newdist.implicit_hazard_integral(
                    record.remaining_exponential_interval, now)
                if record.heap_entry is not None:
                    if when_fire<record.heap_entry._item[0]:
                        self.priority.adjust_key(record.heap_entry,
                            (when_fire, transition))
                    else:
                        self.priority.delete(record.heap_entry)
                        record.heap_entry=self.priority.insert(
                            (when_fire, transition))
                else:
                    record.heap_entry=self.priority.insert(
                        (when_fire, transition))
                record.last_modification_time=now
            else:
                interval=-math.log(self.rng.uniform(0, 1))
                firing_time=newdist.implicit_hazard_integral(
                    interval, now)
                transition._nr=NextReactionRecord()
                transition._nr.remaining_exponential_interval=interval
                transition._nr.last_modification_time=now
                transition._nr.heap_entry=self.priority.insert(
                    (firing_time, transition))
        else:
            record=transition._nr
            self.priority.delete(transition._nr.heap_entry)
            transition._nr.heap_entry=None
            time_penalty=olddist.hazard_integral(
                transition._nr.last_modification_time, now)
            transition._nr.remaining_exponential_interval-=time_penalty
            transition._nr.last_modification_time=now


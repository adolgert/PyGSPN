import logging
import sample

logger=logging.getLogger(__file__)

# This design makes mixin versions of the Place and Transition
# meaning that their base classes make them intrusive entries in
# an adjacency list to define a graph. The subclasses
# define whatever the particular process needs.
class Place:
    def __init__(self):
        self._adjacency=list()

class Transition:
    def __init__(self):
        self.place=dict()
        self.dep=dict()
        self._distribution=None

class LLCP:
    """
    Long-lived competing processes.
    """
    def __init__(self):
        self._current_time=0.0
        self.p=list()
        self.t=list()

    def add_place(self, place):
        self.p.append(place)

    def add_transition(self, transition):
        self.t.append(transition)
        for p in transition.place.values():
            p._adjacency.append(transition)
        for d in transition.dep.values():
            d._adjacency.append(transition)

    def init(self):
        self._initial_enable()

    def current_time(self):
        return self._current_time

    def fire(self, transition, when, report=None):
        self._current_time=when
        transition.fire()
        transition._distribution=None
        #self._initial_enable()
        self._incremental_update(transition, report)

    def enabled_transitions(self, functor):
        for t in self.t:
            if t.distribution is not None:
                functor(t, t.distribution)

    def _initial_enable(self):
        for t in self.t:
            enabled, dist=t.enable(self._current_time)
            if enabled:
                t.distribution=dist
            else:
                t.distribution=None

    def _incremental_update(self, fired_transition, report):
        affected_transitions=set()
        for p in fired_transition.place.values():
            affected_transitions.update(p._adjacency)
        for t in affected_transitions:
            was_enabled=t.distribution is not None
            enabled, dist=t.enable(self._current_time)
            if enabled:
                t.distribution=dist
            else:
                t.distribution=None
            if was_enabled or enabled:
                report(t, was_enabled, enabled)



class CountPlace(Place):
    def __init__(self, id):
        super(CountPlace, self).__init__()
        self.id=id
        self.count=0

class RecoverTransition(Transition):
    def __init__(self):
        super(RecoverTransition, self).__init__()
        self._nr=sample.NextReactionRecord()

    def enable(self, now):
        if self.place['i'].count>0:
            return (True, lambda x, y: now+y.exponential(1.0))
        else:
            return (False, None)

    def fire(self):
        self.place['i'].count=0
        self.place['r'].count=1


class InfectTransition(Transition):
    def __init__(self):
        super(InfectTransition, self).__init__()
        self._nr=sample.NextReactionRecord()

    def enable(self, now):
        if self.place['i'].count>0 and self.place['s'].count>0:
            return (True, lambda x, y: now+y.exponential(0.5))
        else:
            return (False, None)

    def fire(self):
        self.place['s'].count=0
        self.place['n'].count=1


def BuildSIR(individual_cnt):
    net=LLCP()
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
            t.place['i']=places[(internal_idx, 'i')]
            t.place['r']=places[(internal_idx, 'r')]
            net.add_transition(t)

    for source_idx in range(individual_cnt):
        for target_idx in range(individual_cnt):
            if source_idx!=target_idx:
                t=InfectTransition()
                t.place['s']=places[(target_idx, 's')]
                t.place['i']=places[(source_idx, 'i')]
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

import logging
import sample

logger=logging.getLogger(__file__)

# This design makes mixin versions of the Place and Transition
# meaning that their base classes make them intrusive entries in
# an adjacency list to define a graph. The subclasses
# define whatever the particular process needs.

class LLCP:
    """
    Long-lived competing processes.
    """
    def __init__(self):
        self._current_time=0.0
        self.p=list()
        self.t=list()

    def add_place(self, place):
        place._adjacency=list() # inject
        self.p.append(place)

    def add_transition(self, transition):
        """
        A transition must have two dictionary members,
        place and dep, which are places used in firing
        and dependencies for determining hazard rates.
        """
        transition._distribution=None # inject
        self.t.append(transition)
        for d in transition.dep.values():
            d._adjacency.append(transition)

    def init(self, report=None):
        self._initial_enable(report)

    def current_time(self):
        return self._current_time

    def fire(self, transition, when, report=None):
        self._current_time=when
        transition.fire()
        if report is not None:
            report(transition, transition._distribution, None,
                self._current_time)
        transition._distribution=None
        self._incremental_update(transition, report)

    def enabled_transitions(self, functor):
        for t in self.t:
            if t._distribution is not None:
                functor(t, t._distribution, self._current_time)

    def _initial_enable(self, report):
        for t in self.t:
            enabled, dist=t.enable(self._current_time)
            if enabled:
                if report is not None:
                    report(t, None, dist, self._current_time)
                t._distribution=dist
            else:
                t._distribution=None

    def _incremental_update(self, fired_transition, report):
        affected_transitions=set()
        for p in fired_transition.place.values():
            affected_transitions.update(p._adjacency)
        for t in affected_transitions:
            was_enabled=t._distribution is not None
            enabled, dist=t.enable(self._current_time)
            if report is not None and (was_enabled or enabled):
                report(t, t._distribution, dist, self._current_time)
            t._distribution=dist


import logging

logger=logging.getLogger(__file__)

# This design makes mixin versions of the Place and Transition
# meaning that their base classes make them intrusive entries in
# an adjacency list to define a graph. The subclasses
# define whatever the particular process needs.

class LLCPTransition(object):
    def depends(self):
        """
        Upon which places does this transition depend, for
        whether it is enabled, for its hazard rate, and for any
        choices it makes when it fires.
        """
        return set()
    def affected(self):
        """
        When this transition fires, which places does it change.
        """
        return set()

class LLCP:
    """
    Long-lived competing processes.
    """
    def __init__(self):
        self._current_time=0.0
        self.t=list()

    def add_place(self, place):
        place._adjacency=list() # inject

    def add_transition(self, transition):
        """
        A transition must have two dictionary members,
        place and dep, which are places used in firing
        and dependencies for determining hazard rates.
        """
        transition._distribution=None # inject
        self.t.append(transition)
        for d in transition.depends():
            d._adjacency.append(transition)

    def init(self, report=None):
        self._initial_enable(report)

    def current_time(self):
        return self._current_time

    def fire(self, transition, when, rng, report=None):
        self._current_time=when
        transition.fire(when, rng)
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
            enabled, dist=t.enabled(self._current_time)
            if enabled:
                if report is not None:
                    report(t, None, dist, self._current_time)
                t._distribution=dist
            else:
                t._distribution=None

    def _incremental_update(self, fired_transition, report):
        affected_transitions=set()
        for p in fired_transition.affected():
            affected_transitions.update(p._adjacency)
        for t in affected_transitions:
            was_enabled=t._distribution is not None
            enabled, dist=t.enabled(self._current_time)
            if report is not None and (was_enabled or enabled):
                report(t, t._distribution, dist, self._current_time)
            t._distribution=dist


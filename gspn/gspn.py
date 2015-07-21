"""
How do we implement an API for making modules that build
a GSPN and then run the transitions within?
"""
import logging


logger=logging.getLogger(__file__)


class Place:
    """
    A place has a unique key and holds tokens.
    """
    def __init__(self, key):
        self.key=key
        self.tokens=list()


class TokenFlow:
    """
    A TokenFlow represents movement of tokens among places to
    and from a transition. It encapsulates
    1. stoichiometry
    2. policy of which token is taken (first, last, random, specific)
    3. policy of how and whether new tokens are created
    """
    def __init__(self, take_edge, take_cnt, give_edge, give_cnt):
        self.take_edge=take_edge
        self.take_cnt=take_cnt
        self.give_edge=give_edge
        self.give_cnt=give_cnt
    def fire(self, localstate):
        take=list()
        for i in range(f.take_cnt):
            take.append(localstate[f.take_edge].pop())
        for j in range(f.give_cnt):
            if length(take)>0:
                localstate[f.give_edge].append(take.pop())

class TransitionModifier:
    """
    A TransitionModifier is part of a transition used to change
    its enabling or hazard rate.
    """
    def __init__(self, token_flows):
        self.flows=token_flows

    def places(self):
        return None

    def enabled(self, globalstate, localstate, offset, te, t0, rng):
        """
        The offset is an integer offset into which integer-indexed
        part of the local state corresponds to places required by
        this TransitionModifier.
        """
        return False

    def hazard_modifier(self):
        return 1.0


class Transition:
    """
    Transition keys aren't unique.
    """
    def __init__(self, key, token_flows, stochastic_variable):
        self.key=key
        self.flows=token_flows
        self.stochastic=stochastic_variable

    def enabled(self, globalstate, localstate, t0, rng):
        """
        globalstate would have the scenario.
        localstate is an integer-indexed set of edges to local marking.
        t0 is current time.
        rng is random number generator
        """
        return False

    def fire(self, globalstate, localstate, t0, rng):
        """
        globalstate would have the scenario.
        localstate is an integer-indexed set of edges to local marking.
        t0 is current time.
        rng is random number generator
        """
        pass



class StoichiometricTransition:
    """
    This transition enforces stoichiometric coefficients and no more.
    But there may be more than one token type, each with different
    stoichiometry, expressed in a different flow.
    When enabled() is called, if the transition was already enabled,
    it compares the current invariant with the previous value to 
    determine whether to set a new te.
    """
    def __init__(self, key, token_flows, stochastic_variable):
        """
        key is the non-unique id for this transition.
        token_flows is a list of TokenFlows.
        stochastic_variable is a stochastic variable.
        """
        self.key=key
        self.flows=token_flows
        self.stochastic=stochastic_variable

    def enabled(self, globalstate, localstate, t0, rng):
        """
        globalstate would have the scenario.
        localstate is an integer-indexed set of edges to local marking.
        te is enabling time.
        t0 is current time.
        rng is random number generator
        returns distribution, and an invariant.
        """
        for f in self.flows:
            if length(localstate[f.take_edge])<f.take_cnt:
                return None

        return (self.stochastic.build(te), [])

    def fire(self, globalstate, localstate, t0, rng):
        """
        globalstate would have the scenario.
        localstate is an integer-indexed set of edges to local marking.
        t0 is current time.
        rng is random number generator
        """
        for f in self.flows:
            f.fire(localstate)


class ModifiedTransition(Transition):
    """
    This is a transition that has been modified by TransitionModifier objects.
    """
    def __init__(self, base_transition):
        self.base_transition=base_transition
        self.modifiers=list()
        self.key=self.base_transition.key
        self.flows=base_transition.flows

    def enabled(self, globalstate, localstate, te, t0, rng):
        b_enabled=self.base_transition.enabled()
        for m in self.modifiers:
            if not m.enabled(globalstate, localstate, te, t0, rng):
                b_enabled=False
        modifier=1.0
        for m in self.modifiers:
            modifier*=m.hazard_modifier(globalstate, localstate, t0, rng)
        hazard=base_transition.hazard(globalstate, localstate, t0, rng, modifier)
        return hazard, b_enabled

    def fire(self, globalstate, localstate, t0, rng):
        pass


class GSPN:
    """
    A GSPN is a bipartite graph between places and transitions
    and a directed dependency graph for transitions.
    Places have unique keys. Transitions have non-unique keys.
    Transitions and places both have integer ids which are indices
    into the adjacency list.
    """
    def __init__(self):
        # places will be a list of lists, an adjacency list
        self.p=list()
        # transitions will be a list of lists, an adjacency list
        self.t=list()
        self.p_key_to_id=dict()

    #### Construction
    def add_place(self, pkey):
        pid=len(self.p)
        self.p.append([Place(pkey), list()])
        self.p_key_to_id[pkey]=pid
        return pid

    def add_transition(self, transition, place_keys, dep_keys):
        """
        XXX dep_keys not used. Looks like the place should just
        point to those transitions which depend on it.
        """
        tid=len(self.t)
        # each entry is the transition, the places, the dependencies.
        transition_entry=[transition, list(), list()]
        self.t.append(transition_entry)
        for pkey in place_keys:
            pid=self.p_key_to_id[pkey]
            self.p[pid][1].append(tid)
            transition_entry[1].append(pid)
        for dkey in dep_keys:
            transition_entry[2].append(self.p_key_to_id[dkey])

    #### Access
    def transition_places(self, tid):
        """
        tid is an integer transition id.
        returns list of place ids, integer indices into places.
        """
        return self.t[tid][1]

    def transition_dependency(self, tid):
        return self.t[tid][2]

    def dependent_transitions(self, tid):
        dep_trans=set()
        for pid in self.t[tid][2]:
            dep_trans.update(self.p[pid][1])
        return dep_trans


class LRCPProcess:
    """
    This process uses the GSPN for places and transitions
    but doesn't restrict what enables a transition or what the transition
    does when it fires.
    """
    def __init__(self):
        self.gspn=GSPN()
        self.current_time=0

    # Builder methods pass through to GSPN.
    def add_place(self, pkey):
        self.gspn.add_place(self, pkey)

    def add_transition(self, transition, place_keys):
        self.gspn.add_transition(transition, place_keys)

    # Set initial marking.
    def add_token(self, pkey, token):
        self.gspn.p[self.gspn.p_key_to_id[pkey]][0].tokens.append(token)


    def transition_distribution(self, transition_entry):
        transition=transition_entry[0]
        places=transition_entry[1]
        dist=transition.distribution(self.globalstate, self.places,
            self.current_time, rng)
        return dist


    def enabled_transitions(self, functor):
        enabled_t=list()
        results=[0]*length(self.gspn.t)
        for idx, tentry in enumerate(self.gspn.t):
            results[idx]=functor(self.transition_distribution(tentry))
        return results

    


class GSPNProcess:
    def __init__(self):
        self.gspn=GSPN()
        self.current_time=0

    # Builder methods pass through to GSPN.
    def add_place(self, pkey):
        self.gspn.add_place(self, pkey)

    def add_transition(self, transition, place_keys, dep_keys):
        self.gspn.add_transition(transition, place_keys, dep_keys)

    # Set initial marking.
    def add_token(self, pkey, token):
        self.gspn.p[self.gspn.p_key_to_id[pkey]][0].tokens.append(token)

    # These are for sampling.
    def stoichiometry_satisfied(self, tid):
        """
        Check that input stoichiometry is satisfied.
        """
        satisfied=True
        transition_entry=self.gspn.t[tid]
        places=transition_entry[1]
        for f in transition_entry[0].flows:
            if length(self.gspn.p[f.take_edge][0].tokens)<f.take_cnt:
                satisfied=False
        return satisfied

    def transition_distribution(self, tid, te):


    def enabled_transitions(self):
        enabled_t=list()
        for t in self.gspn.t:
            if self.stoichiometry_satisfied(tid):


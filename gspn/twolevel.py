import logging
import gspn


logger=logging.getLogger(__file__)

# Given the structure above for the GSPN, the Scenario
# will contain models which interact with the GSPN in two
# passes. First the models in the Scenario build the GSPN.
# Then, when transitions run, they call objects within the
# Scenario in order to help decide when they are enabled
# and disabled. For the second pass, there must be no
# mutable state in the Scenario that is modified.

class Scenario:
    def __init__(self):
        self.unit_id=list()
        self.unit_disease=list()
        self.disease_models=list()
        self.spread_models=list()

# These are the rules for how to name places given a model in this scenario.
def ModelPlaceToScenarioPlace(unit_idx, place):
    return (unit_idx, unit_idx, place)

def ModelPlaceToScenarioPlace(unit_a_idx, unit_b_idx, place):
    return (unit_a_idx, unit_b_idx, place)


class Unit:
    """
    A Unit is a flyweight on the scenario.
    """
    def __init__(self, scenario, idx):
        self.scenario=scenario
        self.idx=idx

    def indirect_infectious_ask(self):
        disease_model_idx=self.scenario.unit_disease[self.idx]
        disease=self.scenario.disease_models[disease_model_idx]
        places=disease.build_indirect_infectious(self.idx)

    def indirect_infectious(self):
        return False


class DiseaseModel:
    def __init__(self):
        pass

class HPAI(DiseaseModel):
    def __init__(self):
        self.name="HPAI"
        self.states=["susceptible", "clinical", "recovered"]

class FMD(DiseaseModel):
    def __init__(self):
        self.name="FMD"
        self.states=["susceptible", "subclinical", "clinical", "recovered"]

    def build_indirect_infectious(self, unit_idx):
        pass

    def place_of_state(self, unit_idx, state):
        return (unit_idx, unit_idx, state)

class AirborneSpread:
    def __init__(self):
        pass

class DirectContact:
    def __init__(self):
        pass


def BuildScenario():
    s=Scenario();
    s.unit_id.extend([0, 7, 19])
    s.disease_models.append(HPAI(), FMD())
    s.unit_disease.extend([0, 0, 1])
    u=Unit(s, 1)

    gspn=GSPN()

    for disease_unit in s.GetUnits():
        disease_unit.add_disease(gspn)

    for transport_unit in s.GetUnits():
        transport_unit.add_transport(gspn)

    return (s, gspn)


if __name__ == "__main__":
    BuildScenario()

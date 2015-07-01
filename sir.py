import logging
import gspn
import stochvar

logger=logging.getLogger(__file__)

def BuildSIR(individual_cnt):
    net=gspn.GSPNProcess()
    for add_idx in range(individual_cnt):
        for disease_state in ['s', 'i', 'r']:
            net.add_place((disease_state, add_idx))

    for internal_idx in range(individual_cnt):
        for internal_trans in ['r']:
            flows=[gspn.TokenFlow(0, 1, 1, 1),]
            recover=gspn.StoichiometricTransition(internal_trans, flows,
                stochvar.Weibull(1.0, 2.0, 0.0))
            net.add_transition(recover,
                [('i', internal_idx), ('r', internal_idx)], [])

    for source_idx in range(individual_cnt):
        for target_idx in range(individual_cnt):
            flows=[gspn.TokenFlow(0, 1, 2, 1), gspn.TokenFlow(1, 1, 1, 1)]
            infect=gspn.StoichiometricTransition('i', flows,
                stochvar.Exponential(1.0))
            net.add_transition(infect,
                [('s', target_idx), ('i', source_idx), ('i', target_idx)], [])

    return net


if __name__ == "__main__":
    individual_cnt=10
    net=BuildSIR(individual_cnt)

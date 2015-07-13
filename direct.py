"""
Compare a transliteration of NAADSM 3.2's direct contact spread
with an algorithm which produces a stochastic, continuous-time
direct spread. This runs both on the same scenario and produces
plots to compare.
"""
import logging
import numpy as np
import scipy.spatial as spatial
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import point_process

logger=logging.getLogger(__file__)

def farm_locations(farm_cnt):
    # Put farms in space with some non-uniform distribution.
    kappa=5
    mu=farm_cnt/kappa
    pts=point_process.thomas_point_process_2D(kappa, 0.2, mu, (0, 1, 0, 1))
    cnt=len(pts)
    logger.debug("There are {0} points".format(cnt))

    # Create a zone somewhere in the bounds and zone all inside.
    zone_center=np.random.rand(2).reshape(1,2)
    zone_distance=spatial.distance.cdist(pts, zone_center).flatten()
    logger.debug("zone distance {0}".format(zone_distance))
    zoned=(zone_distance<0.3).flatten()

    # Quarantine everything near the zone center, including outside it.
    p=np.zeros(cnt)
    p[zone_distance<0.6]=0.4
    quarantined=np.random.binomial(1, p)==1
    logger.debug("zoned {0}".format(zoned))
    logger.debug("quarantined {0}".format(quarantined))
    possible_starts=np.where(np.logical_and(
            np.logical_not(zoned), np.logical_not(quarantined))==True)[0]
    start_idx=np.random.choice(possible_starts)
    return pts, zoned, quarantined, start_idx


def naadsm_direct(scenario, run_cnt=1):
    """
    This tries to replicate A4.1 Direct contact spread from
    NAADSM Model Specification 1.2.1.
    """
    pts, zoned, quarantined, start_idx=scenario
    uniform_kernel_max=0.4
    start_loc=pts[start_idx].reshape(1,2)

    #3a Make a list of those that can be a recipient.
    recipients_idx=np.where(quarantined==False)[0]
    recipients_idx=recipients_idx[recipients_idx!=start_idx]
    recipients=pts[recipients_idx]
    zoned_recipients=zoned[recipients_idx]

    epsilon=0.01
    distances=spatial.distance.cdist(recipients, start_loc).flatten()

    hits=np.zeros(pts.shape[0], dtype=np.int)
    for run_idx in range(run_cnt):
        #3c Sample a number, distance, from the movement distance distribution.
        # Using uniform here.
        travel_distance=np.random.uniform(0, uniform_kernel_max)

        #3d Choose the closest
        closest_idx=np.argmin(np.abs(distances-travel_distance))
        closest_distance=distances[closest_idx]
        possibles=np.where(distances-closest_distance<epsilon)[0]
        logger.debug("possibles {0}".format(possibles))

        #3e If forbidden by zoning, drop it.
        unzoned=possibles[zoned_recipients[possibles]==False]
        if np.any(unzoned):
            hits[recipients_idx[np.random.choice(unzoned)]]+=1
        else:
            pass # nothing sent
    return hits


def faithful(scenario, run_cnt=1):
    """
    This is an attempt at an equivalent continuous-time model.
    """
    pts, zoned, quarantined, start_idx=scenario
    uniform_kernel_max=0.4
    start_loc=pts[start_idx].reshape(1,2)

    #3a Make a list of those that can be a recipient.
    recipients_idx=np.where(quarantined==False)[0]
    recipients_idx=recipients_idx[recipients_idx!=start_idx]
    recipients=pts[recipients_idx]
    zoned_recipients=zoned[recipients_idx]

    #3d For each recipient, find whether it will really receive.
    epsilon=0.01
    distances=spatial.distance.cdist(recipients, start_loc).flatten()
    sort_idx=np.argsort(distances)
    # What is the probability of sending a truck, ignoring zones?
    prob_basket=np.zeros(len(distances), dtype=np.double)
    inner=0.0
    for didx in range(1, len(distances)-1):
        ptidx=sort_idx[didx]
        if inner<uniform_kernel_max:
            outer=0.5*(distances[ptidx]+distances[sort_idx[didx+1])
            prob_basket[ptidx]=outer-inner
            inner=outer
    prob_basket/=np.sum(prob_basket)

    # Now take trucks that would have gone to zones and either
    # don't send them or give them to the nearby available farms.
    adj_prob=np.zeros(len(distances), dtype=np.double)
    for ridx in range(1, len(distances)-1):
        possibles=np.where(distances-distances[ridx]<epsilon)[0]
        unzoned=possibles[zoned_recipients[possibles=False]]
        if len(unzoned)>0:
            adj_prob[unzoned]+=prob_basket[ridx]/len(unzoned)

    # Losing trucks makes the probability less than one. Readjust
    # and take that probability out by reducing the number of times
    # we draw a value.
    total_prob=np.sum(adj_prob)
    adj_prob /= total_prob
    indices=np.array(range(len(adj_prob)))

    hits=np.zeros(pts.shape[0], dtype=np.int)
    bins=np.array(range(len(adj_prob)+1))
    hits=np.histogram(np.random.choice(indices, p=adj_prob,
        size=run_cnt*total_prob), bins=bins)[0]
    return hits    


def plot_comparison(scenario, hits0, hits1):
    pts, zoned, quarantined, start_idx=scenario
    max_size=200
    min_size=10

    hitpts=pts[hits0>0]
    hithits=hits0[hits0>0]
    hithits=hithits*max_size/np.max(hithits) + min_size

    pts2=np.delete(pts, start_idx, 0)
    hits2=np.delete(hits0, start_idx)
    zoned2=np.delete(zoned, start_idx)
    quarantined2=np.delete(quarantined, start_idx)
    mispts=pts2[hits2==0]
    miszone=zoned2[hits2==0]
    misquar=quarantined2[hits2==0]

    colors=np.zeros(len(mispts))
    colors[:]=0.9
    colors[miszone]=0.3
    colors[misquar]=0.5
    #logger.debug("plotsizes {0}".format(sizes))

    plt.scatter(hitpts[:,0], hitpts[:,1], s=hithits, c="blue", marker='o',
        alpha=0.5)
    plt.scatter(mispts[:,0], mispts[:,1], s=min_size, c=colors, marker='s',
        alpha=0.5)
    plt.scatter(pts[start_idx,0], pts[start_idx,1], s=min_size, c="orange", marker='*',
        alpha=1)
    plt.savefig("scatter.pdf", format="pdf")


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    cnt=100
    scenario=farm_locations(cnt)
    run_cnt=10000
    hits=naadsm_direct(scenario, run_cnt)
    print("hits {0}".format(hits))
    plot_comparison(scenario, hits, hits)


import logging
import numpy as np
import scipy.stats

logger=logging.getLogger(__file__)

def poisson_point_process_2D(lam, bounds):
    """
    lam is the intensity.
    bounds are (xlow, xhigh, ylow, yhigh).
    Used this article.
    http://connor-johnson.com/2014/02/25/spatial-point-processes/
    """
    N=scipy.stats.poisson(lam*(bounds[1]-bounds[0])*(bounds[3]-bounds[2])).rvs()
    logger.debug("{0} points".format(N))
    x=scipy.stats.uniform.rvs(loc=bounds[0], scale=(bounds[1]-bounds[0]),
        size=((N, 1)))
    y=scipy.stats.uniform.rvs(loc=bounds[2], scale=(bounds[3]-bounds[2]),
        size=((N, 1)))
    return np.hstack((x, y))


def thomas_point_process_2D(kappa, sigma, mu, bounds):
    """
    kappa is the intensity of the high-level process.
    sigma is the Gaussian variance.
    mu is the intensity of the low-level process.
    bounds are (xlow, xhigh, ylow, yhigh).
    Used this article.
    http://connor-johnson.com/2014/02/25/spatial-point-processes/
    """
    parents=poisson_point_process_2D(kappa, bounds)
    parent_cnt=parents.shape[0]
    list_of_arrays=list()
    sub_scale=np.array([sigma, sigma])
    for i in range(parent_cnt):
        children_cnt=scipy.stats.poisson(mu).rvs()
        pdf=scipy.stats.norm(loc=np.array(parents[i, :2]), scale=sub_scale)
        list_of_arrays.append(pdf.rvs((children_cnt, 2)))
    return np.vstack(list_of_arrays)



import numpy as np

def compute_msd(r, interval):
    result = np.zeros( (r.shape[0],) )
    count = np.zeros( (r.shape[0],) )
    n_total = r.shape[0]
    n_interval = n_total/interval
    for i in range(n_interval):
        rsqr = np.sum( (r[i*interval:]-r[i*interval])**2, axis=2 )
        print rsqr.shape, i*interval
        if i>0:
            result[:-i*interval] += rsqr.mean(axis=1)
            count[:-i*interval] += 1
        else:
            result += rsqr.mean(axis=1)
            count += 1
    mask = count > 0
    result[mask] /= count[mask]
    return result, count


import numpy as np

def random_walk(steps, N=10, dim=3):
    """Returns a dataset for N particles in dimensions dim that perform a random
    walk."""
    return np.cumsum(np.random.random_integers(-1,1,(steps,N,dim)),axis=0)



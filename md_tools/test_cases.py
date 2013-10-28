# -*- coding: utf-8 -*-
# Copyright 2013 Pierre de Buyl
#
# This file is part of md_tools
#
# md_tools is free software and is licensed under the modified BSD license (see
# LICENSE file).

import numpy as np

def random_walk(steps, N=10, dim=3):
    """Returns a dataset for N particles in dimensions dim that perform a random
    walk."""
    return np.cumsum(np.random.random_integers(-1,1,(steps,N,dim)),axis=0)

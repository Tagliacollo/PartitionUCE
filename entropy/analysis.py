from functions import *
import numpy as np

# test the entropy function
p1 = np.array([0.25, 0.25, 0.25, 0.25]) # should be an entropy 2
entropy_calc(p1)

p2 = np.array([0.0, 0.0, 0.5, 0.5]) # should be an entropy of 1
entropy_calc(p2)
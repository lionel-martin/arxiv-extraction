"""  Pickles HOWTO

To extract information from the pickles in this directory do the following
"""

import pickle
import pygsp

# pik is a dict containing 'W' and 'active_authors' the indices of the connected authors
pik = pickle.load(open(filename, 'rb'))

# You can then do:
auth_idx = pik['active_authors']
W = pik['W'][np.ix_(auth_idx, auth_idx)]
G = pygsp.graphs.Graph(W)

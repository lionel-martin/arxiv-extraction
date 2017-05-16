""" CSV2PyData
Transform CSV into large python objects.

This is step 1 for the transformation of the dataset.

"""

import pickle
import os
from collections import defaultdict


def def_dict_rec():
    return defaultdict(def_dict_rec)

temporal_links = def_dict_rec()
papers_auth = {}
authors = []

os.chdir('..')
if not os.path.isdir('pickles'):
    os.mkdir('pickles')

for line in open('pscp-data_1991-2013.csv', 'r'):
    if line[0] != '#':
        values = line.split(';')
        this_pauth = []
        for auth in values[5].split(','):
            auth_inlist = filter(lambda x: x[1] == auth, authors)
            if not len(auth_inlist):
                this_pauth.append(len(authors))
                authors.append((len(authors), auth))
            else:
                this_pauth.append(auth_inlist[0][0])

        papers_auth[values[0]] = this_pauth

pickle.dump(authors, open(os.path.join('pickles', 'authors'), 'wb'))
pickle.dump(papers_auth, open(os.path.join('pickles', 'papers_auth'), 'wb'))

for line in open('pscp-data_1991-2013.csv', 'r'):
    if line[0] != '#':
        values = line.split(';')
        date = values[0].split('/')[1][:4] if '/' in values[0] else values[0].split('.')[0]
        if isinstance(temporal_links[date[:2]][date[2:]], defaultdict):
            temporal_links[date[:2]][date[2:]] = []
        auth_refs = reduce(lambda x, y: x + y, [papers_auth[pp] for pp in values[4].split(',')]) if int(values[2]) > 0 else None
        if auth_refs is not None:
            temporal_links[date[:2]][date[2:]].append((papers_auth[values[0]], auth_refs))

pickle.dump(temporal_links, open(os.path.join('pickles', 'temporal_links'), 'wb'))

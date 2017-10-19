""" CSV2PyData
Transform CSV into large python objects.

This is step 1 for the transformation of the dataset.

"""

import pickle
import os
from collections import defaultdict, Counter


def def_dict_rec():
    return defaultdict(def_dict_rec)

temporal_links = def_dict_rec()
papers_auth = {}
authors = []
topics = []
authors_top = defaultdict(list)

os.chdir('..')
if not os.path.isdir('pickles'):
    os.mkdir('pickles')

for line in open('pscp-data_1991-2013.csv', 'r'):
    if line[0] != '#':
        values = line.split(';')

        main_topic = values[1].split(',')[0]
        top_inlist = list(filter(lambda x: x[1] == main_topic, topics))

        if len(top_inlist):
            top_id = top_inlist[0][0]
        else:
            top_id = len(topics)
            topics.append((top_id, main_topic))

        this_pauth = []
        for auth in values[5].split(','):
            auth_inlist = list(filter(lambda x: x[1] == auth, authors))
            if len(auth_inlist):
                auth_id = auth_inlist[0][0]
            else:
                auth_id = len(authors)
                authors.append((auth_id, auth))

            this_pauth.append(auth_id)
            authors_top[auth_id].append(top_id)

        papers_auth[values[0]] = this_pauth


for i in xrange(len(authors_top)):
    authors_top[i] = Counter(authors_top[i]).most_common(1)[0][0]

pickle.dump(authors, open(os.path.join('pickles', 'authors'), 'wb'))
pickle.dump(papers_auth, open(os.path.join('pickles', 'papers_auth'), 'wb'))
pickle.dump(topics, open(os.path.join('pickles', 'topics'), 'wb'))
pickle.dump(authors_top, open(os.path.join('pickles', 'authors_top'), 'wb'))

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

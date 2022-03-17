import json
from scipy.stats import binom as binomrv
import numpy as np

ps = np.atleast_2d([55, 65, 75, 85, 95])
qs = np.array([[5, 25, 75, 95]]).T

lookup_npq = dict()
for n in range(1, 41):
  block = binomrv.ppf(qs / 100, n, ps / 100)

  lookup_npq[n] = dict()
  for pi, p in enumerate(ps.ravel()):
    p = int(p)
    lookup_npq[n][p] = dict()
    for qi, q in enumerate(qs.ravel()):
      q = int(q)
      lookup_npq[n][p][q] = int(block[qi, pi])

with open('src/npq2.json', 'w') as fid:
  json.dump(lookup_npq, fid)

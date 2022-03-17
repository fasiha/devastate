import numpy as np
from scipy.special import binom as binomcoef, betaln
import pylab as plt

from confidenceIntervals import pmfNormalizingDenominator, pmfUnnormalized

plt.style.use('ggplot')
plt.ion()

pmf1 = lambda k, n, p: float(pmfUnnormalized(k, n, p) / pmfNormalizingDenominator(n, p))
pmf = np.vectorize(pmf1)
ns = [2, 5, 10]
ps = [0.55, 0.85]

fig, ax = plt.subplots(len(ns))

for pi, p in enumerate(ps):
  for ni, n in enumerate(ns):
    ks = np.linspace(0, n + 1)
    ax[ni].plot(ks, pmf(ks, n, p), label=f'{p=}', ls='--' if pi == 1 else '-')

maxy = max([max(a.get_ylim()) for a in ax])
maxx = max([max(a.get_xlim()) for a in ax])
minx = min([min(a.get_xlim()) for a in ax])
for ai, a in enumerate(ax):
  a.set_ylim([0, maxy])
  a.set_xlim([minx, maxx])
  a.set_xlabel('k')
  a.set_ylabel('PMF(k; n,p)')
  a.set_title(f'n={ns[ai]}')
  a.legend()

fig.set_figheight(6)
fig.set_figwidth(4)

fig.tight_layout()
fig.savefig('continuous_binomial_pmf.png', dpi=300)
fig.savefig('continuous_binomial_pmf.svg')

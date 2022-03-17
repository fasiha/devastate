import numpy as np
from scipy.special import binom as binomcoef, betaln
import pylab as plt

plt.style.use('ggplot')
plt.ion()

pmf = lambda k, n, p: (p**k * (1 - p)**(n - k) * binomcoef(n, k))
ns = [5, 10, 20]
ps = [0.55, 0.85]

fig, ax = plt.subplots(2)

for pi, p in enumerate(ps):
  for n in ns:
    ks = np.linspace(0, n + 1)
    ax[pi].plot(ks, pmf(ks, n, p), label=f'{n=}')

ax[0].legend()
maxy = max([max(a.get_ylim()) for a in ax])
ax[0].set_ylim([0, maxy])
ax[1].set_ylim([0, maxy])
ax[0].set_xlabel('k')
ax[1].set_xlabel('k')
ax[0].set_ylabel('PMF(k;n,p)')
ax[1].set_ylabel('PMF(k;n,p)')
ax[0].set_title(f'p={ps[0]}')
ax[1].set_title(f'p={ps[1]}')

fig.tight_layout()
fig.savefig('continuous_binomial_pmf.png', dpi=300)
fig.savefig('continuous_binomial_pmf.svg')

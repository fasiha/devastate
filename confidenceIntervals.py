# semi-related but different distribution? https://rdrr.io/cran/cbinom/man/cbinom.html

# You'll need scipy and numpy to run this
from scipy.stats import binom as binomrv
from scipy.special import binom as binomcoef, betaln
from scipy.interpolate import interp1d
import mpmath as mp
import numpy as np
import pylab as plt
import tqdm
import json

plt.ion()


def _binomln(n, k):
  "Log of scipy.special.binom calculated entirely in the log domain"
  return -betaln(1 + n - k, 1 + k) - np.log(n + 1)


pmf = lambda k, n, p: (p**k * (1 - p)**(n - k) * binomcoef(n, k))


def pmf2(k, n, p):
  logp = np.log(p)
  logq = np.log(1 - p)
  return np.exp(k * logp + (n - k) * logq + _binomln(n, k))


mpbinomcoef = lambda x, y: 1 / ((x + 1) * mp.beta(x - y + 1, y + 1))
mppmf = lambda k, n, p: (p**k * (1 - p)**(n - k) * mpbinomcoef(n, k))

n, p = 5, 0.55
n, p = 4, 0.55
n, p = 1, 0.55
n, p = 2, 0.55
n, p = 15, 0.55
viz = False

ks = np.linspace(0, n + 1, 20)

if viz:
  plt.figure()
  [plt.plot(np.linspace(0, n + 1), pmf(np.linspace(0, n + 1), n, p)) for n in range(2, 10)]

tops = np.linspace(0, n + 1, 20)
denominator = mp.quad(lambda k: mppmf(k, n, p), [0.0, n + 1])
cdfs = np.array([mp.quad(lambda k: mppmf(k, n, p), [0.0, top]) for top in tops]) / denominator
pmfs = [mppmf(top, n, p) / denominator for top in tops]
([f'k/n={x}, pmf={w}, cdf={y}, 1-cdf={z}' for x, y, z, w in zip(*[tops / n, cdfs, 1 - cdfs, pmfs])])

lookup_n_p = dict()
for n in tqdm.tqdm(range(1, 41)):
  tops = np.linspace(0, n + 1, 20)
  for p in [0.55, 0.65, 0.75, 0.85, 0.95]:
    if (n, p) not in lookup_n_p:
      denominator = mp.quad(lambda k: mppmf(k, n, p), [0.0, n + 1])
      cdfs = np.array([mp.quad(lambda k: mppmf(k, n, p), [0.0, top]) for top in tops]) / denominator
      lookup_n_p[(n, p)] = (interp1d(cdfs, tops), denominator, tops, cdfs)

f = lambda top: mp.quad(lambda k: mppmf(k, n, p), [0.0, top])
mp.findroot(lambda x: f(x) - 0.95, (n + 1) * 0.95)


def ppf(q, n, p):
  interpolator, den, *_ = lookup_n_p[(n, p)]
  f = lambda top: mp.quad(lambda k: mppmf(k, n, p) / den, [0.0, top])
  return mp.findroot(lambda x: f(x) - q, float(interpolator(q)), tol=1e-4)


lookup_npq = dict()
for n in tqdm.tqdm(range(1, 40 + 1)):
  for p in [0.55, 0.65, 0.75, 0.85, 0.95]:
    for q in [0.05, 0.95, 0.25, 0.75, 0.45, 0.55]:
      if (n, p, q) not in lookup_npq:
        lookup_npq[(n, p, q)] = float(ppf(q, n, p))

with open('npq.pydat', 'w') as fid:
  fid.write(str(lookup_npq))

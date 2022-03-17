# semi-related but different distribution? https://rdrr.io/cran/cbinom/man/cbinom.html

# You'll need scipy and numpy to run this
from fractions import Fraction
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

if viz:
  ks = np.linspace(0, n + 1, 20)
  plt.figure()
  [plt.plot(np.linspace(0, n + 1), pmf(np.linspace(0, n + 1), n, p)) for n in range(2, 10)]


def ppf(q, n, p, den, init, tol=1e-4):
  f = lambda top: mp.quad(lambda k: mppmf(k, n, p), [0.0, top])
  return mp.findroot(lambda x: f(x) - q * den, float(init), tol=tol)


def maybeFloatInt(x: float | int | Fraction) -> int | float:
  "Returns input if it's truly float, else an int"
  return int(x) if int(x) == x else float(x)


def confIntervalsToQuantiles(qs: list[int]) -> list[float | int]:
  "Given confidence intervals, e.g., [0.99, 0.9], generate min/max quantiles"
  ret: list[int | float] = []
  for q in qs:
    assert 0 < q < 100
    ret.append(maybeFloatInt(Fraction(100 + q, 2)))
    ret.append(maybeFloatInt(Fraction(100 - q, 2)))
  return ret


lookup_npq = dict()
# try to rehydrate this object from disk
try:

  def object_hook(d):
    return {maybeFloatInt(float(key)): d[key] for key in d}

  with open('npq.json', 'r') as fid:
    lookup_npq = json.load(fid, object_hook=object_hook)
except FileNotFoundError:
  pass

for n in tqdm.tqdm(range(1, 41)):
  # this is the number of questions in a given confidence
  if n not in lookup_npq:
    lookup_npq[n] = dict()

  kgrid = np.linspace(0, n + 1, 20)
  for pPct in [55, 65, 75, 85, 95]:
    p = pPct / 100
    # this is the actual confidence of all `n` questions
    # (Per Galef's book, these are the only inputs acceptable)
    if p not in lookup_npq[n]:
      lookup_npq[n][pPct] = dict()

    denominator = mp.quad(lambda k: mppmf(k, n, p), [0.0, n + 1])
    cdfs = [mp.quad(lambda k: mppmf(k, n, p), [0.0, k]) / denominator for k in kgrid]
    interp = interp1d(cdfs, kgrid)

    for qPct in confIntervalsToQuantiles([90, 50, 10]):
      q = qPct / 100
      # these are quantiles
      lookup_npq[n][pPct][qPct] = float(ppf(q, n, p, denominator, float(interp(q))))

with open('npq.json', 'w') as fid:
  json.dump(lookup_npq, fid)

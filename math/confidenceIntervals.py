# semi-related but different distribution? https://rdrr.io/cran/cbinom/man/cbinom.html

# You'll need scipy and numpy to run this
from fractions import Fraction
from scipy.interpolate import interp1d
import mpmath as mp
import numpy as np
import tqdm
import json
from functools import cache

JSON_FILE = '../src/npq.json'


def binomcoef(x, y):
  """The binomial coefficient (x choose y) for real numbers"""
  # https://en.wikipedia.org/w/index.php?title=Binomial_coefficient&oldid=1077044362#Two_real_or_complex_valued_arguments
  return 1 / ((x + 1) * mp.beta(x - y + 1, y + 1))


def pmfUnnormalized(k, n, p):
  """
  The pseudo-probability mass function (PMF) for the continuous binomial

  "Pseudo" because it needs to be normalized (see `pmfNormalizingDenominator`).

  This is just the ordinary binomial's PMF except extended to the reals and
  `0 <= k <= n + 1` is valid. (Don't worry, the PMF at `n+1` is 0!) 
  
  - `k` = number of successes
  - `n` = total number of trials
  - `p` = probability of each success
  """
  return p**k * (1 - p)**(n - k) * binomcoef(n, k)


@cache
def pmfNormalizingDenominator(n, p):
  ""
  return mp.quad(lambda k: pmfUnnormalized(k, n, p), [0.0, n + 1])


def ppf(q, n, p, init, tol=1e-4):
  """
  Almost the percent point function (inverse of the CDF)
  
  Following scipy.stats, this is the inverse of the cumulative 
  distribution function (CDF), except it needs an initial guess, `init`.

  Returns `k1` such that, for given `n` number of trials and `p` underlying 
  binomial probability:

  `Probability(k <= k1; n) = q`

  to within `tol` tolerance. Since I cannot be bothered to actually do math,
  we use numerical integration and root-finding via mpmath.
  """
  unnormCdf = lambda top: mp.quad(lambda k: pmfUnnormalized(k, n, p), [0.0, top])
  den = pmfNormalizingDenominator(n, p)
  return mp.findroot(lambda x: unnormCdf(x) - q * den, float(init), tol=tol)


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

  with open(JSON_FILE, 'r') as fid:
    lookup_npq = json.load(fid, object_hook=object_hook)
except FileNotFoundError:
  pass

# Run!
qs = confIntervalsToQuantiles([90, 50, 10])
for n in tqdm.tqdm(range(1, 41)):
  # this is the number of questions in a given confidence
  if n not in lookup_npq:
    lookup_npq[n] = dict()

  kgrid = np.linspace(0, n + 1, 20)
  for pPct in [55, 65, 75, 85, 95]:
    p = pPct / 100
    # this is the actual confidence of all `n` questions
    # (Per Galef's book, these are the only inputs acceptable)
    if pPct not in lookup_npq[n]:
      lookup_npq[n][pPct] = dict()
    if all(q in lookup_npq[n][pPct] for q in qs):
      continue

    denominator = pmfNormalizingDenominator(n, p)
    cdfs = [mp.quad(lambda k: pmfUnnormalized(k, n, p), [0.0, k]) / denominator for k in kgrid]
    interp = interp1d(cdfs, kgrid)

    for qPct in qs:
      q = qPct / 100
      # these are quantiles
      if qPct not in lookup_npq[n][pPct]:
        lookup_npq[n][pPct][qPct] = float(ppf(q, n, p, float(interp(q))))

# persist to disk
with open(JSON_FILE, 'w') as fid:
  json.dump(lookup_npq, fid)

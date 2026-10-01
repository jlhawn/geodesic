"""The sweep's response surface:  python3 scripts/sweep/fit.py [runs/sweep]

Fits a full quadratic (linear, square and pairwise interaction terms) in
the parameters coded to [-1, 1] over their ranges to the score of every
row of results.csv, by ridge regression with the penalty chosen by
leave-one-out error (78 coefficients from 41 runs need one). Writes
fit.txt (the coefficients, the leave-one-out error, each parameter's
sensitivities, each score term's linear sensitivity to each parameter)
and candidates.json: the four best design points, the surface's minimum
over the box (projected gradient from many starts; 'qmin'), with the
stationary point clipped to the box reported beside it, and the minimum
of the score composed from each term's own fit ('cmin'), which cannot
fall below zero where the quadratic of the score can, and the defaults
('base') for reference.
"""
import csv, json, sys
import numpy as np

SWEEP = sys.argv[1] if len(sys.argv) > 1 else 'runs/sweep'
PARAMETERS = [
    ('varianceScale', 2, 10), ('mixingLength', 150, 600), ('stratiformHours', 1, 6), ('cloudHours', 0.5, 2),
    ('plumeEntrainment', 0.05, 0.2), ('plumeCape', 40, 200), ('minimumInversion', 2, 6), ('criticalHumidity', 0.7, 0.9),
    ('seaDrag', 1.0e-3, 1.5e-3), ('stableMixingLength', 10, 60), ('cumulusCeiling', 1500, 2500),
]
KEYS = [p[0] for p in PARAMETERS]
WEIGHTS = {'balance': 4, 'albedo': 2, 'rain': 1, 'sepLow': 1, 'peruLow': 1, 'sepLwp': 1, 'peruLwp': 1, 'sepRain': 1, 'peruRain': 1, 'itczRain': 1, 'itczPeak': 0.5, 'zonalPeak': 0.5, 'stress': 1, 'arctic': 2}
LOW = np.array([p[1] for p in PARAMETERS]); HIGH = np.array([p[2] for p in PARAMETERS])
D = len(KEYS)

rows = list(csv.DictReader(open(f'{SWEEP}/results.csv')))
rows.sort(key=lambda r: int(r['point']))
X = np.array([[float(r[k]) for k in KEYS] for r in rows])
y = np.array([float(r['score']) for r in rows])
coded = lambda x: 2 * (x - LOW) / (HIGH - LOW) - 1
uncoded = lambda c: LOW + (c + 1) * (HIGH - LOW) / 2
PAIRS = [(j, k) for j in range(D) for k in range(j + 1, D)]
NAMES = ['1'] + KEYS + [f'{k}^2' for k in KEYS] + [f'{KEYS[j]}*{KEYS[k]}' for j, k in PAIRS]


def features(c):
    c = np.atleast_2d(c)
    return np.hstack([np.ones((len(c), 1)), c, c ** 2, np.stack([c[:, j] * c[:, k] for j, k in PAIRS], axis=1)])


def ridge(F, y, lam):
    P = np.eye(F.shape[1]) * lam
    P[0, 0] = 0
    A = F.T @ F + P
    beta = np.linalg.solve(A, F.T @ y)
    H = F @ np.linalg.solve(A, F.T)
    loo = (y - F @ beta) / (1 - np.diag(H))
    return beta, np.sqrt(np.mean(loo ** 2))


C = coded(X)
F = features(C)
grid = np.logspace(-4, 3, 71)
fits = [(lam, *ridge(F, y, lam)) for lam in grid]
lam, beta, loo = min(fits, key=lambda f: f[2])
fitted = F @ beta
r2 = 1 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2)
q2 = 1 - loo ** 2 * len(y) / np.sum((y - y.mean()) ** 2)
predict = lambda c: features(c) @ beta


def gradient(c):
    g = beta[1:1 + D] + 2 * beta[1 + D:1 + 2 * D] * c
    for n, (j, k) in enumerate(PAIRS):
        b = beta[1 + 2 * D + n]
        g[j] += b * c[k]
        g[k] += b * c[j]
    return g


Q = np.diag(beta[1 + D:1 + 2 * D]).astype(float)
for n, (j, k) in enumerate(PAIRS):
    Q[j, k] = Q[k, j] = beta[1 + 2 * D + n] / 2
eig = np.linalg.eigvalsh(Q)
try:
    stationary = np.linalg.solve(2 * Q, -beta[1:1 + D])
except np.linalg.LinAlgError:
    stationary = np.zeros(D)
clipped = np.clip(stationary, -1, 1)

rng = np.random.default_rng(1)
starts = np.vstack([C, rng.uniform(-1, 1, (400, D)), clipped[None, :]])
best, bestValue = None, np.inf
step = 1 / (2 * max(abs(eig).max(), 1e-9))
for c in starts:
    c = c.copy()
    for _ in range(3000):
        n = np.clip(c - step * gradient(c), -1, 1)
        if np.max(abs(n - c)) < 1e-10:
            break
        c = n
    v = predict(c)[0]
    if v < bestValue:
        best, bestValue = c, v

base = coded(np.array([float(rows[0][k]) for k in KEYS])) if rows[0]['point'] == '0' else np.zeros(D)
out = []
say = out.append
say(f'quadratic response surface of the score over {len(y)} runs (coded to [-1, 1] over each range); ridge penalty {lam:.3g} by leave-one-out')
say(f'fit R2 {r2:.3f}, leave-one-out RMSE {loo:.2f} (score spread {y.std():.2f}, mean {y.mean():.2f}, Q2 {q2:.3f}); quadratic part eigenvalues {", ".join(f"{e:.2f}" for e in eig)}')
say('coefficients (coded units): ' + ', '.join(f'{n} {b:+.2f}' for n, b in zip(NAMES, beta)))
say('')
say('sensitivity of the score at the defaults (point 0): slope over the half range, curvature, and the surface along each range with the others at the defaults (low / default / high)')
for k, key in enumerate(KEYS):
    line = []
    for v in (-1, base[k], 1):
        c = base.copy(); c[k] = v
        line.append(predict(c)[0])
    say(f'  {key:18s} slope {gradient(base)[k]:+7.2f}  curvature {2 * beta[1 + D + k]:+7.2f}  surface {line[0]:7.2f} / {line[1]:7.2f} / {line[2]:7.2f}')
say('largest interactions: ' + ', '.join(f'{NAMES[1 + 2 * D + n]} {beta[1 + 2 * D + n]:+.2f}' for n in np.argsort(-abs(beta[1 + 2 * D:]))[:8]))
say('')
terms = [k[2:] for k in rows[0].keys() if k.startswith('e_')]
say('linear sensitivity of each normalized error to each parameter over its full range (least squares, coded units x2), and R2:')
say('  ' + ' ' * 10 + ''.join(f'{k[:9]:>10s}' for k in KEYS) + '        R2')
L = np.hstack([np.ones((len(C), 1)), C])
for t in terms:
    e = np.array([float(r[f'e_{t}']) for r in rows])
    b, *_ = np.linalg.lstsq(L, e, rcond=None)
    spread = np.sum((e - e.mean()) ** 2)
    rr = 1 - np.sum((e - L @ b) ** 2) / spread if spread > 0 else float('nan')
    say(f'  {t:10s}' + ''.join(f'{2 * x:+10.2f}' for x in b[1:]) + f'{rr:10.2f}')
say('')
weights = WEIGHTS
G = np.hstack([np.ones((len(C), 1)), C, C ** 2])
termFits, looTerms = {}, np.zeros((len(C), len(terms)))
for n, t in enumerate(terms):
    e = np.array([float(r[f'e_{t}']) for r in rows])
    lamT, b, _ = min(((l, *ridge(G, e, l)) for l in grid), key=lambda f: f[2])
    P = np.eye(G.shape[1]) * lamT; P[0, 0] = 0
    H = G @ np.linalg.solve(G.T @ G + P, G.T)
    looTerms[:, n] = e - (e - G @ b) / (1 - np.diag(H))
    termFits[t] = b
W = np.array([weights[t] for t in terms])
B = np.array([termFits[t] for t in terms])


def composedAll(c):
    e = B[:, 0][None, :] + c @ B[:, 1:1 + D].T + (c ** 2) @ B[:, 1 + D:].T
    value = (e ** 2) @ W
    grad = 2 * ((W * e) @ B[:, 1:1 + D]) + 4 * c * ((W * e) @ B[:, 1 + D:])
    return value, grad


composed = lambda c: composedAll(np.atleast_2d(c))[0][0]
looComposed = looTerms ** 2 @ W
rmseComposed = np.sqrt(np.mean((looComposed - y) ** 2))
c = starts.copy()
for _ in range(20000):
    _, g = composedAll(c)
    c = np.clip(c - 0.001 * g, -1, 1)
values, _ = composedAll(c)
cbest, cvalue = c[np.argmin(values)], values.min()
say(f'cross-check: the score composed from each normalized error fitted as linear plus square terms (ridge by leave-one-out), Σ w (ê)²: leave-one-out RMSE {rmseComposed:.2f}')
say('its minimum over the box: ' + ', '.join(f'{k} {v:.4g}' for k, v in zip(KEYS, uncoded(cbest))) + f'; composed {cvalue:.2f}, quadratic surface there {predict(cbest)[0]:.2f}')
say('')
order = np.argsort(y)
say('runs by score: ' + ', '.join(f'{rows[i]["point"]} {y[i]:.1f}' for i in order))
say('surface minimum over the box: ' + ', '.join(f'{k} {v:.4g}' for k, v in zip(KEYS, uncoded(best))) + f'; predicted {bestValue:.2f}')
say('stationary point clipped to the box: ' + ', '.join(f'{k} {v:.4g}' for k, v in zip(KEYS, uncoded(clipped))) + f'; predicted {predict(clipped)[0]:.2f}')
open(f'{SWEEP}/fit.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))

candidates = [{'name': f'p{int(rows[i]["point"]):02d}', 'point': int(rows[i]['point']), 'screenScore': y[i], 'predicted': predict(C[i])[0], **{k: float(rows[i][k]) for k in KEYS}} for i in order[:4]]
minimum = {k: float(f'{v:.4g}') for k, v in zip(KEYS, uncoded(best))}
candidates.append({'name': 'qmin', 'point': None, 'screenScore': None, 'predicted': float(bestValue), **minimum})
candidates.append({'name': 'cmin', 'point': None, 'screenScore': None, 'predicted': float(predict(cbest)[0]), 'composed': float(cvalue), **{k: float(f'{v:.4g}') for k, v in zip(KEYS, uncoded(cbest))}})
zero = next(i for i, r in enumerate(rows) if r['point'] == '0')
candidates.append({'name': 'base', 'point': 0, 'screenScore': y[zero], 'predicted': predict(C[zero])[0], **{k: float(rows[zero][k]) for k in KEYS}})
json.dump({'surfaceMinimum': minimum, 'composedMinimum': {k: float(f'{v:.4g}') for k, v in zip(KEYS, uncoded(cbest))}, 'composedValue': float(cvalue), 'composedLoo': float(rmseComposed), 'lambda': float(lam), 'loo': float(loo), 'r2': float(r2)}, open(f'{SWEEP}/fit.json', 'w'), indent=1)
json.dump(candidates, open(f'{SWEEP}/candidates.json', 'w'), indent=1)

# Quarter, year and rolling-year means of a spin-up log's daily lines (365-day
# model year, day 0 at the March equinox, quarter ends at round(k*365/4)), with
# the regional rain of each segment's 'regions after' line against Earth's:
#   python3 scripts/rounds.py runs/verda-eleven/eleven128.log [more logs]
# A bare tag names runs/<tag>.log. Environment: PER_YEAR (4 quarters), SEGMENTS
# (1; 0 leaves out each segment's own regional line).
import re, sys, os
from datetime import date, timedelta
RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'runs')
YEAR, PER_YEAR, SEGMENTS = 365, int(os.environ.get('PER_YEAR', 4)), os.environ.get('SEGMENTS', '1') != '0'
cal = lambda d: (date(2001, 3, 20) + timedelta(days=d)).strftime('%b %-d')
qend = lambda k: int(k * YEAR / PER_YEAR + 0.5)
EARTH = {'ts': 14.85, 'albedo': 0.30, 'rain': 2.7}
EARTH_RAIN = {'sahara': 0.2, 'arabia': 0.2, 'sahel': 1.5, 'india': 3.0, 'congo': 4.5, 'amazon': 6.0, 'seAsia': 5.0, 'borneo': 8.0, 'europe': 2.0, 'eastUS': 3.0, 'siberia': 1.2, 'ausInterior': 0.7, 'kalahari': 1.0, 'gobi': 0.5, 'usSouthwest': 0.8, 'cerrado': 4.0, 'ausNorth': 2.7, 'ausEast': 2.5, 'ausSoutheast': 1.8, 'ausWest': 0.7, 'newGuinea': 8.0}
DAY = re.compile(r'day (\d+) .*Ts (-?[\d.]+) °C, ASR ([\d.]+)(?: \(atmosphere [\d.]+\))? OLR ([\d.]+).*precip ([\d.]+) mm/d, ice ([\d.]+)% \(N ([\d.]+) S ([\d.]+) Mkm²\), albedo ([\d.]+)')

def load(tag):
    rows, segments, last = [], [], 0
    for line in open(tag if tag.endswith('.log') else f'{RUNS}/{tag}.log'):
        m = DAY.match(line)
        if m: rows.append([float(x) for x in m.groups()]); last = int(rows[-1][0])
        elif line.startswith('regions after'):
            days = int(re.match(r'regions after (\d+) days', line).group(1))
            segments.append({'end': last, 'days': days, 'regions': {k: tuple(float(x) for x in v.split('/')) for k, v in re.findall(r'(\w+) ([\d.]+/[\d.]+/-?[\d.]+)', line)}})
    return rows, segments

def summary(rows, a, b):
    win = [r for r in rows if a < r[0] <= b]
    if len(win) < (b - a) - 2: return None
    n = len(win); mean = lambda k: sum(r[k] for r in win) / n
    return f"Ts {mean(1):.2f} °C, TOA {mean(2)-mean(3):+.1f} W/m², albedo {mean(8):.3f}, rain {mean(4):.2f} mm/d, ice N {min(r[6] for r in win):.1f}–{max(r[6] for r in win):.1f} S {min(r[7] for r in win):.1f}–{max(r[7] for r in win):.1f} Mkm² (mean N {mean(6):.1f} S {mean(7):.1f})"

def regional(segments, a, b):
    inside = [s for s in segments if a < s['end'] <= b]
    if not inside: return None
    names = inside[-1]['regions'].keys(); total = sum(s['days'] for s in inside)
    return ', '.join(f"{k} {sum(s['regions'][k][0] * s['days'] for s in inside if k in s['regions']) / total:.1f}/{EARTH_RAIN.get(k, float('nan')):.1f} v{inside[-1]['regions'][k][1]:.2f}" for k in names)

for arg in sys.argv[1:]:
    rows, segments = load(arg)
    tag = os.path.basename(arg)[:-4] if arg.endswith('.log') else arg
    last = int(max(r[0] for r in rows)) if rows else 0
    k = 1
    while qend(k - 1) < last:
        a, b = qend(k - 1), min(qend(k), last)
        s = summary(rows, a, b)
        if s: print(f"{tag} Q{k} days {a+1}-{b} ({cal(a+1)}–{cal(b)}): {s}")
        for seg in segments if SEGMENTS else []:
            if a < seg['end'] <= b: print(f"  regions to day {seg['end']} ({seg['days']} d, rain mm/d / vegetation / surface °C): " + ', '.join(f"{n} {v[0]:.1f}/{v[1]:.2f}/{v[2]:.0f}" for n, v in seg['regions'].items()))
        k += 1
    for y in range(last // YEAR):
        print(f"{tag} YEAR {y+1}: {summary(rows, YEAR * y, YEAR * (y + 1))}  | Earth: Ts {EARTH['ts']}, albedo {EARTH['albedo']}, rain {EARTH['rain']}, ice N 15/6, S 18/3")
    for k in range(PER_YEAR, last // (YEAR // PER_YEAR) + 2):
        b = qend(k)
        if b > last or b % YEAR == 0: continue
        s = summary(rows, b - YEAR, b)
        if not s: continue
        print(f"{tag} ROLLING YEAR to day {b} ({cal(b)}): {s}")
        r = regional(segments, b - YEAR, b)
        if r: print(f"  rolling-year regional rain / Earth (mm/d), latest vegetation: {r}")
    if last >= YEAR:
        r = regional(segments, (last // YEAR - 1) * YEAR, (last // YEAR) * YEAR)
        if r: print(f"{tag} YEAR {last // YEAR} regional rain / Earth (mm/d), latest vegetation: {r}")

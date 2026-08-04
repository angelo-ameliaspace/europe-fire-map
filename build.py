#!/usr/bin/env python3
"""
Rebuild the "Where Europe is burning right now" page from live NASA FIRMS data.

Standard library only — no numpy, geopandas or scikit-learn. Runs in any Python 3.9+.

  python3 build.py [--out europe_fire_activity.html]

Reads  : boundaries.json, template.html  (alongside this script)
Writes : europe_fire_activity.html
Exits non-zero on any fetch or sanity-check failure, so a scheduled run fails loudly
rather than publishing a broken or empty page.
"""

import argparse, csv, io, json, math, os, sys, time, urllib.error, urllib.request
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
FIRMS = "https://firms.modaps.eosdis.nasa.gov/data/active_fire"
SATS = [
    ("suomi-npp-viirs-c2", "SUOMI_VIIRS_C2", "Suomi-NPP"),
    ("noaa-20-viirs-c2",   "J1_VIIRS_C2",    "NOAA-20"),
    ("noaa-21-viirs-c2",   "J2_VIIRS_C2",    "NOAA-21"),
]

# map frame (drawn) and projection constants
BB = (-25.0, 33.0, 40.0, 64.0)
LON0, PH1, PH2, PH0 = 13.0, 40.0, 60.0, 50.0
W, H = 1160.0, 760.0
GRID_DEG = 0.03          # density-layer cell size
LINK_KM = 3.0            # single-linkage distance for complexes
MIN_CLUSTER = 5          # detections required to be called a complex
FRP_BREAKS = [5, 20, 60, 150]
R_EARTH_KM = 6371.0

# minimum plausible detection count; guards against a truncated or empty feed
SANITY_MIN_7D = 500


def log(*a):
    print("[build]", *a, flush=True)


# ---------------------------------------------------------------- fetch
def fetch(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "europe-fire-map/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read().decode("utf-8", "replace")
            if not data.lstrip().lower().startswith("latitude"):
                raise ValueError("unexpected payload (not a FIRMS CSV): " + data[:120])
            return data
        except Exception as e:  # noqa: BLE001
            last = e
            log(f"  attempt {i+1}/{tries} failed: {e}")
            time.sleep(3 * (i + 1))
    raise SystemExit(f"FATAL: could not fetch {url}: {last}")


def load_window(win):
    rows = []
    for slug, stem, plat in SATS:
        url = f"{FIRMS}/{slug}/csv/{stem}_Europe_{win}.csv"
        log(f"fetching {stem} {win}")
        rd = csv.DictReader(io.StringIO(fetch(url)))
        n = 0
        for r in rd:
            if r.get("confidence") == "low":
                continue
            try:
                lat = float(r["latitude"]); lon = float(r["longitude"])
                frp = float(r["frp"] or 0.0)
                t = datetime.strptime(r["acq_date"] + r["acq_time"].zfill(4),
                                      "%Y-%m-%d%H%M").replace(tzinfo=timezone.utc)
            except (KeyError, ValueError):
                continue
            rows.append({"lat": lat, "lon": lon, "frp": frp, "t": t, "plat": plat})
            n += 1
        log(f"  {n:,} usable rows")
    return rows


# ---------------------------------------------------------------- projection
_n = math.log(math.cos(math.radians(PH1)) / math.cos(math.radians(PH2))) / math.log(
    math.tan(math.pi / 4 + math.radians(PH2) / 2) / math.tan(math.pi / 4 + math.radians(PH1) / 2))
_F = math.cos(math.radians(PH1)) * math.tan(math.pi / 4 + math.radians(PH1) / 2) ** _n / _n
_rho0 = _F / math.tan(math.pi / 4 + math.radians(PH0) / 2) ** _n


def lcc(lon, lat):
    """Lambert conformal conic, unit sphere."""
    lat = max(min(lat, 89.0), -89.0)
    rho = _F / math.tan(math.pi / 4 + math.radians(lat) / 2) ** _n
    th = _n * math.radians(lon - LON0)
    return rho * math.sin(th), _rho0 - rho * math.cos(th)


# fit the frame
_xs, _ys = [], []
for i in range(41):
    for j in range(41):
        x, y = lcc(BB[0] + (BB[2] - BB[0]) * i / 40, BB[1] + (BB[3] - BB[1]) * j / 40)
        _xs.append(x); _ys.append(y)
X0, X1, Y0, Y1 = min(_xs), max(_xs), min(_ys), max(_ys)
S = min(W / (X1 - X0), H / (Y1 - Y0))
OX, OY = (W - (X1 - X0) * S) / 2, (H - (Y1 - Y0) * S) / 2


def to_px(lon, lat):
    x, y = lcc(lon, lat)
    return (x - X0) * S + OX, (Y1 - y) * S + OY      # flip y for SVG


def to_km(lon, lat):
    x, y = lcc(lon, lat)
    return x * R_EARTH_KM, y * R_EARTH_KM


# ---------------------------------------------------------------- countries
def load_boundaries():
    with open(os.path.join(HERE, "boundaries.json")) as f:
        return json.load(f)


def in_ring(x, y, ring):
    inside = False
    n = len(ring)
    j = n - 1
    for i in range(n):
        xi, yi = ring[i]; xj, yj = ring[j]
        if (yi > y) != (yj > y):
            if x < (xj - xi) * (y - yi) / (yj - yi) + xi:
                inside = not inside
        j = i
    return inside


def attribute(rows, bnd):
    names, bbs, polys = bnd["names"], bnd["bbox"], bnd["polys"]
    order = sorted(range(len(names)), key=lambda i: -sum(len(r) for r in polys[i]))
    hit = 0
    for r in rows:
        x, y = r["lon"], r["lat"]
        r["country"] = None
        for i in order:
            b = bbs[i]
            if not (b[0] <= x <= b[2] and b[1] <= y <= b[3]):
                continue
            for ring in polys[i]:
                if in_ring(x, y, ring):
                    r["country"] = names[i]
                    break
            if r["country"]:
                break
        if r["country"]:
            hit += 1
    log(f"country attribution: {hit:,}/{len(rows):,} ({100*hit/max(len(rows),1):.1f}%)")


# ---------------------------------------------------------------- clustering
def cluster(rows):
    """Single-linkage at LINK_KM via grid hashing + union-find."""
    parent = list(range(len(rows)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    cells = {}
    pts = []
    for i, r in enumerate(rows):
        kx, ky = to_km(r["lon"], r["lat"])
        pts.append((kx, ky))
        cells.setdefault((int(kx // LINK_KM), int(ky // LINK_KM)), []).append(i)

    lim = LINK_KM * LINK_KM
    for (cx, cy), idx in cells.items():
        near = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                near.extend(cells.get((cx + dx, cy + dy), ()))
        for a in idx:
            ax, ay = pts[a]
            for b in near:
                if b <= a:
                    continue
                bx, by = pts[b]
                if (ax - bx) ** 2 + (ay - by) ** 2 <= lim:
                    union(a, b)

    groups = {}
    for i in range(len(rows)):
        groups.setdefault(find(i), []).append(i)
    return [g for g in groups.values() if len(g) >= MIN_CLUSTER]


def hull_area_ha(pts):
    """Convex hull footprint, geodesic-ish via the km plane. Not a burned area."""
    p = sorted(set(pts))
    if len(p) < 3:
        return 0.0

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lo = []
    for q in p:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], q) <= 0:
            lo.pop()
        lo.append(q)
    up = []
    for q in reversed(p):
        while len(up) >= 2 and cross(up[-2], up[-1], q) <= 0:
            up.pop()
        up.append(q)
    h = lo[:-1] + up[:-1]
    if len(h) < 3:
        return 0.0
    a = 0.0
    for i in range(len(h)):
        x1, y1 = h[i]; x2, y2 = h[(i + 1) % len(h)]
        a += x1 * y2 - x2 * y1
    return abs(a) / 2 * 100.0        # km² -> ha


def bucket(v):
    b = 0
    for t in FRP_BREAKS:
        if v >= t:
            b += 1
    return b


# ---------------------------------------------------------------- build
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "europe_fire_activity.html"))
    args = ap.parse_args()

    r7 = load_window("7d")
    r24 = load_window("24h")
    if len(r7) < SANITY_MIN_7D:
        raise SystemExit(f"FATAL: only {len(r7)} 7-day detections — feed looks broken, refusing to publish")

    bnd = load_boundaries()
    attribute(r7, bnd)
    attribute(r24, bnd)

    cut = min(r["t"] for r in r24)
    inframe = lambda r: BB[0] <= r["lon"] <= BB[2] and BB[1] <= r["lat"] <= BB[3]
    f7 = [r for r in r7 if inframe(r)]
    log(f"{len(r7)-len(f7):,} of {len(r7):,} detections outside map frame "
        f"({100*(len(r7)-len(f7))/len(r7):.1f}%)")

    # density grid
    cellmax = {}
    for r in f7:
        k = (round(r["lon"] / GRID_DEG), round(r["lat"] / GRID_DEG))
        if r["frp"] > cellmax.get(k, -1):
            cellmax[k] = r["frp"]
    dots = []
    for (gx, gy), frp in cellmax.items():
        x, y = to_px(gx * GRID_DEG, gy * GRID_DEG)
        dots.append([round(x, 1), round(y, 1), bucket(frp)])

    # complexes
    comp = []
    for g in cluster(f7):
        rs = [f7[i] for i in g]
        cs = {}
        for r in rs:
            if r["country"]:
                cs[r["country"]] = cs.get(r["country"], 0) + 1
        name = max(cs, key=cs.get) if cs else "offshore / unattributed"
        lat = sum(r["lat"] for r in rs) / len(rs)
        lon = sum(r["lon"] for r in rs) / len(rs)
        x, y = to_px(lon, lat)
        n24 = sum(1 for r in rs if r["t"] >= cut)
        comp.append({
            "x": round(x, 1), "y": round(y, 1), "n": len(rs),
            "fs": round(sum(r["frp"] for r in rs)), "fm": round(max(r["frp"] for r in rs), 1),
            "c": name, "act": n24 > 0, "n24": n24,
            "fp": round(hull_area_ha([to_km(r["lon"], r["lat"]) for r in rs])),
            "b": bucket(max(r["frp"] for r in rs)),
            "lat": round(lat, 3), "lon": round(lon, 3),
            "t0": min(r["t"] for r in rs).strftime("%d %b %H:%M"),
            "t1": max(r["t"] for r in rs).strftime("%d %b %H:%M"),
        })
    comp.sort(key=lambda d: -d["fs"])
    log(f"{len(comp):,} complexes, {sum(1 for c in comp if c['act']):,} active in last 24 h")

    # country table
    agg = {}
    for r in r7:
        if not r["country"]:
            continue
        a = agg.setdefault(r["country"], {"det": 0, "frp": 0.0, "det24": 0})
        a["det"] += 1; a["frp"] += r["frp"]
    for r in r24:
        if r["country"] in agg:
            agg[r["country"]]["det24"] += 1
    countries = sorted(
        [{"name": k, "det": v["det"], "frp": round(v["frp"]),
          "det24": v["det24"], "fmean": round(v["frp"] / v["det"], 1)}
         for k, v in agg.items()], key=lambda d: -d["frp"])[:14]

    # country outlines + neatline
    paths = []
    for rings in bnd["polys"]:
        d = ""
        for ring in rings:
            pts = [to_px(a, b) for a, b in ring]
            d += "M" + "L".join(f"{a:.1f} {b:.1f}" for a, b in pts) + "Z"
        if d:
            paths.append(d)
    edge = ([(BB[0] + (BB[2] - BB[0]) * i / 119, BB[1]) for i in range(120)]
            + [(BB[2], BB[1] + (BB[3] - BB[1]) * i / 119) for i in range(120)]
            + [(BB[0] + (BB[2] - BB[0]) * (119 - i) / 119, BB[3]) for i in range(120)]
            + [(BB[0], BB[1] + (BB[3] - BB[1]) * (119 - i) / 119) for i in range(120)])
    fp = [to_px(a, b) for a, b in edge]
    frame = "M" + "L".join(f"{a:.1f} {b:.1f}" for a, b in fp) + "Z"

    now = datetime.now(timezone.utc)
    payload = {
        "meta": {
            "win7": [min(r["t"] for r in r7).strftime("%d %b %Y %H:%M"),
                     max(r["t"] for r in r7).strftime("%d %b %Y %H:%M")],
            "win24": [min(r["t"] for r in r24).strftime("%d %b %H:%M"),
                      max(r["t"] for r in r24).strftime("%d %b %H:%M")],
            "n7": len(r7), "n24": len(r24), "ncomp": len(comp),
            "nact": sum(1 for c in comp if c["act"]), "ncell": len(dots),
            "outpct": round(100 * (len(r7) - len(f7)) / len(r7), 1),
            "retrieved": now.strftime("%-d %B %Y, %H:%M UTC") if os.name != "nt"
                         else now.strftime("%d %B %Y, %H:%M UTC"),
            "W": W, "H": H, "breaks": FRP_BREAKS,
        },
        "paths": paths, "frame": frame, "dots": dots, "comp": comp, "countries": countries,
    }

    with open(os.path.join(HERE, "template.html")) as f:
        tpl = f.read()
    if "/*PAYLOAD*/" not in tpl:
        raise SystemExit("FATAL: template.html has no /*PAYLOAD*/ marker")
    html = tpl.replace("/*PAYLOAD*/", "const FIRE=" + json.dumps(payload, separators=(",", ":")) + ";")
    with open(args.out, "w") as f:
        f.write(html)

    top = next(c for c in comp if c["act"])
    live = sum(c["fs"] for c in comp if c["act"])
    log(f"wrote {args.out} ({os.path.getsize(args.out)/1024:.0f} KB)")
    log(f"headline: {top['c']} {top['lat']}/{top['lon']} {top['fs']:,} MW "
        f"= {round(100*top['fs']/max(live,1))}% of active radiative power")
    log(f"as of {payload['meta']['retrieved']}")


if __name__ == "__main__":
    main()

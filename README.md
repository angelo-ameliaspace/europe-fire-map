# europe-fire-map

Rebuilds the **Where Europe is burning right now** page from live NASA FIRMS VIIRS
active-fire detections, and republishes it to a fixed Claude artifact URL.

**Artifact (do not change):**
`https://claude.ai/code/artifact/200730f4-270f-4758-b887-70fcef594c1e`

Refreshed hourly by a Claude Code cloud routine. Publishing to that exact URL is
what keeps the link stable for everyone it has been shared with — publishing without it
mints a new artifact and viewers silently keep seeing stale data on the old link.

## Run it

```bash
python3 build.py                     # -> europe_fire_activity.html
python3 build.py --out /tmp/page.html
```

No dependencies. Standard library only, Python 3.9+. Takes about 30 seconds, almost all
of it downloading.

## Files

| File | Purpose |
|---|---|
| `build.py` | Fetch, cluster, project, render. The whole pipeline. |
| `template.html` | Page shell — CSS, markup, client JS. Contains a `/*PAYLOAD*/` marker that `build.py` replaces with the data. |
| `brand.css` | Amelia Space design tokens: Inter (variable, SIL OFL 1.1) as taken from ameliaspace.com, the brand primary `#008cff`, and the Amelia wordmark in both polarities. All base64-inlined — the artifact CSP forbids remote fonts and images. |
| `places.json` | Natural Earth 1:10 m populated places inside the frame (1,282 settlements) used to name fire locations by nearest settlement, with distance and bearing. |
| `boundaries.json` | Simplified Natural Earth 1:50 m country outlines, clipped to 30°W–45°E / 30–72°N. Used both for point-in-polygon country attribution and for drawing the basemap. |

## What the pipeline does

1. Downloads eight CSVs from FIRMS — 24-hour and 7-day windows for VIIRS on Suomi-NPP,
   NOAA-20 and NOAA-21, plus MODIS C6.1 (Terra and Aqua). Open feeds; no API key needed,
   but they cap at 7 days.
2. Drops `confidence == "low"` detections.
3. Attributes each detection to a country by ray-casting against `boundaries.json`,
   with a bounding-box prefilter.
4. Groups detections into **complexes** by single-linkage clustering at 3 km, requiring
   at least 5 detections. Distances are computed in a Lambert conformal conic plane
   scaled to kilometres.
5. Projects everything to Lambert conformal conic (standard parallels 40°/60°N, central
   meridian 13°E), flipping y for SVG.
6. Emits a JSON payload — basemap paths, a 3 km density grid, complexes, country
   aggregates — and injects it into `template.html`.

## Editing the page

Copy changes, colours, layout and interaction all live in `template.html`. Every
data-dependent statement is generated at build time from the payload — the headline
claim, the country comparison in the bar-table lede, the map callout, the as-of stamp,
and the out-of-frame percentage. **Do not hardcode a country name, a percentage or a
date into the prose**; it will go stale on the next refresh and start lying.

## Two sensors, one analysis

`SATS` in `build.py` carries an `in_analysis` flag per platform. VIIRS 375 m drives
clustering, classification, FRP totals and the country table. MODIS 1 km is fetched for
coverage only — it appears in the density layer and the freshness figures and nowhere else.
Its footprint is seven times coarser, its detection threshold higher, and the classification
thresholds (`WF_PEAK` and friends) are calibrated on VIIRS radiative power; letting MODIS
into the totals would move every headline number silently. Note MODIS reports confidence as
0–100 rather than low/nominal/high, handled by `MODIS_MIN_CONF`.

## Scope

`EXCLUDE_CONTINENTS` in `build.py` drops detections attributed to those continents from
every figure. It is set to `{"Africa"}`: NASA's "Europe" feed reaches into the Maghreb and
those fires are out of scope. Add `"Asia"` for a strict continental filter, but note that
removes Turkey, Cyprus and the Caucasus, and Turkey is consistently among the largest
contributors of radiative power in the basin. African coastlines are still drawn for
cartographic context.

## Freshness

The page reports its own data age rather than warning about staleness in prose. `build.py`
measures the gap between the newest detection and build time, the per-platform lag, and the
longest unobserved gap in the window, and marks a platform STALE past
`STALE_PLATFORM_H` hours (default 12; normal end-to-end lag is around 4).

Hourly rebuilds cannot beat the upstream floor: NASA processes near-real-time detections
roughly three hours after observation. Five platforms across two sensors give about ten
passes a day, which held the longest unobserved gap to 3.4 h when this was written — adding
MODIS more than halved it from 8.8 h. There is no way to make the page live — the
artifact runtime grants no network-fetch capability, so a published page cannot call FIRMS
itself.

## Times

Every timestamp rendered on the page is **Europe/London**, labelled BST or GMT, because
the audience is UK-based. FIRMS source timestamps are UTC and are converted at build
time only. `uk()` uses `zoneinfo` when the container has tzdata and falls back to the EU
DST rule otherwise, so it never crashes on a minimal image.

## Failure behaviour

`build.py` exits non-zero and writes nothing if a feed fails after 4 retries, if a
response isn't a FIRMS CSV, if fewer than 500 7-day detections come back, or if the
template has lost its payload marker. A scheduled run that fails leaves the previously
published page untouched, which is the desired outcome — a stale page carrying an honest
older timestamp beats a broken or empty one.

## Caveats carried on the page itself

FRP is not burned area. Detections are not fires. Cloud and smoke block detection, and
each satellite passes roughly twice daily, so absence of detections is not evidence of
absence of fire. VIIRS near-real-time data is unvalidated and gets revised.

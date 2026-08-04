# europe-fire-map

Rebuilds the **Where Europe is burning right now** page from live NASA FIRMS VIIRS
active-fire detections, and republishes it to a fixed Claude artifact URL.

**Artifact (do not change):**
`https://claude.ai/code/artifact/200730f4-270f-4758-b887-70fcef594c1e`

Refreshed twice daily by a Claude Code cloud routine. Publishing to that exact URL is
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
| `boundaries.json` | Simplified Natural Earth 1:50 m country outlines, clipped to 30°W–45°E / 30–72°N. Used both for point-in-polygon country attribution and for drawing the basemap. |

## What the pipeline does

1. Downloads six CSVs from FIRMS — 24-hour and 7-day windows for Suomi-NPP, NOAA-20 and
   NOAA-21. These are the open feeds; no API key needed, but they cap at 7 days.
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

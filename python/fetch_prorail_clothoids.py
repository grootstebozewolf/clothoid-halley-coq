"""
Fetch all 'Overgangsboog' (clothoid transition) records from the public
ProRail Spoorgeometrie ArcGIS FeatureServer and write a normalised
snapshot to data/prorail_clothoids.json.gz.

Each output record carries the fields needed by the chord-length solver:

    {
      "objectid":      int,                # ArcGIS OBJECTID for traceability
      "P0":            [x_rd, y_rd],       # start point, RD New (EPSG:28992), metres
      "P1":            [x_rd, y_rd],       # end point
      "k0":            float,              # start curvature (1/m, signed)
      "k1":            float,              # end curvature (1/m, signed)
      "L_design":      float,              # design length, ELEMENT_LENGTE (m)
      "straal_begin":  float | None,       # raw radius at start (m); 0 / None -> straight
      "straal_eind":   float | None,       # raw radius at end
      "rotatie_begin": "CCW" | "CW" | None,
      "rotatie_eind":  "CCW" | "CW" | None
    }

Data is licensed CC BY 4.0 by ProRail (Spoorgeometrie). When this snapshot
is redistributed, that attribution applies; see data/LICENSE_DATA.txt.
"""
from __future__ import annotations

import gzip
import json
import math
import sys
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen, Request

URL = "https://maps.prorail.nl/arcgis/rest/services/Spoorgeometrie/FeatureServer/11/query"
PAGE_SIZE = 2000
OUT = Path(__file__).resolve().parent.parent / "data" / "prorail_clothoids.json.gz"


def _curv(straal: float | None, rotatie: str | None) -> float:
    """Convert (radius, rotation) -> signed curvature kappa = sign/straal.

    Straight endpoint (None, 0, or non-positive radius) -> 0.
    CCW = +, CW = -.
    """
    if straal is None or straal == 0 or not math.isfinite(straal) or straal <= 0:
        return 0.0
    sign = +1.0 if rotatie == "CCW" else (-1.0 if rotatie == "CW" else 0.0)
    return sign / straal


def _fetch_page(offset: int) -> dict:
    params = {
        "where": "ELEMENT_TYPE='Overgangsboog'",
        "outFields": "OBJECTID,STRAAL_BEGIN,STRAAL_EIND,ROTATIE_BEGIN,ROTATIE_EIND,ELEMENT_LENGTE",
        "returnGeometry": "true",
        "outSR": "28992",
        "f": "json",
        "resultOffset": offset,
        "resultRecordCount": PAGE_SIZE,
    }
    req = Request(f"{URL}?{urlencode(params)}", headers={"Accept": "application/json"})
    with urlopen(req, timeout=120) as r:
        return json.load(r)


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    records: list[dict] = []
    offset = 0
    page = 0
    t0 = time.time()
    while True:
        page += 1
        data = _fetch_page(offset)
        feats = data.get("features") or []
        if not feats:
            break
        for f in feats:
            attrs = f.get("attributes", {})
            geom = f.get("geometry", {})
            paths = geom.get("paths") or []
            if not paths or not paths[0]:
                continue
            poly = paths[0]
            P0 = poly[0]
            P1 = poly[-1]
            if len(P0) < 2 or len(P1) < 2:
                continue
            sb = attrs.get("STRAAL_BEGIN")
            se = attrs.get("STRAAL_EIND")
            rb = attrs.get("ROTATIE_BEGIN")
            re_ = attrs.get("ROTATIE_EIND")
            L_des = attrs.get("ELEMENT_LENGTE")
            if L_des is None or L_des <= 0 or not math.isfinite(L_des):
                continue
            records.append({
                "objectid":      attrs.get("OBJECTID"),
                "P0":            [float(P0[0]), float(P0[1])],
                "P1":            [float(P1[0]), float(P1[1])],
                "k0":            _curv(sb, rb),
                "k1":            _curv(se, re_),
                "L_design":      float(L_des),
                "straal_begin":  None if sb in (None, 0) else float(sb),
                "straal_eind":   None if se in (None, 0) else float(se),
                "rotatie_begin": rb if rb in ("CCW", "CW") else None,
                "rotatie_eind":  re_ if re_ in ("CCW", "CW") else None,
            })
        got = len(feats)
        sys.stderr.write(f"  page {page:3d}  +{got:4d}  total={len(records):6d}\n")
        if got < PAGE_SIZE:
            break
        offset += got
        time.sleep(0.05)
    elapsed = time.time() - t0
    sys.stderr.write(f"\nfetched {len(records)} clothoid records in {elapsed:.1f}s\n")

    payload = {
        "source":   URL,
        "layer":    11,
        "filter":   "ELEMENT_TYPE='Overgangsboog'",
        "license":  "CC BY 4.0 ProRail Spoorgeometrie",
        "fetched":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "count":    len(records),
        "records":  records,
    }
    with gzip.open(OUT, "wt", encoding="utf-8") as g:
        json.dump(payload, g, separators=(",", ":"))
    sys.stderr.write(f"wrote {OUT.relative_to(OUT.parent.parent)}  ({OUT.stat().st_size/1024:.1f} KB)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

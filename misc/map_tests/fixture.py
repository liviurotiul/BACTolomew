"""Shared sample-data fixture so every map backend renders the SAME thing.

Approximates the ST492 geographic distribution: per-country sample counts and
one point per sample (deterministically jittered around the country centroid).
"""
from __future__ import annotations
from typing import Dict, List

# country (display name) -> (sample_count, lat, lon)
_COUNTRY_DATA = {
    "Serbia": (24, 44.79, 20.45),
    "Romania": (22, 44.43, 26.10),
    "USA": (30, 38.90, -77.04),
    "China": (28, 39.90, 116.40),
    "Switzerland": (8, 46.95, 7.45),
    "Montenegro": (4, 42.44, 19.26),
    "Croatia": (2, 45.81, 15.98),
    "Hungary": (5, 47.50, 19.04),
    "Germany": (2, 52.52, 13.40),
    "Saudi Arabia": (4, 24.71, 46.68),
    "Slovenia": (3, 46.05, 14.51),
    "Kosovo": (1, 42.67, 21.17),
    "France": (2, 48.85, 2.35),
    "Lithuania": (1, 54.69, 25.28),
    "UK (Scotland)": (7, 55.95, -3.19),
    "Lebanon": (2, 33.89, 35.50),
    "Thailand": (1, 13.75, 100.50),
    "Italy": (1, 41.90, 12.50),
    "Canada": (2, 45.42, -75.70),
    "Japan": (2, 35.68, 139.69),
    "Sudan": (1, 15.50, 32.56),
    "India": (2, 28.61, 77.21),
    "Singapore": (1, 1.35, 103.82),
    "Brazil": (1, -15.79, -47.88),
    "Fiji": (1, -18.14, 178.44),
}


def _jitter(i: int) -> tuple[float, float]:
    """Deterministic small offset so co-located points fan out a little."""
    import math
    ang = (i * 2.399963)  # golden angle in radians
    r = 0.35 * (1 + (i % 5))
    return (r * math.cos(ang), r * math.sin(ang))


def make_fixture() -> Dict[str, object]:
    points: List[dict] = []
    countries: Dict[str, int] = {}
    for country, (count, lat, lon) in _COUNTRY_DATA.items():
        countries[country] = count
        for j in range(count):
            dlat, dlon = _jitter(j)
            points.append(
                {
                    "sample": f"{country[:3].upper()}_{j+1:02d}",
                    "country": country,
                    "lat": lat + dlat,
                    "lon": lon + dlon,
                }
            )
    max_count = max(countries.values())
    return {"points": points, "countries": countries, "max_count": max_count}


FIXTURE = make_fixture()
POINTS = FIXTURE["points"]
COUNTRIES = FIXTURE["countries"]
MAX_COUNT = FIXTURE["max_count"]

# Output geometry every backend targets.
WIDTH = 1100
HEIGHT = 620
TITLE = "ST492 isolate geographic distribution"

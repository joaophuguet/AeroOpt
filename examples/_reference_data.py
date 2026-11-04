"""Explicit, opt-in retrieval of pinned UIUC originals to ignored outputs only."""

import hashlib
import json
from datetime import datetime
from pathlib import Path
from urllib.request import urlopen
from zoneinfo import ZoneInfo

CATALOG_PATH = Path(__file__).with_name("reference_catalog.json")
POLICY_PAGES = ["https://m-selig.ae.illinois.edu/ads/coord_database.html",
                "https://m-selig.ae.illinois.edu/ads_faq.html",
                "https://m-selig.ae.illinois.edu/pd.html"]


def published_sources(directory: Path, *, fetch: bool = False) -> list[dict]:
    """Verify preserved raw bytes against the reviewed catalog, optionally fetch.

    No license permission is inferred from download access. Performance-data
    conditions on pd.html are not assumed to license coordinate redistribution.
    """
    catalog = json.loads(CATALOG_PATH.read_text(encoding="utf-8"))
    prior = directory / "provenance.json"
    dates = {row["file"]: row.get("obtained_on") for row in
             json.loads(prior.read_text(encoding="utf-8"))} if prior.exists() else {}
    records = []
    for entry in catalog:
        target = directory / entry["file"]
        if fetch:
            if target.exists():
                raise FileExistsError(f"Refusing to overwrite {target}")
            with urlopen(entry["url"], timeout=30) as response:
                raw = response.read()
            directory.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as stream:
                stream.write(raw)  # Preserve even an upstream revision for diagnosis.
        else:
            raw = target.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != entry["sha256"]:
            raise ValueError(f"Source hash changed for {target}: expected {entry['sha256']}, got {digest}; "
                             "review upstream revision explicitly")
        records.append({**entry, "sha256": digest, "byte_count": len(raw),
                        "obtained_on": datetime.now(ZoneInfo("America/Sao_Paulo")).date().isoformat()
                        if fetch else dates.get(entry["file"]),
                        "origin": "published coordinates, not asserted to be experimental measurements",
                        "coordinate_units": "x/c,y/c nominal chord units",
                        "ordering": "upper TE -> LE -> lower TE",
                        "license_pages_checked": POLICY_PAGES,
                        "redistribution": "local use only; bytes excluded from version control"})
    return records

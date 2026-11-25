"""The reviewed API removals recorded with the published fixture snapshot."""

import json
from pathlib import Path
from typing import cast


def removed_apis() -> frozenset[str]:
    environment = Path(__file__).resolve().parents[1] / "results/environment.json"
    parsed = cast(object, json.loads(environment.read_text(encoding="utf-8")))
    if not isinstance(parsed, dict):
        raise TypeError("Invalid published environment metadata")
    metadata = cast(dict[str, object], parsed)
    history = metadata.get("prune_history", {})
    if not isinstance(history, dict):
        raise TypeError("Invalid prune history")
    removed = cast(dict[str, object], history).get("removed_apis", [])
    if not isinstance(removed, list):
        raise TypeError("Invalid reviewed API removals")
    removed_items = cast(list[object], removed)
    if not all(isinstance(api, str) for api in removed_items):
        raise TypeError("Invalid reviewed API removals")
    return frozenset(cast(list[str], removed_items))

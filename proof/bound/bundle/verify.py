#!/usr/bin/env python3
"""Replay the fixed natural-bound arithmetic bundle; never generate a tree.

All paths resolve relative to this file, independently of the caller's cwd.
Metadata checks are deliberately distinguished from fresh arithmetic replay.
"""
from __future__ import annotations

import argparse
from collections import deque
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import threading
import time

BUNDLE = Path(__file__).resolve().parent
MANIFEST = BUNDLE / "proof-manifest.json"


def resolve(path: str) -> Path:
    return (BUNDLE / path).resolve()


def values_at(value, parts):
    if not parts:
        return [value]
    part, tail = parts[0], parts[1:]
    if part == "*":
        children = value.values() if isinstance(value, dict) else value
        assert isinstance(value, (dict, list)), f"wildcard on {type(value).__name__}"
        values = []
        for child in children:
            values.extend(values_at(child, tail))
        assert values, "empty wildcard cannot prove a universal report predicate"
        return values
    assert isinstance(value, dict) and part in value, f"missing report field {part}"
    return values_at(value[part], tail)


def check_report(spec):
    path = resolve(spec["path"])
    data = json.loads(path.read_text())
    for name, expected in spec["checks"].items():
        values = values_at(data, name.split("."))
        assert all(value == expected and type(value) is type(expected) for value in values), (
            f"{path.name}: {name} != {expected!r}"
        )
    return {"path": spec["path"], "report_predicates_passed": len(spec["checks"])}


def check_inputs(manifest):
    for artifact in manifest["input_artifacts"]:
        path = resolve(artifact["path"])
        raw = path.read_bytes()
        assert len(raw) == artifact["bytes"], f"input size changed: {path}"
        assert hashlib.sha256(raw).hexdigest() == artifact["sha256"], f"input hash changed: {path}"
    own_bytes = sum(p.stat().st_size for p in BUNDLE.rglob("*") if p.is_file())
    assert own_bytes <= manifest["budget"]["bundle_additional_limit_bytes"], "bundle exceeds 1MB"
    assert not manifest["budget"]["historical_031_tree_included"]
    assert not manifest["budget"]["new_search_tree_generated_by_entry"]
    return {"input_files": len(manifest["input_artifacts"]), "bundle_bytes": own_bytes}


def run_operation(unit):
    operation = unit["operation"]
    assert operation["kind"] in {"independent_replay", "verified_fixed_algebra_checker"}
    script = resolve(operation["script"])
    arguments = [str(resolve(arg["path"])) if isinstance(arg, dict) else arg
                 for arg in operation["args"]]
    command = [sys.executable, str(script), *arguments]
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=resolve(operation["cwd"]),
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, bufsize=1)
    # Retain bounded diagnostics and drain the pipe even when a legacy checker
    # prints a full JSON report. The entry itself emits only compact progress.
    tail = deque(maxlen=32)
    def drain():
        assert process.stdout is not None
        for line in process.stdout:
            tail.append(line[-600:])
    reader = threading.Thread(target=drain, daemon=True)
    reader.start()
    while True:
        try:
            exit_code = process.wait(timeout=30)
            break
        except subprocess.TimeoutExpired:
            print(f"  running {unit['id']} ({int(time.monotonic()-started)}s)", flush=True)
    reader.join()
    assert exit_code == 0, f"{unit['id']} exited {exit_code}: {''.join(tail)[-5000:]}"
    result = {"id": unit["id"], "exit_code": exit_code,
              "seconds": round(time.monotonic() - started, 3)}
    if "replay_report" in unit:
        result.update(check_report(unit["replay_report"]))
    return result


def selected_units(manifest, selected):
    units = {unit["id"]: unit for unit in manifest["units"]}
    if not selected:
        return manifest["units"]
    required = set()
    def include(name):
        assert name in units, f"unknown unit {name}"
        if name in required:
            return
        required.add(name)
        for dependency in units[name]["requires"]:
            include(dependency)
    for name in selected:
        include(name)
    return [unit for unit in manifest["units"] if unit["id"] in required]


def main():
    if not __debug__:
        raise RuntimeError("Verification requires assertions: disable -O and PYTHONOPTIMIZE.")
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--list", action="store_true", help="list fixed arithmetic units")
    mode.add_argument("--check-metadata", action="store_true", help="check hashes and existing evidence; default")
    mode.add_argument("--replay", action="store_true", help="fresh replay; canonical helper uses --deep")
    parser.add_argument("--unit", action="append", default=[], help="select a unit plus its dependencies")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    units = selected_units(manifest, args.unit)
    if args.list:
        for unit in units:
            print(f"{unit['id']}: {unit['status']} ({unit['group']})")
        return 0
    result = {"mode": "fresh_arithmetic_replay" if args.replay else "metadata_only",
              "success": False, "algorithm_version": manifest["algorithm_version"],
              "semantic_pending_items": manifest["semantic_pending_items"],
              "missing_evidence_reports": [],
              "units": []}
    output = BUNDLE / ("last-replay-run.json" if args.replay else "last-metadata-run.json")
    try:
        result.update(check_inputs(manifest))
        for index, unit in enumerate(units, 1):
            print(f"[{index}/{len(units)}] {unit['id']}", flush=True)
            if args.replay:
                checked = run_operation(unit)
            elif "current_evidence" in unit:
                if resolve(unit["current_evidence"]["path"]).is_file():
                    checked = {"id": unit["id"], **check_report(unit["current_evidence"])}
                else:
                    result["missing_evidence_reports"].append(unit["id"])
                    checked = {"id": unit["id"],
                               "metadata": "input hashes checked; generate evidence with --replay"}
            else:
                checked = {"id": unit["id"], "metadata": "source hash checked; identities require --replay"}
            result["units"].append(checked)
        result["success"] = True
        result["entire_registered_arithmetic_freshly_replayed"] = (
            args.replay and len(units) == len(manifest["units"])
        )
        result["manifest_has_no_pending_items"] = not manifest["semantic_pending_items"]
        result["all_registered_checks_passed_in_this_run"] = (
            result["entire_registered_arithmetic_freshly_replayed"]
        )
    except Exception as error:
        result["error"] = str(error)
        print(f"FAIL: {error}", file=sys.stderr, flush=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"success": result["success"], "mode": result["mode"],
                      "units_checked": len(result["units"]),
                      "missing_evidence_reports": len(result["missing_evidence_reports"]),
                      "entire_registered_arithmetic_freshly_replayed": result.get(
                          "entire_registered_arithmetic_freshly_replayed", False),
                      "pending_items": len(result["semantic_pending_items"]),
                      "report": str(output)}, ensure_ascii=False), flush=True)
    return 0 if result["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

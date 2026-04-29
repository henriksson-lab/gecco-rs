#!/usr/bin/env python3
"""Export GECCO's fitted sklearn type classifier for Rust inference."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


def _repo_root() -> pathlib.Path:
    return pathlib.Path(__file__).resolve().parents[1]


def _as_floats(values):
    return [float(v) for v in values]


def _export_tree(estimator):
    tree = estimator.tree_
    return {
        "children_left": [int(v) for v in tree.children_left.tolist()],
        "children_right": [int(v) for v in tree.children_right.tolist()],
        "feature": [int(v) for v in tree.feature.tolist()],
        "threshold": _as_floats(tree.threshold.tolist()),
        "value": tree.value.astype(float).tolist(),
    }


def export(output: pathlib.Path) -> None:
    root = _repo_root()
    sys.path.insert(0, str(root / "GECCO"))

    from gecco.types import TypeClassifier

    classifier = TypeClassifier.trained()
    model = classifier.model

    payload = {
        "n_features_in": int(model.n_features_in_),
        "outputs": [
            {"classes": _as_floats(classes.tolist())} for classes in model.classes_
        ],
        "estimators": [_export_tree(estimator) for estimator in model.estimators_],
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="\n") as fh:
        json.dump(payload, fh, separators=(",", ":"))
        fh.write("\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        default=_repo_root() / "GECCO/gecco/types/type_classifier.rf.json",
        help="Output JSON path",
    )
    args = parser.parse_args(argv)
    export(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

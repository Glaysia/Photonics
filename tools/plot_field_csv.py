#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class FieldMeta:
    nx: int
    ny: int
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    z: float
    component: str
    part: str
    time: float


def parse_header(line: str) -> FieldMeta:
    if not line.startswith("#"):
        raise ValueError("missing header comment")
    items: dict[str, str] = {}
    for token in line[1:].strip().split():
        if "=" not in token:
            continue
        k, v = token.split("=", 1)
        items[k.strip()] = v.strip()

    def req(name: str) -> str:
        if name not in items:
            raise ValueError(f"header missing {name}")
        return items[name]

    return FieldMeta(
        nx=int(req("nx")),
        ny=int(req("ny")),
        xmin=float(req("xmin")),
        xmax=float(req("xmax")),
        ymin=float(req("ymin")),
        ymax=float(req("ymax")),
        z=float(req("z")),
        component=req("component"),
        part=req("part"),
        time=float(req("time")),
    )


def load_field_csv(path: Path) -> tuple[FieldMeta, np.ndarray]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline()
    meta = parse_header(header)
    data = np.loadtxt(path, delimiter=",", comments="#")
    if data.shape != (meta.ny, meta.nx):
        raise ValueError(f"unexpected array shape {data.shape}, expected {(meta.ny, meta.nx)}")
    return meta, data


def load_holes_json(path: Path) -> list[tuple[float, float, float]]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    holes = []
    for h in obj.get("holes", []):
        holes.append((float(h["x"]), float(h["y"]), float(h["r"])))
    return holes


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot a Meep-exported field slice CSV with hole overlay.")
    ap.add_argument("--field", type=Path, required=True, help="CSV from write_field_slice_csv")
    ap.add_argument("--geometry", type=Path, required=False, help="Geometry JSON from --export-geometry (optional)")
    ap.add_argument("--out", type=Path, required=True, help="Output PNG path")
    ap.add_argument("--title", type=str, default="", help="Plot title")
    ap.add_argument("--cmap", type=str, default="RdBu_r", help="Matplotlib colormap")
    ap.add_argument("--percentile", type=float, default=99.0, help="Color limit percentile (symmetric)")
    ap.add_argument("--no-holes", action="store_true", help="Do not overlay hole circles")
    ap.add_argument("--alpha-holes", type=float, default=0.9, help="Hole outline alpha")
    ap.add_argument("--lw-holes", type=float, default=1.2, help="Hole outline line width")
    args = ap.parse_args()

    meta, data = load_field_csv(args.field)
    vmax = np.percentile(np.abs(data), float(args.percentile))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    fig, ax = plt.subplots(figsize=(8.0, 3.2), constrained_layout=True)
    im = ax.imshow(
        data,
        origin="lower",
        extent=(meta.xmin, meta.xmax, meta.ymin, meta.ymax),
        cmap=args.cmap,
        vmin=-vmax,
        vmax=vmax,
        interpolation="bilinear",
    )
    fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)

    if args.title:
        ax.set_title(args.title)
    else:
        ax.set_title(f"{meta.component} ({meta.part}) @ t={meta.time:.1f}, z={meta.z:g}")

    ax.set_xlabel("x (a.u.)")
    ax.set_ylabel("y (a.u.)")
    ax.set_aspect("equal", adjustable="box")

    if args.geometry and not args.no_holes:
        holes = load_holes_json(args.geometry)
        for x, y, r in holes:
            circ = plt.Circle((x, y), r, fill=False, color="k", alpha=float(args.alpha_holes), lw=float(args.lw_holes))
            ax.add_patch(circ)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


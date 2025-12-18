#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
from dataclasses import dataclass


@dataclass(frozen=True)
class Hole:
    x: float
    y: float
    r: float


def read_csv(path: str) -> list[Hole]:
    holes: list[Hole] = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            holes.append(Hole(float(row["x"]), float(row["y"]), float(row["r"])))
    return holes


def read_json(path: str) -> tuple[str, list[Hole]]:
    with open(path, "r") as f:
        data = json.load(f)
    units = str(data.get("units", "unknown"))
    holes: list[Hole] = []
    for h in data.get("holes", []):
        holes.append(Hole(float(h["x"]), float(h["y"]), float(h["r"])))
    return units, holes


def nice_step(span: float, target_ticks: int) -> float:
    if not math.isfinite(span) or span <= 0:
        return 1.0
    target_ticks = max(2, target_ticks)
    raw = span / (target_ticks - 1)
    exp = 10 ** math.floor(math.log10(raw))
    frac = raw / exp
    if frac <= 1:
        nice = 1
    elif frac <= 2:
        nice = 2
    elif frac <= 5:
        nice = 5
    else:
        nice = 10
    return nice * exp


def ticks(min_v: float, max_v: float, target_ticks: int) -> list[float]:
    if min_v > max_v:
        min_v, max_v = max_v, min_v
    span = max_v - min_v
    step = nice_step(span, target_ticks)
    if step <= 0:
        return [min_v, max_v]

    start = math.ceil(min_v / step) * step
    end = math.floor(max_v / step) * step

    values: list[float] = []
    v = start
    # Guard against infinite loops due to FP issues.
    for _ in range(10000):
        if v > end + 0.5 * step:
            break
        values.append(v)
        v += step

    if len(values) < 2:
        return [min_v, max_v]
    return values


def write_svg(
    path: str,
    holes: list[Hole],
    width: int,
    height: int,
    margin: int,
    invert_y: bool,
    units: str,
    tick_count: int,
    grid: bool,
) -> None:
    if not holes:
        raise SystemExit("No holes to plot")

    min_x = min(h.x - h.r for h in holes)
    max_x = max(h.x + h.r for h in holes)
    min_y = min(h.y - h.r for h in holes)
    max_y = max(h.y + h.r for h in holes)

    w = max(1e-9, max_x - min_x)
    h = max(1e-9, max_y - min_y)
    usable_w = max(1.0, width - 2 * margin)
    usable_h = max(1.0, height - 2 * margin)
    scale = min(usable_w / w, usable_h / h)

    def map_x(x: float) -> float:
        return margin + (x - min_x) * scale

    def map_y(y: float) -> float:
        v = margin + (y - min_y) * scale
        return (height - v) if invert_y else v

    def f3(v: float) -> str:
        if not math.isfinite(v):
            return "0"
        return f"{v:.3f}"

    def fmt_tick(v: float) -> str:
        av = abs(v)
        if av != 0 and (av < 1e-3 or av >= 1e4):
            return f"{v:.3e}"
        if abs(v - round(v)) < 1e-10:
            return str(int(round(v)))
        return f"{v:.3f}".rstrip("0").rstrip(".")

    def median(values: list[float]) -> float:
        vals = [v for v in values if math.isfinite(v)]
        if not vals:
            return float("nan")
        vals.sort()
        n = len(vals)
        mid = n // 2
        if n % 2 == 1:
            return vals[mid]
        return 0.5 * (vals[mid - 1] + vals[mid])

    def estimate_spacing_and_pair() -> tuple[float, tuple[Hole, Hole] | None]:
        if len(holes) < 2:
            return float("nan"), None
        nearest_dists: list[float] = []
        nearest_neighbor: list[int] = [-1] * len(holes)

        for i, a in enumerate(holes):
            best_d2 = float("inf")
            best_j = -1
            for j, b in enumerate(holes):
                if i == j:
                    continue
                dx = a.x - b.x
                dy = a.y - b.y
                d2 = dx * dx + dy * dy
                if d2 < best_d2:
                    best_d2 = d2
                    best_j = j
            if best_j >= 0 and math.isfinite(best_d2) and best_d2 > 0:
                nearest_neighbor[i] = best_j
                nearest_dists.append(math.sqrt(best_d2))

        a_est = median(nearest_dists)
        if not math.isfinite(a_est):
            return float("nan"), None

        # Choose a representative pair whose nearest-neighbor distance is closest to the median.
        best_i = -1
        best_err = float("inf")
        for i, j in enumerate(nearest_neighbor):
            if j < 0:
                continue
            dx = holes[i].x - holes[j].x
            dy = holes[i].y - holes[j].y
            d = math.hypot(dx, dy)
            err = abs(d - a_est)
            if err < best_err:
                best_err = err
                best_i = i

        if best_i < 0:
            return a_est, None
        j = nearest_neighbor[best_i]
        if j < 0:
            return a_est, None
        return a_est, (holes[best_i], holes[j])

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as out:
        out.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        out.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n')
        out.write(f'  <rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>\n')

        x_ticks = ticks(min_x, max_x, tick_count)
        y_ticks = ticks(min_y, max_y, tick_count)
        a_est, pair = estimate_spacing_and_pair()

        out.write("  <defs>\n")
        # Use auto-start-reverse so marker-start points outward too.
        out.write('    <marker id="arrow" markerWidth="10" markerHeight="10" refX="6" refY="3" orient="auto-start-reverse" markerUnits="strokeWidth">\n')
        out.write('      <path d="M0,0 L6,3 L0,6 Z" fill="#000000"/>\n')
        out.write("    </marker>\n")
        out.write("  </defs>\n")

        if grid:
            out.write('  <g fill="none" stroke="#e6e6e6" stroke-width="1">\n')
            for xv in x_ticks:
                xpx = map_x(xv)
                out.write(f'    <line x1="{f3(xpx)}" y1="{margin}" x2="{f3(xpx)}" y2="{height - margin}"/>\n')
            for yv in y_ticks:
                ypx = map_y(yv)
                out.write(f'    <line x1="{margin}" y1="{f3(ypx)}" x2="{width - margin}" y2="{f3(ypx)}"/>\n')
            out.write("  </g>\n")

        # Axes at plot borders (paper-style): left + bottom.
        out.write('  <g fill="none" stroke="#000000" stroke-width="2">\n')
        out.write(f'    <rect x="{margin}" y="{margin}" width="{width - 2 * margin}" height="{height - 2 * margin}" fill="none"/>\n')
        out.write(f'    <line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}"/>\n')
        out.write(f'    <line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}"/>\n')
        out.write("  </g>\n")

        # Ticks + labels.
        out.write('  <g fill="none" stroke="#000000" stroke-width="1">\n')
        tick_len = 8
        for xv in x_ticks:
            xpx = map_x(xv)
            out.write(f'    <line x1="{f3(xpx)}" y1="{height - margin}" x2="{f3(xpx)}" y2="{height - margin + tick_len}"/>\n')
        for yv in y_ticks:
            ypx = map_y(yv)
            out.write(f'    <line x1="{margin}" y1="{f3(ypx)}" x2="{margin - tick_len}" y2="{f3(ypx)}"/>\n')
        out.write("  </g>\n")

        out.write('  <g fill="#000000" font-family="Times New Roman, Times, serif" font-size="18">\n')
        for xv in x_ticks:
            xpx = map_x(xv)
            out.write(f'    <text x="{f3(xpx)}" y="{height - margin + 28}" text-anchor="middle">{fmt_tick(xv)}</text>\n')
        for yv in y_ticks:
            ypx = map_y(yv)
            out.write(f'    <text x="{margin - 14}" y="{f3(ypx + 6)}" text-anchor="end">{fmt_tick(yv)}</text>\n')
        x_label = f"x ({units})" if units else "x"
        y_label = f"y ({units})" if units else "y"
        out.write(f'    <text x="{(width // 2)}" y="{height - 10}" text-anchor="middle">{x_label}</text>\n')
        out.write(f'    <text x="18" y="{(height // 2)}" text-anchor="middle" transform="rotate(-90 18 {(height // 2)})">{y_label}</text>\n')
        out.write("  </g>\n")

        if pair is not None and math.isfinite(a_est) and a_est > 0:
            h1, h2 = pair
            x1 = map_x(h1.x)
            y1 = map_y(h1.y)
            x2 = map_x(h2.x)
            y2 = map_y(h2.y)

            vx = x2 - x1
            vy = y2 - y1
            vlen = math.hypot(vx, vy)
            if vlen > 1e-6:
                nx = -vy / vlen
                ny = vx / vlen
                offset = 22.0
                lx1 = x1 + nx * offset
                ly1 = y1 + ny * offset
                lx2 = x2 + nx * offset
                ly2 = y2 + ny * offset
                mx = 0.5 * (lx1 + lx2)
                my = 0.5 * (ly1 + ly2)

                label = f"a ≈ {fmt_tick(a_est)}"
                if units:
                    label += f" {units}"

                out.write('  <g fill="none" stroke="#000000" stroke-width="2" marker-start="url(#arrow)" marker-end="url(#arrow)">\n')
                out.write(f'    <line x1="{f3(lx1)}" y1="{f3(ly1)}" x2="{f3(lx2)}" y2="{f3(ly2)}"/>\n')
                out.write("  </g>\n")
                out.write('  <g fill="#000000" font-family="Times New Roman, Times, serif" font-size="18">\n')
                out.write(f'    <text x="{f3(mx)}" y="{f3(my - 10)}" text-anchor="middle">{label}</text>\n')
                out.write("  </g>\n")

        out.write('  <g fill="none" stroke="#000000" stroke-width="1">\n')
        for hole in holes:
            cx = map_x(hole.x)
            cy = map_y(hole.y)
            r = max(0.0, hole.r * scale)
            out.write(f'    <circle cx="{f3(cx)}" cy="{f3(cy)}" r="{f3(r)}"/>\n')
        out.write("  </g>\n")
        out.write("</svg>\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot photonic crystal holes into an SVG (no dependencies).")
    ap.add_argument("input", help="Input .csv or .json (exported by the C++ app)")
    ap.add_argument("-o", "--output", default="out/holes_py.svg", help="Output SVG path")
    ap.add_argument("--width", type=int, default=1200)
    ap.add_argument("--height", type=int, default=900)
    ap.add_argument("--margin", type=int, default=60)
    ap.add_argument("--no-invert-y", action="store_true")
    ap.add_argument("--units", default="", help="Units label for axes when reading CSV (e.g. um)")
    ap.add_argument("--ticks", type=int, default=6, help="Target tick count per axis")
    ap.add_argument("--no-grid", action="store_true", help="Disable light grid lines")
    args = ap.parse_args()

    units = args.units
    if args.input.endswith(".csv"):
        holes = read_csv(args.input)
    elif args.input.endswith(".json"):
        units, holes = read_json(args.input)
    else:
        raise SystemExit("Input must be .csv or .json")

    write_svg(
        args.output,
        holes,
        args.width,
        args.height,
        args.margin,
        invert_y=not args.no_invert_y,
        units=units,
        tick_count=args.ticks,
        grid=not args.no_grid,
    )
    print(f"Wrote {args.output} ({len(holes)} holes)")


if __name__ == "__main__":
    main()

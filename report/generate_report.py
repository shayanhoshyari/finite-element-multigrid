#!/usr/bin/env python3
"""Generate SVG figures plus Markdown and PDF benchmark reports."""

import csv
import math
from array import array
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.graphics import renderSVG
from reportlab.graphics.shapes import Circle, Drawing, Line, Path as ShapePath, Polygon, Rect, String
from reportlab.lib import colors


HERE = Path(__file__).resolve().parent
DATA = HERE / "generated"
METHOD_LABELS = {
    "mg-v": "V-cycle MG",
    "mg-w": "W-cycle MG",
    "gs": "single-grid GS",
    "pcg-gs": "PCG + GS",
    "pcg-mg-v": "PCG + V-cycle MG",
}
METHOD_COLORS = {
    "mg-v": colors.HexColor("#2166ac"),
    "mg-w": colors.HexColor("#b2182b"),
    "gs": colors.HexColor("#4d4d4d"),
    "pcg-gs": colors.HexColor("#1b7837"),
    "pcg-mg-v": colors.HexColor("#762a83"),
}


def mesh_sequence_figure():
    mesh_directory = HERE.parent / "meshes" / "dirt"
    panel_width, panel_height = 320, 330
    figure = Image.new("RGB", (4 * panel_width, panel_height), "white")
    font = ImageFont.load_default()

    for level in range(1, 5):
        node_path = mesh_directory / f"squaremg.{level}.node"
        element_path = mesh_directory / f"squaremg.{level}.ele"
        with node_path.open() as stream:
            vertex_count = int(stream.readline().split()[0])
            xx = array("d", [0.0]) * vertex_count
            yy = array("d", [0.0]) * vertex_count
            for line in stream:
                fields = line.split()
                if not fields or fields[0].startswith("#"):
                    continue
                vertex = int(fields[0])
                xx[vertex] = float(fields[1])
                yy[vertex] = float(fields[2])

        # All levels use the same complete unit-square view so refinement is
        # visually meaningful.
        crop_width = 1.0
        crop_min = 0.0
        plot_size = 280
        plot_top = 45
        panel = Image.new("RGB", (panel_width, panel_height), "white")
        draw = ImageDraw.Draw(panel)
        mesh_image = Image.new("RGB", (plot_size, plot_size), "white")
        mesh_draw = ImageDraw.Draw(mesh_image)
        draw.text((12, 8), f"Level {level}: {vertex_count:,} vertices", fill="black", font=font)
        draw.text((12, 23), "full unit square; element edges", fill="#555555", font=font)

        def point(vertex):
            return ((xx[vertex] - crop_min) * plot_size / crop_width,
                    plot_size - (yy[vertex] - crop_min) * plot_size / crop_width)

        with element_path.open() as stream:
            stream.readline()
            for line in stream:
                fields = line.split()
                if not fields or fields[0].startswith("#"):
                    continue
                vertices = (int(fields[1]), int(fields[2]), int(fields[3]))
                mesh_draw.line([point(vertices[0]), point(vertices[1]), point(vertices[2]), point(vertices[0])],
                               fill="#626262", width=1)
        mesh_draw.rectangle((0, 0, plot_size - 1, plot_size - 1), outline="black", width=1)
        panel.paste(mesh_image, (0, plot_top))
        figure.paste(panel.crop((0, 0, plot_size, panel_height)),
                     (((level - 1) % 4) * panel_width + 20, ((level - 1) // 4) * panel_height))

    output = DATA / "mesh_sequence.png"
    figure.save(output, optimize=True)


def field_heatmaps():
    panel_width, panel_height = 320, 360
    plot_size, plot_top = 280, 58
    figure = Image.new("RGB", (4 * panel_width, 2 * panel_height), "white")
    title_font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 14)
    subtitle_font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 11)

    def color(value, scale):
        normalized = max(-1.0, min(1.0, value / scale if scale else 0.0))
        if normalized < 0:
            fraction = normalized + 1
            return (round(35 + 220 * fraction), round(90 + 165 * fraction), 210)
        return (210, round(255 * (1 - normalized)), round(255 * (1 - normalized)))

    for level in range(1, 5):
        fields = read_csv(f"fields_level{level}.csv")
        xx = [float(row["x"]) for row in fields]
        yy = [float(row["y"]) for row in fields]
        solution = [float(row["solution"]) for row in fields]
        residual = [float(row["residual"]) for row in fields]
        residual_scale = max(abs(value) for value in residual)
        residual_norm = math.sqrt(sum(value * value for value in residual))
        triangles = []
        with (HERE.parent / "meshes" / "dirt" / f"squaremg.{level}.ele").open() as stream:
            stream.readline()
            for line in stream:
                parts = line.split()
                if parts and not parts[0].startswith("#"):
                    triangles.append((int(parts[1]), int(parts[2]), int(parts[3])))

        for row_index, (name, values, scale) in enumerate((
            ("final solution", solution, 1.0),
            ("final residual", residual, residual_scale),
        )):
            panel = Image.new("RGB", (panel_width, panel_height), "white")
            panel_draw = ImageDraw.Draw(panel)
            heatmap = Image.new("RGB", (plot_size, plot_size), "white")
            heatmap_draw = ImageDraw.Draw(heatmap)
            panel_draw.text((10, 7), f"Level {level}: {name}", fill="black", font=title_font)
            if row_index == 0:
                subtitle = "common scale: -1 (blue) to +1 (red)"
            else:
                subtitle = f"symmetric scale: +/-{scale:.2e}; ||r||2={residual_norm:.2e}"
            panel_draw.text((10, 29), subtitle, fill="#444444", font=subtitle_font)
            for triangle in triangles:
                points = [(xx[v] * (plot_size - 1), (1 - yy[v]) * (plot_size - 1)) for v in triangle]
                mean_value = sum(values[v] for v in triangle) / 3
                heatmap_draw.polygon(points, fill=color(mean_value, scale))
            heatmap_draw.rectangle((0, 0, plot_size - 1, plot_size - 1), outline="black", width=1)
            panel.paste(heatmap, (0, plot_top))
            x_position = (level - 1) * panel_width + 20
            y_position = row_index * panel_height
            figure.paste(panel.crop((0, 0, plot_size, panel_height)), (x_position, y_position))

    figure.save(DATA / "solution_residual_heatmaps.png", optimize=True)


def multigrid_cycle_diagram():
    width, height = 900, 330
    drawing = Drawing(width, height)
    drawing.add(Rect(0, 0, width, height, fillColor=colors.white, strokeColor=None))
    drawing.add(String(width / 2, height - 28, "Multigrid cycle shapes",
                       textAnchor="middle", fontName="Helvetica-Bold", fontSize=18))

    def arrow(x1, y1, x2, y2, color):
        drawing.add(Line(x1, y1, x2, y2, strokeColor=color, strokeWidth=3))
        angle = math.atan2(y2 - y1, x2 - x1)
        size = 8
        points = [x2, y2]
        for offset in (2.55, -2.55):
            points.extend([x2 + size * math.cos(angle + offset),
                           y2 + size * math.sin(angle + offset)])
        drawing.add(Polygon(points, fillColor=color, strokeColor=color))

    def panel(x_offset, title, levels):
        panel_width = 390
        level_y = {0: 245, 1: 175, 2: 105}
        drawing.add(String(x_offset + panel_width / 2, 275, title,
                           textAnchor="middle", fontName="Helvetica-Bold", fontSize=16))
        for level, label in ((0, "fine grid"), (1, "intermediate grid"), (2, "coarse grid")):
            guide = Line(x_offset + 30, level_y[level], x_offset + panel_width - 10, level_y[level],
                         strokeColor=colors.HexColor("#dddddd"), strokeWidth=1)
            guide.strokeDashArray = [4, 5]
            drawing.add(guide)
            drawing.add(String(x_offset + 30, level_y[level] + 8, label,
                               fontSize=10, fillColor=colors.HexColor("#555555")))

        spacing = (panel_width - 110) / (len(levels) - 1)
        points = [(x_offset + 90 + index * spacing, level_y[level])
                  for index, level in enumerate(levels)]
        for index, ((x1, y1), (x2, y2)) in enumerate(zip(points, points[1:])):
            color = colors.HexColor("#2166ac") if y2 < y1 else colors.HexColor("#b2182b")
            arrow(x1, y1, x2, y2, color)
        for x, y in points:
            drawing.add(Circle(x, y, 5, fillColor=colors.white,
                               strokeColor=colors.HexColor("#333333"), strokeWidth=2))

    panel(20, "V-cycle: one recursive coarse visit", [0, 1, 2, 1, 0])
    panel(470, "W-cycle: two recursive coarse visits", [0, 1, 2, 1, 2, 1, 0])
    drawing.add(String(90, 52, "blue: restrict the residual", fontSize=12,
                       fillColor=colors.HexColor("#2166ac")))
    drawing.add(String(300, 52, "red: prolongate and correct", fontSize=12,
                       fillColor=colors.HexColor("#b2182b")))
    drawing.add(String(640, 52, "smoothing is applied around coarse-grid correction",
                       textAnchor="middle", fontSize=12))
    renderSVG.drawToFile(drawing, str(DATA / "multigrid_cycles.svg"))


def read_csv(name):
    with (DATA / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def transform(value, logarithmic):
    return math.log10(max(value, 1e-300)) if logarithmic else value


def scientific_latex(value, digits=3):
    mantissa, exponent = f"{value:.{digits}e}".split("e")
    return f"{mantissa} \\times 10^{{{int(exponent)}}}"


def log_log_slope(rows):
    points = rows[-4:]
    xx = [math.log(row["vertices"]) for row in points]
    yy = [math.log(row["solve_seconds"]) for row in points]
    x_mean = sum(xx) / len(xx)
    y_mean = sum(yy) / len(yy)
    return sum((x - x_mean) * (y - y_mean) for x, y in zip(xx, yy)) / sum(
        (x - x_mean) ** 2 for x in xx
    )


def chart(series, title, x_label, y_label, x_log=False, y_log=False, legend_position="top-left"):
    width, height = 900, 500
    left, right, bottom, top = 95, 30, 70, 60
    plot_width, plot_height = width - left - right, height - bottom - top
    all_x = [x for item in series for x, _ in item[1]]
    all_y = [y for item in series for _, y in item[1] if y > 0]
    tx = [transform(x, x_log) for x in all_x]
    ty = [transform(y, y_log) for y in all_y]
    xmin, xmax = min(tx), max(tx)
    ymin, ymax = min(ty), max(ty)
    if xmin == xmax: xmax += 1
    if ymin == ymax: ymax += 1

    def sx(x): return left + (transform(x, x_log) - xmin) * plot_width / (xmax - xmin)
    def sy(y): return bottom + (transform(y, y_log) - ymin) * plot_height / (ymax - ymin)

    drawing = Drawing(width, height)
    drawing.add(Rect(0, 0, width, height, fillColor=colors.white, strokeColor=None))
    drawing.add(String(width / 2, height - 25, title, textAnchor="middle", fontName="Helvetica-Bold", fontSize=18))
    drawing.add(Line(left, bottom, left, bottom + plot_height, strokeColor=colors.black))
    drawing.add(Line(left, bottom, left + plot_width, bottom, strokeColor=colors.black))

    def tick_values(minimum, maximum, logarithmic):
        if not logarithmic:
            return [minimum + i * (maximum - minimum) / 5 for i in range(6)]
        first = math.ceil(minimum)
        last = math.floor(maximum)
        step = max(1, math.ceil(max(1, last - first + 1) / 6))
        values = list(range(first, last + 1, step))
        if values and values[-1] != last:
            values.append(last)
        return values or [minimum, maximum]

    for x_value in tick_values(xmin, xmax, x_log):
        xp = left + (x_value - xmin) * plot_width / (xmax - xmin)
        drawing.add(Line(xp, bottom, xp, bottom + plot_height, strokeColor=colors.HexColor("#e5e5e5")))
        x_text = f"10^{int(x_value)}" if x_log else f"{x_value:.2g}"
        drawing.add(String(xp, bottom - 20, x_text, textAnchor="middle", fontSize=11))

    for y_value in tick_values(ymin, ymax, y_log):
        yp = bottom + (y_value - ymin) * plot_height / (ymax - ymin)
        drawing.add(Line(left, yp, left + plot_width, yp, strokeColor=colors.HexColor("#e5e5e5")))
        y_text = f"10^{int(y_value)}" if y_log else f"{y_value:.2g}"
        drawing.add(String(left - 10, yp - 3, y_text, textAnchor="end", fontSize=11))

    drawing.add(String(left + plot_width / 2, 20, x_label, textAnchor="middle", fontSize=14))
    drawing.add(String(left, bottom + plot_height + 10, y_label, fontSize=14))

    legend_x = left + 8
    if legend_position == "bottom-left":
        legend_y = bottom + 12 + 20 * (len(series) - 1)
    else:
        legend_y = bottom + plot_height - 12
    for index, item in enumerate(series):
        name, values, color = item[:3]
        dash = item[3] if len(item) > 3 else None
        path = ShapePath()
        for point_index, (x, y) in enumerate(values):
            if point_index == 0: path.moveTo(sx(x), sy(y))
            else: path.lineTo(sx(x), sy(y))
        path.strokeColor = color
        path.strokeWidth = 3
        path.strokeDashArray = dash
        path.fillColor = None
        drawing.add(path)
        if len(values) < 20 and not dash:
            for x, y in values:
                drawing.add(Circle(sx(x), sy(y), 3, fillColor=color, strokeColor=color))
        y_position = legend_y - 20 * index
        legend_line = Line(legend_x, y_position, legend_x + 28, y_position,
                           strokeColor=color, strokeWidth=3)
        legend_line.strokeDashArray = dash
        drawing.add(legend_line)
        drawing.add(String(legend_x + 35, y_position - 4, name, fontSize=12))
    return drawing


def main():
    results = read_csv("results.csv")
    history = read_csv("history.csv")
    mesh_sequence_figure()
    field_heatmaps()
    multigrid_cycle_diagram()
    for row in results:
        for key in ("fine", "levels", "vertices", "triangles", "success", "iterations"):
            row[key] = int(row[key])
        for key in ("residual", "setup_seconds", "solve_seconds", "l1", "l2", "linf"):
            row[key] = float(row[key])
    grouped = defaultdict(list)
    for row in results: grouped[row["method"]].append(row)
    for rows in grouped.values(): rows.sort(key=lambda item: item["fine"])

    finest = max(row["fine"] for row in results)
    history_grouped = defaultdict(list)
    for row in history:
        if int(row["fine"]) == finest:
            history_grouped[row["method"]].append((int(row["iteration"]) + 1, float(row["residual"])))

    convergence = chart(
        [(METHOD_LABELS[m], history_grouped[m], METHOD_COLORS[m]) for m in METHOD_LABELS],
        f"Residual convergence on mesh level {finest}", "iteration + 1", "residual L2 norm",
        x_log=True, y_log=True, legend_position="bottom-left",
    )
    convergence_time_series = []
    for method in METHOD_LABELS:
        samples = history_grouped[method]
        total_time = next(row["solve_seconds"] for row in grouped[method] if row["fine"] == finest)
        timed_samples = [
            (total_time * sample_number / len(samples), residual)
            for sample_number, (_, residual) in enumerate(samples, start=1)
        ]
        convergence_time_series.append((METHOD_LABELS[method], timed_samples, METHOD_COLORS[method]))
    convergence_time = chart(
        convergence_time_series,
        f"Residual convergence versus solve time on mesh level {finest}",
        "estimated cumulative solve time (s)", "residual L2 norm",
        x_log=True, y_log=True, legend_position="bottom-left",
    )
    timing_series = []
    for method in METHOD_LABELS:
        converged_rows = [row for row in grouped[method] if row["success"]]
        label = "single-grid GS (converged)" if method == "gs" else METHOD_LABELS[method]
        timing_series.append(
            (label, [(row["vertices"], row["solve_seconds"]) for row in converged_rows], METHOD_COLORS[method])
        )
    reference_vertices = sorted({row["vertices"] for row in results})
    reference_n = 30873
    reference_time = 0.1
    timing_series.extend([
        ("O(N) reference",
         [(n, reference_time * (n / reference_n)) for n in reference_vertices],
         colors.HexColor("#222222"), [10, 7]),
        ("O(N^1.5) reference",
         [(n, reference_time * (n / reference_n) ** 1.5) for n in reference_vertices],
         colors.HexColor("#777777"), [4, 6]),
    ])
    timing = chart(
        timing_series,
        "Solve time versus problem size", "vertices", "solve seconds", x_log=True, y_log=True,
    )
    accuracy = chart(
        [("measured L2 error", [(r["vertices"], r["l2"]) for r in grouped["mg-v"]], METHOD_COLORS["mg-v"])],
        "Finite-element discretization error", "vertices", "nodal L2 error", x_log=True, y_log=True,
    )
    renderSVG.drawToFile(convergence, str(DATA / "convergence.svg"))
    renderSVG.drawToFile(convergence_time, str(DATA / "convergence_time.svg"))
    renderSVG.drawToFile(timing, str(DATA / "timing.svg"))
    renderSVG.drawToFile(accuracy, str(DATA / "accuracy.svg"))

    summary_rows = []
    for level in sorted({row["fine"] for row in grouped["mg-v"]}):
        level_rows = {row["method"]: row for row in results if row["fine"] == level}
        summary_rows.append([
            level, level_rows["mg-v"]["vertices"],
            level_rows["mg-v"]["iterations"], level_rows["mg-w"]["iterations"],
            level_rows["gs"]["iterations"], level_rows["pcg-gs"]["iterations"],
            level_rows["pcg-mg-v"]["iterations"],
        ])

    finest_rows = {row["method"]: row for row in results if row["fine"] == finest}
    measured_slopes = {method: log_log_slope(grouped[method]) for method in METHOD_LABELS}
    markdown = [
        "# Multigrid Convergence for a Finite-Element Poisson Problem",
        "",
        "This repository implements and benchmarks geometric multigrid for a two-dimensional finite-element Poisson problem on non-nested triangular meshes.",
        "",
        "## Contents",
        "",
        "- [1. Overview](#1-overview)",
        "- [2. Problem and numerical method](#2-problem-and-numerical-method)",
        "- [3. How multigrid works](#3-how-multigrid-works)",
        "- [4. Mesh sequence](#4-mesh-sequence)",
        "- [5. Results](#5-results)",
        "- [6. Interpretation](#6-interpretation)",
        "",
        "## 1. Overview",
        "",
        "The experiment compares V-cycle and W-cycle multigrid with single-grid Gauss–Seidel and Preconditioned Conjugate Gradient methods. It measures residual convergence, solve time, scaling with mesh size, and finite-element discretization error.",
        "",
        "## 2. Problem and numerical method",
        "",
        "The code solves the following equation on the unit square with linear triangular finite elements:",
        "",
        "$$",
        "-\\nabla \\cdot (K \\nabla u) = f.",
        "$$",
        "",
        "The diffusion tensor is the identity matrix:",
        "",
        "$$",
        "K = I.",
        "$$",
        "",
        "The manufactured exact solution is:",
        "",
        "$$",
        "u(x,y) = \\sin(\\pi x)\\sin(2\\pi y).",
        "$$",
        "",
        "The source term is chosen so that this function is the exact solution:",
        "",
        "$$",
        "f(x,y) = 5\\pi^2 \\sin(\\pi x)\\sin(2\\pi y).",
        "$$",
        "",
        "The code approximates the Dirichlet boundary value $u=g$ with a large-coefficient Robin condition:",
        "",
        "$$",
        "\\nabla u \\cdot n + \\beta u = \\beta g,",
        "\\qquad \\beta = 10^{10}.",
        "$$",
        "",
        "Here, $n$ is the outward unit normal. Equivalently, the condition is $\\nabla u \\cdot n + \\beta(u-g)=0$. As $\\beta \\to \\infty$, it approaches the Dirichlet condition $u=g$.",
        "",
        "The report compares five solvers:",
        "",
        "- geometric multigrid (MG) with a V-cycle;",
        "- geometric multigrid with a W-cycle;",
        "- single-grid Gauss–Seidel (GS), with one residual-correction sweep per iteration;",
        "- Preconditioned Conjugate Gradient (PCG) with a GS preconditioner; and",
        "- PCG with one V-cycle as the preconditioner.",
        "",
        "Each multigrid cycle uses three pre-smoothing sweeps and three post-smoothing sweeps per level. The coarsest level uses 500 sweeps. Every solver uses the following stopping criterion:",
        "",
        "$$",
        "\\lVert r \\rVert_2 < 10^{-10}.",
        "$$",
        "",
        "## 3. How multigrid works",
        "",
        "### 3.1. Residual correction",
        "",
        "After finite-element discretization, the problem is a linear system:",
        "",
        "$$",
        "Ax=b.",
        "$$",
        "",
        "For a current approximation $x$, the residual is:",
        "",
        "$$",
        "r=b-Ax.",
        "$$",
        "",
        "A residual-correction method updates the solution by:",
        "",
        "$$",
        "\\delta x=A^*r,",
        "$$",
        "",
        "$$",
        "x \\leftarrow x+\\delta x,",
        "$$",
        "",
        "where $A^*$ is intended to approximate $A^{-1}$. If $A^*=A^{-1}$, one update gives the exact solution. Useful iterative methods choose an inexpensive $A^*$ that captures enough of the structure of $A$ to reduce the error quickly.",
        "",
        "**Richardson iteration:** The simplest choice is $A^*=I$, which adds the residual directly to the current solution. A relaxation factor is often included to control the size of the correction.",
        "",
        "**Gauss–Seidel (GS):** Write the matrix as $A=D+L+U$, where $D$, $L$, and $U$ are its diagonal, lower-triangular, and upper-triangular parts. A forward GS correction uses the inexpensive triangular inverse $(D+L)^{-1}$. This code uses symmetric GS, whose correction operator is:",
        "",
        "$$",
        "A^*_{\\mathrm{SGS}}=(D+U)^{-1}D(D+L)^{-1}.",
        "$$",
        "",
        "**Preconditioned Conjugate Gradient (PCG):** PCG applies a preconditioner $M^{-1}\\approx A^{-1}$ to the residual, then combines the resulting vectors into mutually $A$-conjugate search directions. This avoids repeatedly undoing progress made in earlier directions. In this report, $M^{-1}$ is either symmetric GS or one multigrid V-cycle.",
        "",
        "### 3.2. Two-level multigrid correction",
        "",
        "GS quickly reduces error that oscillates from one fine-grid vertex to the next, but it reduces smooth, long-wavelength error very slowly. Multigrid moves that smooth error to a coarse mesh, where it appears more oscillatory and is cheaper to correct.",
        "",
        "Starting from the fine-grid residual:",
        "",
        "$$",
        "r_h=b_h-A_hx_h,",
        "$$",
        "",
        "restriction transfers it to the coarse mesh:",
        "",
        "$$",
        "r_H=Rr_h.",
        "$$",
        "",
        "The coarse-grid error equation is:",
        "",
        "$$",
        "A_He_H=r_H.",
        "$$",
        "",
        "Here, this code obtains $A_H$ by rediscretizing the differential equation on the independently generated coarse mesh. The coarse equation is solved approximately—using many GS sweeps at the coarsest level or another multigrid cycle above it. The correction is interpolated back to the fine mesh and added to the solution:",
        "",
        "$$",
        "x_h \\leftarrow x_h+Pe_H.",
        "$$",
        "",
        "The meshes in this experiment are generated independently and are not nested. The prolongation operator $P$ uses barycentric interpolation to transfer a coarse-grid correction to fine-grid vertices. The restriction operator is the transpose of prolongation (the standard variational choice in finite-element multigrid):",
        "",
        "$$",
        "R=P^T.",
        "$$",
        "",
        "A complete two-level cycle therefore performs fine-grid pre-smoothing, restricts the residual, approximately solves the coarse error equation, prolongates the correction, and finally performs post-smoothing. Restriction is the fine-to-coarse operation $R$; prolongation is the coarse-to-fine operation $P$.",
        "",
        "### 3.3. V-cycles and W-cycles",
        "",
        "With more than two meshes, the coarse solve is performed recursively. A V-cycle visits each coarser level once before returning to the fine mesh. A W-cycle revisits coarse levels, spending more work on the coarse correction. It can be more robust, but each cycle is more expensive.",
        "",
        "<img src=\"report/generated/multigrid_cycles.svg\" alt=\"V-cycle and W-cycle multigrid diagrams\" width=\"100%\">",
        "",
        "## 4. Mesh sequence",
        "",
        "The meshes were regenerated with Jonathan Shewchuk's Triangle 1.6 and the original area constraints. Levels 1–8 contain 150, 525, 1,989, 7,813, 30,873, 127,786, 511,082, and 2,184,494 vertices, respectively.",
        "",
        "![The first four regenerated triangular meshes at a common full-domain scale](report/generated/mesh_sequence.png)",
        "",
        "Each panel shows the complete unit square at the same scale. The element edges are read directly from the corresponding `.ele` file. Levels 5–8 are omitted because their elements cannot be distinguished at the report's display resolution. These finer meshes are still included in all numerical results.",
        "",
        "## 5. Results",
        "",
        "### 5.1. Solver convergence",
        "",
        "The table shows how many iterations or multigrid cycles each solver required. A value of 2,000 in the GS column means that Gauss–Seidel reached the iteration limit before satisfying the residual tolerance.",
        "",
        "| Level | Vertices | V cycles | W cycles | GS sweeps | PCG+GS | PCG+V |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        markdown.append("| " + " | ".join(f"{value:,}" if isinstance(value, int) else str(value) for value in row) + " |")
    markdown += [
        "",
        "The first plot shows residual reduction as a function of iteration count.",
        "",
        "**Level-8 single-grid GS does not achieve any net residual reduction during the measured run.** Its residual starts at $1.708 \\times 10^{-2}$ and ends at $2.032 \\times 10^{-2}$ after 2,000 sweeps. It rises initially and then decreases too slowly even to recover its starting value.",
        "",
        "<img src=\"report/generated/convergence.svg\" alt=\"Residual convergence\" width=\"100%\">",
        "",
        "The next plot relates residual reduction to solve time. The original benchmark recorded total solve time, not a timestamp for each residual sample. The horizontal positions therefore use each method's measured average time per iteration. Setup time is excluded.",
        "",
        "<img src=\"report/generated/convergence_time.svg\" alt=\"Residual convergence versus estimated solve time\" width=\"100%\">",
        "",
        "### 5.2. Final solution and residual fields",
        "",
        "The first row shows the converged finite-element solution on mesh levels 1–4. All four solution plots use the same color scale. The second row shows the final algebraic residual:",
        "",
        "$$",
        "r = b - Au.",
        "$$",
        "",
        "Residual magnitudes differ between levels, so each residual plot uses its own symmetric color scale. The scale and $\\lVert r \\rVert_2$ are printed above each plot.",
        "",
        "<img src=\"report/generated/solution_residual_heatmaps.png\" alt=\"Final solution and residual heatmaps for mesh levels 1 through 4\" width=\"100%\">",
        "",
        "### 5.3. Runtime scaling",
        "",
        "This plot compares time to convergence as the number of mesh vertices, $N$, increases. It excludes setup work such as matrix assembly and construction of the inter-grid transfer operators. Single-grid GS is shown only for levels 1 and 2, where it converged. On levels 3–8, it reached the 2,000-sweep limit before satisfying the residual tolerance, so those invalid time-to-convergence points are omitted. The dashed $O(N)$ and $O(N^{3/2})$ lines are slope guides normalized at $N=30{,}873$; they are not fitted timing models.",
        "",
        "<img src=\"report/generated/timing.svg\" alt=\"Solve time versus problem size\" width=\"100%\">",
        "",
        f"A log–log fit over the four finest meshes gives exponents of {measured_slopes['mg-v']:.2f} for V-cycle MG, {measured_slopes['mg-w']:.2f} for W-cycle MG, {measured_slopes['pcg-mg-v']:.2f} for PCG with V-cycle MG, and {measured_slopes['pcg-gs']:.2f} for PCG with GS. The multigrid methods are therefore close to $O(N)$ over this range, while PCG with GS is close to $O(N^{{3/2}})$.",
        "",
        "### 5.4. Discretization accuracy",
        "",
        "This plot shows the nodal $L_2$ error of the V-cycle solution. All converged solvers produce the same discrete solution, so the plot measures finite-element discretization error rather than incomplete solver convergence.",
        "",
        "<img src=\"report/generated/accuracy.svg\" alt=\"Finite-element discretization error\" width=\"100%\">",
        "",
        "## 6. Interpretation",
        "",
        f"The V-cycle requires only 8–9 cycles on every mesh from level 2 through level 8. Its iteration count is therefore essentially independent of mesh size. The W-cycle requires 7–8 cycles, but each cycle does more work. On level 8, the W-cycle takes {finest_rows['mg-w']['solve_seconds']:.2f} s, compared with {finest_rows['mg-v']['solve_seconds']:.2f} s for the V-cycle.",
        "",
        f"Single-grid Gauss–Seidel reaches the 2,000-sweep limit on levels 3–8. On level 8, it takes {finest_rows['gs']['solve_seconds']:.1f} s and stops with $\\lVert r \\rVert_2 = {scientific_latex(finest_rows['gs']['residual'], 2)}$. GS-preconditioned CG converges, but its iteration count grows from 40 on level 2 to {finest_rows['pcg-gs']['iterations']:,} on level 8. In contrast, V-cycle-preconditioned CG requires only 6–8 iterations.",
        "",
        f"The regenerated level-8 nodal $L_2$ error is ${scientific_latex(finest_rows['mg-v']['l2'])}$. The archived 2018 value is $3.409150602 \\times 10^{{-7}}$, a difference of about 0.013%. The error decreases by a factor of about four each time the characteristic mesh spacing is halved. This is consistent with the expected $O(h^2)$ convergence of linear finite elements in the nodal error measure used here.",
        "",
        "Development and reproduction instructions are in [`DEVELOPMENT.md`](DEVELOPMENT.md), and the important `src/` interfaces are summarized in [`API.md`](API.md). Raw results and per-iteration residuals are retained in `report/generated/results.csv` and `report/generated/history.csv`.",
    ]
    (HERE.parent / "README.md").write_text("\n".join(markdown) + "\n")

if __name__ == "__main__":
    main()

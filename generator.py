"""
foam_generator.py
=================
Parametric generator for 1-D thermal simulations in OpenFOAM.
Supports: conduction, convection, planar / cylindrical / spherical geometry,
N material layers with individual thermal properties.

Usage:
    python foam_generator.py                        # interactive wizard
    python foam_generator.py --config config.yaml   # direct YAML mode
"""

import os
import sys
import math
import argparse
import yaml
from jinja2 import Environment, FileSystemLoader


# ─────────────────────────────────────────────────────────────
#  Utilities
# ─────────────────────────────────────────────────────────────

def compute_DT(layer: dict) -> float:
    """Thermal diffusivity  α = λ / (ρ·cp)  [m²/s]"""
    return layer["lambda"] / (layer["rho"] * layer["cp"])


def fmt(v: float, decimals: int = 6) -> str:
    return f"{v:.{decimals}f}"


# ─────────────────────────────────────────────────────────────
#  Sphere topology constants  (cube → sphere projection)
# ─────────────────────────────────────────────────────────────

_SQRT3_INV = 1.0 / math.sqrt(3)   # corner coordinate of unit cube inscribed in unit sphere
_SQRT2_INV = 1.0 / math.sqrt(2)   # midpoint coordinate for arc midpoints

# Signs (±1) of the 8 corners of the cube inscribed in the sphere.
# Vertex k at radius r  →  position = r * SQRT3_INV * (sx, sy, sz)
_SPH_SIGNS = [
    (+1, +1, +1),  # 0: +X+Y+Z
    (-1, +1, +1),  # 1: -X+Y+Z
    (-1, -1, +1),  # 2: -X-Y+Z
    (+1, -1, +1),  # 3: +X-Y+Z
    (+1, +1, -1),  # 4: +X+Y-Z
    (-1, +1, -1),  # 5: -X+Y-Z
    (-1, -1, -1),  # 6: -X-Y-Z
    (+1, -1, -1),  # 7: +X-Y-Z
]

# The 6 hex blocks: (a0,a1,a2,a3) = inner face vertex indices.
# Radial direction = d3 (v0→v4),  grading = (1 1 n_radial).
# Verified right-handed: cross(v0→v1, v0→v3) · (v0→v4_radial) > 0
_SPH_BLOCKS = [
    (0, 1, 2, 3),  # +Z face
    (0, 3, 7, 4),  # +X face
    (1, 5, 6, 2),  # -X face
    (0, 4, 5, 1),  # +Y face
    (2, 6, 7, 3),  # -Y face
    (7, 6, 5, 4),  # -Z face
]

# 12 cube edges: (v_from, v_to, (dx, dy, dz))
# Arc midpoint = r * SQRT2_INV * (dx, dy, dz)
_SPH_ARCS = [
    (0, 1, ( 0,  1,  1)), (1, 2, (-1,  0,  1)),
    (2, 3, ( 0, -1,  1)), (3, 0, ( 1,  0,  1)),
    (4, 5, ( 0,  1, -1)), (5, 6, (-1,  0, -1)),
    (6, 7, ( 0, -1, -1)), (7, 4, ( 1,  0, -1)),
    (0, 4, ( 1,  1,  0)), (1, 5, (-1,  1,  0)),
    (2, 6, (-1, -1,  0)), (3, 7, ( 1, -1,  0)),
]


# ─────────────────────────────────────────────────────────────
#  Planar geometry
# ─────────────────────────────────────────────────────────────

def build_planar_mesh(config: dict) -> dict:
    """
    Build vertices and blocks for a 1-D extruded planar geometry.
    Variation direction: X.  Y and Z are the transverse (extrusion) dimensions.
    """
    w = config["dimensions"]["width"]
    h = config["dimensions"]["height"]

    x_coords = [0.0]
    for layer in config["layers"]:
        x_coords.append(round(x_coords[-1] + layer["thickness"], 8))

    vertices = []
    for x in x_coords:
        vertices += [
            f"({fmt(x)} 0 0)",
            f"({fmt(x)} {fmt(w)} 0)",
            f"({fmt(x)} {fmt(w)} {fmt(h)})",
            f"({fmt(x)} 0 {fmt(h)})",
        ]

    blocks = []
    for i, layer in enumerate(config["layers"]):
        v = i * 4
        indices = f"{v} {v+4} {v+5} {v+1} {v+3} {v+7} {v+6} {v+2}"
        n = layer["cells"]
        blocks.append(
            f'hex ({indices}) {layer["name"]} ({n} 1 1) simpleGrading (1 1 1)'
        )

    last_vi = (len(x_coords) - 1) * 4
    return {
        "vertices":       vertices,
        "blocks":         blocks,
        "arcs":           [],
        "interior_faces": ["0 3 2 1"],
        "exterior_faces": [f"{last_vi} {last_vi+3} {last_vi+2} {last_vi+1}"],
        "radii":          [],
        "project_faces":  [],
    }


# ─────────────────────────────────────────────────────────────
#  Cylindrical geometry  (full 360°, 4 × 90° sectors per layer)
# ─────────────────────────────────────────────────────────────

def build_cylindrical_mesh(config: dict) -> dict:
    """
    Full 360° cylinder — 4 sectors of 90° per radial layer.
    Each layer generates 4 hex blocks with circular arcs → proper cylinder in ParaView.

    Vertex indexing:  j * 8 + k
      j = radius index  (0 = r_inner, n_layers = r_outer)
      k = 0..3  at z=0  (angles 0°, 90°, 180°, 270°)
      k = 4..7  at z=h  (same angles, upper level)

    Hex per sector k, layer i:
      (ib+k, ob+k, ob+k1, ib+k1,  ib+k+4, ob+k+4, ob+k1+4, ib+k1+4)
      radial direction = d1 (v0→v1),  grading (n_radial n_angular 1)
    """
    r_inner = config.get("r_inner", 0.01)
    h       = config["dimensions"]["height"]

    radii = [r_inner]
    for layer in config["layers"]:
        radii.append(round(radii[-1] + layer["thickness"], 8))

    # 4 angles × 2 Z levels = 8 vertices per radius
    a_main = [math.radians(a) for a in [0, 90, 180, 270]]

    vertices = []
    for r in radii:
        for a in a_main:
            vertices.append(f"({fmt(r * math.cos(a))} {fmt(r * math.sin(a))} 0)")
        for a in a_main:
            vertices.append(f"({fmt(r * math.cos(a))} {fmt(r * math.sin(a))} {fmt(h)})")

    # n_angular: cells in the circumferential direction per 90° sector
    # Higher = smoother cylinder in ParaView  (default 12 = 48 cells total around circle)
    n_ang = config.get("n_angular", 12)

    # Hex blocks: 4 sectors per radial layer
    # Block dimensions: (n_radial, n_angular, 1)
    blocks = []
    for i, layer in enumerate(config["layers"]):
        ib = i * 8
        ob = (i + 1) * 8
        n  = layer["cells"]
        for k in range(4):
            k1  = (k + 1) % 4
            idx = (f"{ib+k} {ob+k} {ob+k1} {ib+k1} "
                   f"{ib+k+4} {ob+k+4} {ob+k1+4} {ib+k1+4}")
            blocks.append(f'hex ({idx}) {layer["name"]} ({n} {n_ang} 1) simpleGrading (1 1 1)')

    # Circular arcs: midpoints at 45°, 135°, 225°, 315° for each radius
    a_mid = [math.radians(a) for a in [45, 135, 225, 315]]
    arcs  = []
    for j, r in enumerate(radii):
        base = j * 8
        for k in range(4):
            k1 = (k + 1) % 4
            xm = r * math.cos(a_mid[k])
            ym = r * math.sin(a_mid[k])
            arcs.append(f"arc {base+k}   {base+k1}   ({fmt(xm)} {fmt(ym)} 0)")
            arcs.append(f"arc {base+k+4} {base+k1+4} ({fmt(xm)} {fmt(ym)} {fmt(h)})")

    # Boundary faces
    # interior (r_inner): normal toward axis  → order (k, k+4, k1+4, k1)
    # exterior (r_outer): normal outward      → order (k, k1, k1+4, k+4)
    ib0     = 0
    ob_last = (len(radii) - 1) * 8

    interior_faces = []
    exterior_faces = []
    for k in range(4):
        k1 = (k + 1) % 4
        interior_faces.append(f"{ib0+k}     {ib0+k+4}     {ib0+k1+4}     {ib0+k1}")
        exterior_faces.append(f"{ob_last+k} {ob_last+k1} {ob_last+k1+4} {ob_last+k+4}")

    # Project faces onto cylinder surfaces — one entry per radius (not per layer!)
    # Iterating per layer would generate duplicates at shared intermediate radii.
    project_faces = []
    for j, _r in enumerate(radii):
        base = j * 8
        for k in range(4):
            k1 = (k + 1) % 4
            project_faces.append((f"cylinder_{j}",
                                   f"{base+k} {base+k1} {base+k1+4} {base+k+4}"))

    return {
        "vertices":       vertices,
        "blocks":         blocks,
        "arcs":           arcs,
        "interior_faces": interior_faces,
        "exterior_faces": exterior_faces,
        "radii":          radii,
        "project_faces":  project_faces,
        "height":         h,
    }


# ─────────────────────────────────────────────────────────────
#  Spherical geometry  (full sphere, 6-block cube→sphere projection)
# ─────────────────────────────────────────────────────────────

def build_spherical_mesh(config: dict) -> dict:
    """
    Full sphere via cube→sphere projection (6 blocks per radial layer).
    Each layer generates 6 hex blocks with circular arcs → proper sphere in ParaView.

    The idea: take a cube, project its 8 corners onto the sphere (r/√3 × (±1,±1,±1)),
    define 6 hex blocks (one per cube face), and let blockMesh project all interior
    face points onto the exact sphere surface via the 'faces { project }' section.

    Vertex indexing:  j * 8 + k
      j = radius index  (0 = r_inner, n_layers = r_outer)
      k = 0..7  = the 8 cube corners (see _SPH_SIGNS)

    Hex per block (a0,a1,a2,a3), layer i:
      (ib+a0, ib+a1, ib+a2, ib+a3,  ob+a0, ob+a1, ob+a2, ob+a3)
      radial direction = d3 (v0→v4),  grading (n_angular n_angular n_radial)
    """
    r_inner = config.get("r_inner", 0.01)

    radii = [r_inner]
    for layer in config["layers"]:
        radii.append(round(radii[-1] + layer["thickness"], 8))

    # 8 cube corners per sphere shell
    vertices = []
    for r in radii:
        c = r * _SQRT3_INV
        for sx, sy, sz in _SPH_SIGNS:
            vertices.append(f"({fmt(sx * c)} {fmt(sy * c)} {fmt(sz * c)})")

    # n_angular: cells per cube-face tangential direction
    # Higher = smoother sphere  (default 16: each face is 16×16 = 256 cells)
    n_ang = config.get("n_angular", 16)

    # 6 hex blocks per radial layer
    # Block dimensions: (n_angular, n_angular, n_radial)
    blocks = []
    for i, layer in enumerate(config["layers"]):
        ib = i * 8
        ob = (i + 1) * 8
        n  = layer["cells"]
        for a0, a1, a2, a3 in _SPH_BLOCKS:
            idx = (f"{ib+a0} {ib+a1} {ib+a2} {ib+a3} "
                   f"{ob+a0} {ob+a1} {ob+a2} {ob+a3}")
            blocks.append(f'hex ({idx}) {layer["name"]} ({n_ang} {n_ang} {n}) simpleGrading (1 1 1)')

    # 12 curved edges per sphere shell (arc midpoints lie exactly on the sphere)
    arcs = []
    for j, r in enumerate(radii):
        base = j * 8
        m    = r * _SQRT2_INV
        for va, vb, (dx, dy, dz) in _SPH_ARCS:
            arcs.append(f"arc {base+va} {base+vb} ({fmt(m*dx)} {fmt(m*dy)} {fmt(m*dz)})")

    # Boundary faces
    # interior (r_inner): reversed inner face → normal toward origin
    # exterior (r_outer): direct outer face   → normal away from origin
    ib0     = 0
    ob_last = (len(radii) - 1) * 8

    interior_faces = []
    exterior_faces = []
    for a0, a1, a2, a3 in _SPH_BLOCKS:
        interior_faces.append(f"{ib0+a3}     {ib0+a2}     {ib0+a1}     {ib0+a0}")
        exterior_faces.append(f"{ob_last+a0} {ob_last+a1} {ob_last+a2} {ob_last+a3}")

    # Project faces onto sphere surfaces — one entry per radius (not per layer!)
    # Iterating per layer would generate duplicates at shared intermediate radii.
    project_faces = []
    for j, _r in enumerate(radii):
        base = j * 8
        for a0, a1, a2, a3 in _SPH_BLOCKS:
            project_faces.append((f"sphere_{j}",
                                   f"{base+a0} {base+a1} {base+a2} {base+a3}"))

    return {
        "vertices":       vertices,
        "blocks":         blocks,
        "arcs":           arcs,
        "interior_faces": interior_faces,
        "exterior_faces": exterior_faces,
        "radii":          radii,
        "project_faces":  project_faces,
        "sph_blocks":     [list(b) for b in _SPH_BLOCKS],
    }


# ─────────────────────────────────────────────────────────────
#  Layer processing & validation
# ─────────────────────────────────────────────────────────────

def process_layers(config: dict) -> list:
    """
    Add computed DT (thermal diffusivity), x_start, and x_end to each layer.
    For radial geometries (cylindrical/spherical), x represents radius and
    starts from r_inner, not from 0.
    """
    geometry = config.get("geometry", "planar")
    x = config.get("r_inner", 0.0) if geometry != "planar" else 0.0
    processed = []
    for layer in config["layers"]:
        dt = compute_DT(layer)
        x_end = round(x + layer["thickness"], 8)
        processed.append({
            **layer,
            "DT":      dt,
            "x_start": x,
            "x_end":   x_end,
        })
        x = x_end
    return processed


def validate_config(config: dict):
    required_fields = ["lambda", "rho", "cp", "thickness", "cells", "name"]
    for i, layer in enumerate(config["layers"]):
        for field in required_fields:
            if field not in layer:
                raise ValueError(
                    f"Layer {i} ('{layer.get('name', '?')}') is missing field '{field}' in config.yaml"
                )

    geom = config.get("geometry", "planar")
    if geom not in ("planar", "cylindrical", "spherical"):
        raise ValueError(
            f"Unknown geometry '{geom}'. Supported values: planar, cylindrical, spherical."
        )

    bc = config.get("boundary_conditions", {})
    for side_new, side_old in (("inner", "left"), ("outer", "right")):
        t = bc.get(side_new, bc.get(side_old, {})).get("type")
        if t not in ("fixedValue", "zeroGradient", "convection"):
            raise ValueError(
                f"Unknown boundary condition type for '{side_new}': '{t}'. "
                "Accepted values: fixedValue, zeroGradient, convection"
            )


# ─────────────────────────────────────────────────────────────
#  Case generation
# ─────────────────────────────────────────────────────────────

def create_case(config: dict, env: Environment):
    case_path = config["case_name"]
    geometry  = config.get("geometry", "planar")

    # Create directory structure
    for folder in ("0", "constant", "system"):
        os.makedirs(f"{case_path}/{folder}", exist_ok=True)

    # Process layers
    layers_computed = process_layers(config)
    # Default DT: innermost layer for planar; outermost layer for radial.
    # setFieldsDict applies layers from outside in, overwriting each time,
    # so the outermost layer must be the default field value.
    if geometry == "planar":
        DT_default = layers_computed[0]["DT"]
    else:
        DT_default = layers_computed[-1]["DT"]

    # Build mesh
    mesh_builders = {
        "planar":      build_planar_mesh,
        "cylindrical": build_cylindrical_mesh,
        "spherical":   build_spherical_mesh,
    }
    mesh = mesh_builders[geometry](config)

    # Boundary conditions (support both new keys and legacy left/right)
    bc       = config["boundary_conditions"]
    bc_inner = bc.get("inner", bc.get("left"))
    bc_outer = bc.get("outer", bc.get("right"))

    # Initial temperature = arithmetic mean of the two boundary values
    T_inner = bc_inner.get("value", bc_inner.get("Tf", 300.0))
    T_outer = bc_outer.get("value", bc_outer.get("Tf", 300.0))
    T_init  = round((T_inner + T_outer) / 2, 2)

    # Template context shared by all templates
    ctx = {
        **config,
        **mesh,
        "geometry":        geometry,
        "layers_computed": layers_computed,
        "DT_default":      DT_default,
        "T_init":          T_init,
        "bc_inner":        bc_inner,
        "bc_outer":        bc_outer,
    }

    def render(template_name: str, dest_path: str):
        content = env.get_template(template_name).render(ctx)
        with open(dest_path, "w") as f:
            f.write(content)
        print(f"  [OK] {dest_path}")

    print(f"\nGenerating case: {case_path}  |  geometry: {geometry}")
    print(f"  Layers: {len(layers_computed)}")
    for layer in layers_computed:
        print(f"    • {layer['name']:15s}  λ={layer['lambda']:.3f} W/m·K  "
              f"ρ={layer['rho']:.1f} kg/m³  cp={layer['cp']:.1f} J/kg·K  "
              f"α={layer['DT']:.3e} m²/s")
    print()

    # system/
    for tmpl, fname in [
        ("controlDict.j2",        "controlDict"),
        ("blockMeshDict.j2",      "blockMeshDict"),
        ("fvSchemes.j2",          "fvSchemes"),
        ("fvSolution.j2",         "fvSolution"),
        ("setFieldsDict.j2",      "setFieldsDict"),
    ]:
        render(tmpl, f"{case_path}/system/{fname}")

    # constant/
    render("physicalProperties.j2", f"{case_path}/constant/physicalProperties")

    # 0/
    render("T.j2",  f"{case_path}/0/T")
    render("DT.j2", f"{case_path}/0/DT")

    # Run script
    run_script = f"""#!/bin/bash
# Auto-generated by foam_generator.py
cd "$(dirname "$0")"

echo ">>> blockMesh  (building the mesh)"
blockMesh

echo ">>> setFields  (assigning thermal diffusivity per layer)"
setFields

echo ">>> {config['simulation']['solver']}  (running the solver)"
{config['simulation']['solver']}

echo ">>> Done! Visualise results with:"
echo "    paraFoam"
"""
    script_path = f"{case_path}/Allrun"
    with open(script_path, "w") as f:
        f.write(run_script)
    os.chmod(script_path, 0o755)
    print(f"  [OK] {script_path}")

    print("\n" + "─" * 55)
    print("Case generated successfully!")
    print(f"\nRun it:\n  cd {case_path}\n  ./Allrun")
    print("─" * 55)


# ─────────────────────────────────────────────────────────────
#  Interactive wizard
# ─────────────────────────────────────────────────────────────

MATERIAL_PRESETS = {
    "steel":       {"lambda": 50.0,  "rho": 7800.0, "cp": 500.0},
    "copper":      {"lambda": 385.0, "rho": 8900.0, "cp": 385.0},
    "aluminium":   {"lambda": 205.0, "rho": 2700.0, "cp": 900.0},
    "concrete":    {"lambda": 1.7,   "rho": 2300.0, "cp": 880.0},
    "wood":        {"lambda": 0.15,  "rho": 600.0,  "cp": 1700.0},
    "glass_wool":  {"lambda": 0.04,  "rho": 30.0,   "cp": 1000.0},
    "polystyrene": {"lambda": 0.035, "rho": 20.0,   "cp": 1450.0},
    "air":         {"lambda": 0.026, "rho": 1.2,    "cp": 1005.0},
    "water":       {"lambda": 0.6,   "rho": 1000.0, "cp": 4182.0},
    "glass":       {"lambda": 1.0,   "rho": 2500.0, "cp": 840.0},
}


def clr(code: str, text: str) -> str:
    """ANSI terminal colour helper."""
    codes = {"cyan": "\033[96m", "yellow": "\033[93m", "green": "\033[92m",
             "red": "\033[91m", "bold": "\033[1m", "reset": "\033[0m", "dim": "\033[2m"}
    return f"{codes.get(code, '')}{text}{codes['reset']}"


def ask(prompt: str, default=None, cast=str, validator=None):
    """Read a value from stdin with optional default and validation."""
    default_str = f"  [{clr('dim', str(default))}]" if default is not None else ""
    while True:
        raw = input(f"  {prompt}{default_str}: ").strip()
        if raw == "" and default is not None:
            return default
        try:
            val = cast(raw)
            if validator and not validator(val):
                raise ValueError
            return val
        except (ValueError, TypeError):
            print(clr("red", "    ✗ Invalid value. Please try again."))


def ask_choice(prompt: str, choices: list, default: int = 0) -> str:
    """Pick from a numbered list."""
    print(f"\n  {prompt}")
    for i, c in enumerate(choices):
        marker = clr("green", "▶") if i == default else " "
        print(f"    {marker} {i+1}. {c}")
    while True:
        raw = input(f"  Choose [1-{len(choices)}]  [{clr('dim', str(default+1))}]: ").strip()
        if raw == "":
            return choices[default]
        try:
            idx = int(raw) - 1
            if 0 <= idx < len(choices):
                return choices[idx]
        except ValueError:
            pass
        print(clr("red", f"    ✗ Please enter a number between 1 and {len(choices)}."))


def ask_bc(surface_label: str) -> dict:
    """Wizard for a single boundary condition."""
    print(f"\n  {clr('yellow', f'Boundary condition — {surface_label}:')}")
    bc_type = ask_choice(
        "Type:",
        ["fixedValue   — constant imposed temperature",
         "zeroGradient — insulated surface, zero heat flux",
         "convection   — heat exchange with a fluid (h + Tf)"],
    )
    if "fixedValue" in bc_type:
        val = ask("Temperature [K]", default=300.0, cast=float, validator=lambda x: x > 0)
        return {"type": "fixedValue", "value": val}
    elif "zeroGradient" in bc_type:
        return {"type": "zeroGradient"}
    else:
        h  = ask("Convection coefficient h [W/m²·K]", default=10.0,  cast=float, validator=lambda x: x > 0)
        Tf = ask("Fluid temperature Tf [K]",           default=300.0, cast=float, validator=lambda x: x > 0)
        return {"type": "convection", "h": h, "Tf": Tf}


def ask_layer(idx: int) -> dict:
    """Wizard for a single material layer."""
    print(f"\n  {clr('cyan', f'── Layer {idx} ──')}")

    preset_names = list(MATERIAL_PRESETS.keys()) + ["custom"]
    choice = ask_choice("Material (preset or custom):", preset_names, default=len(preset_names) - 1)

    name = ask("Layer name", default=choice if choice != "custom" else f"layer_{idx}", cast=str)

    if choice in MATERIAL_PRESETS:
        p   = MATERIAL_PRESETS[choice]
        lam = ask("λ thermal conductivity [W/m·K]", default=p["lambda"], cast=float, validator=lambda x: x > 0)
        rho = ask("ρ density [kg/m³]",              default=p["rho"],    cast=float, validator=lambda x: x > 0)
        cp  = ask("cp specific heat [J/kg·K]",      default=p["cp"],     cast=float, validator=lambda x: x > 0)
    else:
        lam = ask("λ thermal conductivity [W/m·K]", cast=float, validator=lambda x: x > 0)
        rho = ask("ρ density [kg/m³]",              cast=float, validator=lambda x: x > 0)
        cp  = ask("cp specific heat [J/kg·K]",      cast=float, validator=lambda x: x > 0)

    thickness = ask("Thickness [m]", default=0.01, cast=float, validator=lambda x: x > 0)
    cells     = ask("Cells across thickness", default=max(5, int(thickness / 0.001)),
                    cast=int, validator=lambda x: x >= 1)

    dt = lam / (rho * cp)
    print(clr("dim", f"    → α = λ/(ρ·cp) = {dt:.3e} m²/s"))

    return {"name": name, "lambda": lam, "rho": rho, "cp": cp,
            "thickness": thickness, "cells": cells, "convection": False}


def run_wizard() -> dict:
    """Full interactive wizard — returns a ready-to-use config dict."""
    print("\n" + clr("bold", "=" * 62))
    print(clr("bold", "   1-D THERMAL SIMULATION GENERATOR — OpenFOAM"))
    print(clr("bold", "=" * 62))
    print(clr("dim", "  Press Enter to accept the default value shown in [...]"))
    print(clr("dim", "  Temperatures in Kelvin  (300 K = 27°C  |  373 K = 100°C  |  800 K = 527°C)"))

    # ── Case name ─────────────────────────────────────────────
    print()
    case_name = ask("Case name (a folder with this name will be created)",
                    default="Thermal_Simulation", cast=str)

    # ── Geometry ──────────────────────────────────────────────
    geom_choice = ask_choice(
        "Geometry:",
        ["planar       — flat wall, temperature varies along X",
         "cylindrical  — tube / pipe, temperature varies radially",
         "spherical    — spherical shell, temperature varies radially"],
    )
    geometry = geom_choice.split()[0]   # extract "planar" / "cylindrical" / "spherical"

    # ── Dimensions ────────────────────────────────────────────
    n_ang = 8
    if geometry == "planar":
        print(f"\n  {clr('yellow', 'Extrusion dimensions (do not affect the physics, only visualisation):')}")
        width  = ask("Width Y [m]",  default=0.05, cast=float, validator=lambda x: x > 0)
        height = ask("Height Z [m]", default=0.05, cast=float, validator=lambda x: x > 0)
        dimensions = {"width": width, "height": height}
        r_inner = 0.0
    else:
        geom_label = "tube / pipe" if geometry == "cylindrical" else "spherical shell"
        print(f"\n  {clr('yellow', f'Dimensions — {geom_label}:')}")
        r_inner = ask(
            "INNER radius r_inner [m]  (inner surface, measured from the centre)",
            default=0.05, cast=float, validator=lambda x: x > 0,
        )
        if geometry == "cylindrical":
            height = ask(
                "Axial length Z [m]  (how long is the tube)",
                default=0.20, cast=float, validator=lambda x: x > 0,
            )
        else:
            height = r_inner * 0.01   # sphere: symbolic height, does not affect physics
        dimensions = {"width": r_inner, "height": height}
        print(clr("dim", "  ℹ  Outer radius = r_inner + sum of all layer thicknesses below"))
        ang_default = 16 if geometry == "spherical" else 12
        ang_hint    = "sphere: 16+ recommended" if geometry == "spherical" else "cylinder: 12+ recommended"
        n_ang = ask(
            f"Angular cells per sector  (more = smoother shape  |  {ang_hint})",
            default=ang_default, cast=int, validator=lambda x: x >= 2,
        )

    # ── Material layers ───────────────────────────────────────
    print(f"\n  {clr('yellow', 'Material layers  (order: INNER → OUTER):')}")
    print(clr("dim", "  ℹ  Layer 1 = closest to the centre / inner surface"))
    n_layers = ask("Number of layers", default=2, cast=int, validator=lambda x: x >= 1)
    layers = []
    for i in range(1, n_layers + 1):
        layers.append(ask_layer(i))

    # ── Boundary conditions ───────────────────────────────────
    if geometry == "planar":
        lbl_inner = "LEFT surface   (x = 0)"
        lbl_outer = "RIGHT surface  (x = total thickness)"
    elif geometry == "cylindrical":
        lbl_inner = f"INNER surface  (inner pipe wall, r = {r_inner} m)"
        lbl_outer =  "OUTER surface  (outer pipe wall, r = r_inner + total thickness)"
    else:
        lbl_inner = f"INNER surface  (inner shell, r = {r_inner} m)"
        lbl_outer =  "OUTER surface  (outer shell, r = r_inner + total thickness)"

    print(f"\n  {clr('yellow', 'Boundary conditions:')}")
    print(clr("dim", "  ℹ  fixedValue   — constant temperature imposed on the surface (e.g. 800 K)"))
    print(clr("dim", "  ℹ  zeroGradient — insulated surface, no heat flows through it"))
    print(clr("dim", "  ℹ  convection   — surface cooled/heated by a surrounding fluid"))
    bc_inner = ask_bc(lbl_inner)
    bc_outer = ask_bc(lbl_outer)

    # ── Simulation parameters ─────────────────────────────────
    print(f"\n  {clr('yellow', 'Simulation parameters:')}")
    print(clr("dim", "  ℹ  Run long enough for the temperature to reach a steady state"))
    end_time  = ask("End time [s]",                        default=500,  cast=float, validator=lambda x: x > 0)
    delta_t   = ask("Time step deltaT [s]",                default=0.1,  cast=float, validator=lambda x: x > 0)
    write_int = ask("Save results every [s]",              default=50,   cast=float, validator=lambda x: x > 0)

    # ── Summary ───────────────────────────────────────────────
    def bc_str(bc):
        if bc["type"] == "fixedValue":
            return f"T = {bc['value']} K  ({bc['value'] - 273.15:.1f}°C)"
        if bc["type"] == "zeroGradient":
            return "insulated (zero flux)"
        return f"convection  h = {bc['h']} W/m²K  Tf = {bc['Tf']} K"

    total_thick = sum(l["thickness"] for l in layers)
    print("\n" + clr("bold", "─" * 62))
    print(clr("green", "  SUMMARY:"))
    print(f"    Case:     {case_name}")
    print(f"    Geometry: {geometry.upper()}")
    if geometry == "planar":
        print(f"    Total thickness: {total_thick * 1000:.1f} mm  |  {n_layers} layers")
    else:
        r_out = r_inner + total_thick
        print(f"    r_inner = {r_inner*1000:.1f} mm  →  r_outer = {r_out*1000:.1f} mm  |  {n_layers} layers")
        print(f"    Angular cells per sector: {n_ang}")
    for layer in layers:
        dt = layer["lambda"] / (layer["rho"] * layer["cp"])
        print(f"      • {layer['name']:12s}  {layer['thickness']*1000:.1f} mm  "
              f"λ = {layer['lambda']} W/m·K  α = {dt:.2e} m²/s")
    print(f"    Inner BC: {bc_str(bc_inner)}")
    print(f"    Outer BC: {bc_str(bc_outer)}")
    print(f"    t_end = {end_time}s  |  dt = {delta_t}s  |  write every {write_int}s")
    print(clr("bold", "─" * 62))

    confirm = input(f"\n  {clr('green', 'Generate the case? [Y/n]')}: ").strip().lower()
    if confirm in ("n", "no"):
        print("Cancelled.")
        sys.exit(0)

    return {
        "case_name":  case_name,
        "geometry":   geometry,
        "dimensions": dimensions,
        "r_inner":    r_inner,
        "n_angular":  n_ang,
        "layers":     layers,
        "boundary_conditions": {"inner": bc_inner, "outer": bc_outer},
        "simulation": {
            "solver":        "laplacianFoam",
            "endTime":       end_time,
            "deltaT":        delta_t,
            "writeInterval": write_int,
            "maxCo":         0.5,
        },
    }


# ─────────────────────────────────────────────────────────────
#  Entry point
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Parametric OpenFOAM thermal simulation case generator"
    )
    parser.add_argument("--config", default=None,
                        help="YAML config file (omit for interactive wizard)")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    if args.config:
        with open(args.config, "r") as f:
            config = yaml.safe_load(f)
    else:
        config = run_wizard()
        # Save the generated config for reuse with --config
        out_yaml = config["case_name"] + "_config.yaml"
        os.makedirs(config["case_name"], exist_ok=True)
        with open(out_yaml, "w") as f:
            yaml.dump(config, f, allow_unicode=True, default_flow_style=False)
        print(clr("dim", f"\n  (Config saved to {out_yaml} — reuse with: --config {out_yaml})"))

    try:
        validate_config(config)
    except ValueError as e:
        print(clr("red", f"\n[ERROR] {e}"))
        sys.exit(1)

    env = Environment(loader=FileSystemLoader(os.path.join(script_dir, "templates")))
    create_case(config, env)


if __name__ == "__main__":
    main()

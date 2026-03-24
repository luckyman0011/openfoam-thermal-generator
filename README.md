
# OpenFOAM Thermal Generator

> A parametric 1-D thermal simulation case generator for OpenFOAM.

---

## Why does this exist?

My heat transfer professor doesn't allow internet access during exams.

That means every time I need to analyse heat conduction through a layered wall, a pipe, or a spherical shell, I have to work through the full analytical solution by hand — which is fine conceptually, but slow and error-prone under exam pressure when you're simultaneously trying to remember four boundary condition formulas and not spill coffee on your exam sheet.

So I built this tool. The idea is simple: describe your geometry and materials in a config file (or answer a few prompts in the terminal), and the generator produces a complete, ready-to-run OpenFOAM case. Run it, open ParaView, and you have a visual simulation of the temperature profile in seconds.

Will my professor actually let me use it during the exam? Genuinely no idea. Probably not — not because he's against it (he's actually one of the more tech-forward professors I've had, the kind that gets genuinely excited about numerical methods at 8am), but because the university regulation is what it is. So this might be the most elaborate study tool that never gets used during an actual exam.

Either way, building it taught me more about heat transfer, mesh topology, and numerical methods than any problem set had — which I suspect is exactly the kind of irony my professor would appreciate.

---

## What it does

Given a geometry and a stack of material layers, the generator produces:

- `blockMeshDict` — a structured hex mesh (flat wall, full cylinder, or full sphere)
- `T` — temperature field with correct boundary conditions
- `setFieldsDict` — thermal diffusivity assigned per layer
- All supporting OpenFOAM dictionaries (`fvSchemes`, `fvSolution`, `controlDict`, etc.)
- An `Allrun` script that builds and runs everything in one command

Three geometry modes are supported:

| Geometry | Description | Physics |
|----------|-------------|---------|
| `planar` | Flat wall, extruded slab | Temperature varies along X |
| `cylindrical` | Full 360° tube or pipe | Temperature varies radially |
| `spherical` | Full spherical shell | Temperature varies radially |

For cylindrical and spherical cases the mesh uses `searchableCylinder` / `searchableSphere` projection, so the geometry is mathematically exact — not a faceted approximation.

**What the solver actually does:**

The simulation runs `laplacianFoam`, which solves the transient heat conduction equation at every cell, at every time step:

$$\frac{\partial T}{\partial t} = \nabla \cdot (\alpha \, \nabla T)$$

where $\alpha = \frac{\lambda}{\rho \, c_p}$ is the thermal diffusivity $[\text{m}^2/\text{s}]$.

In plain English: the rate of change of temperature at any point depends on how much heat is flowing in minus how much is flowing out — weighted by how fast that material conducts heat. The solver marches this forward in time until you tell it to stop (or until the solution stops changing, which means steady state).

Each material layer gets its own $\alpha$ assigned by `setFields` before the solver starts, so the discontinuities at layer interfaces are handled correctly.

---

## Quick start

**Requirements:** OpenFOAM 13, Python 3.8+, `pyyaml`, `jinja2`


**Step 1 — Run the generator and answer the prompts:**

```bash
python3 generator.py
```

The wizard asks for geometry, material layers, boundary conditions, and simulation time.
At the end it creates a folder with the complete OpenFOAM case inside.

**Step 2 — Run the simulation:**

```bash
cd <your_case_name>
./Allrun
```

**Step 3 — Visualise:**

```bash
paraFoam
```

That's it. No manual file editing required.

---

## Config file overview

```yaml
geometry: "cylindrical"   # planar | cylindrical | spherical

r_inner: 0.05             # [m] inner radius (cylindrical/spherical only)
n_angular: 12             # angular cells per sector — higher = smoother shape

layers:
  - name: "steel"
    thickness: 0.01       # [m]
    cells: 10
    lambda: 50.0          # [W/m·K] thermal conductivity
    rho: 7800.0           # [kg/m³] density
    cp: 500.0             # [J/kg·K] specific heat

boundary_conditions:
  inner:
    type: "fixedValue"
    value: 800.0          # [K]
  outer:
    type: "fixedValue"
    value: 300.0          # [K]
```

Three boundary condition types: `fixedValue` (imposed temperature), `zeroGradient` (insulated surface), `convection` (heat exchange with a surrounding fluid).

---

## What the output should look like

Cylindrical case — three material layers, 800 K inner wall, 300 K outer wall:

```
Generating case: Thermal_Simulation  |  geometry: cylindrical
  Layers: 3
    • steel            λ=50.000 W/m·K  ρ=7800.0 kg/m³  cp=500.0 J/kg·K  α=1.282e-05 m²/s
    • insulation       λ=0.040 W/m·K   ρ=30.0 kg/m³    cp=1000.0 J/kg·K  α=1.333e-06 m²/s
    • copper           λ=385.000 W/m·K ρ=8900.0 kg/m³  cp=385.0 J/kg·K  α=1.124e-04 m²/s

  [OK] Thermal_Simulation/system/blockMeshDict
  [OK] Thermal_Simulation/system/setFieldsDict
  [OK] Thermal_Simulation/0/T
  ...

Case generated successfully!
Run it:
  cd Thermal_Simulation
  ./Allrun
```

---

## Technical notes

**Mesh topology:**
- Planar: standard hex extrusion along X
- Cylindrical: 4 × 90° sectors with circular arcs, full 360°
- Spherical: 6-block cube→sphere projection — each of the 6 cube faces becomes one curved block; `faces { project }` forces every interior face point onto the exact sphere surface

**Thermal diffusivity assignment:**
`setFields` applies layers from outside in, each overwriting the previous region. This handles the multi-layer case correctly regardless of geometry.

**Solver:** `laplacianFoam` (transient heat conduction). The field being solved is `T`; `DT` (thermal diffusivity α) is set per cell by `setFields`.

---

## License

MIT — use it, modify it, cite it if it helps you pass an exam.

# Sun Analemma – Reflection Simulation
Based on the following [research paper](https://drive.google.com/file/d/18Ge1nyZpSrt2elqEW2IscvNKPeJti0cU/view) authored by me.

This repository contains a Python implementation of the **Sun Analemma reflection model**, which simulates how sunlight reflects off the **concave glass façade** of 20 Fenchurch Street (the “Walkie-Talkie”) at various times of day. The program uses the Sun’s azimuth and altitude angles to calculate where the incoming solar rays strike the building and how they reflect onto the ground.

---

## Overview

Given:
- The **Sun’s azimuth** and **altitude** (in degrees)
- A **starting seed point** near the building façade  
the script:
1. Determines the direction vector of the incident sunlight.
2. Solves for the intersection between the incident ray and the curved façade.
3. Computes the façade’s surface normal at that point.
4. Calculates the reflected ray direction using the vector law of reflection.
5. Optionally finds where that reflected ray hits the ground (`y = 0`).

---
## Mathematical Model

### 1. Facade Geometry

The façade is modeled as a 3D paraboloid:

```
x = 0.0022(y² + z²) - 6.22
```

This represents the building’s concave surface with curvature coefficients derived from real façade measurements.

---

### 2. Incident Ray

Given the Sun’s azimuth **β** and altitude **α**:

```
dₓ = cos(β) * cos(α)
d_y = cos(β) * sin(α)
d_z = sin(β)
```

The incident ray is parameterized as:

```
(x, y, z) = (x₀, y₀, z₀) + λ(dₓ, d_y, d_z)
```

---

### 3. Intersection

Plugging this line into the façade equation yields a **quadratic in λ**:

```
aλ² + bλ + c = 0
```

where:

```
a = 0.0022(d_y² + d_z²)
b = 0.0044(y₀d_y + z₀d_z) - dₓ
c = 0.0022(y₀² + z₀²) - x₀ - 6.22
```

The script solves this quadratic and selects the smallest **positive real root** (the physical intersection).

---

### 4. Surface Normal

The surface normal is derived from the gradient of the façade function:

```
∇F = ⟨1, -0.0044y, -0.0044z⟩
```

---

### 5. Reflection

The reflection vector is computed using the **vector law of reflection**:

```
R = D - 2 * ((D · N) / |N|²) * N
```

where:
- `D` = incident direction vector  
- `N` = surface normal vector  
- `R` = reflected direction vector
```
---

## Usage

Run in any Python 3 environment:

```bash
python sun_analemma_reflection.py
```

Input values when prompted:
- **Azimuth (°)** – Sun’s compass bearing  
- **Altitude (°)** – Sun’s elevation angle  
- **Seed coordinates** – A point near the building façade (x₀, y₀, z₀)

Example:
```
Azimuth: 174.74
Altitude: 15.35
Seed X: -3.48
Seed Y: 33
Seed Z: 0
```

Output:
```
Intersection point on facade: x=-6.1172, y=32.7869, z=0.0378
Surface normal: Nx=1.0000, Ny=-0.1442, Nz=-0.0002
Reflected direction vector: Rx=0.8693, Ry=-0.4931, Rz=0.0081
Ground hit (y=0): x=51.06, y=0, z=0.84
```

---

## Physical Interpretation

- The **intersection point** is where the Sun’s ray hits the curved façade.
- The **reflected direction** indicates where the sunlight is directed after reflection.
- The **ground hit** shows where that reflection would appear on the street.
- 1 unit = 4 m (scaling factor used in the original analysis).
- See research paper for a visualization of these cordinates

---

## Full Source Code

```python
# sun_analemma_reflection.py
# Compute reflected solar ray directions from a concave facade model.

import math as mt

def deg_to_rad(angle_deg): return mt.radians(angle_deg)

def incident_dir_from_az_alt(azimuth_deg, altitude_deg):
    az = deg_to_rad(azimuth_deg)
    alt = deg_to_rad(altitude_deg)
    dx = mt.cos(az) * mt.cos(alt)
    dy = mt.cos(az) * mt.sin(alt)
    dz = mt.sin(az)
    return dx, dy, dz

def solve_lambda_for_intersection(x0, y0, z0, dx, dy, dz):
    a = 0.0022 * (dy**2 + dz**2)
    b = 0.0044 * (y0 * dy + z0 * dz) - dx
    c = 0.0022 * (y0**2 + z0**2) - x0 - 6.22
    disc = b*b - 4*a*c
    if disc < 0:
        raise ValueError("No real intersection (discriminant < 0).")
    sqrt_disc = mt.sqrt(disc)
    l1 = (-b + sqrt_disc) / (2*a) if a != 0 else (-c / b)
    l2 = (-b - sqrt_disc) / (2*a) if a != 0 else (-c / b)
    positives = [l for l in (l1, l2) if l > 0]
    return min(positives) if positives else min((l1, l2), key=abs)

def normal_at(x, y, z):
    return (1.0, -0.0044 * y, -0.0044 * z)

def reflect(D, N):
    dx, dy, dz = D
    nx, ny, nz = N
    dn = dx*nx + dy*ny + dz*nz
    nn = nx*nx + ny*ny + nz*nz
    rx = dx - 2.0 * dn / nn * nx
    ry = dy - 2.0 * dn / nn * ny
    rz = dz - 2.0 * dn / nn * nz
    return rx, ry, rz

def main():
    azimuth = float(input("Azimuth (deg): "))
    altitude = float(input("Altitude (deg): "))
    x0 = float(input("Seed X: "))
    y0 = float(input("Seed Y: "))
    z0 = float(input("Seed Z: "))

    dx, dy, dz = incident_dir_from_az_alt(azimuth, altitude)
    lam = solve_lambda_for_intersection(x0, y0, z0, dx, dy, dz)
    xi, yi, zi = x0 + lam*dx, y0 + lam*dy, z0 + lam*dz
    Nx, Ny, Nz = normal_at(xi, yi, zi)
    Rx, Ry, Rz = reflect((dx, dy, dz), (Nx, Ny, Nz))

    print(f"\nIntersection: ({xi:.4f}, {yi:.4f}, {zi:.4f})")
    print(f"Normal: ({Nx:.4f}, {Ny:.4f}, {Nz:.4f})")
    print(f"Reflected dir: ({Rx:.4f}, {Ry:.4f}, {Rz:.4f})")

    if abs(Ry) > 1e-12:
        t = -yi / Ry
        xg, yg, zg = xi + t*Rx, yi + t*Ry, zi + t*Rz
        print(f"Ground hit (y=0): ({xg:.4f}, {yg:.4f}, {zg:.4f})")
        print("Note: multiply by 4 for meters (1 unit = 4 m).")
    else:
        print("Reflected ray is parallel to ground; no intersection.")

if __name__ == "__main__":
    main()
```

---

## Example Application

Use this script to:
- Generate daily reflection coordinates for a range of azimuth/altitude pairs.
- Plot the **analemma pattern** of sunlight reflections on the ground.
- Visualize how façade curvature concentrates light (as in the Walkie-Talkie glare phenomenon).

---

## Notes
- Uses purely geometric optics (no atmospheric effects).
- The building’s curvature parameters (0.0022, 6.22) were fitted from façade data.
- Scaling factor: 1 unit = 4 m.
- Adjust the coefficients if modeling other parabolic façades.

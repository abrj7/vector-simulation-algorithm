# sun_analemma_reflection.py
# Compute reflected solar ray directions from a concave facade model.
# Based on the derivation and parameters in the associated write-up.

import math as mt

# ---------------------------
# Angle helpers
# ---------------------------
def deg_to_rad(angle_deg: float) -> float:
    """Convert degrees to radians."""
    return mt.radians(angle_deg)

# ---------------------------
# Spherical (SunCalc) -> direction
# ---------------------------
def incident_dir_from_az_alt(azimuth_deg: float, altitude_deg: float):
    """
    Map SunCalc-style azimuth (clockwise from North) and altitude (elevation)
    to a 3D direction vector (dx, dy, dz), following the paper's convention.

    Returns:
        (dx, dy, dz)
    """
    az = deg_to_rad(azimuth_deg)
    alt = deg_to_rad(altitude_deg)
    dx = mt.cos(az) * mt.cos(alt)
    dy = mt.cos(az) * mt.sin(alt)
    dz = mt.sin(az)
    return dx, dy, dz

# ---------------------------
# Intersect line with facade: x = 0.0022*(y^2 + z^2) - 6.22
# by solving a*λ^2 + b*λ + c = 0
# ---------------------------
def solve_lambda_for_intersection(x0, y0, z0, dx, dy, dz):
    """
    Solve quadratic for lambda where the parametric line
        (x, y, z) = (x0, y0, z0) + λ (dx, dy, dz)
    intersects the implicit facade:
        x - 0.0022*(y^2 + z^2) + 6.22 = 0  <=>  x = 0.0022*(y^2 + z^2) - 6.22

    Returns:
        best_lambda (float) or raises ValueError if no real root.
    """
    a = 0.0022 * (dy**2 + dz**2)
    b = 0.0044 * (y0 * dy + z0 * dz) - dx
    c = 0.0022 * (y0**2 + z0**2) - x0 - 6.22

    disc = b*b - 4*a*c
    if disc < 0:
        # Complex: physically meaningless in this setting
        raise ValueError("No real intersection (discriminant < 0).")

    sqrt_disc = mt.sqrt(disc)
    l1 = (-b + sqrt_disc) / (2*a) if a != 0 else (-c / b)
    l2 = (-b - sqrt_disc) / (2*a) if a != 0 else (-c / b)

    # Prefer the smallest positive root; otherwise fall back to the
    # root with the smallest absolute value (matches original paper's workflow).
    candidates = [l1, l2]
    positives = [l for l in candidates if l > 0]
    if positives:
        return min(positives)
    # fallback
    return min(candidates, key=lambda v: abs(v))

# ---------------------------
# Surface normal from implicit F(x,y,z) = x - 0.0022*(y^2+z^2) + 6.22
# ∇F = <1, -0.0044*y, -0.0044*z>
# ---------------------------
def normal_at(x, y, z):
    return (1.0, -0.0044 * y, -0.0044 * z)

# ---------------------------
# Vector reflection
# R = D - 2 * (D·N)/(N·N) * N            (valid even if N is not unit length)
# ---------------------------
def reflect(D, N):
    dx, dy, dz = D
    nx, ny, nz = N
    dn = dx*nx + dy*ny + dz*nz
    nn = nx*nx + ny*ny + nz*nz
    # Guard against degenerate normal (shouldn't happen for this surface)
    if nn == 0:
        raise ZeroDivisionError("Degenerate normal encountered.")
    rx = dx - 2.0 * dn / nn * nx
    ry = dy - 2.0 * dn / nn * ny
    rz = dz - 2.0 * dn / nn * nz
    return (rx, ry, rz)

# ---------------------------
# Main I/O flow
# ---------------------------
def main():
    print("Enter Sun angles (degrees, SunCalc convention) and seed point near the facade.\n")
    azimuth = float(input("Azimuth (clockwise from North, degrees): ").strip())
    altitude = float(input("Altitude (elevation above horizon, degrees): ").strip())

    x0 = float(input("Seed X (units): ").strip())
    y0 = float(input("Seed Y (units): ").strip())
    z0 = float(input("Seed Z (units): ").strip())

    # Incident direction from angles
    dx, dy, dz = incident_dir_from_az_alt(azimuth, altitude)

    # Solve for intersection lambda
    lam = solve_lambda_for_intersection(x0, y0, z0, dx, dy, dz)

    # Intersection point on the facade
    xi, yi, zi = (x0 + lam*dx, y0 + lam*dy, z0 + lam*dz)

    # Surface normal (unnormalized is fine)
    Nx, Ny, Nz = normal_at(xi, yi, zi)

    # Reflected direction
    Rx, Ry, Rz = reflect((dx, dy, dz), (Nx, Ny, Nz))

    # Results
    print("\n--- Results ---")
    print(f"Intersection point on facade: x={xi:.6f}, y={yi:.6f}, z={zi:.6f}")
    print(f"Surface normal (unnormalized): Nx={Nx:.6f}, Ny={Ny:.6f}, Nz={Nz:.6f}")
    print(f"Reflected direction vector:   Rx={Rx:.6f}, Ry={Ry:.6f}, Rz={Rz:.6f}")

    # (Optional) Example of intersecting the reflected ray with the ground plane y=0:
    # If yi + t*Ry = 0 => t = -yi / Ry (provided Ry != 0)
    if abs(Ry) > 1e-12:
        t_ground = -yi / Ry
        xg = xi + t_ground * Rx
        yg = yi + t_ground * Ry   # ~0
        zg = zi + t_ground * Rz
        print(f"\nGround hit (y=0): x={xg:.6f}, y={yg:.6f}, z={zg:.6f}")
        print("Note: multiply by 4 to convert units -> meters (1 unit = 4 m).")
    else:
        print("\nReflected ray is parallel to the ground plane; no y=0 intersection.")

if __name__ == "__main__":
    main()

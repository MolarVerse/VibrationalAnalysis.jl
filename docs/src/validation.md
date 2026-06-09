# Validation

This page summarizes the current regression coverage behind the package examples and release surface.

## Reference H₂O Values

The fixture-backed H₂O calculation in the test suite currently reproduces:

- vibrational wavenumbers: `1492.535`, `3669.024`, `3784.327 cm^-1`
- infrared intensities: `150.521`, `86.433`, `164.086 km mol^-1`
- force constants: `1.417`, `8.312`, `9.171 mdyn Å^-1`
- reduced masses: `1.080`, `1.048`, `1.087 amu`

Those values come from the restart file, Hessian, and moldescriptor fixtures shipped in `test/data`.

## Normal-mode regression

Current tests compare computed normal modes against reference data for:

- H₂O
- MeOH

The H₂O comparison is checked to `1e-5` absolute tolerance on the vibrational subspace. The MeOH comparison currently allows a maximum absolute mode difference of `1.1e-4`.

## Linear and near-linear molecules

The package has explicit regression coverage for linear and near-linear systems, including CO₂ geometries. The current implementation drops rotational axes whose norm is below a relative tolerance of `1e-6`, which keeps exactly linear systems numerically stable while still allowing clearly bent geometries to use the nonlinear path.

In the current tests:

- exactly linear CO₂ is handled as a linear molecule
- `1e-6` out-of-plane perturbations are still treated on the linear path
- a `1e-5` out-of-plane perturbation is treated as bent

That behavior matters because real structures are often not perfectly planar or perfectly linear after optimization or finite-temperature sampling.

## Near-linear CO₂ Stability Example

The example below uses a synthetic identity Hessian. It is not a physical CO₂ spectrum; it is only a numerical stability check for a near-linear geometry.

```@example validation_co2
using LinearAlgebra
using VibrationalAnalysis

atom_masses = [15.9994, 12.0107, 15.9994]
atom_coords = [-1.16 0.0 0.0; 0.0 0.0 1e-5; 1.16 0.0 0.0]
hessian = Matrix{Float64}(I, 9, 9)

wavenumbers, force_constants, reduced_masses, modes =
    calculate(atom_masses, atom_coords, hessian, hessian_sign = :positive)

(
    all(isfinite, wavenumbers),
    length(wavenumbers),
    round.(wavenumbers[6:end], digits = 3),
)
```

## Hessian sign handling

Regression coverage also checks these cases:

- automatic sign selection through `hessian_sign = :auto`
- explicit positive and negative sign conventions
- numeric sign input with `1` or `-1`

If your upstream program uses the opposite Hessian sign convention, prefer setting `hessian_sign` explicitly in user-facing pipelines so the behavior is obvious.

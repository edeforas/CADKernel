# NURBS Surface-Surface Intersection Solver Specification

## Overview
Define a robust, production-grade surface-surface intersection solver for the CADKernel NURBS subsystem. The solver must reliably find intersection curves between two `NurbsSurface` objects, refine them to high geometric accuracy in 3D and their parametric UV domains, and provide outputs suitable for trimming and boolean operations.

## Goals
- Find all intersection curves (open and closed) between two surfaces.
- Produce ordered sample chains with parametric coordinates on both surfaces and ability to convert these samples into smooth NURBS curve representations.
- Robustly handle partial overlaps, nearly-tangent contacts, and degenerate configurations with clear diagnostics.
- Be integrable with existing `NurbsIntersectionResult` / `NurbsBoolean` code paths.

## API (proposal)
- compute_surface_intersection(const NurbsSurface& A, const NurbsSurface& B, NurbsIntersectionResult& out, const IntersectionOptions& opt)

IntersectionOptions (suggested fields):
- `double geomTol` — absolute 3D tolerance for equality (default: scale-based, e.g., diag*1e-8)
- `double paramTol` — UV tolerance for parameter refinement
- `double step` — predictor step length in parameter-space arclength
- `int maxNewtonIterations`
- `int maxChains`
- `int seedSamplingRes` — initial grid resolution for seed detection
- `bool fitToNurbs` — whether to fit the sampled curve to NURBS patches
- `int fitDegreeU` / `fitDegreeV` — degree suggestions for fitted curves
- logging/diagnostics flags

Return: `NurbsIntersectionResult` containing chains of `NurbsIntersectionSample` with 3D points and paired (u,v) coordinates for both surfaces, plus flags for partial overlap diagnostics.

## High-level algorithm
1. Pre-checks / culling:
   - Compute bounding boxes of both surfaces. If disjoint (with tolerance), return no intersections.
   - Quick rejection if distance between bboxes larger than threshold.

2. Seed detection (coarse):
   - Sample each surface on a coarse adaptive grid (option `seedSamplingRes`).
   - For each sample on A, find nearest sample on B. If distance <= seedThreshold, emit a seed (uA,vA,uB,vB).
   - Optionally perform surface pair bounding-box subdivision to focus refinement in overlapping regions.

3. Local refinement of seeds (Newton/least-squares):
   - For each seed, run a constrained Newton/Levenberg–Marquardt solver to refine parameters so that S1(u1,v1) == S2(u2,v2) within `geomTol`.
   - Math formulation:
     - Define F(u1,v1,u2,v2) = S1(u1,v1) - S2(u2,v2) (3 residuals).
     - Jacobian J is 3x4: [dS1/du1 dS1/dv1 -dS2/du2 -dS2/dv2].
     - The solution manifold is 1D (curve) — J is rank 3 in regular intersection cases and has a 1D nullspace corresponding to curve tangent.
   - Use Gauss-Newton on the least-squares objective ||F||^2 with damping (Levenberg–Marquardt) and robust line search. Stop when ||F|| < geomTol.
   - If Gauss-Newton is underdetermined/stiff, fallback to solving an augmented system (see continuation corrector below).

4. Curve tracing / predictor-corrector continuation:
   - From each refined seed, trace the intersection curve in both directions using predictor-corrector:
     - Compute tangent in parameter-space: find nullspace vector t of J (4-vector) such that J * t = 0. Normalize t by a chosen parameter metric (e.g., scale by typical parameter ranges).
     - Predictor: p_pred = p_current + step * t.
     - Corrector: solve for delta p using augmented system that enforces both F(p)=0 (3 eqns) and an orthogonality condition to the predictor direction, e.g. (p - p_pred)·t = 0 (1 eqn). This yields a 4x4 linear solve for the correction.
     - Use adaptive step control: reduce step on failed corrector, increase on successive successes, clamp to UV domain boundaries.
   - Collect samples along the traced path until stopping criteria: final correction convergence fails repeatedly, parameters reach domain boundary, or loop closes (distance to start < closureTol and parameter direction consistent).
   - Detect branching: if nullspace dimension changes or bifurcation indicators occur, mark branch points and optionally start new traces.

5. Post-processing chains:
   - Merge nearby chains (within mergeTol) and remove degenerate short chains.
   - Resample chains to approximately equal arclength spacing using `resample_curve_equal_arclength` (configurable target spacing).
   - Final refinement: run a few Newton iterations per sample to reproject onto both surfaces for increased accuracy.

6. Curve fitting to NURBS (optional):
   - Fit the 3D sample chain to a rational NURBS curve or series of Bezier patches (piecewise NURBS). Options:
     - Global parametric fit: choose degree and number of spans, solve weighted least-squares for control points and weights (nonlinear: alternate linear solves for control points with weight updates, or use constrained optimization).
     - Local Bezier extraction from surfaces: where possible, convert intersection segments that lie inside single surface spans to exact/Bezier representation by intersecting surface Bezier patches with the other surface.
   - Preserve parametric correspondence: compute / store associated UV params for each fitted control point or provide mapping functions from curve parameter to (u,v) on both surfaces.

7. Output and diagnostics:
   - `NurbsIntersectionResult` should include:
     - `curves`: list of `NurbsIntersectionCurve`, each with ordered `NurbsIntersectionSample` entries (point, uA,vA,uB,vB);
     - `closed` flag per curve;
     - `hasPartialOverlap` flag (true if algorithm detected large overlapping areas, not simple 1D intersections);
     - Debug info: number of seeds, failed seeds, branching points, convergence stats.

## Numerical details & linear algebra
- Jacobian evaluation: use `evaluate_derivatives(u,v,du,dv,duu,... )` to obtain partial derivatives `dS/du` and `dS/dv` for each surface.
- Nullspace computation: compute SVD (3x4) or use QR with column pivoting to find the nullspace vector corresponding to smallest singular value. For small matrices this is inexpensive and robust.
- Corrector linear system: solve 4x4 symmetric normal/EQ augmented system with stable pivoting (use small direct solver, avoid naive inversion).
- Damping strategies: Levenberg–Marquardt with adaptive damping factor; fallback to small-step secant steps when Jacobian is ill-conditioned.
- Parameterization scaling: scale parameter increments by typical knot domain sizes (knots normalized to [0,1] internally) and optionally by physical surface scale.

## Special cases & fallbacks
- Nearly tangent surfaces: intersection degenerates to a tangency curve; detection via small singular values and small normal angle; handle using adaptive sampling, mark as special (may require trimming or partial overlap handling).
- Partial overlaps / coincident patches: detect using bounding-box overlap and dense sampling; if large-area overlap present, set `hasPartialOverlap=true` and optionally return a set of trimmed surfaces approximating the overlap boundary (see existing `build_partial_overlap_*` stubs).
- Missed seeds: if coarse sampling misses narrow intersections, provide progressive subdivision of param domains (quadtree subdivision) to locate seeds.
- Robust domain handling: clamp parameter corrections to knot intervals; support periodic/closed directions.

## Diagnostics & testing
- Produce tests for:
  - Intersections of simple primitives (plane-plane, plane-sphere, cylinder-plane).
  - Known analytic intersections (sphere-sphere circle, cylinder-cylinder curves) to validate geometric error.
  - Trimming and boolean integration flow: fit -> trim -> boolean and verify topology.
- Logging levels: minimal (counts), verbose (per-seed convergence traces), debug (Jacobian/SVD values per step).

## Integration notes
- Replace current `src/NurbsIntersection.cpp` placeholder algorithm with the new implementation in incremental stages:
  1. Implement seed-detection + Newton local refinement and unit tests.
  2. Add predictor-corrector tracing; test closed/open curves.
  3. Add fitting to NURBS/Bezier patches and trimming conversion.
  4. Integrate with `NurbsBoolean.cpp` to use fitted loops when available and fall back to sampled-loops otherwise.
- Keep the old coarse-path as a fallback mode selectable via `IntersectionOptions` for quick diagnostics.

## Performance & parallelization
- Seed detection and local refinements are embarrassingly parallel — process seeds concurrently.
- Tracing is inherently sequential per chain but multiple chains can be traced in parallel.
- Reuse basis evaluations and derivative computations across nearby param evaluations where possible.

## Deliverables for the implementation phase
1. `NurbsIntersectionSpec` document (this file).
2. `IntersectionOptions` struct and updated `compute_surface_intersection` signature.
3. Unit tests under `tests/` for seed detection and local refinement.
4. Incremental implementation: coarse seeds + Newton, then predictor-corrector, then fitting and trimming conversion.


---

Authors: suggested implementer
Date: 2026-06-29


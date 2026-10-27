# DJL solver

MATLAB code for computing steady internal solitary waves in a continuously stratified fluid between rigid lids, by solving the DJL (Dubreil-Jacotin–Long) equation with Newton's method and numerical continuation. The baroclinic vorticity term can be taken with or without the Boussinesq approximation (`fe.Boussinesq`).

## Quick start

Run one of the examples in root. Each loads a converged solution and continues it:

| Script | Wave | Stratification | 
|---|---|---|
| `exampleConvexESW.m` | convex embedded solitary wave | double pycnocline |
| `exampleConvexGSW.m` |convex generalised solitary wave | double pycnocline |
| `exampleConcave.m` | concave solitary wave, not embedded | double pycnocline |
| `exampleMSL.m` | embedded solitary wave  | MSL (mixed surface layer)  |

As shipped, `co.vals` holds only the stored value, so each example re-converges its stored solution. To follow a branch, step `co.vals` (e.g. `co.vals = pa.H3 + (0:0.001:0.01)`). Solutions in this problem lie on complicated bifurcation branches, so new waves are best found by continuing from a known solution rather than from the WNL (weakly nonlinear) initial guess.

`stratified_solver.m` is the general template: it can start from the WNL initial guess (`fe.Loaddata=false`, settings in `core/buildInitialData.m`) or from a saved file.

## Structures

| Struct | Holds |
|---|---|
| `pa` | parameters, including target values of fixed quantities (`pa.amp`, `pa.Q`, `pa.K`) |
| `un` | unknowns and solution diagnostics (`un.psi`, `un.Q`, measured tail curvature `un.K`, energies) |
| `fe` | features / solver options, including the continuation settings `fe.co` |
| `me` | mesh and differentiation matrices |
| `de` | derivatives |
| `li` | lists of solutions along a branch |
| `ja` | Jacobian bookkeeping |

## Continuation options

Each script sets a struct `co` before running `core/solveBranch.m`:

```matlab
co.fix  = "K";          % quantities held at their pa values (one equation each)
co.free = [];           % pa fields solved for, in addition to un.Q
co.step = "z0";         % pa field looped over
co.vals = pa.z0+0.01;   % its values
```

Since `un.Q` is always free, `numel(co.fix) = 1 + numel(co.free)`.

| `co.fix` entry | Equation |
|---|---|
| `"amp"` | amplitude at `fe.Amploc` equals `pa.amp` |
| `"Q"` | `un.Q` equals `pa.Q` |
| `"area"`, `"area2"` | wave area (whole domain, or the column at `fe.Embedloc`) equals `pa.area` / `pa.area2` |
| `"K"` | tail curvature `d²ψ/ds²` at `fe.Embedloc` equals `pa.K` (e.g. `pa.K=0` for embedded solitary waves) |
| `"arclength_z0_area2"`, `"arclength_H3_Q"`, ... | pseudo-arclength steps in a pair of quantities, using `pa.arclength_ds` and `pa.arclength_factor` |

Examples:

| Aim | `co.fix` | `co.free` |
|---|---|---|
| branch at fixed amplitude | `"amp"` | `[]` |
| branch at fixed speed | `"Q"` | `[]` |
| zero tail curvature, speed and amplitude free | `"K"` | `[]` |
| zero tail curvature at fixed amplitude, solve for H3 | `["K","amp"]` | `"H3"` |

`fe.Solitary` additionally fixes the streamline heights in the far field, and is ony adviced to use when computing mode-1 ISWs. `co` is stored in saved files as `fe.co`, so the examples reuse it with `co = fe.co`.

## Stratification models

Set in `theory/stratificationModels.m` and evaluated in `core/evalStratification.m`: `single_pycnocline`, `MSL`.

## Saved files

`saveSolutionData` stores a branch; `loadSolutionData` reads it and applies `core/loadUpdateLegacy.m`, which maps old field names and options to current ones (including the old `fe.Fixparameter`, `fe.Embed` and `fe.Embedvary`, which become `fe.co`). `utils/fixLegacy(folder)` converts a folder of old files, and `utils/makeExample(infile,outfile)` makes a one-solution starting file.

## Credits

`core/gauss.m` is from Nick Trefethen. `core/djles_compute_apedens.m` is adapted from the DJLES package (Dunphy, Subich and Stastna 2011).

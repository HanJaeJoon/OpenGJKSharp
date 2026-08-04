# OpenGJKSharp round 2 dogfooding findings

Target: local prerelease package `OpenGJKSharp 0.6.0-dogfood.2/.3` (built from `main` after the
post-round-1 work), consumed from a local NuGet feed by `Dogfood2/` (net8.0 console app).
Library source was not opened while writing the scenarios - only the packaged README and XML docs.
Round 1 findings: `../OpenGJKSharp-dogfood/findings.md`.

Purpose: the API surface changed a lot since 0.5.0 (`GjkWorkspace`, `TransformedShape`,
double-precision overloads, `TryGetPenetration`, `GkPolytope`/`GkSimplex` removed), so the
"1.0 is reachable" verdict from round 1 no longer covers what consumers actually see.

## Scenarios

- [x] 1. Character controller frame loop (workspace + TransformedShape + TryGetPenetration)
- [x] 2. Moving bodies (TransformedShape vs rebuilding a hull, cross-check against oriented Box)
- [x] 3. Depenetration (deep overlap recovery, sign convention, 20-sphere pile)
- [x] 4. Large coordinates through the double-precision span overloads
- [x] 5. Migration from 0.5.0
- [x] 6. Public surface audit
- [x] 7. Threading (workspace per thread, and what misuse looks like)
- [x] 8. Edge cases and the two contracts that changed

---

# Verdict: 1.0 is reachable. Two bugs were found and fixed during this round; three gaps remain, none blocking.

The round-1 P0/P1/P2 work holds up under real use: the physics loop is allocation-free, the
pose wrapper is exact and faster than rebuilding, the double path recovers precision the
`Vector3` path cannot, and every documented contract checked out. The one serious defect found
here predates all of it and would have shipped in 1.0.

## Bugs found (both fixed in this round)

### 1. [BUG, was critical] No penetration depth when the two bodies share an interior point

Crossed capsules, axes through the origin, half-length 2 - analytic depth is `2 * radius`:

| radius | 0.5.0 | 0.6.0 before fix | 0.6.0 after fix |
|---|---|---|---|
| 0.25-0.75 | exact | exact | exact |
| 1.00 | 0.00000 (-100%) | 3.02523 (+51%) | 2.00000 (0.00%) |
| 1.50 | 0.00000 (-100%) | 3.97316 (+32%) | 3.00000 (0.00%) |

Root cause: GJK seeded simplex vertex 0 with the difference of the two `GetCenter()` results.
That is an interior point of the Minkowski difference - exactly the origin when both bodies
share a center - and it stayed in the simplex. EPA builds its polytope from those vertices and
needs boundary points; an interior vertex puts a face through the origin, so the closest-face
search reports distance 0 forever. In 0.5.0 that surfaced as depth 0 (a solver would never push
the bodies apart); the round-1 EPA change turned it into a large overestimate instead.

Affects any pair sharing an interior point: concentric bodies, one body inside another, and
long bodies crossing at their common center. Fixed by using the center difference only as the
initial search direction and taking the first support point along it, which is what the
algorithm calls for (`fix: seed the GJK simplex with a support point`). Side effect: a separated
support-mapped query dropped from about 610 ns to 470 ns.

Note this bug hid behind round 1's own measurements: round 1 tested concentric *spheres* (which
happen to recover) but never crossed bodies of comparable thickness.

### 2. [polish] `TryGetPenetration` returned negative zero on exact touch

`depth` came back as `-0.0` for bodies that touch exactly, which formats as `-0` and flips the
sign of anything divided by it. Fixed.

## What was verified as working

**Physics loop (scenario 1).** 600-frame capsule character controller, ground + wall, up to 4
resolution passes per frame: resting y = 0.7000 and blocked x = 2.7000 (both analytic), 0 NaN
frames, 0 residual overlap after resolution, and **0-4 B per frame** total allocation (run to
run variation, i.e. nothing but JIT noise). Round 1's "every frame allocates a few KB" friction
is gone; 109 us/frame for ~10 queries, dominated by the EPA path on contact frames.

**TransformedShape (scenario 2).** Distance difference against the built-in oriented `Box`
across 200 random poses: **0.00E+000** - bit-identical, not just close. Moving a 64-point hull
over 2000 frames: rebuilding the world-space array takes 27.4 ms and 4417 B/frame; the wrapper
takes 11.6 ms and **2 B/frame** (2.4x faster, ~2000x less garbage) with an identical checksum.

**Depenetration (scenario 3).** Every tested overlap (0.5 to 1.99 of a 2.0 diameter) is resolved
in a **single step** - the first correction equals the analytic depth to 3 decimals. A 20-sphere
pile settles after 600 frames with worst remaining sphere-sphere overlap 5.7e-5 and
sphere-ground overlap 8.0e-5, and the lowest sphere rests at y = 0.4999 (0.5 exact).
`TryGetPenetration` and `ComputeCollisionInformation` agree exactly.

**Double precision (scenario 4).** A 1e-3 gap survives at any coordinate the `Vector3` path
destroys:

| origin | double path | Vector3 path |
|---|---|---|
| 0 | 1.000000E-003 | 1.000000E-003 |
| 1e3 | 1.000000E-003 | 9.765625E-004 |
| 1e6 | 1.000000E-003 | 0 |
| 1e9 | 1.000047E-003 | 0 |

Witness points and contact normal come back through caller `Span<double>` buffers, 0 B/query
with a workspace (2808 B without), and results match the obsolete `CsFunction` bit for bit.

**Migration (scenario 5).** `CsFunction`'s `[Obsolete]` message names the replacement including
the coordinate layout change, and all three entry points return the identical distance.

**Surface (scenario 6).** 8 public types, no internals leaked, `GjkWorkspace` is an opaque handle
with a single constructor and no public members, `TransformedShape` is the only mutable shape
(`SetPose`, no public setters). Workspace overloads exist for every query family that exists.

**Threading (scenario 7).** `ThreadLocal<GjkWorkspace>` over 4000 pairs: 0 result mismatches,
24.0 ms -> 15.0 ms on 12 cores.

**Contracts (scenario 8).** "Colliding bodies report exactly 0" holds literally (`== 0`) for
700/700 overlapping pairs spanning scales 1e-3 to 1e6 - round 1 measured 2.2e-16 here. EPA depth
error is at most 0.24% for spheres and exactly 0 for boxes and capsules (round 1: -5.62%).
Degenerate-shape policy matches the README exactly. NaN input rejected.

## Open findings (not fixed)

### A. [gap] `ComputeCollisionInformation` does not exist for `Vector2[]`

2D polygons can be tested for collision and distance, but there is no way to get penetration
depth or contact normal for them - a 2D consumer has to lift everything to 3D. This predates
0.5.0 but is more visible now that every other family covers all four input kinds.
`TryGetPenetration` has the same 2D gap.

### B. [gap] `TryGetPenetration` is missing for the double-precision path

Deliberate (an `out` pair reads badly next to `Span<double>` output buffers), but it means the
convenience API is not available to large-coordinate consumers, who arguably need protection
from the sign convention just as much.

### C. [friction] Sharing a workspace across threads fails confusingly

The docs say one workspace per thread. Breaking that rule produced `IndexOutOfRangeException`
and `ArgumentOutOfRangeException` from inside the library plus **3618 of 4000 silently wrong
results**. A single `Interlocked` ownership flag on the workspace would turn that into an
`InvalidOperationException` naming the actual mistake, for roughly the cost of one atomic per
query.

### D. [friction] Nothing points a 0.5.0 consumer from the compiler error to `GjkWorkspace`

Code using `GkPolytope`/`GkSimplex` now fails with CS0122/CS0246. The type is gone, so no
`[Obsolete]` message can carry the migration hint; only the CHANGELOG has it. A short
"Migrating from 0.5.0" section in the README would be where a stuck consumer actually looks.

### E. [tooling, contributor-facing] `dotnet pack` alone silently packs stale output

Because the csproj sets `GeneratePackageOnBuild=true`, running `dotnet pack` on its own does not
rebuild: it packaged a DLL from an earlier build (missing `TryGetPenetration` entirely), and
after cleaning `bin` it failed with NU5026 instead of building. The CI workflow is safe - it runs
`dotnet build` then `dotnet pack --no-build` - but a contributor packing locally can produce a
package that does not match the source. Worth a line in CLAUDE.md.

## Assets

`Dogfood2/` holds the scenarios (`dotnet run -- 1..8`), `Compare050/` and `CompareNoFix/` are the
same probe against published 0.5.0 and against the commit before the round-1 EPA change, which is
how bug 1 was attributed. `feed/` holds the local prerelease packages.

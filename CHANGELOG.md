# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `TryGetPenetration` overloads that also return the contact points (`pointA`, `pointB`): the
  deepest point of each body along the contact normal, each on its own body's surface and
  separated along the normal by the depth (`dot(pointA - pointB, contactNormal) == depth`).
  When the bodies are apart they are the closest points, as `ComputeMinimumDistance` reports.
  Previously the only witness points came from `ComputeMinimumDistance`, which collapse onto an
  interior point once the bodies overlap. Available for `ISupportMappable` shapes, `Vector3[]`
  and `Vector2[]` vertex arrays, and `ReadOnlySpan<double>` coordinates, each with a
  `GjkWorkspace` variant; flat and collinear pairs are covered too.
- `TryCastRay`, a ray cast against an `ISupportMappable` shape: it returns whether the segment
  from an origin along a translation hits the shape, the fraction of the translation at the
  first hit, and the outward surface normal there. A ray that starts inside the shape reports a
  hit at fraction 0 with a zero normal. Available with a `GjkWorkspace` variant.
- `TryCastShape`, a linear shape cast for continuous collision detection: it sweeps one
  `ISupportMappable` shape along a translation against another and returns the fraction of the
  translation at first contact, the contact normal (pointing from the moving shape toward the
  static one), and the contact point. Shapes that already overlap report a hit at fraction 0
  with a zero normal. Available with a `GjkWorkspace` variant.
  Both casts use the GJK ray cast (van den Bergen 2004) and never report a fraction past the
  true one. Shapes count as touching within about 1e-6 of the shape size, or a few float ulps
  of the coordinates far from the origin; a cast that starts within that distance reports
  fraction 0 with a zero normal. Normals are exact on flat faces; on curved surfaces they are
  accurate to about 1e-3 near the origin (up to 3e-3 on a capsule side), and to about 1e-2 for
  a unit sphere at 1e3, where float resolution limits them.

### Changed

- `HasCollision` stops GJK as soon as a lower bound on the distance clears the precision (or
  an upper bound falls within it) instead of running to convergence, and returns the same
  answer as `ComputeMinimumDistance(a, b) <= precision`. On separated pairs with a reused
  workspace it went from 888 to 480 ns for 32-vertex clouds and from 1509 to 159 ns for a
  sphere and a box (about 1.9x and 9.5x faster).

## [0.6.0] - 2026-08-03

### Breaking

Read this section before upgrading from 0.5.0. The README has a
[Migrating from 0.5.0](README.md#migrating-from-050) section with before/after code for each item.

- `GkPolytope` and `GkSimplex` are now `internal`, along with the `OpenGjk` overloads that took
  them. They exposed internal algorithm state - an undocumented `double[][]` coordinate layout and
  settable counters that could break internal invariants - so code using them no longer compiles
  (CS0122/CS0246). Vertex-array callers should use the `Vector3[]`/`Vector2[]` overloads; callers
  who used them to avoid per-query allocation should use `GjkWorkspace`. This is the only
  compile-time break; an API diff against 0.5.0 shows no other removal or signature change.
- Four results change value without any compile error, all of them fixes to answers that were
  wrong before. Check any code that depended on the old values (each has its own entry below):
  - colliding bodies now report a minimum distance of exactly 0, where a small residual was left,
  - the penetration depth of curved shapes is larger and no longer an underestimate,
  - bodies whose interior points coincide report a depth instead of 0,
  - flat pairs report a depth measured in the plane they share, instead of 0.
- `netstandard2.0` consumers gain a `System.Memory` package dependency, needed by the new
  `ReadOnlySpan<double>` overloads. This affects .NET Framework only; the other targets are
  unchanged.

### Added

- `TryGetPenetration`, an unsigned convenience overload for depenetration: it returns whether
  the bodies collide, a positive depth, and the contact normal, so callers no longer have to
  interpret the sign convention of `ComputeCollisionInformation` (negative means colliding).
  Available for `ISupportMappable` shapes, vertex arrays, and `ReadOnlySpan<double>` coordinates,
  each with a `GjkWorkspace` variant.
- Double-precision query overloads taking coordinates as a `ReadOnlySpan<double>` of
  consecutive (x, y, z) triples (`HasCollision`, `ComputeMinimumDistance` with and without
  witness points, `ComputeCollisionInformation`, each with a `GjkWorkspace` variant).
  The obsolete `CsFunction` used to be the only double-precision entry point; these
  overloads replace it with the full result set. Adds a `System.Memory` dependency
  on `netstandard2.0` only.
- Penetration depth for flat pairs, computed in the plane the two bodies share. EPA needs a
  tetrahedron, which a flat Minkowski difference has no room for, so coplanar bodies reported a
  depth of 0 - two unit squares overlapping by 0.25 in the z = 0 plane came back as if they were
  merely touching. Coplanar bodies now run a two-dimensional expanding-polygon pass instead, and
  bodies on a common line the one-dimensional equivalent. This covers `Vector2[]` polygons, and
  also `Vector3[]` vertices and flat `ConvexHull` shapes that happen to share a plane.
- `ComputeCollisionInformation` and `TryGetPenetration` overloads for `Vector2[]` polygons, each
  with a `GjkWorkspace` variant, returning the depth and a `Vector2` contact normal in the plane.
  Previously the 2D API stopped at `HasCollision` and `ComputeMinimumDistance`.
- `TransformedShape`: places any `ISupportMappable` at a pose (rotation + translation),
  so a shape can be defined once in local space and moved with `SetPose` instead of being
  rebuilt every frame. Previously only `Box` accepted an orientation, and moving a
  `ConvexHull` meant rebuilding its vertex array.

- `GjkWorkspace`: reusable scratch memory for queries. Every `HasCollision`,
  `ComputeMinimumDistance`, and `ComputeCollisionInformation` overload (vertex arrays,
  2D polygons, and `ISupportMappable` shapes) now has a counterpart taking a workspace.
  Reusing one makes repeated queries allocation-free once its buffers have grown -
  previously a support-mapped query allocated about 4 KB, and a colliding query
  about 150 KB for the EPA pass.

- A `GjkWorkspace` now detects being used by two queries at once and throws
  `InvalidOperationException` naming the mistake, instead of letting the queries overwrite each
  other's buffers. Sharing one workspace across threads used to either throw
  `IndexOutOfRangeException` from deep inside the library or, more often, return silently wrong
  results (3618 of 4000 queries in a parallel probe). Re-entering a query with the workspace it
  was called from - for example querying from a custom `GetSupportPoint` - is caught by the same
  check. A failing query releases its workspace, so it stays usable afterwards. Allocation-free
  reuse is unaffected (still 0 B/query) and the added cost is below measurement noise.

### Fixed

- Penetration depth is now reported for bodies whose interior points coincide - concentric
  bodies, one body containing another, or two long bodies crossing at their common center.
  GJK seeded its simplex with the difference of the two `GetCenter` results, an interior point
  of the Minkowski difference (the origin itself when the centers match); EPA inherited it as
  a polytope vertex, which puts a face through the origin and stalls the closest-face search
  at distance zero. The simplex now starts from a support point along that direction instead,
  as the algorithm intends. Crossed capsules of radius 1 reported a depth of 0 (0.5.0) and are
  now exact; separated queries also got slightly faster.
- EPA no longer underestimates the penetration depth of deeply overlapping curved shapes.
  A polytope face is a lower bound on the depth and can never resolve a curved surface, which
  cost up to 5.6% (two concentric unit spheres reported 1.888 instead of 2.0); the depth is
  now taken from the tightest support plane found during the search, which is an upper bound
  and exact when the contact normal is. Worst-case error over the sphere overlap sweep drops
  from 5.62% to 0.24%, polyhedra stay exact, and the cost is unchanged.

### Changed

- A flat pair of bodies now reports the penetration depth measured in the plane they share,
  where it used to report 0. Measured in three dimensions 0 is the honest answer - a flat pair
  separates when either body leaves the plane by any amount - but it is useless to a caller
  working with flat geometry and indistinguishable from "not colliding". The consequence is a
  discontinuity in thickness: a pair 1e-6 thick is genuinely three-dimensional and still reports
  its out-of-plane depth of 1e-6, while the same pair at thickness 0 reports the in-plane depth.
- Colliding bodies now report a minimum distance of exactly 0. The iteration used to leave a
  residual (about 1e-16 at unit scale, larger for bigger bodies), so `distance == 0` missed
  collisions; the residual is now measured against the size of the bodies and flushed to zero.
- Documented the policy for degenerate shapes: a thickness parameter (`Sphere`/`Capsule`
  radius, `Box` half extents) must be positive, while geometry given as explicit points may be
  degenerate - a `ConvexHull` can be a single point, a segment, or a flat polygon. A `Capsule`
  with coincident endpoints is a sphere, not an exception to the rule.
- Documented the result contract that dogfooding found missing: the direction passed to
  `GetSupportPoint` is never zero but its length is unspecified; `GetCenter` must return an
  interior point and only affects convergence speed; witness points are meaningless while the
  bodies overlap (both collapse to a common interior point); the contact normal points from
  the first body toward the second and is arbitrary for concentric bodies.
- GJK/EPA internals no longer allocate temporaries per iteration (face normals, horizon
  edges, witness computation and the distance sub-algorithm now work on reused buffers).

### Deprecated

- The `[Obsolete]` message on `OpenGjk.CsFunction` no longer points at the now-internal
  `GkPolytope` overload.

## [0.5.0] - 2026-08-01

### Added

- `ISupportMappable` support-mapping interface: GJK/EPA queries on implicit convex shapes without tessellation
- Built-in shapes in `OpenGJKSharp.Shapes`: `Sphere`, `Capsule`, `Box` (oriented), `ConvexHull`
- `OpenGjk` overloads for `ISupportMappable`: `HasCollision`, `ComputeMinimumDistance` (with and without witness points), and `ComputeCollisionInformation`

### Changed

- Internal GJK/EPA loops now run on a support-function abstraction; the vertex-array and `GkPolytope` APIs are unchanged and keep their exact behavior

## [0.4.0] - 2026-07-23

### Added

- `netstandard2.1` target: the package can now be used from .NET Core 3.0+, Unity 2021.2+, and Mono in addition to .NET 8
- `netstandard2.0` target: the package can now be used from .NET Framework 4.6.2+ (via `System.Numerics.Vectors`)
- `net10.0` target for the current LTS runtime; tests now run on both .NET 8 and .NET 10
- SourceLink and symbol package (`.snupkg`) for source-level debugging
- Package icon and richer NuGet metadata (repository type, copyright, release notes link)
- License section in README clarifying GPL-3.0 copyleft implications

### Changed

- NuGet publishing now happens only on `v*` release tag push (previously on every push to `main`), with a tag-vs-csproj version check

### Deprecated

- `OpenGjk.CsFunction` (legacy openGJK C interface, already `[Obsolete]`) will be removed in the next major version - use `ComputeMinimumDistance` instead

## [0.3.0] - 2026-07-10

### Changed

- **Breaking:** rename the entry-point class `OpenGJKSharp` to `OpenGjk`. The old name collided with the namespace, so consumers could not reference it without a namespace alias
- Replace `Console.WriteLine` noise in the library with `Debug.WriteLine`
- Ship IntelliSense XML documentation with the package

### Added

- Input validation with clear exceptions: null/empty vertex arrays, non-finite (NaN/Infinity) coordinates, negative precision, and `GkPolytope` shape mismatches
- `ComputeMinimumDistance` overloads for `Vector3[]`/`Vector2[]` that return the distance and the closest point on each body (witness points)

### Fixed

- EPA reporting zero penetration depth when both bodies share their first vertex
- EPA mis-orienting face normals when a polytope face passes through the origin, which stalled EPA at ~0 penetration depth for sphere-like hulls
- Missing returns in degenerate witness branches that could produce NaN results

## [0.2.0] - 2026-07-09

### Changed

- Apply modern C# conventions to the model types (public fields to properties)
- Fix the `_presicion` typo to `_precision`

## [0.1.0] - 2026-07-09

### Added

- EPA (Expanding Polytope Algorithm) ported from upstream openGJK: `ComputeCollisionInformation` returns the penetration depth and contact normal
- `Vector3` convenience overload and a `ComputeMinimumDistance` overload accepting a `GkSimplex`

### Changed

- NuGet publishing switched from an API key to Trusted Publishing (OIDC)

## [0.0.5] - 2025-03-08

### Added

- Initial release: GJK minimum distance and `HasCollision` for convex polyhedra (`Vector3[]`) and flat polygons (`Vector2[]`)

[Unreleased]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/HanJaeJoon/OpenGJKSharp/compare/v0.0.5...v0.1.0
[0.0.5]: https://github.com/HanJaeJoon/OpenGJKSharp/releases/tag/v0.0.5

# Support Function API (ISupportMappable) - Design

Date: 2026-08-01
Target version: 0.5.0

## Motivation

The current public API only accepts vertex clouds (`Vector3[]` / `Vector2[]`).
GJK is fundamentally a support-function-based algorithm, so shapes like spheres
and capsules can be handled exactly - without tessellation - by exposing a
support mapping interface. This makes the library usable for physics-style
workloads (sphere/capsule/box queries) and lets users plug in custom convex
shapes.

## Public API

### Interface

```csharp
namespace OpenGJKSharp;

public interface ISupportMappable
{
    /// <summary>Returns the furthest point of the shape in the given direction (world space).
    /// The direction is not required to be normalized and is never zero.</summary>
    Vector3 GetSupportPoint(Vector3 direction);

    /// <summary>Returns an interior point of the shape (typically the center).
    /// Used as the initial guess for GJK and as the seed for EPA.</summary>
    Vector3 GetCenter();
}
```

- No default interface members (netstandard2.0 compatibility).
- Contract notes documented on the interface: implementations must describe a
  convex shape; `GetSupportPoint` must be exact for the shape (GJK/EPA accuracy
  depends on it).

### Built-in shapes

Namespace `OpenGJKSharp.Shapes`, all immutable classes taking world-space
parameters. No generic transform wrapper in v1 (YAGNI); only `Box` carries an
orientation.

| Shape | Constructor | Validation (ArgumentException / ArgumentOutOfRangeException) |
|---|---|---|
| `Sphere` | `(Vector3 center, float radius)` | `radius > 0`, finite inputs |
| `Capsule` | `(Vector3 p0, Vector3 p1, float radius)` | `radius > 0`, finite inputs (`p0 == p1` is allowed and degenerates to a sphere) |
| `Box` | `(Vector3 center, Vector3 halfExtents, Quaternion orientation)` + convenience `(center, halfExtents)` = identity | all half extents `> 0`, finite inputs, non-zero-length quaternion (normalized internally) |
| `ConvexHull` | `(Vector3[] points)` | non-null, at least 1 point, finite inputs |

Support mappings:

- `Sphere`: `center + radius * normalize(direction)`
- `Capsule`: pick the endpoint with the larger dot product, then add
  `radius * normalize(direction)`
- `Box`: rotate direction into local space, pick `sign(dir) * halfExtents`
  per axis, rotate back and translate
- `ConvexHull`: linear vertex scan (same as the existing internal support)

### OpenGjk overloads

Mirroring the existing `Vector3[]` surface:

```csharp
public static bool HasCollision(ISupportMappable a, ISupportMappable b, double precision = 1e-6); // same default as the array overloads
public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b);
public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b, out Vector3 closestA, out Vector3 closestB);
public static double ComputeCollisionInformation(ISupportMappable a, ISupportMappable b, out Vector3 contactNormal);
```

Existing `Vector3[]` / `Vector2[]` overloads keep their exact signatures and
behavior. Internally the `Vector3[]` overloads are free to keep their current
`GkPolytope` fast path; result parity with the `ConvexHull` route is enforced
by regression tests.

## Internal design

### Support abstraction

Today the GJK loop calls `Support(GkPolytope body, double[] v)` (linear vertex
scan tracking `SIdx`), and EPA calls `SupportEpa(...)` which relies on vertex
indices for duplicate-support detection.

Change:

- Introduce an internal support contract (e.g. `IGkBody` or an internal
  abstract wrapper) with two implementations:
  - vertex-cloud body: current `GkPolytope` behavior, keeps `SIdx`;
  - support-mapped body: adapts a user `ISupportMappable`, converting the
    double-precision search direction to `Vector3` and the returned point back
    to doubles.
- The GJK distance loop and the EPA loop call the abstraction instead of the
  concrete vertex scan.

### EPA duplicate-support detection

EPA currently detects "no new support point" by comparing vertex indices.
Implicit shapes have no indices, so generalize to a coordinate comparison:
a new support point is a duplicate when it is within an epsilon (relative to
shape scale, reusing the existing tolerance constants) of an existing polytope
vertex. The vertex-cloud path may keep the index shortcut as an optimization,
but the coordinate check is the shared correctness mechanism.

### Witness points

Witness points are reconstructed from per-body support point coordinates
stored on the simplex, so the mechanism carries over to support-mapped bodies.
Implementation must verify that per-body support coordinates (not only
Minkowski-difference points) are tracked on every code path that feeds
`Witnesses`.

### Numerics

- Internal math stays double-based, as today.
- The public interface uses `Vector3` (float), so custom-shape support points
  are float-quantized. The existing array API has the same input quantization;
  document this limitation in the README instead of introducing a double-based
  public interface.

### Termination on curved shapes

For polytopes GJK terminates in finitely many steps; for curved shapes
(sphere/capsule) it converges asymptotically. Existing iteration caps and
convergence tolerances apply unchanged; on hitting the cap the current best
result is returned (same policy as today). No new configuration knobs in v1.

## Error handling

- Shape constructors validate their parameters (table above) so that query
  methods never see invalid shapes.
- `OpenGjk` overloads validate `precision >= 0` (matching existing behavior)
  and null-check both shapes (`ArgumentNullException`).

## Testing

xunit, `tests/OpenGJKSharp.Tests`:

1. **Analytic cases**
   - Sphere-sphere: separated distance `= |c1-c2| - r1 - r2`; overlapping
     penetration depth `= r1 + r2 - |c1-c2|`, contact normal along the center
     line; witness points on the surface along that line.
   - Capsule-capsule (parallel and skew segments) against hand-computed
     segment-segment distances; capsule with `p0 == p1` behaves as a sphere.
   - Box: rotated box vs the equivalent 8-point hull - results match within
     tolerance.
2. **Parity / regression**
   - `ConvexHull` route vs existing `Vector3[]` route: identical results on
     the existing test fixtures (cubes, random hulls).
   - Full existing test suite stays green (no behavior change for array APIs).
3. **Validation**
   - Constructor argument checks; null shape arguments; negative precision.

## Documentation and release

- README: new "Shapes and custom support functions" section with a
  sphere-vs-hull example and a custom `ISupportMappable` example.
- CHANGELOG: Added - `ISupportMappable` + built-in shapes; note on EPA
  duplicate detection generalization.
- Version bump to 0.5.0 (minor: additive API).

## Out of scope (deliberate)

- Generic transform wrapper shape (add later if requested).
- Struct-generic devirtualization / BenchmarkDotNet performance work.
- 2D (`Vector2`) support-mapped shapes.
- Minkowski-sum composite shapes (e.g. rounded boxes via shape + radius).

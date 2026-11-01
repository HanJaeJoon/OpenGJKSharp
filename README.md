## OpenGJKSharp

[![NuGet](https://img.shields.io/nuget/v/OpenGJKSharp)](https://www.nuget.org/packages/OpenGJKSharp)
[![NuGet Downloads](https://img.shields.io/nuget/dt/OpenGJKSharp)](https://www.nuget.org/packages/OpenGJKSharp)
[![Build](https://github.com/HanJaeJoon/OpenGJKSharp/actions/workflows/nuget-publish.yml/badge.svg)](https://github.com/HanJaeJoon/OpenGJKSharp/actions/workflows/nuget-publish.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

OpenGJKSharp is a C# native implementation of the GJK (Gilbert-Johnson-Keerthi) algorithm, designed for efficient collision detection between convex polyhedra in 3D space.
This project is inspired by the original [openGJK](https://github.com/MattiaMontanari/openGJK) (written in C) and reimagined for the .NET ecosystem.

Targets `netstandard2.0`, `netstandard2.1`, `net8.0`, and `net10.0`, so it can be used from .NET Framework 4.6.2+, .NET Core 2.0+, .NET 5+, Unity 2021.2+, and Mono.

## Installation
You can install `OpenGJKSharp` via NuGet:

``` bash
dotnet add package OpenGJKSharp
```

Or using the Package Manager:

``` bash
Install-Package OpenGJKSharp
```

> **Note:** Starting with 0.3.0, the entry-point class is named `OpenGjk`
> (previously `OpenGJKSharp`, which collided with the namespace and could not be
> referenced without a namespace alias).

## Usage

Here is an example of detecting a collision between two overlapping cubes:

``` c#
using System.Numerics;
using OpenGJKSharp;

// Cube 1
var a = new Vector3[]
{
    new(0, 0, 0),
    new(1, 0, 0),
    new(0, 1, 0),
    new(1, 1, 0),
    new(0, 0, 1),
    new(1, 0, 1),
    new(0, 1, 1),
    new(1, 1, 1),
};

// Cube 2
var b = new Vector3[]
{
    new(0.5f, 0.5f, 0),
    new(1.5f, 0.5f, 0),
    new(0.5f, 1.5f, 0),
    new(1.5f, 1.5f, 0),
    new(0.5f, 0.5f, 1),
    new(1.5f, 0.5f, 1),
    new(0.5f, 1.5f, 1),
    new(1.5f, 1.5f, 1),
};

bool hasCollision = OpenGjk.HasCollision(a, b);

Console.WriteLine($"Collision detected: {hasCollision}"); // Outputs: true
```

The following image illustrates the collision between the two cubes:

![Collision Example](https://raw.githubusercontent.com/HanJaeJoon/OpenGJKSharp/refs/heads/main/example.png)

### Penetration depth and contact normal (EPA)

`ComputeCollisionInformation` runs the GJK algorithm and, when a collision is detected,
the EPA (Expanding Polytope Algorithm) to compute the penetration depth and contact normal:

``` c#
double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

if (distance < 0)
{
    // Colliding: -distance is the penetration depth, and contactNormal is a unit vector
    // pointing from the first body toward the second - move b along it (or a against it)
    // by the depth to separate them.
    Console.WriteLine($"Penetration depth: {-distance}, normal: {contactNormal}");
}
else
{
    // Separated: distance is the minimum distance between the two bodies
    Console.WriteLine($"Distance: {distance}");
}
```

For depenetration, `TryGetPenetration` reports the same result without a sign convention to
get backwards - the depth is always positive and the normal always points from the first body
toward the second:

``` c#
if (OpenGjk.TryGetPenetration(body, ground, out double depth, out Vector3 normal))
{
    position -= normal * (float)depth;   // move the first body out of the second
}
```

Add `pointA`/`pointB` to also get the contact points - the deepest point of each body along
the normal, each on its own body's surface, so that `dot(pointA - pointB, normal) == depth`:

``` c#
if (OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal,
        out Vector3 pointA, out Vector3 pointB))
{
    // pointA is on a's surface, pointB on b's, depth apart along the normal.
    Vector3 contact = (pointA + pointB) / 2;
}
```

Where flat faces meet, the two points face each other inside the contact region; on a curved
surface the point is the body's support point along the normal. When the bodies are apart
(`false`), they are the closest points, the same as `ComputeMinimumDistance` reports. Every
`TryGetPenetration` overload has this variant (`Vector2[]` polygons return `Vector2` points,
`ReadOnlySpan<double>` coordinates write them to 3-element spans).

For exactly concentric bodies every direction is a valid minimum translation, so the
returned normal direction is arbitrary (still unit length).

The depth is exact for polyhedra. Curved shapes (spheres, capsules, and custom support
mappings) cannot be resolved exactly by a finite polytope, so their depth is reported as a
slight overestimate - within about 0.3% for spheres, even when one body is fully inside the
other - rather than an underestimate, so pushing the bodies apart by it always separates them.

### Minimum distance and closest points

`ComputeMinimumDistance` returns the minimum distance between two bodies
(0 when they collide) and can also report the closest point on each body:

``` c#
double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

Console.WriteLine($"Distance: {distance}");
Console.WriteLine($"Closest points: {closestA} <-> {closestB}");
```

The closest points are only meaningful while the bodies are separated, where they lie on
the respective surfaces. Once the bodies overlap there is no closest pair: the distance is
0 and both witness points collapse onto a common point inside the intersection, which is
not on either surface. Use `ComputeCollisionInformation` for that case, or `TryGetPenetration`
with `pointA`/`pointB` for contact points on the surfaces.

### 2D polygons and other flat geometry

Every query also accepts `Vector2[]` polygons, treated as flat shapes on the z = 0 plane, and
reports the depth and normal in that plane:

``` c#
Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
Vector2[] b = [new(0.75f, 0), new(1.75f, 0), new(1.75f, 1), new(0.75f, 1)];

bool hit = OpenGjk.HasCollision(a, b);
double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector2 onA, out Vector2 onB);
double info = OpenGjk.ComputeCollisionInformation(a, b, out Vector2 contactNormal);

if (OpenGjk.TryGetPenetration(a, b, out double depth, out Vector2 normal))
{
    Console.WriteLine($"{depth} along {normal}");   // 0.25 along <1, 0>
}
```

The same holds for a coplanar pair given in 3D - `Vector3[]` vertices or flat `ConvexHull`
shapes that share a plane - and for two bodies on a common line: the depth is the minimum
translation within the plane (or along the line).

Note that this is a deliberate reading of a degenerate case. In three dimensions the honest
depth of a flat pair is 0, because moving either body out of the shared plane by any amount at
all separates them. Reporting that is useless to a caller working with flat geometry, and
indistinguishable from "not colliding", so the in-plane answer is reported instead. The
consequence is a discontinuity in the thickness of the bodies: a pair 1e-6 thick is genuinely
three-dimensional and gets its out-of-plane depth of 1e-6, while the same pair at thickness 0
gets the in-plane depth.

With `Vector2[]` input the pair is flat by construction. With single-precision `Vector3`
vertices it has to be flat *as quantized*: rotating a coplanar pair by an arbitrary angle moves
each vertex to the nearest `float` and generally leaves the two bodies in slightly different
planes, which is a three-dimensional configuration where the two flat bodies may not intersect
at all. Use the `ReadOnlySpan<double>` overloads, or keep flat geometry on a plane whose
coordinates are exact, when this matters.

### Double-precision input

`Vector3` is single precision, so vertex coordinates are quantized on input - at a
coordinate of 1e6 the spacing of `float` is about 6e-2, and smaller gaps are lost.
Internal math always runs in double precision, and every query also accepts coordinates
as a `ReadOnlySpan<double>` of consecutive (x, y, z) triples when the input itself needs
that range:

``` c#
double[] a = [1e6, 0, 0];
double[] b = [1e6 + 1e-3, 0, 0];

double distance = OpenGjk.ComputeMinimumDistance(a, b);   // 0.001, not 0

Span<double> normal = stackalloc double[3];
double info = OpenGjk.ComputeCollisionInformation(a, b, normal);

if (OpenGjk.TryGetPenetration(a, b, normal, out double depth))
{
    // Same unsigned result as the Vector3 overload: depth is positive, normal points
    // from the first body toward the second.
}
```

These overloads replace the legacy `CsFunction` (the obsolete openGJK C interface,
which took a `[3, n]` array); they are the same algorithm with a fuller result set.

### Reusing buffers in a query loop (`GjkWorkspace`)

Every query needs working buffers (simplex, support points, EPA polytope). The overloads
shown above allocate them per call, which is fine occasionally but becomes steady garbage
in a physics loop. Pass a `GjkWorkspace` to reuse the buffers instead - after the first few
queries have grown them, further queries of the same kind allocate nothing:

``` c#
var workspace = new GjkWorkspace();

foreach (var (a, b) in pairs)
{
    double distance = OpenGjk.ComputeMinimumDistance(a, b, workspace);
    double info = OpenGjk.ComputeCollisionInformation(a, b, workspace, out Vector3 normal);
}
```

Every `HasCollision`, `ComputeMinimumDistance`, and `ComputeCollisionInformation` overload
has a workspace counterpart, for vertex arrays, 2D polygons, and `ISupportMappable` shapes
alike. A workspace holds no result state, so it can be reused as soon as a query returns.

It is **not thread-safe**: use one workspace per thread (a `[ThreadStatic]` field, or one per
worker). Sharing one between concurrent queries would have them overwrite each other's buffers,
so a query that finds its workspace already occupied throws `InvalidOperationException` instead
of returning a wrong result:

``` c#
var shared = new GjkWorkspace();

Parallel.ForEach(pairs, pair =>
{
    // Throws InvalidOperationException as soon as two queries overlap.
    OpenGjk.ComputeMinimumDistance(pair.A, pair.B, shared);
});

// Correct: one workspace per worker.
var perThread = new ThreadLocal<GjkWorkspace>(() => new GjkWorkspace());

Parallel.ForEach(pairs, pair =>
{
    OpenGjk.ComputeMinimumDistance(pair.A, pair.B, perThread.Value!);
});
```

The same check catches re-entering a query with the workspace it was called from - for example
querying from inside a custom `GetSupportPoint`. A failed query releases its workspace, so the
workspace stays usable after the exception.

### Ray casts and shape casts

`TryCastRay` finds where a ray - the segment from `origin` to `origin + translation` - first
hits an `ISupportMappable` shape, and `TryCastShape` sweeps one shape along a translation and
finds where it first touches another, for continuous collision detection. Both report the hit
as a fraction of the translation in [0, 1]:

``` c#
var box = new Box(new Vector3(5, 0, 0), halfExtents: new Vector3(1, 1, 1));

if (OpenGjk.TryCastRay(box, origin: Vector3.Zero, translation: new Vector3(10, 0, 0),
        out double fraction, out Vector3 normal))
{
    // fraction = 0.4: the ray hits at (4, 0, 0), normal = (-1, 0, 0)
}

var ball = new Sphere(new Vector3(0, 0, 0), radius: 1f);

if (OpenGjk.TryCastShape(ball, box, translation: new Vector3(10, 0, 0),
        out double toi, out Vector3 contactNormal, out Vector3 contactPoint))
{
    // toi = 0.3: moving the ball by 0.3 * translation makes it touch the box,
    // contactNormal = (1, 0, 0) points from the ball toward the box,
    // contactPoint = about (4, 0, 0) - on a curved surface it lags slightly
}
```

The reported fraction never lies past the true time of impact, so moving by it never leaves
the shapes overlapping. Shapes count as touching when they come closer than the touch
distance: about 1e-6 times the size of the larger shape, or a few float ulps of the
coordinates far from the origin, where float support points are coarser than that (one ulp is
about 6e-5 at 1e3). A ray that starts inside the shape or within the touch distance of it,
or a sweep whose shapes overlap or touch at the start, reports a hit with fraction 0 and a
zero normal - use `TryGetPenetration` to separate such a pair.

The normal is exact on flat faces (except within the touch distance of an edge, where it may
belong to the adjacent face). On curved surfaces it lags the true normal by about 1e-3 near
the origin (up to 3e-3 on the side of a capsule); far from the origin float resolution
limits it to about 1e-2 for a unit sphere at 1e3. Both queries have `GjkWorkspace` overloads.

### Shapes and custom support functions

GJK is a support-function-based algorithm, so shapes like spheres and capsules can be
handled exactly - without tessellation - through the `ISupportMappable` interface.
Built-in shapes live in `OpenGJKSharp.Shapes`:

``` c#
using OpenGJKSharp;
using OpenGJKSharp.Shapes;
using System.Numerics;

var sphere = new Sphere(new Vector3(0, 0, 0), radius: 1f);
var box = new Box(new Vector3(3, 0, 0), halfExtents: new Vector3(1, 1, 1));

bool hit = OpenGjk.HasCollision(sphere, box);
double distance = OpenGjk.ComputeMinimumDistance(sphere, box, out Vector3 onSphere, out Vector3 onBox);
double info = OpenGjk.ComputeCollisionInformation(sphere, box, out Vector3 contactNormal);
```

`Sphere`, `Capsule`, `Box` (oriented, via a `Quaternion`), and `ConvexHull` are provided.

`TransformedShape` places any of them - or your own shape - at a pose, so a shape can be
defined once in local space and moved without rebuilding it. It is the only mutable shape:
`SetPose` updates the pose in place, which keeps a physics loop free of per-frame rebuilds:

``` c#
var hull = new ConvexHull(localPoints);   // built once, in local space
var body = new TransformedShape(hull);
var workspace = new GjkWorkspace();

foreach (var frame in frames)
{
    body.SetPose(frame.Rotation, frame.Position);
    double distance = OpenGjk.ComputeMinimumDistance(body, ground, workspace);
}
```

You can plug in any custom convex shape by implementing the interface:

``` c#
public sealed class Point : ISupportMappable
{
    private readonly Vector3 _position;

    public Point(Vector3 position) => _position = position;

    // Furthest point of the shape in the given direction (world space).
    public Vector3 GetSupportPoint(Vector3 direction) => _position;

    // An interior point, used as the initial guess for GJK.
    public Vector3 GetCenter() => _position;
}
```

#### Degenerate shapes

Shapes defined by a thickness parameter require it to be positive: `Sphere` and `Capsule`
reject a zero or negative radius, and `Box` rejects a zero or negative half extent, since
those are almost always a bug rather than an intent. A `Capsule` whose endpoints coincide is
not one of those cases - it is exactly a sphere and is handled as such.

Zero-thickness geometry is supported, through the shapes that take explicit points:
`ConvexHull` accepts a single point, two points (a segment), or coplanar points (a flat
polygon), and a custom `ISupportMappable` may describe any convex set. All of them, including
two coincident points, are queried without special cases.

Their penetration depth follows the rule above - the minimum translation within the plane or
along the line the pair shares - which is worth a moment's thought for a degenerate pair:

- A point inside a segment reports the distance to the *nearer end*, not 0: that is how far it
  has to travel to leave the segment. Same for a point inside a flat polygon, which reports the
  distance to the nearest edge.
- A segment lying exactly on the edge of a flat polygon reports 0. The two do touch, but the
  overlap has no area, so no translation is needed to separate them.
- Two coincident points report 0, which is the same reasoning one dimension lower.

#### Custom shapes

Implementations must describe a convex shape, and `GetSupportPoint` must be exact for
the shape - GJK/EPA accuracy depends on it. Note that the public interface uses
single-precision `Vector3`, so support points are float-quantized on input (the
vertex-array API has the same input quantization); internal math runs in double
precision.

## Migrating from 0.5.0

0.6.0 removes two types from the public surface and fixes three results that used to be wrong.
Everything else in the 0.5.0 API is unchanged.

### `GkPolytope` / `GkSimplex` are now internal

They exposed internal algorithm state - an undocumented `double[][]` coordinate layout and
settable counters that could break internal invariants - so they are `internal`, along with the
`OpenGjk` overloads that took them. Code using them fails to compile (CS0122/CS0246). Pass the
vertices directly instead:

``` c#
// 0.5.0
var bd1 = new GkPolytope { NumPoints = 3, Coord = [[0, 0, 0], [1, 0, 0], [0, 1, 0]] };
var bd2 = new GkPolytope { NumPoints = 3, Coord = [[2, 0, 0], [3, 0, 0], [2, 1, 0]] };
double distance = OpenGjk.ComputeMinimumDistance(bd1, bd2);

// 0.6.0
Vector3[] a = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0)];
Vector3[] b = [new(2, 0, 0), new(3, 0, 0), new(2, 1, 0)];
double distance = OpenGjk.ComputeMinimumDistance(a, b);
```

Note that `GkPolytope.Coord` was a `[point][axis]` array, while the double-precision overloads
take one flat span of consecutive (x, y, z) triples.

If you reused a `GkPolytope`/`GkSimplex` pair across queries to avoid per-query allocation, that
is what `GjkWorkspace` is for now - and it covers every query kind, not just the vertex arrays:

``` c#
// 0.5.0
var simplex = new GkSimplex();
foreach (var (bd1, bd2) in pairs)
{
    double distance = OpenGjk.ComputeMinimumDistance(bd1, bd2, simplex);
    OpenGjk.ComputeCollisionInformation(bd1, bd2, simplex, ref distance, normal);
}

// 0.6.0
var workspace = new GjkWorkspace();
foreach (var (a, b) in pairs)
{
    double distance = OpenGjk.ComputeMinimumDistance(a, b, workspace);
    double info = OpenGjk.ComputeCollisionInformation(a, b, workspace, out Vector3 normal);
}
```

A workspace must not be shared between threads; since 0.6.0 a query that finds its workspace
already in use throws `InvalidOperationException` instead of quietly returning a wrong result (see
[Reusing buffers in a query loop](#reusing-buffers-in-a-query-loop-gjkworkspace)).

### `CsFunction` callers

`CsFunction` (the legacy openGJK C interface, `[Obsolete]` since 0.4.0) still exists but is
superseded by the `ReadOnlySpan<double>` overloads, which are the same algorithm with the full
result set - witness points and penetration depth, not just the distance. It took a `[axis, point]`
array; the new overloads take consecutive (x, y, z) triples:

``` c#
// 0.5.0
double[,] a = { { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 0 } };   // [axis, point]
double[,] b = { { 2, 3, 2 }, { 0, 0, 1 }, { 0, 0, 0 } };
double distance = OpenGjk.CsFunction(3, a, 3, b);

// 0.6.0
double[] a = [0, 0, 0, 1, 0, 0, 0, 1, 0];                 // (x, y, z) triples
double[] b = [2, 0, 0, 3, 0, 0, 2, 1, 0];
double distance = OpenGjk.ComputeMinimumDistance(a, b);
```

### Results that changed

These are bug fixes, but if your code depended on the old values, check these three:

- **Colliding bodies now report a minimum distance of exactly 0.** `ComputeMinimumDistance` used
  to leave a residual (about 1e-16 at unit scale, larger for bigger bodies), so a `distance == 0`
  test missed collisions. If you compensated with an epsilon, the epsilon is no longer needed.
- **Penetration depth is larger for curved shapes.** A polytope face is a lower bound on the depth
  and can never resolve a curved surface, which cost up to 5.6% (two concentric unit spheres
  reported 1.888 instead of 2.0). The depth now comes from the tightest support plane, an upper
  bound, so worst-case error over a sphere overlap sweep is 0.24%. Depenetration that used to
  leave a sliver of overlap now separates the bodies; polyhedra are unaffected (still exact).
- **Flat pairs report an in-plane depth instead of 0.** Two coplanar bodies, or two bodies on a
  common line, used to come back with depth 0 - see
  [2D polygons and other flat geometry](#2d-polygons-and-other-flat-geometry) for what is
  reported now and why.
- **Bodies whose interior points coincide now report a depth instead of 0.** Concentric bodies, one
  body containing another, and two long bodies crossing at their common center used to come back
  with depth 0 - indistinguishable from "not colliding" for callers that tested the depth.

### `netstandard2.0` consumers gain a `System.Memory` dependency

The `ReadOnlySpan<double>` overloads need `Span<T>`, which is not in the box before
`netstandard2.1`, so the `netstandard2.0` target adds a `System.Memory` package reference (4.6.3).
This affects .NET Framework consumers only; `netstandard2.1`, `net8.0`, and `net10.0` are
unchanged. On .NET Framework, prefer a PackageReference project (SDK-style csproj or
`packages.config` migrated to PackageReference) - transitive dependencies do not resolve
automatically under classic `packages.config`.

## License

OpenGJKSharp is licensed under [GPL-3.0](LICENSE), following the upstream
[openGJK](https://github.com/MattiaMontanari/openGJK) project. Note that GPL-3.0
has copyleft requirements — if you plan to use this library in closed-source or
commercial software, please review the license terms carefully.
# Support Function API (ISupportMappable) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a public support-mapping interface (`ISupportMappable`) with built-in Sphere/Capsule/Box/ConvexHull shapes and `OpenGjk` overloads, backed by an internal body abstraction so GJK/EPA work on both vertex clouds and implicit shapes.

**Architecture:** Introduce an internal `GkBody` abstraction with two implementations: `GkPolytopeBody` (wraps the existing `GkPolytope`, preserving exact current behavior including stateful support and `S`/`SIdx` mutation) and `GkMappedBody` (adapts a user `ISupportMappable`; every support point it returns is appended to an internal point store and identified by its store index). Because points are index-addressable in both implementations, ALL existing index plumbing (`GkSimplex.VrtxIdx`, EPA `VIdx`) keeps working unchanged; witness reconstruction just switches from `bd.Coord[idx]` to `bd.GetPoint(idx)`. Note: EPA duplicate-support detection is ALREADY coordinate-based in the current code (`epsSq` comparisons), so the spec section about generalizing it requires no code change.

**Tech Stack:** C# multi-target `netstandard2.0;netstandard2.1;net8.0;net10.0`, System.Numerics (`Vector3`, `Quaternion`), xunit.

**Spec:** `docs/superpowers/specs/2026-08-01-support-function-api-design.md`

## Global Constraints

- netstandard2.0 compatibility: no `Math.Clamp`, `float.IsFinite`, `MathF` in library code (tests target net8.0 and may use them). Use internal helpers.
- Commit messages in English, conventional format `<type>: <subject>`, imperative, <= 50 chars subject.
- Existing public API signatures and behavior must not change; full existing test suite stays green after every task.
- Test project: `tests/OpenGJKSharp.Tests` (xunit). Run with `dotnet test tests/OpenGJKSharp.Tests`.
- New public types get XML doc comments (`GenerateDocumentationFile` is on; missing docs cause warnings).
- Target version 0.5.0; CHANGELOG follows Keep a Changelog.

## File Structure

- Create: `src/OpenGJKSharp/ISupportMappable.cs` - public interface (namespace `OpenGJKSharp`).
- Create: `src/OpenGJKSharp/Shapes/ShapeArgument.cs` - internal validation helpers shared by shapes.
- Create: `src/OpenGJKSharp/Shapes/Sphere.cs`, `Capsule.cs`, `Box.cs`, `ConvexHull.cs` - built-in shapes (namespace `OpenGJKSharp.Shapes`).
- Create: `src/OpenGJKSharp/Models/GkBody.cs` - internal abstract body contract.
- Create: `src/OpenGJKSharp/Models/GkPolytopeBody.cs` - vertex-cloud implementation.
- Create: `src/OpenGJKSharp/Models/GkMappedBody.cs` - support-mapped implementation.
- Modify: `src/OpenGJKSharp/OpenGjk.cs` - core loops take `GkBody`; public `GkPolytope` overloads become thin wrappers; new `ISupportMappable` overloads.
- Create: `tests/OpenGJKSharp.Tests/ShapeTests.cs` - shape constructor validation + support-mapping unit tests.
- Create: `tests/OpenGJKSharp.Tests/SupportMappableQueryTests.cs` - analytic, parity, and validation tests for the new overloads.
- Modify: `src/OpenGJKSharp/OpenGJKSharp.csproj` (version), `CHANGELOG.md`, `README.md`.

---

### Task 1: ISupportMappable interface + Sphere shape

**Files:**
- Create: `src/OpenGJKSharp/ISupportMappable.cs`
- Create: `src/OpenGJKSharp/Shapes/ShapeArgument.cs`
- Create: `src/OpenGJKSharp/Shapes/Sphere.cs`
- Create: `tests/OpenGJKSharp.Tests/ShapeTests.cs`

**Interfaces:**
- Consumes: nothing.
- Produces: `OpenGJKSharp.ISupportMappable` with `Vector3 GetSupportPoint(Vector3 direction)` and `Vector3 GetCenter()`; `OpenGJKSharp.Shapes.Sphere(Vector3 center, float radius)` with `Vector3 Center { get; }`, `float Radius { get; }`; internal `OpenGJKSharp.Shapes.ShapeArgument` with `static bool IsFinite(float)`, `static void EnsureFinite(Vector3 value, string paramName, string description)`, `static void EnsurePositiveFinite(float value, string paramName, string description)`.

- [ ] **Step 1: Write the failing tests**

Create `tests/OpenGJKSharp.Tests/ShapeTests.cs`:

```csharp
using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class ShapeTests
{
    [Fact]
    public void Should_ReturnSurfacePoint_When_SphereSupportIsQueried()
    {
        var sphere = new Sphere(new Vector3(1, 2, 3), 2f);

        // Direction does not need to be normalized
        Vector3 support = sphere.GetSupportPoint(new Vector3(10, 0, 0));

        Assert.Equal(3.0, support.X, 5);
        Assert.Equal(2.0, support.Y, 5);
        Assert.Equal(3.0, support.Z, 5);
        Assert.Equal(new Vector3(1, 2, 3), sphere.GetCenter());
    }

    [Theory]
    [InlineData(0f)]
    [InlineData(-1f)]
    [InlineData(float.NaN)]
    [InlineData(float.PositiveInfinity)]
    public void Should_ThrowArgumentOutOfRangeException_When_SphereRadiusIsInvalid(float radius)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new Sphere(Vector3.Zero, radius));
    }

    [Fact]
    public void Should_ThrowArgumentException_When_SphereCenterIsNotFinite()
    {
        Assert.Throws<ArgumentException>(() => new Sphere(new Vector3(float.NaN, 0, 0), 1f));
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: compile FAILURE ("Sphere" / namespace not found).

- [ ] **Step 3: Write the implementation**

Create `src/OpenGJKSharp/ISupportMappable.cs`:

```csharp
using System.Numerics;

namespace OpenGJKSharp;

/// <summary>
/// A convex shape described by its support mapping.
/// Implementations must describe a convex shape, and <see cref="GetSupportPoint"/>
/// must be exact for the shape - GJK/EPA accuracy depends on it.
/// </summary>
public interface ISupportMappable
{
    /// <summary>
    /// Returns the furthest point of the shape in the given direction (world space).
    /// The direction is not required to be normalized and is never zero.
    /// </summary>
    /// <param name="direction">The (non-zero) search direction.</param>
    /// <returns>The furthest point of the shape along <paramref name="direction"/>.</returns>
    Vector3 GetSupportPoint(Vector3 direction);

    /// <summary>
    /// Returns an interior point of the shape (typically the center).
    /// Used as the initial guess for GJK and as the seed for EPA.
    /// </summary>
    /// <returns>An interior point of the shape.</returns>
    Vector3 GetCenter();
}
```

Create `src/OpenGJKSharp/Shapes/ShapeArgument.cs`:

```csharp
using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// Shared constructor-argument validation for the built-in shapes.
/// </summary>
internal static class ShapeArgument
{
    // float.IsFinite is unavailable on netstandard2.0
    public static bool IsFinite(float value)
    {
        return !float.IsNaN(value) && !float.IsInfinity(value);
    }

    public static void EnsureFinite(Vector3 value, string paramName, string description)
    {
        if (!IsFinite(value.X) || !IsFinite(value.Y) || !IsFinite(value.Z))
        {
            throw new ArgumentException($"{description} contains a non-finite coordinate: {value}.", paramName);
        }
    }

    public static void EnsurePositiveFinite(float value, string paramName, string description)
    {
        if (!IsFinite(value) || value <= 0)
        {
            throw new ArgumentOutOfRangeException(paramName, value, $"{description} must be positive and finite.");
        }
    }
}
```

Create `src/OpenGJKSharp/Shapes/Sphere.cs`:

```csharp
using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// A sphere defined by a world-space center and radius.
/// </summary>
public sealed class Sphere : ISupportMappable
{
    /// <summary>The center of the sphere.</summary>
    public Vector3 Center { get; }

    /// <summary>The radius of the sphere.</summary>
    public float Radius { get; }

    /// <summary>
    /// Creates a sphere.
    /// </summary>
    /// <param name="center">The center of the sphere. Must be finite.</param>
    /// <param name="radius">The radius of the sphere. Must be positive and finite.</param>
    /// <exception cref="ArgumentException"><paramref name="center"/> is not finite.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="radius"/> is not positive and finite.</exception>
    public Sphere(Vector3 center, float radius)
    {
        ShapeArgument.EnsureFinite(center, nameof(center), "Center");
        ShapeArgument.EnsurePositiveFinite(radius, nameof(radius), "Radius");

        Center = center;
        Radius = radius;
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        return Center + Radius * Vector3.Normalize(direction);
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return Center;
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: PASS (all ShapeTests).

- [ ] **Step 5: Commit**

```bash
git add src/OpenGJKSharp/ISupportMappable.cs src/OpenGJKSharp/Shapes/ShapeArgument.cs src/OpenGJKSharp/Shapes/Sphere.cs tests/OpenGJKSharp.Tests/ShapeTests.cs
git commit -m "feat: add ISupportMappable interface and Sphere shape"
```

---

### Task 2: Capsule shape

**Files:**
- Create: `src/OpenGJKSharp/Shapes/Capsule.cs`
- Modify: `tests/OpenGJKSharp.Tests/ShapeTests.cs` (add tests)

**Interfaces:**
- Consumes: `ISupportMappable`, `ShapeArgument` (Task 1).
- Produces: `OpenGJKSharp.Shapes.Capsule(Vector3 p0, Vector3 p1, float radius)` with `Vector3 P0 { get; }`, `Vector3 P1 { get; }`, `float Radius { get; }`. `p0 == p1` is allowed (degenerates to a sphere).

- [ ] **Step 1: Write the failing tests**

Add to `ShapeTests.cs`:

```csharp
[Fact]
public void Should_ReturnEndCapPoint_When_CapsuleSupportIsQueried()
{
    var capsule = new Capsule(new Vector3(0, 0, 0), new Vector3(0, 2, 0), 0.5f);

    Vector3 top = capsule.GetSupportPoint(new Vector3(0, 3, 0));
    Vector3 bottom = capsule.GetSupportPoint(new Vector3(0, -1, 0));

    Assert.Equal(2.5, top.Y, 5);
    Assert.Equal(-0.5, bottom.Y, 5);
    Assert.Equal(1.0, capsule.GetCenter().Y, 5);
}

[Fact]
public void Should_AllowDegenerateCapsule_When_EndpointsCoincide()
{
    var capsule = new Capsule(new Vector3(1, 1, 1), new Vector3(1, 1, 1), 0.5f);

    Vector3 support = capsule.GetSupportPoint(new Vector3(1, 0, 0));

    Assert.Equal(1.5, support.X, 5);
}

[Theory]
[InlineData(0f)]
[InlineData(-0.5f)]
[InlineData(float.NaN)]
public void Should_ThrowArgumentOutOfRangeException_When_CapsuleRadiusIsInvalid(float radius)
{
    Assert.Throws<ArgumentOutOfRangeException>(() => new Capsule(Vector3.Zero, Vector3.UnitY, radius));
}

[Fact]
public void Should_ThrowArgumentException_When_CapsuleEndpointIsNotFinite()
{
    Assert.Throws<ArgumentException>(() => new Capsule(new Vector3(float.NaN, 0, 0), Vector3.UnitY, 1f));
    Assert.Throws<ArgumentException>(() => new Capsule(Vector3.Zero, new Vector3(0, float.PositiveInfinity, 0), 1f));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: compile FAILURE ("Capsule" not found).

- [ ] **Step 3: Write the implementation**

Create `src/OpenGJKSharp/Shapes/Capsule.cs`:

```csharp
using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// A capsule: the set of points within <see cref="Radius"/> of the world-space
/// segment from <see cref="P0"/> to <see cref="P1"/>.
/// </summary>
public sealed class Capsule : ISupportMappable
{
    /// <summary>The first endpoint of the capsule axis.</summary>
    public Vector3 P0 { get; }

    /// <summary>The second endpoint of the capsule axis.</summary>
    public Vector3 P1 { get; }

    /// <summary>The radius of the capsule.</summary>
    public float Radius { get; }

    /// <summary>
    /// Creates a capsule.
    /// </summary>
    /// <param name="p0">The first endpoint of the axis. Must be finite.</param>
    /// <param name="p1">The second endpoint of the axis. Must be finite.
    /// May equal <paramref name="p0"/>, in which case the capsule degenerates to a sphere.</param>
    /// <param name="radius">The radius of the capsule. Must be positive and finite.</param>
    /// <exception cref="ArgumentException"><paramref name="p0"/> or <paramref name="p1"/> is not finite.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="radius"/> is not positive and finite.</exception>
    public Capsule(Vector3 p0, Vector3 p1, float radius)
    {
        ShapeArgument.EnsureFinite(p0, nameof(p0), "P0");
        ShapeArgument.EnsureFinite(p1, nameof(p1), "P1");
        ShapeArgument.EnsurePositiveFinite(radius, nameof(radius), "Radius");

        P0 = p0;
        P1 = p1;
        Radius = radius;
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        Vector3 endpoint = Vector3.Dot(P0, direction) >= Vector3.Dot(P1, direction) ? P0 : P1;

        return endpoint + Radius * Vector3.Normalize(direction);
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return (P0 + P1) * 0.5f;
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/OpenGJKSharp/Shapes/Capsule.cs tests/OpenGJKSharp.Tests/ShapeTests.cs
git commit -m "feat: add Capsule shape"
```

---

### Task 3: Box shape

**Files:**
- Create: `src/OpenGJKSharp/Shapes/Box.cs`
- Modify: `tests/OpenGJKSharp.Tests/ShapeTests.cs` (add tests)

**Interfaces:**
- Consumes: `ISupportMappable`, `ShapeArgument` (Task 1).
- Produces: `OpenGJKSharp.Shapes.Box(Vector3 center, Vector3 halfExtents, Quaternion orientation)` and convenience `Box(Vector3 center, Vector3 halfExtents)` (identity orientation), with `Vector3 Center { get; }`, `Vector3 HalfExtents { get; }`, `Quaternion Orientation { get; }` (normalized).

- [ ] **Step 1: Write the failing tests**

Add to `ShapeTests.cs` (tests target net8.0, so `MathF` is fine here):

```csharp
[Fact]
public void Should_ReturnCorner_When_AxisAlignedBoxSupportIsQueried()
{
    var box = new Box(new Vector3(0, 0, 0), new Vector3(1, 2, 3));

    Vector3 support = box.GetSupportPoint(new Vector3(1, -1, 1));

    Assert.Equal(1.0, support.X, 5);
    Assert.Equal(-2.0, support.Y, 5);
    Assert.Equal(3.0, support.Z, 5);
}

[Fact]
public void Should_ReturnRotatedCorner_When_RotatedBoxSupportIsQueried()
{
    // 90 degrees about +Z maps local +X to world +Y.
    var orientation = Quaternion.CreateFromAxisAngle(Vector3.UnitZ, MathF.PI / 2);
    var box = new Box(new Vector3(10, 0, 0), new Vector3(1, 2, 3), orientation);

    // Local direction = inverse-rotated (1, 0.2, 0.5) = (0.2, -1, 0.5)
    // -> local corner (+1, -2, +3) -> world (2, 1, 3) + center.
    Vector3 support = box.GetSupportPoint(new Vector3(1, 0.2f, 0.5f));

    Assert.Equal(12.0, support.X, 4);
    Assert.Equal(1.0, support.Y, 4);
    Assert.Equal(3.0, support.Z, 4);
}

[Theory]
[InlineData(0f, 1f, 1f)]
[InlineData(1f, -1f, 1f)]
[InlineData(1f, 1f, float.NaN)]
public void Should_ThrowArgumentException_When_BoxHalfExtentsAreInvalid(float x, float y, float z)
{
    Assert.Throws<ArgumentException>(() => new Box(Vector3.Zero, new Vector3(x, y, z)));
}

[Fact]
public void Should_ThrowArgumentException_When_BoxOrientationIsInvalid()
{
    var zero = new Quaternion(0, 0, 0, 0);
    var nan = new Quaternion(float.NaN, 0, 0, 1);

    Assert.Throws<ArgumentException>(() => new Box(Vector3.Zero, Vector3.One, zero));
    Assert.Throws<ArgumentException>(() => new Box(Vector3.Zero, Vector3.One, nan));
}

[Fact]
public void Should_NormalizeOrientation_When_BoxIsCreated()
{
    var doubled = new Quaternion(0, 0, 0, 2f); // identity scaled by 2
    var box = new Box(Vector3.Zero, Vector3.One, doubled);

    Assert.Equal(1.0, box.Orientation.Length(), 5);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: compile FAILURE ("Box" not found).

- [ ] **Step 3: Write the implementation**

Create `src/OpenGJKSharp/Shapes/Box.cs`:

```csharp
using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// An oriented box defined by a world-space center, half extents, and orientation.
/// </summary>
public sealed class Box : ISupportMappable
{
    /// <summary>The center of the box.</summary>
    public Vector3 Center { get; }

    /// <summary>The half extents of the box along its local axes.</summary>
    public Vector3 HalfExtents { get; }

    /// <summary>The orientation of the box (normalized).</summary>
    public Quaternion Orientation { get; }

    /// <summary>
    /// Creates an axis-aligned box (identity orientation).
    /// </summary>
    /// <param name="center">The center of the box. Must be finite.</param>
    /// <param name="halfExtents">The half extents of the box. Every component must be positive and finite.</param>
    /// <exception cref="ArgumentException">An argument is not finite or a half extent is not positive.</exception>
    public Box(Vector3 center, Vector3 halfExtents)
        : this(center, halfExtents, Quaternion.Identity)
    {
    }

    /// <summary>
    /// Creates an oriented box.
    /// </summary>
    /// <param name="center">The center of the box. Must be finite.</param>
    /// <param name="halfExtents">The half extents of the box. Every component must be positive and finite.</param>
    /// <param name="orientation">The orientation of the box. Must be finite with non-zero length; normalized internally.</param>
    /// <exception cref="ArgumentException">An argument is not finite, a half extent is not positive,
    /// or <paramref name="orientation"/> has zero length.</exception>
    public Box(Vector3 center, Vector3 halfExtents, Quaternion orientation)
    {
        ShapeArgument.EnsureFinite(center, nameof(center), "Center");

        if (!ShapeArgument.IsFinite(halfExtents.X) || !ShapeArgument.IsFinite(halfExtents.Y) || !ShapeArgument.IsFinite(halfExtents.Z)
            || halfExtents.X <= 0 || halfExtents.Y <= 0 || halfExtents.Z <= 0)
        {
            throw new ArgumentException($"Half extents must be positive and finite: {halfExtents}.", nameof(halfExtents));
        }

        if (!ShapeArgument.IsFinite(orientation.X) || !ShapeArgument.IsFinite(orientation.Y)
            || !ShapeArgument.IsFinite(orientation.Z) || !ShapeArgument.IsFinite(orientation.W))
        {
            throw new ArgumentException($"Orientation contains a non-finite component: {orientation}.", nameof(orientation));
        }

        float lengthSquared = orientation.LengthSquared();
        if (lengthSquared <= 0)
        {
            throw new ArgumentException("Orientation must have a non-zero length.", nameof(orientation));
        }

        Center = center;
        HalfExtents = halfExtents;
        Orientation = Quaternion.Normalize(orientation);
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        // Rotate the direction into local space, pick the corner per axis, rotate back.
        Vector3 localDirection = Vector3.Transform(direction, Quaternion.Conjugate(Orientation));
        var localSupport = new Vector3(
            localDirection.X >= 0 ? HalfExtents.X : -HalfExtents.X,
            localDirection.Y >= 0 ? HalfExtents.Y : -HalfExtents.Y,
            localDirection.Z >= 0 ? HalfExtents.Z : -HalfExtents.Z);

        return Center + Vector3.Transform(localSupport, Orientation);
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return Center;
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/OpenGJKSharp/Shapes/Box.cs tests/OpenGJKSharp.Tests/ShapeTests.cs
git commit -m "feat: add Box shape"
```

---

### Task 4: ConvexHull shape

**Files:**
- Create: `src/OpenGJKSharp/Shapes/ConvexHull.cs`
- Modify: `tests/OpenGJKSharp.Tests/ShapeTests.cs` (add tests)

**Interfaces:**
- Consumes: `ISupportMappable`, `ShapeArgument` (Task 1).
- Produces: `OpenGJKSharp.Shapes.ConvexHull(Vector3[] points)`. Copies the input array (immutability). `GetCenter()` returns the vertex centroid (mean), which is interior for a convex hull.

- [ ] **Step 1: Write the failing tests**

Add to `ShapeTests.cs`:

```csharp
[Fact]
public void Should_ReturnFurthestVertex_When_ConvexHullSupportIsQueried()
{
    var hull = new ConvexHull(
    [
        new Vector3(0, 0, 0), new Vector3(1, 0, 0), new Vector3(0, 1, 0), new Vector3(1, 1, 0),
        new Vector3(0, 0, 1), new Vector3(1, 0, 1), new Vector3(0, 1, 1), new Vector3(1, 1, 1),
    ]);

    Vector3 support = hull.GetSupportPoint(new Vector3(1, 1, 1));

    Assert.Equal(new Vector3(1, 1, 1), support);
    Assert.Equal(0.5, hull.GetCenter().X, 5);
}

[Fact]
public void Should_NotReflectCallerMutations_When_ConvexHullInputArrayChanges()
{
    Vector3[] points = [new(0, 0, 0), new(1, 0, 0)];
    var hull = new ConvexHull(points);

    points[1] = new Vector3(100, 0, 0);

    Assert.Equal(1.0, hull.GetSupportPoint(new Vector3(1, 0, 0)).X, 5);
}

[Fact]
public void Should_ThrowOnInvalidInput_When_ConvexHullIsCreated()
{
    Assert.Throws<ArgumentNullException>(() => new ConvexHull(null!));
    Assert.Throws<ArgumentException>(() => new ConvexHull([]));
    Assert.Throws<ArgumentException>(() => new ConvexHull([new Vector3(float.NaN, 0, 0)]));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: compile FAILURE ("ConvexHull" not found).

- [ ] **Step 3: Write the implementation**

Create `src/OpenGJKSharp/Shapes/ConvexHull.cs`:

```csharp
using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// The convex hull of a set of world-space points.
/// The hull is never computed explicitly; the support mapping scans the vertices.
/// </summary>
public sealed class ConvexHull : ISupportMappable
{
    private readonly Vector3[] _points;
    private readonly Vector3 _centroid;

    /// <summary>
    /// Creates a convex hull from a set of points. The array is copied.
    /// </summary>
    /// <param name="points">The points spanning the hull. Must contain at least one finite point.</param>
    /// <exception cref="ArgumentNullException"><paramref name="points"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="points"/> is empty or contains a non-finite point.</exception>
    public ConvexHull(Vector3[] points)
    {
        if (points is null)
        {
            throw new ArgumentNullException(nameof(points));
        }

        if (points.Length == 0)
        {
            throw new ArgumentException("At least one point is required.", nameof(points));
        }

        var copy = new Vector3[points.Length];
        var sum = Vector3.Zero;
        for (int i = 0; i < points.Length; i++)
        {
            ShapeArgument.EnsureFinite(points[i], nameof(points), $"Point at index {i}");
            copy[i] = points[i];
            sum += points[i];
        }

        _points = copy;
        _centroid = sum / points.Length;
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        Vector3 best = _points[0];
        float bestDot = Vector3.Dot(best, direction);

        for (int i = 1; i < _points.Length; i++)
        {
            float dot = Vector3.Dot(_points[i], direction);
            if (dot > bestDot)
            {
                bestDot = dot;
                best = _points[i];
            }
        }

        return best;
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return _centroid;
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~ShapeTests"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/OpenGJKSharp/Shapes/ConvexHull.cs tests/OpenGJKSharp.Tests/ShapeTests.cs
git commit -m "feat: add ConvexHull shape"
```

---

### Task 5: Internal GkBody abstraction (behavior-preserving refactor)

**Files:**
- Create: `src/OpenGJKSharp/Models/GkBody.cs`
- Create: `src/OpenGJKSharp/Models/GkPolytopeBody.cs`
- Modify: `src/OpenGJKSharp/OpenGjk.cs`

This task is a pure refactor: NO public behavior change. The test gate is the existing suite staying green. Parity rules that must hold exactly:

- `GkPolytopeBody` must mutate the wrapped polytope's public `S`/`SIdx` (users can observe them).
- The stateful GJK support keeps its current semantics: `maxS` starts from the CURRENT `S`, and only a STRICTLY greater dot product replaces it (ties keep the old point).
- The EPA support scan keeps first-strict-max tie-breaking.
- `dot(c, -d)` equals `-dot(c, d)` bitwise in IEEE, so negating the direction for body 2 (instead of negating each dot) preserves results exactly.

**Interfaces:**
- Consumes: `GkPolytope`, `GkSimplex` (existing).
- Produces: internal `abstract class GkBody` in `OpenGJKSharp.Models` with:
  - `abstract double[] S { get; }` - current support point storage
  - `abstract int SIdx { get; }` - id of the current support point
  - `abstract void InitSupport()` - seeds S/SIdx (first vertex / shape center)
  - `abstract void Support(double[] v)` - stateful GJK support
  - `abstract void SupportInto(double[] direction, double[] result, out int idx)` - stateless support (EPA)
  - `abstract double[] GetPoint(int idx)` - resolves an id back to coordinates
  - `internal sealed class GkPolytopeBody : GkBody` with constructor `GkPolytopeBody(GkPolytope polytope)`
  - internal `static double ComputeMinimumDistance(GkBody bd1, GkBody bd2, GkSimplex s)` and internal `static void ComputeCollisionInformation(GkBody bd1, GkBody bd2, GkSimplex simplex, ref double distance, double[] contactNormal)` on `OpenGjk` (Task 6 calls these).

- [ ] **Step 1: Verify green baseline**

Run: `dotnet test tests/OpenGJKSharp.Tests`
Expected: PASS (all existing tests).

- [ ] **Step 2: Create the abstraction**

Create `src/OpenGJKSharp/Models/GkBody.cs`:

```csharp
namespace OpenGJKSharp.Models;

/// <summary>
/// Internal abstraction over the two body kinds the GJK/EPA loops operate on:
/// vertex clouds (<see cref="GkPolytope"/>) and support-mapped shapes (<see cref="ISupportMappable"/>).
/// Support points are identified by an integer id so the index-based simplex/EPA
/// bookkeeping works for both kinds; <see cref="GetPoint"/> resolves ids back to
/// coordinates for witness reconstruction.
/// </summary>
internal abstract class GkBody
{
    /// <summary>Current support point, seeded by <see cref="InitSupport"/> and updated by <see cref="Support"/>.</summary>
    public abstract double[] S { get; }

    /// <summary>Id of the current support point, resolvable via <see cref="GetPoint"/>.</summary>
    public abstract int SIdx { get; }

    /// <summary>Seeds <see cref="S"/>/<see cref="SIdx"/> with the initial point (first vertex or shape center).</summary>
    public abstract void InitSupport();

    /// <summary>Stateful GJK support: moves <see cref="S"/>/<see cref="SIdx"/> to the furthest point along <paramref name="v"/>.</summary>
    public abstract void Support(double[] v);

    /// <summary>Stateless support used by EPA: writes the furthest point along <paramref name="direction"/> into <paramref name="result"/>.</summary>
    public abstract void SupportInto(double[] direction, double[] result, out int idx);

    /// <summary>Resolves a support point id to its coordinates.</summary>
    public abstract double[] GetPoint(int idx);
}
```

Create `src/OpenGJKSharp/Models/GkPolytopeBody.cs`:

```csharp
namespace OpenGJKSharp.Models;

/// <summary>
/// Vertex-cloud body backed by a <see cref="GkPolytope"/>.
/// Delegates support state to the wrapped polytope's <see cref="GkPolytope.S"/>/<see cref="GkPolytope.SIdx"/>
/// so their publicly observable behavior is unchanged.
/// </summary>
internal sealed class GkPolytopeBody : GkBody
{
    private readonly GkPolytope _polytope;

    public GkPolytopeBody(GkPolytope polytope)
    {
        _polytope = polytope;
    }

    public override double[] S => _polytope.S;

    public override int SIdx => _polytope.SIdx;

    public override void InitSupport()
    {
        for (int t = 0; t < 3; ++t)
        {
            _polytope.S[t] = _polytope.Coord[0][t];
        }
        _polytope.SIdx = 0;
    }

    public override void Support(double[] v)
    {
        // Stateful scan: starts from the current support point and only a strictly
        // better vertex replaces it (ties keep the previous point/index).
        double maxS = Dot(_polytope.S, v);
        int better = -1;

        for (int i = 0; i < _polytope.NumPoints; ++i)
        {
            double s = Dot(_polytope.Coord[i], v);
            if (s > maxS)
            {
                maxS = s;
                better = i;
            }
        }

        if (better != -1)
        {
            _polytope.S[0] = _polytope.Coord[better][0];
            _polytope.S[1] = _polytope.Coord[better][1];
            _polytope.S[2] = _polytope.Coord[better][2];
            _polytope.SIdx = better;
        }
    }

    public override void SupportInto(double[] direction, double[] result, out int idx)
    {
        // Fresh scan, first strict max wins (matches the previous EPA scan).
        int best = 0;
        double bestDot = Dot(_polytope.Coord[0], direction);

        for (int i = 1; i < _polytope.NumPoints; ++i)
        {
            double s = Dot(_polytope.Coord[i], direction);
            if (s > bestDot)
            {
                bestDot = s;
                best = i;
            }
        }

        for (int t = 0; t < 3; ++t)
        {
            result[t] = _polytope.Coord[best][t];
        }
        idx = best;
    }

    public override double[] GetPoint(int idx)
    {
        return _polytope.Coord[idx];
    }

    private static double Dot(double[] a, double[] b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
}
```

- [ ] **Step 3: Refactor OpenGjk.cs onto the abstraction**

All edits inside `src/OpenGJKSharp/OpenGjk.cs`:

1. **Public 3-arg GJK entry becomes a wrapper.** Replace the body of `public static double ComputeMinimumDistance(GkPolytope bd1, GkPolytope bd2, GkSimplex s)` so that after the existing `ValidatePolytope`/null checks it returns:

```csharp
return ComputeMinimumDistance(new GkPolytopeBody(bd1), new GkPolytopeBody(bd2), s);
```

2. **Move the loop into an internal GkBody overload.** The existing loop body moves to:

```csharp
internal static double ComputeMinimumDistance(GkBody bd1, GkBody bd2, GkSimplex s)
```

with these edits inside the moved code:
   - The initialization block (previously `v = Coord[0]A - Coord[0]B`, copy `Coord[0]` into `S`, `SIdx = 0`) becomes:

```csharp
bd1.InitSupport();
bd2.InitSupport();

// Initialize search direction
v[0] = bd1.S[0] - bd2.S[0];
v[1] = bd1.S[1] - bd2.S[1];
v[2] = bd1.S[2] - bd2.S[2];

// Initialize simplex
s.NVrtx = 1;
for (int t = 0; t < 3; ++t)
{
    s.Vrtx[0, t] = v[t];
}
s.VrtxIdx[0, 0] = bd1.SIdx;
s.VrtxIdx[0, 1] = bd2.SIdx;
```

   - `Support(bd1, vMinus);` -> `bd1.Support(vMinus);` and `Support(bd2, v);` -> `bd2.Support(v);`.
   - Everything else (tolerances, sub-algorithm, exit conditions, `ComputeWitnesses` call) stays identical.

3. **Delete the private `static void Support(GkPolytope body, double[] v)`** (logic now lives in `GkPolytopeBody.Support`).

4. **Witness functions take GkBody.** Change parameter types `GkPolytope bd1/bd2` -> `GkBody bd1/bd2` on `ComputeWitnesses`, `W0D`, `W1D`, `W2D`, `W3D`, and replace every `bd1.Coord[<idx>]` with `bd1.GetPoint(<idx>)` and every `bd2.Coord[<idx>]` with `bd2.GetPoint(<idx>)` in those functions (W0D: 2 sites; W1D: 4 sites; W2D: 6 sites; W3D: 8 sites).

5. **EPA support goes through SupportInto.** Replace the private `SupportEpa(GkPolytope body1, GkPolytope body2, ...)` with:

```csharp
// Support function for EPA: fresh (stateless) support of each body,
// body2 searched in the opposite direction, combined into the Minkowski difference point.
private static void SupportEpa(GkBody body1, GkBody body2,
    double[] direction, double[] result, int[] resultIdx)
{
    double[] negDirection = [-direction[0], -direction[1], -direction[2]];
    double[] p1 = new double[3];
    double[] p2 = new double[3];

    body1.SupportInto(direction, p1, out int idx1);
    body2.SupportInto(negDirection, p2, out int idx2);

    for (int i = 0; i < 3; i++)
    {
        result[i] = p1[i] - p2[i];
    }
    resultIdx[0] = idx1;
    resultIdx[1] = idx2;
}
```

6. **Public EPA entry becomes a wrapper.** `public static void ComputeCollisionInformation(GkPolytope bd1, GkPolytope bd2, GkSimplex simplex, ref double distance, double[] contactNormal)` delegates:

```csharp
ComputeCollisionInformation(new GkPolytopeBody(bd1), new GkPolytopeBody(bd2), simplex, ref distance, contactNormal);
```

and the existing body moves to `internal static void ComputeCollisionInformation(GkBody bd1, GkBody bd2, GkSimplex simplex, ref double distance, double[] contactNormal)` with every `bd1.Coord[<idx>]` / `bd2.Coord[<idx>]` replaced by `GetPoint(<idx>)` (the simplex-rebuild fallback witness assignments - 4 blocks - plus the converged/duplicate/max-iteration witness blocks - 3 blocks of 6 lookups each).

7. **Array-based EPA entry builds bodies once.** In `public static double ComputeCollisionInformation(Vector3[] a, Vector3[] b, out Vector3 contactNormal)` replace the `var bd1 = ToGkPolytope(a); var bd2 = ToGkPolytope(b);` usage with:

```csharp
var bd1 = new GkPolytopeBody(ToGkPolytope(a));
var bd2 = new GkPolytopeBody(ToGkPolytope(b));
```

(the subsequent `ComputeMinimumDistance(bd1, bd2, simplex)` and `ComputeCollisionInformation(bd1, bd2, simplex, ref distance, normal)` calls now bind to the internal `GkBody` overloads).

- [ ] **Step 4: Run the full suite to verify parity**

Run: `dotnet test tests/OpenGJKSharp.Tests`
Expected: PASS - every existing test, unchanged. Also run `dotnet build src/OpenGJKSharp -c Release` to confirm all four target frameworks compile without warnings.

- [ ] **Step 5: Commit**

```bash
git add src/OpenGJKSharp/Models/GkBody.cs src/OpenGJKSharp/Models/GkPolytopeBody.cs src/OpenGJKSharp/OpenGjk.cs
git commit -m "refactor: route GJK/EPA through internal GkBody abstraction"
```

---

### Task 6: GkMappedBody + OpenGjk ISupportMappable overloads

**Files:**
- Create: `src/OpenGJKSharp/Models/GkMappedBody.cs`
- Modify: `src/OpenGJKSharp/OpenGjk.cs` (new public overloads)
- Create: `tests/OpenGJKSharp.Tests/SupportMappableQueryTests.cs`

**Interfaces:**
- Consumes: `GkBody`, internal `ComputeMinimumDistance(GkBody, GkBody, GkSimplex)` and `ComputeCollisionInformation(GkBody, GkBody, GkSimplex, ref double, double[])` (Task 5); `ISupportMappable`, `Sphere` (Task 1).
- Produces:
  - internal `GkMappedBody(ISupportMappable shape) : GkBody`
  - `public static bool HasCollision(ISupportMappable a, ISupportMappable b, double precision = 1e-6)`
  - `public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b)`
  - `public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b, out Vector3 closestA, out Vector3 closestB)`
  - `public static double ComputeCollisionInformation(ISupportMappable a, ISupportMappable b, out Vector3 contactNormal)`

**Critical invariant:** one query = one pair of `GkMappedBody` instances shared between the GJK and EPA passes, because simplex/EPA indices refer to each instance's point store.

- [ ] **Step 1: Write the failing tests**

Create `tests/OpenGJKSharp.Tests/SupportMappableQueryTests.cs`:

```csharp
using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class SupportMappableQueryTests
{
    [Fact]
    public void Should_ComputeSeparationDistance_When_SpheresAreSeparated()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(5, 0, 0), 2f);

        // Analytic: |c1 - c2| - r1 - r2 = 5 - 1 - 2 = 2
        double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

        Assert.Equal(2.0, distance, 4);
        Assert.Equal(1.0, closestA.X, 3);
        Assert.Equal(0.0, closestA.Y, 3);
        Assert.Equal(3.0, closestB.X, 3);
        Assert.Equal(0.0, closestB.Y, 3);
    }

    [Fact]
    public void Should_ReturnTrueOnlyForOverlap_When_SphereCollisionIsQueried()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var overlapping = new Sphere(new Vector3(1.5f, 0, 0), 1f);
        var separated = new Sphere(new Vector3(3f, 0, 0), 1f);

        Assert.True(OpenGjk.HasCollision(a, overlapping));
        Assert.False(OpenGjk.HasCollision(a, separated));
    }

    [Fact]
    public void Should_ReturnPenetrationDepthAndNormal_When_SpheresOverlap()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(1.5f, 0, 0), 1f);

        // Analytic: penetration = r1 + r2 - |c1 - c2| = 0.5, normal along +X (from a toward b)
        double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 normal);

        Assert.True(distance < 0, $"expected collision, got {distance}");
        Assert.True(Math.Abs(-distance - 0.5) < 0.01, $"expected depth 0.5, got {-distance}");
        Assert.Equal(1.0, normal.Length(), 3);
        Assert.True(normal.X > 0.99f, $"normal should point from a toward b, got {normal}");
    }

    [Fact]
    public void Should_ReturnPenetrationDepth_When_SpheresAreConcentric()
    {
        // Concentric shapes make the initial GJK direction (centerA - centerB) the
        // zero vector; the adapter must substitute a valid direction.
        var a = new Sphere(Vector3.Zero, 1f);
        var b = new Sphere(Vector3.Zero, 0.5f);

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out _);

        Assert.True(distance < 0, $"expected collision, got {distance}");
        Assert.True(Math.Abs(-distance - 1.5) < 0.08, $"expected depth ~1.5, got {-distance}");
    }

    [Fact]
    public void Should_ThrowArgumentNullException_When_ShapeIsNull()
    {
        var sphere = new Sphere(Vector3.Zero, 1f);

        Assert.Throws<ArgumentNullException>(() => OpenGjk.HasCollision(null!, sphere));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.HasCollision(sphere, (ISupportMappable)null!));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeMinimumDistance(null!, sphere));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeMinimumDistance(sphere, (ISupportMappable)null!, out _, out _));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeCollisionInformation((ISupportMappable)null!, sphere, out _));
    }

    [Fact]
    public void Should_ThrowArgumentOutOfRangeException_When_ShapePrecisionIsNegative()
    {
        var a = new Sphere(Vector3.Zero, 1f);
        var b = new Sphere(new Vector3(3, 0, 0), 1f);

        Assert.Throws<ArgumentOutOfRangeException>(() => OpenGjk.HasCollision(a, b, precision: -1));
    }
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `dotnet test tests/OpenGJKSharp.Tests --filter "FullyQualifiedName~SupportMappableQueryTests"`
Expected: compile FAILURE (no `ISupportMappable` overloads on `OpenGjk`).

- [ ] **Step 3: Write the implementation**

Create `src/OpenGJKSharp/Models/GkMappedBody.cs`:

```csharp
using System.Numerics;

namespace OpenGJKSharp.Models;

/// <summary>
/// Adapts a user-provided <see cref="ISupportMappable"/> shape to the internal body contract.
/// Every point the shape returns is appended to a growing store, and its store position acts
/// as the point id tracked by the simplex/EPA bookkeeping; <see cref="GetPoint"/> resolves ids
/// back to coordinates for witness reconstruction. The same instance must therefore be used
/// for the GJK and EPA passes of a single query.
/// </summary>
internal sealed class GkMappedBody : GkBody
{
    private readonly ISupportMappable _shape;
    private readonly List<double[]> _points = [];
    private int _sIdx;

    public GkMappedBody(ISupportMappable shape)
    {
        _shape = shape;
    }

    public override double[] S { get; } = new double[3];

    public override int SIdx => _sIdx;

    public override void InitSupport()
    {
        _sIdx = AddPoint(_shape.GetCenter(), S);
    }

    public override void Support(double[] v)
    {
        SupportInto(v, S, out _sIdx);
    }

    public override void SupportInto(double[] direction, double[] result, out int idx)
    {
        Vector3 point = _shape.GetSupportPoint(ToUnitDirection(direction));
        idx = AddPoint(point, result);
    }

    public override double[] GetPoint(int idx)
    {
        return _points[idx];
    }

    private int AddPoint(Vector3 point, double[] destination)
    {
        double[] coordinates = [point.X, point.Y, point.Z];
        _points.Add(coordinates);

        destination[0] = coordinates[0];
        destination[1] = coordinates[1];
        destination[2] = coordinates[2];

        return _points.Count - 1;
    }

    // Normalizes in double before converting to float so tiny directions survive the
    // conversion. The interface contract promises implementations a non-zero direction,
    // so the zero vector (e.g. concentric shapes make the initial direction
    // centerA - centerB vanish) is replaced with an arbitrary axis here.
    private static Vector3 ToUnitDirection(double[] direction)
    {
        double normSquared = direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2];
        if (normSquared <= 0)
        {
            return Vector3.UnitX;
        }

        double inverseNorm = 1.0 / Math.Sqrt(normSquared);
        return new Vector3(
            (float)(direction[0] * inverseNorm),
            (float)(direction[1] * inverseNorm),
            (float)(direction[2] * inverseNorm));
    }
}
```

Add to `OpenGjk.cs`, next to the existing array overloads (each mirrors its array counterpart):

```csharp
/// <summary>
/// Determines whether two support-mapped convex shapes collide.
/// </summary>
/// <param name="a">The first shape.</param>
/// <param name="b">The second shape.</param>
/// <param name="precision">Distance below which the shapes are considered colliding. Must not be negative.</param>
/// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
/// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
/// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
public static bool HasCollision(ISupportMappable a, ISupportMappable b, double precision = _precision)
{
    ValidateShapes(a, b);
    if (precision < 0)
    {
        throw new ArgumentOutOfRangeException(nameof(precision), precision, "Precision must not be negative.");
    }

    double distance = ComputeMinimumDistance(new GkMappedBody(a), new GkMappedBody(b), new GkSimplex());

    return distance <= precision;
}

/// <summary>
/// Computes the minimum distance between two support-mapped convex shapes.
/// Returns 0 when the shapes collide.
/// </summary>
/// <param name="a">The first shape.</param>
/// <param name="b">The second shape.</param>
/// <returns>The minimum distance between the two shapes.</returns>
/// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b)
{
    return ComputeMinimumDistance(a, b, out _, out _);
}

/// <summary>
/// Computes the minimum distance between two support-mapped convex shapes,
/// along with the closest point on each shape (witness points).
/// Returns 0 when the shapes collide.
/// </summary>
/// <param name="a">The first shape.</param>
/// <param name="b">The second shape.</param>
/// <param name="closestA">Closest point on the first shape.</param>
/// <param name="closestB">Closest point on the second shape.</param>
/// <returns>The minimum distance between the two shapes.</returns>
/// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b, out Vector3 closestA, out Vector3 closestB)
{
    ValidateShapes(a, b);

    var simplex = new GkSimplex();
    double distance = ComputeMinimumDistance(new GkMappedBody(a), new GkMappedBody(b), simplex);

    closestA = new Vector3(
        (float)simplex.Witnesses[0, 0],
        (float)simplex.Witnesses[0, 1],
        (float)simplex.Witnesses[0, 2]);
    closestB = new Vector3(
        (float)simplex.Witnesses[1, 0],
        (float)simplex.Witnesses[1, 1],
        (float)simplex.Witnesses[1, 2]);

    return distance;
}

/// <summary>
/// Computes collision information between two support-mapped convex shapes.
/// Returns the minimum distance when the shapes are separated (positive value),
/// or the penetration depth computed by the EPA algorithm when they collide (negative value).
/// </summary>
/// <param name="a">The first shape.</param>
/// <param name="b">The second shape.</param>
/// <param name="contactNormal">Contact normal, pointing from the first shape toward the second.</param>
/// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
/// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
public static double ComputeCollisionInformation(ISupportMappable a, ISupportMappable b, out Vector3 contactNormal)
{
    ValidateShapes(a, b);

    // The same body instances must serve both the GJK and EPA passes:
    // the simplex records ids into each body's point store.
    var bd1 = new GkMappedBody(a);
    var bd2 = new GkMappedBody(b);
    var simplex = new GkSimplex();

    double distance = ComputeMinimumDistance(bd1, bd2, simplex);

    double[] normal = new double[3];
    ComputeCollisionInformation(bd1, bd2, simplex, ref distance, normal);

    contactNormal = new Vector3((float)normal[0], (float)normal[1], (float)normal[2]);
    return distance;
}

private static void ValidateShapes(ISupportMappable a, ISupportMappable b)
{
    if (a is null)
    {
        throw new ArgumentNullException(nameof(a));
    }

    if (b is null)
    {
        throw new ArgumentNullException(nameof(b));
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `dotnet test tests/OpenGJKSharp.Tests`
Expected: PASS - the new SupportMappableQueryTests AND the full existing suite.

- [ ] **Step 5: Commit**

```bash
git add src/OpenGJKSharp/Models/GkMappedBody.cs src/OpenGJKSharp/OpenGjk.cs tests/OpenGJKSharp.Tests/SupportMappableQueryTests.cs
git commit -m "feat: add support-mapped shape overloads to OpenGjk"
```

---

### Task 7: Analytic capsule/box tests and parity regression tests

**Files:**
- Modify: `tests/OpenGJKSharp.Tests/SupportMappableQueryTests.cs` (add tests)

**Interfaces:**
- Consumes: everything from Tasks 1-6. No new production code expected; if a test exposes a defect, fix it in the smallest scope possible and note it in the commit body.

- [ ] **Step 1: Add the analytic and parity tests**

Add to `SupportMappableQueryTests.cs`:

```csharp
[Fact]
public void Should_ComputeSeparationDistance_When_ParallelCapsulesAreSeparated()
{
    // Parallel segments 2 apart, radii 0.25 each: distance = 2 - 0.5 = 1.5
    var a = new Capsule(new Vector3(0, 0, 0), new Vector3(1, 0, 0), 0.25f);
    var b = new Capsule(new Vector3(0, 2, 0), new Vector3(1, 2, 0), 0.25f);

    double distance = OpenGjk.ComputeMinimumDistance(a, b);

    Assert.Equal(1.5, distance, 4);
}

[Fact]
public void Should_ComputeSeparationDistance_When_SkewCapsulesAreSeparated()
{
    // Segment A along X through origin, segment B along Y at z = 2:
    // segment-segment distance = 2, minus radii 0.2 + 0.3 -> 1.5
    var a = new Capsule(new Vector3(-1, 0, 0), new Vector3(1, 0, 0), 0.2f);
    var b = new Capsule(new Vector3(0, -1, 2), new Vector3(0, 1, 2), 0.3f);

    double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

    Assert.Equal(1.5, distance, 4);
    Assert.Equal(distance, (closestB - closestA).Length(), 3);
}

[Fact]
public void Should_BehaveLikeSphere_When_CapsuleEndpointsCoincide()
{
    var point = new Vector3(1, 2, 3);
    var degenerate = new Capsule(point, point, 0.5f);
    var sphere = new Sphere(point, 0.5f);
    var probe = new Sphere(new Vector3(5, 2, 3), 1f);

    double capsuleDistance = OpenGjk.ComputeMinimumDistance(degenerate, probe);
    double sphereDistance = OpenGjk.ComputeMinimumDistance(sphere, probe);

    Assert.Equal(sphereDistance, capsuleDistance, 5);
}

[Fact]
public void Should_MatchVertexHull_When_RotatedBoxIsQueried()
{
    var orientation = Quaternion.CreateFromAxisAngle(Vector3.Normalize(new Vector3(1, 1, 0)), 0.7f);
    var boxA = new Box(new Vector3(0, 0, 0), new Vector3(1, 0.5f, 0.75f), orientation);
    var boxB = new Box(new Vector3(4, 0.25f, 0), new Vector3(1, 1, 1));

    double shapeDistance = OpenGjk.ComputeMinimumDistance(boxA, boxB);
    double hullDistance = OpenGjk.ComputeMinimumDistance(Corners(boxA), Corners(boxB));

    Assert.Equal(hullDistance, shapeDistance, 4);

    static Vector3[] Corners(Box box)
    {
        var corners = new List<Vector3>();
        foreach (float sx in (float[])[-1f, 1f])
        {
            foreach (float sy in (float[])[-1f, 1f])
            {
                foreach (float sz in (float[])[-1f, 1f])
                {
                    var local = new Vector3(
                        sx * box.HalfExtents.X,
                        sy * box.HalfExtents.Y,
                        sz * box.HalfExtents.Z);
                    corners.Add(box.Center + Vector3.Transform(local, box.Orientation));
                }
            }
        }

        return [.. corners];
    }
}

[Fact]
public void Should_MatchArrayRoute_When_ConvexHullShapeIsQueried()
{
    Vector3[] cubeA =
    [
        new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(1, 1, 0),
        new(0, 0, 1), new(1, 0, 1), new(0, 1, 1), new(1, 1, 1),
    ];
    Vector3[] cubeB = [.. cubeA.Select(p => p + new Vector3(3f, 0, 0))];

    double arrayDistance = OpenGjk.ComputeMinimumDistance(cubeA, cubeB, out Vector3 arrayClosestA, out Vector3 arrayClosestB);
    double shapeDistance = OpenGjk.ComputeMinimumDistance(new ConvexHull(cubeA), new ConvexHull(cubeB), out Vector3 shapeClosestA, out Vector3 shapeClosestB);

    Assert.Equal(arrayDistance, shapeDistance, 5);
    Assert.Equal(arrayClosestA.X, shapeClosestA.X, 4);
    Assert.Equal(arrayClosestB.X, shapeClosestB.X, 4);
    Assert.Equal(OpenGjk.HasCollision(cubeA, cubeB), OpenGjk.HasCollision(new ConvexHull(cubeA), new ConvexHull(cubeB)));
}

[Fact]
public void Should_MatchArrayRoute_For_RandomSeparatedClouds()
{
    var random = new Random(4242);

    for (int iteration = 0; iteration < 100; iteration++)
    {
        Vector3[] a = RandomCloud(random, new Vector3(0, 0, 0));
        Vector3[] b = RandomCloud(random, new Vector3(5, 0, 0));

        double arrayDistance = OpenGjk.ComputeMinimumDistance(a, b);
        double shapeDistance = OpenGjk.ComputeMinimumDistance(new ConvexHull(a), new ConvexHull(b));

        Assert.True(Math.Abs(arrayDistance - shapeDistance) < 1e-4,
            $"iteration {iteration}: array = {arrayDistance}, shape = {shapeDistance}");
    }

    static Vector3[] RandomCloud(Random random, Vector3 center) =>
        [.. Enumerable.Range(0, 20).Select(_ => center + new Vector3(
            (float)random.NextDouble(),
            (float)random.NextDouble(),
            (float)random.NextDouble()))];
}

[Fact]
public void Should_ReturnPenetrationDepth_When_SphereOverlapsBox()
{
    // Sphere r=1 centered at (2.5, 0, 0) vs box [-1, 1]^3:
    // gap between surfaces = 2.5 - 1 - 1 = 0.5 -> separated; move sphere to 1.5 -> penetration 0.5
    var box = new Box(Vector3.Zero, Vector3.One);
    var separated = new Sphere(new Vector3(2.5f, 0, 0), 1f);
    var overlapping = new Sphere(new Vector3(1.5f, 0, 0), 1f);

    double separation = OpenGjk.ComputeMinimumDistance(box, separated);
    double penetration = OpenGjk.ComputeCollisionInformation(box, overlapping, out Vector3 normal);

    Assert.Equal(0.5, separation, 4);
    Assert.True(penetration < 0, $"expected collision, got {penetration}");
    Assert.True(Math.Abs(-penetration - 0.5) < 0.01, $"expected depth 0.5, got {-penetration}");
    Assert.True(normal.X > 0.99f, $"normal should point from box toward sphere, got {normal}");
}
```

- [ ] **Step 2: Run the tests**

Run: `dotnet test tests/OpenGJKSharp.Tests`
Expected: PASS. If an analytic test fails by a hair (curved-shape EPA accuracy), first verify the analytic expectation by hand; only widen a tolerance when the computed value is geometrically correct, and say so in the commit body. If it fails badly, debug with superpowers:systematic-debugging before touching tolerances.

- [ ] **Step 3: Commit**

```bash
git add tests/OpenGJKSharp.Tests/SupportMappableQueryTests.cs
git commit -m "test: add analytic and parity tests for shapes"
```

---

### Task 8: Documentation and release prep (0.5.0)

**Files:**
- Modify: `README.md` (new "Shapes and custom support functions" section)
- Modify: `CHANGELOG.md` (0.5.0 entry, Keep a Changelog format - read the existing file and match its exact structure/link style)
- Modify: `src/OpenGJKSharp/OpenGJKSharp.csproj` (`<Version>0.4.0</Version>` -> `<Version>0.5.0</Version>`)

**Interfaces:**
- Consumes: final public API from Tasks 1-6.

- [ ] **Step 1: Add the README section**

Read `README.md` first and insert a section matching its heading style, after the existing usage examples:

```markdown
## Shapes and custom support functions

GJK is a support-function-based algorithm, so shapes like spheres and capsules can be
handled exactly - without tessellation - through the `ISupportMappable` interface.
Built-in shapes live in `OpenGJKSharp.Shapes`:

```csharp
using OpenGJKSharp;
using OpenGJKSharp.Shapes;
using System.Numerics;

var sphere = new Sphere(new Vector3(0, 0, 0), radius: 1f);
var box = new Box(new Vector3(3, 0, 0), halfExtents: new Vector3(1, 1, 1));

bool hit = OpenGjk.HasCollision(sphere, box);
double distance = OpenGjk.ComputeMinimumDistance(sphere, box, out Vector3 onSphere, out Vector3 onBox);
double info = OpenGjk.ComputeCollisionInformation(sphere, box, out Vector3 contactNormal);
```

`Sphere`, `Capsule`, `Box`, and `ConvexHull` are provided. You can plug in any custom
convex shape by implementing the interface:

```csharp
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

Implementations must describe a convex shape, and `GetSupportPoint` must be exact for
the shape - GJK/EPA accuracy depends on it. Note that the public interface uses
single-precision `Vector3`, so support points are float-quantized on input (the
vertex-array API has the same input quantization); internal math runs in double
precision.
```

- [ ] **Step 2: Update CHANGELOG and version**

In `CHANGELOG.md`, add a `## [0.5.0] - 2026-08-01` section (matching the file's existing conventions) with:

```markdown
### Added

- `ISupportMappable` support-mapping interface: GJK/EPA queries on implicit convex
  shapes without tessellation.
- Built-in shapes in `OpenGJKSharp.Shapes`: `Sphere`, `Capsule`, `Box` (oriented),
  `ConvexHull`.
- `OpenGjk` overloads for `ISupportMappable`: `HasCollision`,
  `ComputeMinimumDistance` (with and without witness points), and
  `ComputeCollisionInformation`.

### Changed

- Internal GJK/EPA loops now run on a support-function abstraction; the vertex-array
  and `GkPolytope` APIs are unchanged and keep their exact behavior.
```

In `src/OpenGJKSharp/OpenGJKSharp.csproj`, change `<Version>0.4.0</Version>` to `<Version>0.5.0</Version>`.

- [ ] **Step 3: Verify build and full suite**

Run: `dotnet build src/OpenGJKSharp -c Release` (all TFMs, no warnings) and `dotnet test tests/OpenGJKSharp.Tests`.
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add README.md CHANGELOG.md src/OpenGJKSharp/OpenGJKSharp.csproj
git commit -m "docs: document shapes API and prepare 0.5.0"
```

---

## Self-Review Notes

- Spec coverage: interface (Task 1), four shapes + validation table (Tasks 1-4), internal support abstraction + witness generalization (Task 5), overloads + error handling (Task 6), analytic/parity/validation tests (Tasks 1-4, 6, 7), README/CHANGELOG/version (Task 8). EPA duplicate detection: already coordinate-based in current code - verified against `OpenGjk.cs` (epsSq comparisons at the `isDuplicate` check and simplex-rebuild paths); no change needed.
- Termination on curved shapes: existing caps (`mk = 25`, `maxIterations = 64`) and tolerances unchanged, per spec.
- Types cross-checked: `GkBody` members used in Task 6 match Task 5 definitions; shape constructors used in Tasks 6-7 match Tasks 1-4.

using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

/// <summary>
/// Penetration depth for pairs whose Minkowski difference is flat (coplanar bodies) or
/// collinear (bodies on one line), where a tetrahedron cannot be built and EPA used to
/// report a depth of 0.
/// </summary>
public class PlanarPenetrationTests
{
    private static readonly Vector3[] _quadA = [new(0, 0, 0), new(1, 0, 0), new(1, 1, 0), new(0, 1, 0)];
    private static readonly Vector3[] _quadB = [new(0.75f, 0, 0), new(1.75f, 0, 0), new(1.75f, 1, 0), new(0.75f, 1, 0)];

    [Fact]
    public void Should_ReportInPlaneDepth_When_CoplanarQuadsOverlap()
    {
        double info = OpenGjk.ComputeCollisionInformation(_quadA, _quadB, out Vector3 normal);

        Assert.Equal(-0.25, info, 6);

        // The normal lies in the shared plane and points from the first body toward the second.
        Assert.Equal(1.0, normal.X, 6);
        Assert.Equal(0.0, normal.Y, 6);
        Assert.Equal(0.0, normal.Z, 6);
        Assert.Equal(1.0, normal.Length(), 6);
    }

    [Fact]
    public void Should_ReportInPlaneDepth_When_FlatConvexHullsOverlap()
    {
        var a = new ConvexHull(_quadA);
        var b = new ConvexHull(_quadB);

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(0.25, depth, 6);
        Assert.Equal(1.0, normal.X, 6);
        Assert.Equal(0.0, normal.Z, 6);
    }

    [Fact]
    public void Should_ReportInPlaneDepth_When_SharedPlaneIsNotAxisAligned()
    {
        // The same pair in the plane spanned by u = (1, 1, 0) and v = (0, 0, 1). The basis is
        // built from exact float values on purpose: rotating a flat pair by an arbitrary angle
        // quantizes the vertices and generally leaves the two bodies in slightly different
        // planes, which is a three-dimensional configuration, not a coplanar one.
        var u = new Vector3(1, 1, 0);
        var v = new Vector3(0, 0, 1);
        Vector3[] a = [Vector3.Zero, u, u + v, v];
        Vector3[] b = [0.75f * u, 1.75f * u, (1.75f * u) + v, (0.75f * u) + v];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));

        // The overlap is a quarter of the edge, whose world length is |u| = sqrt(2).
        Assert.Equal(0.25 * Math.Sqrt(2), depth, 6);

        Vector3 expected = Vector3.Normalize(u);
        Assert.Equal(expected.X, normal.X, 5);
        Assert.Equal(expected.Y, normal.Y, 5);
        Assert.Equal(expected.Z, normal.Z, 5);
    }

    [Fact]
    public void Should_StayThreeDimensional_When_AFlatPairIsQuantizedOutOfItsPlane()
    {
        // Documents the limit of the planar path with single-precision input: rotating a
        // coplanar pair by an arbitrary angle moves each vertex to the nearest float, so the
        // two bodies no longer share a plane and no longer intersect at all. This is honest
        // geometry rather than a missed collision - two flat shapes in different planes only
        // meet along the line where the planes cross.
        var rotation = Quaternion.CreateFromAxisAngle(Vector3.Normalize(new Vector3(1, 2, 3)), 0.7f);
        Vector3[] a = [.. _quadA.Select(v => Vector3.Transform(v, rotation))];
        Vector3[] b = [.. _quadB.Select(v => Vector3.Transform(v, rotation))];

        double distance = OpenGjk.ComputeMinimumDistance(a, b);

        // Whatever the sign, the pair is resolved as a 3D one and the answer stays tiny.
        Assert.True(Math.Abs(distance) < 1e-5, $"expected a near-zero result, got {distance}");
    }

    [Theory]
    [InlineData(0.1)]
    [InlineData(0.25)]
    [InlineData(0.5)]
    [InlineData(0.75)]
    [InlineData(0.9)]
    public void Should_MatchAnalyticDepth_When_CoplanarSquaresOverlapByAKnownAmount(double offset)
    {
        // Unit squares on z = 0, the second shifted along x: the minimum in-plane
        // translation is the remaining x overlap.
        float t = (float)offset;
        Vector3[] a = [new(0, 0, 0), new(1, 0, 0), new(1, 1, 0), new(0, 1, 0)];
        Vector3[] b = [new(t, 0, 0), new(t + 1, 0, 0), new(t + 1, 1, 0), new(t, 1, 0)];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(1 - offset, depth, 6);
        Assert.Equal(1.0, normal.X, 6);
    }

    [Fact]
    public void Should_PickTheSmallerAxis_When_CoplanarSquaresOverlapInBothDirections()
    {
        // The x overlap is 0.1 and the y overlap 0.05, so the minimum translation is along y.
        // An answer taken from the x axis alone would report the wrong pair of values.
        Vector3[] a = [new(0, 0, 0), new(1, 0, 0), new(1, 1, 0), new(0, 1, 0)];
        Vector3[] b = [new(0.9f, 0.95f, 0), new(1.9f, 0.95f, 0), new(1.9f, 1.95f, 0), new(0.9f, 1.95f, 0)];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(0.05, depth, 6);
        Assert.Equal(0.0, normal.X, 6);
        Assert.Equal(1.0, normal.Y, 6);
    }

    [Fact]
    public void Should_ReportAnEdgeNormal_When_CoplanarTrianglesOverlapDiagonally()
    {
        // The closest edge of the Minkowski difference is the hypotenuse, so the normal is
        // diagonal: it cannot come from the coordinate axes or from a single simplex edge.
        Vector2[] a = [new(0, 0), new(1, 0), new(0, 1)];
        Vector2[] b = [new(0.9f, 0.9f), new(-0.1f, 0.9f), new(0.9f, -0.1f)];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector2 normal));

        // The first triangle is x + y <= 1 and the second x + y >= 0.8, so they overlap in a
        // band 0.2 / sqrt(2) thick, and the second clears it by moving along +(1, 1).
        double sqrt2 = Math.Sqrt(2);
        Assert.Equal(0.2 / sqrt2, depth, 5);
        Assert.Equal(1 / sqrt2, normal.X, 5);
        Assert.Equal(1 / sqrt2, normal.Y, 5);
    }

    [Theory]
    [InlineData(1e-3f)]
    [InlineData(1e-4f)]
    [InlineData(1e-6f)]
    public void Should_KeepThreeDimensionalResult_When_PairIsThinButNotFlat(float thickness)
    {
        // The regression risk of the planar path: a pair that is thin but genuinely
        // three-dimensional must keep its out-of-plane answer, which is the smaller
        // translation and therefore the correct one.
        Vector3[] a = [.. _quadA, .. _quadA.Select(v => v + new Vector3(0, 0, thickness))];
        Vector3[] b = [.. _quadB, .. _quadB.Select(v => v + new Vector3(0, 0, thickness))];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(thickness, depth, 1e-9);
        Assert.Equal(1.0, Math.Abs(normal.Z), 4);
    }

    [Fact]
    public void Should_ReportInLineDepth_When_CollinearSegmentsOverlap()
    {
        Vector3[] a = [new(0, 0, 0), new(2, 0, 0)];
        Vector3[] b = [new(1, 0, 0), new(3, 0, 0)];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(1.0, depth, 6);
        Assert.Equal(1.0, normal.X, 6);
        Assert.Equal(1.0, normal.Length(), 6);
    }

    [Fact]
    public void Should_ReportInLineDepth_When_CollinearSegmentContainsTheOther()
    {
        Vector3[] a = [new(0, 0, 0), new(10, 0, 0)];
        Vector3[] b = [new(4, 0, 0), new(6, 0, 0)];

        // To clear the longer segment the shorter one has to travel its whole length: its left
        // end (4) must reach the right end of the other (10), or its right end (6) the left
        // end (0). Both are 6.
        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(6.0, depth, 6);
        Assert.Equal(1.0, Math.Abs(normal.X), 6);
    }

    [Fact]
    public void Should_ReportZeroDepth_When_BodiesAreCoincidentPoints()
    {
        // A zero-dimensional Minkowski difference: the bodies touch, and 0 is the
        // correct depth. Unchanged by the planar path.
        Vector3[] point = [new(1, 2, 3)];

        Assert.True(OpenGjk.TryGetPenetration(point, point, out double depth, out _));
        Assert.Equal(0.0, depth, 12);
        Assert.False(Numbers.IsNegative(depth));
    }

    [Fact]
    public void Should_SeparateBodies_When_InPlaneDepthIsApplied()
    {
        // The property that matters: applying the reported translation must clear the overlap.
        var random = new Random(20260803);

        for (int i = 0; i < 200; i++)
        {
            Vector2[] a = RandomConvexPolygon(random, new Vector2(0, 0));
            Vector2[] b = RandomConvexPolygon(random, new Vector2(
                (float)(random.NextDouble() * 2 - 1),
                (float)(random.NextDouble() * 2 - 1)));

            if (!OpenGjk.TryGetPenetration(a, b, out double depth, out Vector2 normal))
            {
                continue;
            }

            Vector2[] moved = [.. b.Select(v => v + (normal * (float)depth))];
            double remaining = OpenGjk.ComputeMinimumDistance(a, moved);

            Assert.True(remaining >= -1e-4,
                $"still overlapping by {-remaining} after moving by {depth} along {normal}");
        }
    }

    [Fact]
    public void Should_ReportInPlaneDepth_When_PolygonsAreGivenAsVector2()
    {
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [new(0.75f, 0), new(1.75f, 0), new(1.75f, 1), new(0.75f, 1)];

        double info = OpenGjk.ComputeCollisionInformation(a, b, out Vector2 normal);

        Assert.Equal(-0.25, info, 6);
        Assert.Equal(1.0, normal.X, 6);
        Assert.Equal(0.0, normal.Y, 6);

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector2 penetrationNormal));
        Assert.Equal(0.25, depth, 6);
        Assert.Equal(normal, penetrationNormal);
    }

    [Fact]
    public void Should_ReportSeparation_When_PolygonsGivenAsVector2AreApart()
    {
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [new(3, 0), new(4, 0), new(4, 1), new(3, 1)];

        double info = OpenGjk.ComputeCollisionInformation(a, b, out Vector2 normal);

        Assert.Equal(2.0, info, 6);
        Assert.Equal(1.0, normal.X, 6);

        Assert.False(OpenGjk.TryGetPenetration(a, b, out double depth, out _));
        Assert.Equal(0.0, depth);
    }

    [Fact]
    public void Should_ReportCollision_When_PolygonsGivenAsVector2TouchExactly()
    {
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [new(1, 0), new(2, 0), new(2, 1), new(1, 1)];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out _));
        Assert.Equal(0.0, depth, 6);
        Assert.False(Numbers.IsNegative(depth));
    }

    [Fact]
    public void Should_MatchAllocatingOverload_When_PolygonWorkspaceIsReused()
    {
        var workspace = new GjkWorkspace();
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [new(0.75f, 0), new(1.75f, 0), new(1.75f, 1), new(0.75f, 1)];

        for (int i = 0; i < 5; i++)
        {
            double expected = OpenGjk.ComputeCollisionInformation(a, b, out Vector2 expectedNormal);
            double actual = OpenGjk.ComputeCollisionInformation(a, b, workspace, out Vector2 actualNormal);

            Assert.Equal(expected, actual);
            Assert.Equal(expectedNormal, actualNormal);

            Assert.Equal(
                OpenGjk.TryGetPenetration(a, b, out double expectedDepth, out Vector2 expectedTryNormal),
                OpenGjk.TryGetPenetration(a, b, workspace, out double actualDepth, out Vector2 actualTryNormal));
            Assert.Equal(expectedDepth, actualDepth);
            Assert.Equal(expectedTryNormal, actualTryNormal);
        }
    }

    [Fact]
    public void Should_Throw_When_PolygonArgumentsAreInvalid()
    {
        Vector2[] polygon = [new(0, 0), new(1, 0), new(0, 1)];

        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeCollisionInformation(null!, polygon, out Vector2 _));
        Assert.Throws<ArgumentNullException>(
            () => OpenGjk.ComputeCollisionInformation(polygon, polygon, null!, out Vector2 _));
        Assert.Throws<ArgumentException>(() => OpenGjk.TryGetPenetration(polygon, [], out _, out Vector2 _));
    }

#if NET8_0_OR_GREATER
    [Fact]
    public void Should_StopAllocating_When_PlanarQueriesReuseAWorkspace()
    {
        var workspace = new GjkWorkspace();
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [new(0.75f, 0), new(1.75f, 0), new(1.75f, 1), new(0.75f, 1)];

        for (int i = 0; i < 20; i++)
        {
            OpenGjk.ComputeCollisionInformation(a, b, workspace, out Vector2 _);
            OpenGjk.ComputeCollisionInformation(_quadA, _quadB, workspace, out Vector3 _);
        }

        long before = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < 100; i++)
        {
            OpenGjk.ComputeCollisionInformation(a, b, workspace, out Vector2 _);
            OpenGjk.ComputeCollisionInformation(_quadA, _quadB, workspace, out Vector3 _);
        }
        long bytes = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.True(bytes == 0, $"expected no allocation with a reused workspace, got {bytes} bytes");
    }
#endif

    private static Vector2[] RandomConvexPolygon(Random random, Vector2 center)
    {
        // Points on a circle in angle order are convex by construction.
        int count = 3 + random.Next(5);
        float radius = 0.4f + (float)random.NextDouble();
        var angles = new double[count];
        for (int i = 0; i < count; i++)
        {
            angles[i] = random.NextDouble() * Math.PI * 2;
        }

        Array.Sort(angles);

        var points = new Vector2[count];
        for (int i = 0; i < count; i++)
        {
            points[i] = center + new Vector2(
                radius * (float)Math.Cos(angles[i]),
                radius * (float)Math.Sin(angles[i]));
        }

        return points;
    }
}

using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class TryGetPenetrationTests
{
    [Fact]
    public void Should_ReportPositiveDepth_When_ShapesOverlap()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(1.5f, 0, 0), 1f);

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(0.5, depth, 3);
        Assert.Equal(1.0, normal.X, 3);
    }

    [Fact]
    public void Should_ReportZeroDepth_When_ShapesAreApart()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(5, 0, 0), 1f);

        Assert.False(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal));
        Assert.Equal(0.0, depth);

        // The normal still points from the first shape toward the second.
        Assert.Equal(1.0, normal.X, 3);
    }

    [Fact]
    public void Should_ReportCollision_When_ShapesTouchExactly()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(2, 0, 0), 1f);

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out _));
        Assert.Equal(0.0, depth, 6);

        // Positive zero, not negative zero: the latter formats as "-0" and flips the sign of
        // anything divided by it.
        Assert.False(double.IsNegative(depth));
    }

    [Fact]
    public void Should_SeparateBodies_When_DepthIsAppliedAlongTheNormal()
    {
        // The reason the overload exists: applying the result must resolve the overlap
        // without the caller having to reason about the sign.
        var random = new Random(20260808);
        var ground = new Box(new Vector3(0, -1, 0), new Vector3(10, 1, 10));

        for (int i = 0; i < 100; i++)
        {
            var start = new Vector3(
                (float)(random.NextDouble() * 4 - 2),
                (float)(random.NextDouble() * 0.9),
                (float)(random.NextDouble() * 4 - 2));
            var body = new TransformedShape(new Capsule(new Vector3(0, -0.4f, 0), new Vector3(0, 0.4f, 0), 0.3f),
                Quaternion.Identity, start);

            if (!OpenGjk.TryGetPenetration(body, ground, out double depth, out Vector3 normal))
            {
                continue;
            }

            // Moving the first body against the normal by the depth must clear the overlap.
            body.SetPose(Quaternion.Identity, start - (normal * (float)depth));

            double remaining = OpenGjk.ComputeMinimumDistance(body, ground);
            Assert.True(remaining >= -1e-4, $"still overlapping by {-remaining} after depenetration");
        }
    }

    [Fact]
    public void Should_MatchComputeCollisionInformation_When_BodiesAreVertexArrays()
    {
        Vector3[] a = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(1, 1, 0), new(0, 0, 1), new(1, 0, 1), new(0, 1, 1), new(1, 1, 1)];
        Vector3[] b = [new(0.75f, 0, 0), new(1.75f, 0, 0), new(0.75f, 1, 0), new(1.75f, 1, 0), new(0.75f, 0, 1), new(1.75f, 0, 1), new(0.75f, 1, 1), new(1.75f, 1, 1)];

        double info = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 expectedNormal);
        bool colliding = OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal);

        Assert.True(colliding);
        Assert.Equal(-info, depth);
        Assert.Equal(expectedNormal, normal);
        Assert.Equal(0.25, depth, 5);
    }

    [Fact]
    public void Should_MatchAllocatingOverload_When_WorkspaceIsReused()
    {
        var workspace = new GjkWorkspace();
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(1.2f, 0, 0), 1f);
        Vector3[] cloudA = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(0, 0, 1)];
        Vector3[] cloudB = [new(0.5f, 0, 0), new(1.5f, 0, 0), new(0.5f, 1, 0), new(0.5f, 0, 1)];

        for (int i = 0; i < 5; i++)
        {
            bool expected = OpenGjk.TryGetPenetration(a, b, out double expectedDepth, out Vector3 expectedNormal);
            bool actual = OpenGjk.TryGetPenetration(a, b, workspace, out double actualDepth, out Vector3 actualNormal);

            Assert.Equal(expected, actual);
            Assert.Equal(expectedDepth, actualDepth);
            Assert.Equal(expectedNormal, actualNormal);

            Assert.Equal(
                OpenGjk.TryGetPenetration(cloudA, cloudB, out double cloudDepth, out _),
                OpenGjk.TryGetPenetration(cloudA, cloudB, workspace, out double cloudWorkspaceDepth, out _));
            Assert.Equal(cloudDepth, cloudWorkspaceDepth);
        }
    }

    [Fact]
    public void Should_Throw_When_ArgumentsAreInvalid()
    {
        var sphere = new Sphere(Vector3.Zero, 1f);
        Vector3[] cloud = [new(0, 0, 0)];

        Assert.Throws<ArgumentNullException>(() => OpenGjk.TryGetPenetration(null!, sphere, out _, out _));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.TryGetPenetration(sphere, sphere, null!, out _, out _));
        Assert.Throws<ArgumentException>(() => OpenGjk.TryGetPenetration(cloud, [], out _, out _));
    }

    [Fact]
    public void Should_MatchComputeCollisionInformation_When_BodiesAreDoubleCoordinates()
    {
        double[] a = [0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1];
        double[] b = [0.75, 0, 0, 1.75, 0, 0, 0.75, 1, 0, 1.75, 1, 0, 0.75, 0, 1, 1.75, 0, 1, 0.75, 1, 1, 1.75, 1, 1];

        Span<double> expectedNormal = stackalloc double[3];
        double info = OpenGjk.ComputeCollisionInformation(a, b, expectedNormal);

        Span<double> normal = stackalloc double[3];
        bool colliding = OpenGjk.TryGetPenetration(a, b, normal, out double depth);

        Assert.True(colliding);
        Assert.Equal(-info, depth);
        Assert.Equal(0.25, depth, 12);
        for (int t = 0; t < 3; t++)
        {
            Assert.Equal(expectedNormal[t], normal[t]);
        }
    }

    [Fact]
    public void Should_ReportZeroDepth_When_DoubleCoordinatesAreApart()
    {
        double[] a = [0, 0, 0, 1, 0, 0, 0, 1, 0];
        double[] b = [5, 0, 0, 6, 0, 0, 5, 1, 0];

        Span<double> normal = stackalloc double[3];
        Assert.False(OpenGjk.TryGetPenetration(a, b, normal, out double depth));
        Assert.Equal(0.0, depth);

        // The normal still points from the first body toward the second.
        Assert.Equal(1.0, normal[0], 6);
    }

    [Fact]
    public void Should_NotReportNegativeZero_When_DoubleCoordinatesTouchExactly()
    {
        double[] a = [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1];
        double[] b = [1, 0, 0, 2, 0, 0, 1, 1, 0, 1, 0, 1];

        Span<double> normal = stackalloc double[3];
        Assert.True(OpenGjk.TryGetPenetration(a, b, normal, out double depth));
        Assert.Equal(0.0, depth, 12);
        Assert.False(double.IsNegative(depth));
    }

    [Fact]
    public void Should_MatchAllocatingOverload_When_DoubleWorkspaceIsReused()
    {
        var workspace = new GjkWorkspace();
        double[] a = [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1];
        double[] b = [0.5, 0, 0, 1.5, 0, 0, 0.5, 1, 0, 0.5, 0, 1];

        for (int i = 0; i < 5; i++)
        {
            Span<double> expectedNormal = stackalloc double[3];
            bool expected = OpenGjk.TryGetPenetration(a, b, expectedNormal, out double expectedDepth);

            Span<double> normal = stackalloc double[3];
            bool actual = OpenGjk.TryGetPenetration(a, b, workspace, normal, out double depth);

            Assert.Equal(expected, actual);
            Assert.Equal(expectedDepth, depth);
            for (int t = 0; t < 3; t++)
            {
                Assert.Equal(expectedNormal[t], normal[t]);
            }
        }
    }

    [Fact]
    public void Should_Throw_When_DoubleArgumentsAreInvalid()
    {
        double[] coords = [0, 0, 0];
        double[] normal = new double[3];

        Assert.Throws<ArgumentNullException>(
            () => OpenGjk.TryGetPenetration(coords, coords, null!, normal, out _));
        Assert.Throws<ArgumentException>(() => OpenGjk.TryGetPenetration(coords, [], normal, out _));
        Assert.Throws<ArgumentException>(() => OpenGjk.TryGetPenetration(coords, coords, new double[2], out _));
    }

#if NET8_0_OR_GREATER
    [Fact]
    public void Should_StopAllocating_When_DoubleWorkspaceIsReused()
    {
        var workspace = new GjkWorkspace();
        double[] a = [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1];
        double[] b = [0.5, 0, 0, 1.5, 0, 0, 0.5, 1, 0, 0.5, 0, 1];
        double[] normal = new double[3];

        for (int i = 0; i < 20; i++)
        {
            OpenGjk.TryGetPenetration(a, b, workspace, normal, out _);
        }

        long before = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < 100; i++)
        {
            OpenGjk.TryGetPenetration(a, b, workspace, normal, out _);
        }
        long bytes = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.True(bytes == 0, $"expected no allocation with a reused workspace, got {bytes} bytes");
    }
#endif
}

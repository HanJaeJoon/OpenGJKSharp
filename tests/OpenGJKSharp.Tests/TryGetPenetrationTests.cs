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
}

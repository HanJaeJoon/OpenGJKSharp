using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class EpaAccuracyTests
{
    [Theory]
    [InlineData(0.1)]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(1.5)]
    [InlineData(1.7)]
    [InlineData(1.9)]
    [InlineData(1.99)]
    public void Should_ReportAccurateDepth_When_SpheresOverlapDeeply(double overlap)
    {
        // Two unit spheres whose centers are (2 - overlap) apart: the analytic penetration
        // depth is exactly the overlap. This used to be underestimated by up to 5.6% because
        // the polytope can never resolve a curved surface within its iteration budget.
        var a = new Sphere(Vector3.Zero, 1f);
        var b = new Sphere(new Vector3((float)(2.0 - overlap), 0, 0), 1f);

        double depth = -OpenGjk.ComputeCollisionInformation(a, b, out Vector3 normal);

        Assert.True(Math.Abs(depth - overlap) / overlap < 0.01,
            $"expected depth {overlap}, got {depth} ({(depth - overlap) / overlap * 100:F2}%)");
        Assert.Equal(1.0, normal.Length(), 4);
    }

    [Fact]
    public void Should_ReportExactDepth_When_SpheresAreConcentric()
    {
        // The deepest possible case: identical concentric spheres, depth = 2 * radius.
        Assert.Equal(2.0, -OpenGjk.ComputeCollisionInformation(
            new Sphere(Vector3.Zero, 1f), new Sphere(Vector3.Zero, 1f), out _), 4);

        // A small sphere fully inside a large one: depth = r1 + r2.
        Assert.Equal(3.0, -OpenGjk.ComputeCollisionInformation(
            new Sphere(Vector3.Zero, 1f), new Sphere(Vector3.Zero, 2f), out _), 4);
    }

    [Theory]
    [InlineData(0.1)]
    [InlineData(1.0)]
    [InlineData(1.9)]
    [InlineData(2.0)]
    public void Should_StayExact_When_BoxesOverlap(double overlap)
    {
        // A polytope is resolved exactly by EPA, so this must not regress at all.
        var a = new Box(Vector3.Zero, new Vector3(1, 1, 1));
        var b = new Box(new Vector3((float)(2.0 - overlap), 0, 0), new Vector3(1, 1, 1));

        Assert.Equal(overlap, -OpenGjk.ComputeCollisionInformation(a, b, out _), 5);
    }

    [Fact]
    public void Should_NeverUnderestimate_When_SphereDepthIsKnownAnalytically()
    {
        // The reported depth is an upper bound on the true one: a solver that pushes bodies
        // apart by it always resolves the penetration, which underestimating never guarantees.
        var random = new Random(20260807);

        for (int i = 0; i < 300; i++)
        {
            float radiusA = 0.5f + (float)random.NextDouble() * 2;
            float radiusB = 0.5f + (float)random.NextDouble() * 2;
            var direction = new Vector3(
                (float)(random.NextDouble() * 2 - 1),
                (float)(random.NextDouble() * 2 - 1),
                (float)(random.NextDouble() * 2 - 1));
            if (direction.LengthSquared() <= 0)
            {
                direction = Vector3.UnitX;
            }

            // Center distance somewhere inside the overlapping range.
            float centerDistance = (float)(random.NextDouble() * (radiusA + radiusB));
            Vector3 center = Vector3.Normalize(direction) * centerDistance;

            var a = new Sphere(Vector3.Zero, radiusA);
            var b = new Sphere(center, radiusB);
            double expected = radiusA + radiusB - (b.Center - a.Center).Length();

            double depth = -OpenGjk.ComputeCollisionInformation(a, b, out _);

            Assert.True(depth >= expected - 1e-4,
                $"depth {depth} underestimates the analytic {expected}");
            Assert.True(depth <= expected * 1.02 + 1e-4,
                $"depth {depth} overestimates the analytic {expected}");
        }
    }

    [Fact]
    public void Should_KeepDepthConsistentWithShapes_When_CapsulesCross()
    {
        // Crossed capsules: axes 1 apart, radii 0.5 + 0.5, so the depth is exactly 1.
        var a = new Capsule(new Vector3(-2, 0, 0), new Vector3(2, 0, 0), 0.5f);
        var b = new Capsule(new Vector3(0, -2, 0), new Vector3(0, 2, 0), 0.5f);

        Assert.Equal(1.0, -OpenGjk.ComputeCollisionInformation(a, b, out _), 4);
    }

    [Fact]
    public void Should_MatchSeparateQueries_When_WorkspaceIsReusedForDeepOverlaps()
    {
        // The best-plane bookkeeping lives on the reused EPA buffers, so it must not leak
        // between queries.
        var workspace = new GjkWorkspace();
        var shallow = new Sphere(new Vector3(1.9f, 0, 0), 1f);
        var deep = new Sphere(new Vector3(0.1f, 0, 0), 1f);
        var a = new Sphere(Vector3.Zero, 1f);

        for (int i = 0; i < 5; i++)
        {
            Assert.Equal(
                OpenGjk.ComputeCollisionInformation(a, deep, out _),
                OpenGjk.ComputeCollisionInformation(a, deep, workspace, out _));
            Assert.Equal(
                OpenGjk.ComputeCollisionInformation(a, shallow, out _),
                OpenGjk.ComputeCollisionInformation(a, shallow, workspace, out _));
        }
    }
}

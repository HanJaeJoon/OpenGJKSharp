using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

/// <summary>
/// Penetration depth for bodies whose interior points coincide. GJK used to seed its simplex
/// with the difference of the two <see cref="ISupportMappable.GetCenter"/> results, which is
/// an interior point of the Minkowski difference - the origin itself when the centers match.
/// EPA inherited it as a polytope vertex, so a face passed through the origin and the
/// closest-face search never got off distance zero.
/// </summary>
public class CrossedShapeDepthTests
{
    [Theory]
    [InlineData(0.25)]
    [InlineData(0.5)]
    [InlineData(0.75)]
    [InlineData(1.0)]
    [InlineData(1.5)]
    [InlineData(2.0)]
    public void Should_ReportDepth_When_CapsuleAxesCrossAtTheOrigin(double radius)
    {
        // Two perpendicular capsules through the origin: separating them means moving one by
        // 2 * radius along z, whatever their radius.
        var a = new Capsule(new Vector3(-2, 0, 0), new Vector3(2, 0, 0), (float)radius);
        var b = new Capsule(new Vector3(0, -2, 0), new Vector3(0, 2, 0), (float)radius);

        double depth = -OpenGjk.ComputeCollisionInformation(a, b, out Vector3 normal);

        Assert.Equal(2 * radius, depth, 3);
        Assert.Equal(1.0, Math.Abs(normal.Z), 3);
    }

    [Fact]
    public void Should_ReportDepth_When_OneBodyContainsTheOtherAroundACommonCenter()
    {
        Assert.Equal(3.0, -OpenGjk.ComputeCollisionInformation(
            new Sphere(Vector3.Zero, 1f), new Sphere(Vector3.Zero, 2f), out _), 3);

        Assert.Equal(4.0, -OpenGjk.ComputeCollisionInformation(
            new Box(Vector3.Zero, new Vector3(1, 1, 1)), new Box(Vector3.Zero, new Vector3(3, 3, 3)), out _), 3);

        // Mixed shape kinds, still sharing a center.
        Assert.Equal(2.5, -OpenGjk.ComputeCollisionInformation(
            new Sphere(Vector3.Zero, 0.5f), new Box(Vector3.Zero, new Vector3(2, 2, 2)), out _), 3);
    }

    [Fact]
    public void Should_ReportDepth_When_CentersCoincideForCrossedBoxes()
    {
        // Long thin boxes crossing at the origin: 2 * 0.5 along z separates them.
        var a = new Box(Vector3.Zero, new Vector3(3, 0.5f, 0.5f));
        var b = new Box(Vector3.Zero, new Vector3(0.5f, 3, 0.5f));

        double depth = -OpenGjk.ComputeCollisionInformation(a, b, out Vector3 normal);

        Assert.Equal(1.0, depth, 3);
        Assert.Equal(1.0, Math.Abs(normal.Z), 3);
    }

    [Fact]
    public void Should_ReportDepth_When_CustomShapesShareAnInteriorPoint()
    {
        // A custom shape returning the same center as the other body used to hit the same
        // stall; the fix must not depend on the shape being a built-in.
        var a = new CenteredHull(new Vector3(3, 0.5f, 0.5f));
        var b = new CenteredHull(new Vector3(0.5f, 3, 0.5f));

        double depth = -OpenGjk.ComputeCollisionInformation(a, b, out Vector3 normal);

        Assert.Equal(1.0, depth, 3);
        Assert.Equal(1.0, Math.Abs(normal.Z), 3);
    }

    [Fact]
    public void Should_StayConsistent_When_TheSamePairIsQueriedBothWays()
    {
        var a = new Capsule(new Vector3(-2, 0, 0), new Vector3(2, 0, 0), 1f);
        var b = new Capsule(new Vector3(0, -2, 0), new Vector3(0, 2, 0), 1f);

        double forward = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 forwardNormal);
        double backward = OpenGjk.ComputeCollisionInformation(b, a, out Vector3 backwardNormal);

        Assert.Equal(forward, backward, 3);
        Assert.Equal(Math.Abs(forwardNormal.Z), Math.Abs(backwardNormal.Z), 3);
    }

    // An axis-aligned box built from its 8 corners, centered on the origin.
    private sealed class CenteredHull : ISupportMappable
    {
        private readonly Vector3 _halfExtents;

        public CenteredHull(Vector3 halfExtents) => _halfExtents = halfExtents;

        public Vector3 GetSupportPoint(Vector3 direction)
        {
            return new Vector3(
                direction.X >= 0 ? _halfExtents.X : -_halfExtents.X,
                direction.Y >= 0 ? _halfExtents.Y : -_halfExtents.Y,
                direction.Z >= 0 ? _halfExtents.Z : -_halfExtents.Z);
        }

        public Vector3 GetCenter() => Vector3.Zero;
    }
}

using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

/// <summary>
/// Pins the policy for degenerate shapes: a thickness parameter must be positive, while
/// geometry given as explicit points may be degenerate.
/// </summary>
public class DegenerateShapeTests
{
    [Fact]
    public void Should_RejectZeroThickness_When_ShapeIsDefinedByAThicknessParameter()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new Sphere(Vector3.Zero, 0f));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Capsule(Vector3.Zero, Vector3.UnitX, 0f));
        Assert.Throws<ArgumentException>(() => new Box(Vector3.Zero, new Vector3(1, 0, 1)));

        // Negative values are the same mistake and rejected the same way.
        Assert.Throws<ArgumentOutOfRangeException>(() => new Sphere(Vector3.Zero, -1f));
        Assert.Throws<ArgumentException>(() => new Box(Vector3.Zero, new Vector3(1, -1, 1)));
    }

    [Fact]
    public void Should_BehaveAsSphere_When_CapsuleSegmentHasZeroLength()
    {
        // Not a degenerate shape: the radius is positive, so this is exactly a sphere.
        var capsule = new Capsule(new Vector3(2, 0, 0), new Vector3(2, 0, 0), 0.5f);
        var sphere = new Sphere(new Vector3(2, 0, 0), 0.5f);
        var probe = new Sphere(new Vector3(-3, 0, 0), 1f);

        Assert.Equal(
            OpenGjk.ComputeMinimumDistance(sphere, probe),
            OpenGjk.ComputeMinimumDistance(capsule, probe),
            6);
    }

    [Fact]
    public void Should_QueryExactly_When_HullIsASinglePoint()
    {
        var point = new ConvexHull([new Vector3(4, 0, 0)]);
        var sphere = new Sphere(Vector3.Zero, 1f);

        Assert.Equal(3.0, OpenGjk.ComputeMinimumDistance(point, sphere), 5);
        Assert.Equal(0.0, OpenGjk.ComputeMinimumDistance(point, new ConvexHull([new Vector3(4, 0, 0)])), 5);
    }

    [Fact]
    public void Should_QueryExactly_When_HullIsASegment()
    {
        // Segment along x at y = 3, probed by a segment along y: the gap is 3 - 1 = 2.
        var horizontal = new ConvexHull([new Vector3(-2, 3, 0), new Vector3(2, 3, 0)]);
        var vertical = new ConvexHull([new Vector3(0, 0, 0), new Vector3(0, 1, 0)]);

        Assert.Equal(2.0, OpenGjk.ComputeMinimumDistance(horizontal, vertical, out Vector3 onHorizontal, out Vector3 onVertical), 5);
        Assert.Equal(3.0, onHorizontal.Y, 4);
        Assert.Equal(1.0, onVertical.Y, 4);
    }

    [Fact]
    public void Should_QueryExactly_When_HullIsFlat()
    {
        // A unit square on the z = 0 plane and a sphere sitting above its center.
        var square = new ConvexHull([new Vector3(0, 0, 0), new Vector3(1, 0, 0), new Vector3(1, 1, 0), new Vector3(0, 1, 0)]);
        var sphere = new Sphere(new Vector3(0.5f, 0.5f, 3f), 1f);

        Assert.Equal(2.0, OpenGjk.ComputeMinimumDistance(square, sphere), 5);

        // Sunk halfway into the square: the penetration depth is 0.5.
        var sunk = new Sphere(new Vector3(0.5f, 0.5f, 0.5f), 1f);
        double depth = -OpenGjk.ComputeCollisionInformation(square, sunk, out Vector3 normal);

        Assert.Equal(0.5, depth, 3);
        Assert.Equal(1.0, Math.Abs(normal.Z), 3);
    }

    [Fact]
    public void Should_StayStable_When_DegenerateHullsCoincide()
    {
        // Zero-thickness shapes at the same place: no NaN, no hang, depth is zero-ish.
        var pointA = new ConvexHull([new Vector3(1, 1, 1)]);
        var pointB = new ConvexHull([new Vector3(1, 1, 1)]);

        double info = OpenGjk.ComputeCollisionInformation(pointA, pointB, out Vector3 normal);

        Assert.False(double.IsNaN(info));
        Assert.Equal(0.0, info, 6);
        Assert.False(float.IsNaN(normal.X) || float.IsNaN(normal.Y) || float.IsNaN(normal.Z));
        Assert.True(OpenGjk.HasCollision(pointA, pointB));
    }

    [Fact]
    public void Should_StayStable_When_CustomShapeHasNoThickness()
    {
        // The interface contract allows zero-thickness shapes; this is the documented way to
        // model points and segments.
        var segment = new Segment(new Vector3(0, 0, 0), new Vector3(0, 4, 0));
        var probe = new Sphere(new Vector3(3, 2, 0), 1f);

        Assert.Equal(2.0, OpenGjk.ComputeMinimumDistance(segment, probe), 4);
    }

    private sealed class Segment : ISupportMappable
    {
        private readonly Vector3 _p0;
        private readonly Vector3 _p1;

        public Segment(Vector3 p0, Vector3 p1)
        {
            _p0 = p0;
            _p1 = p1;
        }

        public Vector3 GetSupportPoint(Vector3 direction)
        {
            return Vector3.Dot(_p0, direction) >= Vector3.Dot(_p1, direction) ? _p0 : _p1;
        }

        public Vector3 GetCenter() => (_p0 + _p1) * 0.5f;
    }
}

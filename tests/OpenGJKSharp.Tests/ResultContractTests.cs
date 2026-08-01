using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

/// <summary>
/// Pins the result contract documented on <see cref="OpenGjk"/>.
/// </summary>
public class ResultContractTests
{
    [Fact]
    public void Should_ReturnExactlyZero_When_BodiesOverlap()
    {
        var random = new Random(20260806);

        for (int i = 0; i < 200; i++)
        {
            var center = new Vector3(
                (float)(random.NextDouble() * 2 - 1),
                (float)(random.NextDouble() * 2 - 1),
                (float)(random.NextDouble() * 2 - 1));

            // Radii chosen so the spheres always overlap.
            var a = new Sphere(Vector3.Zero, 2f);
            var b = new Sphere(center, 1f);

            Assert.Equal(0.0, OpenGjk.ComputeMinimumDistance(a, b));
        }

        Vector3[] cubeA = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(1, 1, 0), new(0, 0, 1), new(1, 0, 1), new(0, 1, 1), new(1, 1, 1)];
        Vector3[] cubeB = [new(0.5f, 0.5f, 0.5f), new(1.5f, 0.5f, 0.5f), new(0.5f, 1.5f, 0.5f), new(1.5f, 1.5f, 0.5f)];

        Assert.Equal(0.0, OpenGjk.ComputeMinimumDistance(cubeA, cubeB));
        Assert.Equal(0.0, OpenGjk.ComputeMinimumDistance([0.0, 0, 0], [0.0, 0, 0]));
    }

    [Fact]
    public void Should_ReturnUnitContactNormalFromFirstTowardSecond_When_BodiesOverlap()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(1.5f, 0, 0), 1f);

        double info = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 normal);

        Assert.True(info < 0);
        Assert.Equal(1.0, normal.Length(), 4);
        Assert.Equal(1.0, normal.X, 3);
    }

    [Fact]
    public void Should_ReturnSomeUnitNormal_When_BodiesAreConcentric()
    {
        // Every direction is a valid minimum translation here; only unit length is promised.
        var a = new Sphere(Vector3.Zero, 1f);
        var b = new Sphere(Vector3.Zero, 1f);

        double info = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 normal);

        Assert.True(info < 0, $"expected collision, got {info}");
        Assert.Equal(1.0, normal.Length(), 4);
    }

    [Fact]
    public void Should_PlaceWitnessPointsOnSurfaces_When_BodiesAreSeparated()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(5, 0, 0), 2f);

        double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

        Assert.Equal(2.0, distance, 4);
        Assert.Equal(1.0, (closestA - a.Center).Length(), 3);
        Assert.Equal(2.0, (closestB - b.Center).Length(), 3);
        Assert.Equal(distance, (closestB - closestA).Length(), 3);
    }

    [Fact]
    public void Should_CollapseWitnessPointsIntoTheIntersection_When_BodiesOverlap()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1.5f);
        var b = new Sphere(new Vector3(0.5f, 0, 0), 1f);

        OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

        // Documented behavior: a single common point inside both bodies, not on either surface.
        Assert.Equal(0.0, (closestB - closestA).Length(), 6);
        Assert.True((closestA - a.Center).Length() < 1.5f, "witness point should be inside the first sphere");
        Assert.True((closestA - b.Center).Length() < 1f, "witness point should be inside the second sphere");
    }

    [Fact]
    public void Should_NotDependOnWhichInteriorPointGetCenterReturns()
    {
        // Same shape, different (still interior) center seeds. The seed changes the iteration
        // path, so results agree only up to the convergence tolerance - not bit for bit.
        var probe = new Box(new Vector3(6, 2, -1), new Vector3(1, 1, 1));
        var centered = new OffCenterSphere(new Vector3(1, 1, 1), 2f, Vector3.Zero);
        var offCenter = new OffCenterSphere(new Vector3(1, 1, 1), 2f, new Vector3(1.9f, 0, 0));

        double expected = OpenGjk.ComputeMinimumDistance(centered, probe, out Vector3 expectedA, out Vector3 expectedB);
        double actual = OpenGjk.ComputeMinimumDistance(offCenter, probe, out Vector3 actualA, out Vector3 actualB);

        Assert.Equal(expected, actual, 4);
        Assert.Equal(expectedA.X, actualA.X, 4);
        Assert.Equal(expectedB.X, actualB.X, 4);
    }

    [Fact]
    public void Should_NeverReceiveZeroDirection_When_SupportFunctionIsCalled()
    {
        // Concentric shapes make the initial direction (centerA - centerB) vanish, which is
        // the case the "direction is never zero" guarantee has to survive.
        var a = new DirectionRecordingShape(new Sphere(Vector3.Zero, 1f));
        var b = new DirectionRecordingShape(new Sphere(Vector3.Zero, 1f));

        OpenGjk.ComputeCollisionInformation(a, b, out _);

        Assert.True(a.Calls > 0 && b.Calls > 0);
        Assert.Equal(0, a.ZeroDirections);
        Assert.Equal(0, b.ZeroDirections);
    }

    private sealed class OffCenterSphere : ISupportMappable
    {
        private readonly Sphere _sphere;
        private readonly Vector3 _interiorPoint;

        public OffCenterSphere(Vector3 center, float radius, Vector3 interiorOffset)
        {
            _sphere = new Sphere(center, radius);
            _interiorPoint = center + interiorOffset;
        }

        public Vector3 GetSupportPoint(Vector3 direction) => _sphere.GetSupportPoint(direction);

        public Vector3 GetCenter() => _interiorPoint;
    }

    private sealed class DirectionRecordingShape : ISupportMappable
    {
        private readonly ISupportMappable _inner;

        public DirectionRecordingShape(ISupportMappable inner) => _inner = inner;

        public int Calls { get; private set; }

        public int ZeroDirections { get; private set; }

        public Vector3 GetSupportPoint(Vector3 direction)
        {
            Calls++;
            if (direction == Vector3.Zero)
            {
                ZeroDirections++;
            }

            return _inner.GetSupportPoint(direction);
        }

        public Vector3 GetCenter() => _inner.GetCenter();
    }
}

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
        // (float)Math.PI instead of MathF.PI: MathF is unavailable on net48.
        var orientation = Quaternion.CreateFromAxisAngle(Vector3.UnitZ, (float)Math.PI / 2);
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
}

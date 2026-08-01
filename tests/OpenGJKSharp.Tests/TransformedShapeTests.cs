using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class TransformedShapeTests
{
    [Fact]
    public void Should_MatchOrientedBox_When_LocalBoxIsPosed()
    {
        var random = new Random(20260803);
        var probe = new Sphere(new Vector3(6, 1, -2), 0.75f);
        var halfExtents = new Vector3(0.4f, 1.3f, 0.9f);
        var local = new Box(Vector3.Zero, halfExtents);
        var posed = new TransformedShape(local);

        for (int i = 0; i < 200; i++)
        {
            Quaternion rotation = RandomRotation(random);
            var translation = new Vector3(
                (float)(random.NextDouble() * 4 - 2),
                (float)(random.NextDouble() * 4 - 2),
                (float)(random.NextDouble() * 4 - 2));

            // The built-in oriented box is the reference: same shape, same pose.
            var reference = new Box(translation, halfExtents, rotation);
            posed.SetPose(rotation, translation);

            double expected = OpenGjk.ComputeMinimumDistance(reference, probe);
            double actual = OpenGjk.ComputeMinimumDistance(posed, probe);

            Assert.Equal(expected, actual, 5);
        }
    }

    [Fact]
    public void Should_RotateCapsule_When_PoseIsApplied()
    {
        // A capsule along the local x axis, rotated onto the y axis and lifted by 4.
        var local = new Capsule(new Vector3(-1, 0, 0), new Vector3(1, 0, 0), 0.5f);
        Quaternion rotation = Quaternion.CreateFromAxisAngle(Vector3.UnitZ, (float)(Math.PI / 2));
        var posed = new TransformedShape(local, rotation, new Vector3(0, 4, 0));

        // The capsule now spans y in [3, 5] with radius 0.5, so its lowest point is y = 2.5.
        var ground = new Box(new Vector3(0, 0, 0), new Vector3(5, 1, 5));

        Assert.Equal(1.5, OpenGjk.ComputeMinimumDistance(posed, ground), 4);
    }

    [Fact]
    public void Should_TranslateSupportAndCenter_When_PoseIsApplied()
    {
        var local = new Sphere(new Vector3(1, 0, 0), 1f);
        Quaternion rotation = Quaternion.CreateFromAxisAngle(Vector3.UnitZ, (float)(Math.PI / 2));
        var posed = new TransformedShape(local, rotation, new Vector3(0, 0, 3));

        // The local center (1, 0, 0) rotates to (0, 1, 0) and then moves to (0, 1, 3).
        Vector3 center = posed.GetCenter();
        Assert.Equal(0.0, center.X, 4);
        Assert.Equal(1.0, center.Y, 4);
        Assert.Equal(3.0, center.Z, 4);

        // Support along +y is the far side of the sphere: center + (0, 1, 0).
        Vector3 support = posed.GetSupportPoint(Vector3.UnitY);
        Assert.Equal(0.0, support.X, 4);
        Assert.Equal(2.0, support.Y, 4);
        Assert.Equal(3.0, support.Z, 4);
    }

    [Fact]
    public void Should_UpdateResults_When_PoseChangesBetweenQueries()
    {
        var posed = new TransformedShape(new Sphere(Vector3.Zero, 1f));
        var other = new Sphere(new Vector3(10, 0, 0), 1f);
        var workspace = new GjkWorkspace();

        Assert.Equal(8.0, OpenGjk.ComputeMinimumDistance(posed, other, workspace), 4);

        posed.SetPose(Quaternion.Identity, new Vector3(4, 0, 0));
        Assert.Equal(4.0, OpenGjk.ComputeMinimumDistance(posed, other, workspace), 4);

        posed.SetPose(Quaternion.Identity, new Vector3(8.5f, 0, 0));
        Assert.True(OpenGjk.HasCollision(posed, other, workspace));
    }

    [Fact]
    public void Should_NormalizeRotation_When_PoseIsSet()
    {
        var posed = new TransformedShape(new Sphere(Vector3.Zero, 1f), new Quaternion(0, 0, 2, 2), Vector3.Zero);

        Assert.Equal(1.0, posed.Rotation.Length(), 5);
    }

    [Fact]
    public void Should_Throw_When_PoseArgumentsAreInvalid()
    {
        var sphere = new Sphere(Vector3.Zero, 1f);

        Assert.Throws<ArgumentNullException>(() => new TransformedShape(null!));
        Assert.Throws<ArgumentException>(() => new TransformedShape(sphere, new Quaternion(0, 0, 0, 0), Vector3.Zero));
        Assert.Throws<ArgumentException>(() => new TransformedShape(sphere, Quaternion.Identity, new Vector3(float.NaN, 0, 0)));
        Assert.Throws<ArgumentException>(() => new TransformedShape(sphere).SetPose(Quaternion.Identity, new Vector3(float.PositiveInfinity, 0, 0)));
    }

    [Fact]
    public void Should_ComposeNestedTransforms_When_WrapperIsWrapped()
    {
        // Two 90 degree turns about z compose into a 180 degree turn.
        Quaternion quarter = Quaternion.CreateFromAxisAngle(Vector3.UnitZ, (float)(Math.PI / 2));
        var inner = new TransformedShape(new Sphere(new Vector3(1, 0, 0), 0.5f), quarter, Vector3.Zero);
        var outer = new TransformedShape(inner, quarter, Vector3.Zero);

        // Local center (1, 0, 0) -> (0, 1, 0) -> (-1, 0, 0).
        Vector3 center = outer.GetCenter();
        Assert.Equal(-1.0, center.X, 4);
        Assert.Equal(0.0, center.Y, 4);
    }

    private static Quaternion RandomRotation(Random random)
    {
        var axis = new Vector3(
            (float)(random.NextDouble() * 2 - 1),
            (float)(random.NextDouble() * 2 - 1),
            (float)(random.NextDouble() * 2 - 1));

        if (axis.LengthSquared() <= 0)
        {
            axis = Vector3.UnitZ;
        }

        return Quaternion.CreateFromAxisAngle(Vector3.Normalize(axis), (float)(random.NextDouble() * Math.PI * 2));
    }
}

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

        // Fully symmetric case: every direction is equally minimal, so EPA's polytope
        // refinement spreads over the whole Minkowski sphere and underestimates the
        // depth by a few percent within its iteration cap. Allow 10%.
        Assert.True(distance < 0, $"expected collision, got {distance}");
        Assert.True(Math.Abs(-distance - 1.5) < 0.15, $"expected depth ~1.5, got {-distance}");
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
        // Sphere r=1 vs box [-1, 1]^3: at x=2.5 the surfaces are 0.5 apart (separated);
        // at x=1.5 they overlap by 0.5.
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
}

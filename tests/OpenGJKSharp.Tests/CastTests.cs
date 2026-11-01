using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class CastTests
{
    private const double Tolerance = 1e-5;

    // On a curved surface the normal is that of the last separating plane, which lags the
    // surface by the convergence tolerance; the contact point inherits the same error.
    private const double CurvedTolerance = 2e-3;

    public static TheoryData<string, double, Vector3, double> RayHits()
    {
        double rotatedFace = 10 - (1 / Math.Cos(Math.PI / 6));
        float third = (float)(1 / Math.Sqrt(3));

        return new TheoryData<string, double, Vector3, double>
        {
            { "sphere-center", 0.4, new Vector3(-1, 0, 0), Tolerance },
            { "sphere-offset", 0.42, new Vector3(-0.8f, 0.6f, 0), CurvedTolerance },
            { "box-face", 0.4, new Vector3(0, 1, 0), Tolerance },
            { "box-angled", 0.5, new Vector3(-1, 0, 0), Tolerance },
            { "capsule-side", 2.5 / 6, new Vector3(1, 0, 0), Tolerance },
            { "capsule-cap", 0.35, new Vector3(0, 1, 0), Tolerance },
            { "hull-slanted-face", 1.0 / 3, new Vector3(third, third, third), Tolerance },
            { "transformed-box", rotatedFace / 20, new Vector3((float)-Math.Cos(Math.PI / 6), (float)-Math.Sin(Math.PI / 6), 0), Tolerance },
        };
    }

    [Theory]
    [MemberData(nameof(RayHits))]
    public void Should_ReportAnalyticHit_When_RayHitsShape(string scenario, double expectedFraction, Vector3 expectedNormal, double normalTolerance)
    {
        (ISupportMappable shape, Vector3 origin, Vector3 translation) = RayScenario(scenario);

        bool hit = OpenGjk.TryCastRay(shape, origin, translation, out double fraction, out Vector3 normal);

        Assert.True(hit);
        Assert.Equal(expectedFraction, fraction, Tolerance);
        Assert.True(fraction <= expectedFraction + 1e-7, $"fraction {fraction} lies past the true hit {expectedFraction}");
        AssertVector(expectedNormal, normal, normalTolerance);
    }

    [Theory]
    [InlineData(-5f, 1.5f, 10f)] // passes above
    [InlineData(-5f, 0f, -10f)] // points away
    [InlineData(-5f, 0f, 2f)] // stops short
    [InlineData(-5f, 1.001f, 10f)] // grazes just outside
    public void Should_Miss_When_RayDoesNotReachSphere(float x, float y, float length)
    {
        var sphere = new Sphere(Vector3.Zero, 1);

        bool hit = OpenGjk.TryCastRay(sphere, new Vector3(x, y, 0), new Vector3(length, 0, 0), out double fraction, out Vector3 normal);

        Assert.False(hit);
        Assert.Equal(0, fraction);
        Assert.Equal(Vector3.Zero, normal);
    }

    [Fact]
    public void Should_Hit_When_RayGrazesSphereJustInside()
    {
        var sphere = new Sphere(Vector3.Zero, 1);
        double expected = (5 - Math.Sqrt(1 - (0.999 * 0.999))) / 10;

        bool hit = OpenGjk.TryCastRay(sphere, new Vector3(-5, 0.999f, 0), new Vector3(10, 0, 0), out double fraction, out _);

        Assert.True(hit);
        Assert.Equal(expected, fraction, 1e-4);
    }

    [Fact]
    public void Should_Hit_When_RayIsTangentToSphere()
    {
        var sphere = new Sphere(Vector3.Zero, 1);

        bool hit = OpenGjk.TryCastRay(sphere, new Vector3(-5, 1, 0), new Vector3(10, 0, 0), out double fraction, out Vector3 normal);

        // A tangent ray touches at a single point, which conservative advancement approaches slowly.
        Assert.True(hit);
        Assert.InRange(fraction, 0.49, 0.5 + 1e-7);
        Assert.Equal(1, normal.Length(), Tolerance);
    }

    [Fact]
    public void Should_Hit_When_RaySlidesAlongBoxFace()
    {
        var box = new Box(Vector3.Zero, Vector3.One);

        bool hit = OpenGjk.TryCastRay(box, new Vector3(-5, 1, 0), new Vector3(10, 0, 0), out double fraction, out Vector3 normal);

        // The ray first touches the edge at x = -1; any normal of that edge is valid.
        Assert.True(hit);
        Assert.Equal(0.4, fraction, Tolerance);
        Assert.Equal(1, normal.Length(), Tolerance);
        Assert.True(normal.X <= Tolerance && normal.Y >= -Tolerance, $"normal {normal} is not a normal of the edge");
    }

    [Fact]
    public void Should_ReportZeroFractionAndNormal_When_RayStartsInsideShape()
    {
        var sphere = new Sphere(Vector3.Zero, 1);

        bool hit = OpenGjk.TryCastRay(sphere, new Vector3(0.2f, 0.1f, 0), new Vector3(10, 0, 0), out double fraction, out Vector3 normal);

        Assert.True(hit);
        Assert.Equal(0, fraction);
        Assert.Equal(Vector3.Zero, normal);
    }

    [Fact]
    public void Should_TestContainment_When_RayTranslationIsZero()
    {
        var box = new Box(Vector3.Zero, Vector3.One);

        Assert.True(OpenGjk.TryCastRay(box, new Vector3(0.5f, 0, 0), Vector3.Zero, out double inside, out _));
        Assert.Equal(0, inside);
        Assert.False(OpenGjk.TryCastRay(box, new Vector3(3, 0, 0), Vector3.Zero, out _, out _));
    }

    [Fact]
    public void Should_ReportGapOverTranslation_When_SphereIsCastTowardBox()
    {
        var sphere = new Sphere(Vector3.Zero, 1);
        var box = new Box(new Vector3(5, 0, 0), Vector3.One);

        bool hit = OpenGjk.TryCastShape(sphere, box, new Vector3(10, 0, 0), out double fraction, out Vector3 normal, out Vector3 point);

        Assert.True(hit);
        Assert.Equal(0.3, fraction, Tolerance);
        AssertVector(new Vector3(1, 0, 0), normal, Tolerance);
        AssertVector(new Vector3(4, 0, 0), point, CurvedTolerance);
    }

    [Fact]
    public void Should_ReportContactGeometry_When_SphereIsCastObliquelyTowardSphere()
    {
        var a = new Sphere(Vector3.Zero, 1);
        var b = new Sphere(new Vector3(5, 1, 0), 1);
        double s3 = Math.Sqrt(3);

        bool hit = OpenGjk.TryCastShape(a, b, new Vector3(10, 0, 0), out double fraction, out Vector3 normal, out Vector3 point);

        Assert.True(hit);
        Assert.Equal((5 - s3) / 10, fraction, Tolerance);
        AssertVector(new Vector3((float)(s3 / 2), 0.5f, 0), normal, CurvedTolerance);
        AssertVector(new Vector3((float)(5 - (s3 / 2)), 0.5f, 0), point, CurvedTolerance);
    }

    [Fact]
    public void Should_ReportGroundContact_When_CapsuleIsCastDown()
    {
        var capsule = new Capsule(new Vector3(0, -1, 0), new Vector3(0, 1, 0), 0.5f);
        var ground = new Box(new Vector3(0, -5, 0), new Vector3(5, 1, 5));

        bool hit = OpenGjk.TryCastShape(capsule, ground, new Vector3(0, -10, 0), out double fraction, out Vector3 normal, out Vector3 point);

        Assert.True(hit);
        Assert.Equal(0.25, fraction, Tolerance);
        AssertVector(new Vector3(0, -1, 0), normal, Tolerance);
        AssertVector(new Vector3(0, -4, 0), point, Tolerance);
    }

    [Fact]
    public void Should_Miss_When_BoxSlidesPastBox()
    {
        var moving = new Box(new Vector3(0, 3, 0), Vector3.One);
        var fixedBox = new Box(new Vector3(5, 0, 0), Vector3.One);

        bool hit = OpenGjk.TryCastShape(moving, fixedBox, new Vector3(10, 0, 0), out double fraction, out Vector3 normal, out Vector3 point);

        Assert.False(hit);
        Assert.Equal(0, fraction);
        Assert.Equal(Vector3.Zero, normal);
        Assert.Equal(Vector3.Zero, point);
    }

    [Fact]
    public void Should_ReportCommonPoint_When_CastStartsOverlapping()
    {
        var sphere = new Sphere(Vector3.Zero, 1);
        var box = new Box(new Vector3(1.5f, 0, 0), Vector3.One);

        bool hit = OpenGjk.TryCastShape(sphere, box, new Vector3(10, 0, 0), out double fraction, out Vector3 normal, out Vector3 point);

        Assert.True(hit);
        Assert.Equal(0, fraction);
        Assert.Equal(Vector3.Zero, normal);
        Assert.True(point.Length() <= 1 + Tolerance, $"point {point} is outside the sphere");
        Assert.True(point.X >= 0.5 - Tolerance, $"point {point} is outside the box");
    }

    [Fact]
    public void Should_TestOverlap_When_CastTranslationIsZero()
    {
        var sphere = new Sphere(Vector3.Zero, 1);

        Assert.True(OpenGjk.TryCastShape(sphere, new Box(new Vector3(1.5f, 0, 0), Vector3.One), Vector3.Zero,
            out double fraction, out _, out _));
        Assert.Equal(0, fraction);
        Assert.False(OpenGjk.TryCastShape(sphere, new Box(new Vector3(3, 0, 0), Vector3.One), Vector3.Zero,
            out _, out _, out _));
    }

    [Fact]
    public void Should_MatchPlainOverload_When_CastsReuseWorkspace()
    {
        var workspace = new GjkWorkspace();
        var random = new Random(20261101);

        for (int i = 0; i < 200; i++)
        {
            var a = new Sphere(RandomPoint(random), 0.5f + (float)random.NextDouble());
            var b = new Box(RandomPoint(random), new Vector3(0.4f, 0.7f, 1.1f));
            Vector3 origin = RandomPoint(random);
            Vector3 translation = RandomPoint(random) * 2;

            bool expectedRay = OpenGjk.TryCastRay(b, origin, translation, out double expectedRayFraction, out Vector3 expectedRayNormal);
            bool actualRay = OpenGjk.TryCastRay(b, origin, translation, workspace, out double actualRayFraction, out Vector3 actualRayNormal);

            Assert.Equal(expectedRay, actualRay);
            Assert.Equal(expectedRayFraction, actualRayFraction);
            Assert.Equal(expectedRayNormal, actualRayNormal);

            bool expected = OpenGjk.TryCastShape(a, b, translation, out double expectedFraction, out Vector3 expectedNormal, out Vector3 expectedPoint);
            bool actual = OpenGjk.TryCastShape(a, b, translation, workspace, out double actualFraction, out Vector3 actualNormal, out Vector3 actualPoint);

            Assert.Equal(expected, actual);
            Assert.Equal(expectedFraction, actualFraction);
            Assert.Equal(expectedNormal, actualNormal);
            Assert.Equal(expectedPoint, actualPoint);
        }
    }

    [Fact]
    public void Should_AgreeWithDistanceQuery_When_CastsAreRandom()
    {
        var random = new Random(1101);

        for (int i = 0; i < 300; i++)
        {
            var a = new Sphere(RandomPoint(random), 0.3f + (float)random.NextDouble());
            var b = new Box(RandomPoint(random), new Vector3(0.4f, 0.7f, 1.1f));
            Vector3 translation = RandomPoint(random) * 3;

            bool hit = OpenGjk.TryCastShape(a, b, translation, out double fraction, out _, out _);
            if (!hit || fraction == 0)
            {
                continue;
            }

            // At the reported fraction the moved sphere touches the box and has not passed into it.
            var moved = new Sphere(a.Center + ((float)fraction * translation), a.Radius);
            double distance = OpenGjk.ComputeMinimumDistance(moved, b);
            Assert.InRange(distance, 0, 1e-4);
        }
    }

    [Fact]
    public void Should_Throw_When_CastInputIsInvalid()
    {
        var sphere = new Sphere(Vector3.Zero, 1);
        var workspace = new GjkWorkspace();
        var nan = new Vector3(float.NaN, 0, 0);
        var infinite = new Vector3(0, float.PositiveInfinity, 0);

        Assert.Throws<ArgumentNullException>("shape", () => OpenGjk.TryCastRay(null!, Vector3.Zero, Vector3.UnitX, out _, out _));
        Assert.Throws<ArgumentException>("origin", () => OpenGjk.TryCastRay(sphere, nan, Vector3.UnitX, out _, out _));
        Assert.Throws<ArgumentException>("translation", () => OpenGjk.TryCastRay(sphere, Vector3.Zero, infinite, out _, out _));
        Assert.Throws<ArgumentNullException>("workspace", () => OpenGjk.TryCastRay(sphere, Vector3.Zero, Vector3.UnitX, null!, out _, out _));

        Assert.Throws<ArgumentNullException>("a", () => OpenGjk.TryCastShape(null!, sphere, Vector3.UnitX, out _, out _, out _));
        Assert.Throws<ArgumentNullException>("b", () => OpenGjk.TryCastShape(sphere, null!, Vector3.UnitX, out _, out _, out _));
        Assert.Throws<ArgumentException>("translation", () => OpenGjk.TryCastShape(sphere, sphere, nan, out _, out _, out _));
        Assert.Throws<ArgumentNullException>("workspace", () => OpenGjk.TryCastShape(sphere, sphere, Vector3.UnitX, null!, out _, out _, out _));

        // A rejected query must not leave the workspace occupied.
        Assert.False(OpenGjk.TryCastRay(sphere, new Vector3(5, 0, 0), Vector3.UnitX, workspace, out _, out _));
    }

#if NET8_0_OR_GREATER
    [Fact]
    public void Should_StopAllocating_When_CastsReuseWorkspace()
    {
        var workspace = new GjkWorkspace();
        var sphere = new Sphere(Vector3.Zero, 1);
        var box = new Box(new Vector3(5, 0, 0), Vector3.One);
        var translation = new Vector3(10, 0.5f, 0);

        // Warm up: the first queries grow the buffers the workspace then keeps.
        for (int i = 0; i < 20; i++)
        {
            OpenGjk.TryCastRay(sphere, new Vector3(-5, 0.3f, 0), translation, workspace, out _, out _);
            OpenGjk.TryCastShape(sphere, box, translation, workspace, out _, out _, out _);
        }

        const int queries = 100;
        long before = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < queries; i++)
        {
            OpenGjk.TryCastRay(sphere, new Vector3(-5, 0.3f, 0), translation, workspace, out _, out _);
            OpenGjk.TryCastShape(sphere, box, translation, workspace, out _, out _, out _);
        }
        long reusedBytes = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.True(reusedBytes == 0, $"expected no allocation with a reused workspace, got {reusedBytes} bytes");
    }
#endif

    private static (ISupportMappable Shape, Vector3 Origin, Vector3 Translation) RayScenario(string scenario)
    {
        var capsule = new Capsule(new Vector3(0, -1, 0), new Vector3(0, 1, 0), 0.5f);

        return scenario switch
        {
            "sphere-center" => (new Sphere(Vector3.Zero, 1), new Vector3(-5, 0, 0), new Vector3(10, 0, 0)),
            "sphere-offset" => (new Sphere(Vector3.Zero, 1), new Vector3(-5, 0.6f, 0), new Vector3(10, 0, 0)),
            "box-face" => (new Box(Vector3.Zero, new Vector3(1, 2, 3)), new Vector3(0, 10, 0.5f), new Vector3(0, -20, 0)),
            "box-angled" => (new Box(Vector3.Zero, new Vector3(1, 2, 3)), new Vector3(-5, 0.5f, 0.5f), new Vector3(8, 1, 0)),
            "capsule-side" => (capsule, new Vector3(3, 0, 0), new Vector3(-6, 0, 0)),
            "capsule-cap" => (capsule, new Vector3(0, 5, 0), new Vector3(0, -10, 0)),
            "hull-slanted-face" => (
                new ConvexHull([Vector3.Zero, Vector3.UnitX, Vector3.UnitY, Vector3.UnitZ]),
                Vector3.One,
                new Vector3(-2, -2, -2)),
            "transformed-box" => (
                new TransformedShape(new Box(Vector3.Zero, Vector3.One),
                    Quaternion.CreateFromAxisAngle(Vector3.UnitZ, (float)(Math.PI / 6)), new Vector3(10, 0, 0)),
                Vector3.Zero,
                new Vector3(20, 0, 0)),
            _ => throw new ArgumentOutOfRangeException(nameof(scenario), scenario, null),
        };
    }

    private static void AssertVector(Vector3 expected, Vector3 actual, double tolerance)
    {
        Assert.True(Vector3.Distance(expected, actual) <= tolerance, $"expected {expected}, got {actual}");
    }

    private static Vector3 RandomPoint(Random random)
    {
        return new Vector3(
            (float)(random.NextDouble() * 6 - 3),
            (float)(random.NextDouble() * 6 - 3),
            (float)(random.NextDouble() * 6 - 3));
    }
}

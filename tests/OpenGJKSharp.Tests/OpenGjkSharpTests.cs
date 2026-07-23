using OpenGJKSharp.Models;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class OpenGjkSharpTests
{
    [Fact]
    public void Should_ReturnTrue_When_Two3DCubesOverlapping()
    {
        var a = new Vector3[]
        {
            new(0, 0, 0),
            new(1, 0, 0),
            new(0, 1, 0),
            new(1, 1, 0),
            new(0, 0, 1),
            new(1, 0, 1),
            new(0, 1, 1),
            new(1, 1, 1),
        };
        var b = new Vector3[]
        {
            new(0.5f, 0.5f, 0),
            new(1.5f, 0.5f, 0),
            new(0.5f, 1.5f, 0),
            new(1.5f, 1.5f, 0),
            new(0.5f, 0.5f, 1),
            new(1.5f, 0.5f, 1),
            new(0.5f, 1.5f, 1),
            new(1.5f, 1.5f, 1),
        };

        bool hasCollision = OpenGjk.HasCollision(a, b);

        Assert.True(hasCollision);
    }

    /// <summary>
    /// https://github.com/MattiaMontanari/openGJK/issues/52
    /// </summary>
    [Fact]
    public void Should_ReturnTrue_When_Two3DCubesSharingSameCenter()
    {
        var ia = new double[3, 8] {
            { 90945.08, 91077.46, 91454.92, 91322.54, 90945.08, 91077.46, 91454.92, 91322.54 },
            { 54122.54, 53745.082, 53877.46, 54254.918, 54122.54, 53745.082, 53877.46, 54254.918 },
            { 300, 300, 300, 300, 0, 0, 0, 0, },
        };

        var ib = new double[3, 8] {
            { 90919.64, 91162.62, 91480.36, 91237.38, 90919.64, 91162.62, 91480.36, 91237.38 },
            { 54037.383, 53719.637, 53962.617, 54280.363, 54037.383, 53719.637, 53962.617, 54280.363 },
            { 300, 300, 300, 300, 0, 0, 0, 0, },
        };

#pragma warning disable CS0618 // Intentionally testing the obsolete compatibility API
        var distance = OpenGjk.CsFunction(ia.GetLength(1), ia, ib.GetLength(1), ib);
#pragma warning restore CS0618

        Assert.Equal(0, distance);
    }

    /// <summary>
    /// https://github.com/MattiaMontanari/openGJK/issues/56
    /// </summary>
    [Fact]
    public void Should_ReturnZero_When_2DArcDoesNotCollideWithRectangle()
    {
        Vector2[] a = [
            new(14483f, 45848f),
            new(15133f, 44923f),
            new(15147.594f, 44933.43f),
            new(15162.021f, 44944.09f),
            new(15176.277f, 44954.977f),
            new(15190.358f, 44966.086f),
            new(15204.262f, 44977.42f),
            new(15217.983f, 44988.977f),
            new(15231.5205f, 45000.746f),
            new(15244.868f, 45012.727f),
            new(15258.025f, 45024.918f),
            new(15270.986f, 45037.32f),
            new(15283.75f, 45049.926f),
            new(15296.312f, 45062.73f),
            new(15308.668f, 45075.73f),
            new(15320.816f, 45088.93f),
            new(15332.755f, 45102.316f),
            new(15344.479f, 45115.895f),
            new(15355.985f, 45129.652f),
            new(15367.273f, 45143.594f),
            new(15378.338f, 45157.715f),
            new(15389.177f, 45172.008f),
            new(15399.788f, 45186.47f),
            new(15410.169f, 45201.098f),
            new(15420.315f, 45215.89f),
            new(15430.227f, 45230.84f),
            new(15439.899f, 45245.945f),
            new(15449.331f, 45261.203f),
            new(15458.52f, 45276.61f),
            new(15467.463f, 45292.16f),
            new(15476.157f, 45307.848f),
            new(15484.603f, 45323.676f),
            new(15492.796f, 45339.633f),
            new(15500.734f, 45355.72f),
            new(15508.417f, 45371.926f),
            new(15515.841f, 45388.254f),
            new(15523.005f, 45404.7f),
            new(15529.907f, 45421.258f),
            new(15536.547f, 45437.92f),
            new(15542.92f, 45454.688f),
            new(15549.027f, 45471.555f),
            new(15554.865f, 45488.516f),
            new(15560.435f, 45505.566f),
            new(15565.731f, 45522.703f),
            new(15570.757f, 45539.92f),
            new(15575.508f, 45557.223f),
            new(15579.983f, 45574.59f),
            new(15584.184f, 45592.03f),
            new(15588.105f, 45609.535f),
            new(15591.75f, 45627.098f),
            new(15595.115f, 45644.715f),
            new(15598.201f, 45662.387f),
            new(15601.005f, 45680.105f),
            new(15603.528f, 45697.863f),
            new(15605.77f, 45715.66f),
            new(15607.728f, 45733.492f),
            new(15609.403f, 45751.35f),
            new(15610.795f, 45769.234f),
            new(15611.902f, 45787.137f),
            new(15612.726f, 45805.055f),
            new(15613.266f, 45822.984f),
            new(15613.52f, 45840.92f),
            new(15613.49f, 45858.86f),
            new(15613.175f, 45876.793f),
            new(15612.576f, 45894.723f),
            new(15611.692f, 45912.637f),
            new(15610.525f, 45930.54f),
            new(15609.073f, 45948.418f),
            new(15607.339f, 45966.27f),
            new(15605.32f, 45984.094f),
            new(15603.0205f, 46001.883f),
            new(15600.4375f, 46019.633f),
            new(15597.574f, 46037.344f),
            new(15594.43f, 46055f),
            new(15591.005f, 46072.61f),
            new(15587.302f, 46090.16f),
            new(15583.321f, 46107.652f),
            new(15579.0625f, 46125.08f),
            new(15574.528f, 46142.43f),
            new(15569.72f, 46159.71f),
            new(15564.638f, 46176.914f),
            new(15559.283f, 46194.035f),
            new(15553.657f, 46211.066f),
            new(15547.762f, 46228.008f),
            new(15541.599f, 46244.855f),
            new(15535.169f, 46261.6f),
            new(15528.475f, 46278.242f),
            new(15521.517f, 46294.777f),
            new(15514.297f, 46311.195f),
            new(15506.818f, 46327.5f),
            new(15499.082f, 46343.684f),
            new(15491.09f, 46359.742f),
            new(15482.844f, 46375.67f),
            new(15474.346f, 46391.47f),
            new(15465.598f, 46407.13f),
            new(15456.604f, 46422.65f),
            new(15447.363f, 46438.023f),
            new(15437.881f, 46453.25f),
            new(15428.157f, 46468.324f),
            new(15418.196f, 46483.242f),
            new(15408f, 46498f),
            new(14483f, 45848f),
        ];
        Vector2[] b = [
            new(15612f, 43698f),
            new(17462f, 43698f),
            new(17462f, 48644f),
            new(15612f, 48644f),
        ];

        bool hasCollision = OpenGjk.HasCollision(a, b);

        Assert.True(hasCollision);
    }

    [Fact]
    public void Should_ReturnPenetrationDepth_When_TwoCubesOverlap()
    {
        // Unit cube [0, 1]^3
        var a = new Vector3[]
        {
            new(0, 0, 0),
            new(1, 0, 0),
            new(0, 1, 0),
            new(1, 1, 0),
            new(0, 0, 1),
            new(1, 0, 1),
            new(0, 1, 1),
            new(1, 1, 1),
        };

        // Same cube shifted by +0.5 along x only: minimum penetration is 0.5 along x
        var b = a.Select(p => p + new Vector3(0.5f, 0, 0)).ToArray();

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

        Assert.True(distance < 0);
        Assert.Equal(0.5, -distance, 6);
        Assert.Equal(1.0, Math.Abs(contactNormal.X), 6);
        Assert.Equal(0.0, contactNormal.Y, 6);
        Assert.Equal(0.0, contactNormal.Z, 6);
    }

    [Fact]
    public void Should_ReturnSeparationDistance_When_TwoCubesDoNotCollide()
    {
        var a = new Vector3[]
        {
            new(0, 0, 0),
            new(1, 0, 0),
            new(0, 1, 0),
            new(1, 1, 0),
            new(0, 0, 1),
            new(1, 0, 1),
            new(0, 1, 1),
            new(1, 1, 1),
        };

        // Same cube shifted by +2 along x: separated by 1
        var b = a.Select(p => p + new Vector3(2f, 0, 0)).ToArray();

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

        Assert.Equal(1.0, distance, 6);
        Assert.Equal(1.0, contactNormal.X, 6);
    }

    [Fact]
    public void Should_ReturnDeeperPenetrationDepth_When_CubesOverlapMore()
    {
        var a = new Vector3[]
        {
            new(0, 0, 0),
            new(1, 0, 0),
            new(0, 1, 0),
            new(1, 1, 0),
            new(0, 0, 1),
            new(1, 0, 1),
            new(0, 1, 1),
            new(1, 1, 1),
        };

        // Shifted by +0.2 along x: minimum penetration is 0.8 along x
        var b = a.Select(p => p + new Vector3(0.2f, 0, 0)).ToArray();

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

        Assert.True(distance < 0);
        Assert.Equal(0.8, -distance, 6);
        Assert.Equal(1.0, Math.Abs(contactNormal.X), 6);
    }

    [Fact]
    public void Should_ReturnFalse_When_Two3DCubesDoNotOverlap()
    {
        Vector3[] a = UnitCube();
        Vector3[] b = [.. a.Select(p => p + new Vector3(3f, 0, 0))];

        bool hasCollision = OpenGjk.HasCollision(a, b);

        Assert.False(hasCollision);
    }

    [Fact]
    public void Should_ReturnFalse_When_2DRectanglesDoNotOverlap()
    {
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [new(2, 2), new(3, 2), new(3, 3), new(2, 3)];

        bool hasCollision = OpenGjk.HasCollision(a, b);

        Assert.False(hasCollision);
    }

    [Fact]
    public void Should_ComputeWitnessPoints_When_CubesAreSeparated()
    {
        Vector3[] a = UnitCube();
        Vector3[] b = [.. a.Select(p => p + new Vector3(2f, 0, 0))];

        var simplex = new GkSimplex();
        double distance = OpenGjk.ComputeMinimumDistance(ToGkPolytope(a), ToGkPolytope(b), simplex);

        Assert.Equal(1.0, distance, 6);

        // Closest points lie on the x = 1 face of cube A and the x = 2 face of cube B
        Assert.Equal(1.0, simplex.Witnesses[0, 0], 6);
        Assert.Equal(2.0, simplex.Witnesses[1, 0], 6);
    }

    [Fact]
    public void Should_ReturnFullPenetrationDepth_When_CubesAreIdentical()
    {
        // Identical cubes make the initial GJK search direction
        // (vertex A[0] - vertex B[0]) the zero vector, which degenerates the simplex.
        // EPA must still recover the true penetration depth of 1.
        Vector3[] a = UnitCube();
        Vector3[] b = UnitCube();

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

        Assert.Equal(1.0, -distance, 6);
        Assert.Equal(1.0, contactNormal.Length(), 5);
        Assert.True(OpenGjk.HasCollision(a, b));
    }

    [Fact]
    public void Should_ReturnPenetrationDepth_When_CubeIsNestedSharingACorner()
    {
        // a = [0, 1]^3 fully inside b = [0, 2]^3, sharing the corner at the origin.
        // The first vertices coincide, so the initial GJK direction is the zero vector.
        Vector3[] a = UnitCube();
        Vector3[] b = [.. UnitCube().Select(p => p * 2f)];

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

        Assert.Equal(1.0, -distance, 6);
        Assert.Equal(1.0, contactNormal.Length(), 5);
    }

    [Fact]
    public void Should_ReturnZeroPenetration_When_CubesTouchFaceToFace()
    {
        Vector3[] a = UnitCube();
        Vector3[] b = [.. a.Select(p => p + new Vector3(1f, 0, 0))];

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out _);

        Assert.True(Math.Abs(distance) < 1e-6);
    }

    [Fact]
    public void Should_ReturnPenetrationDepth_When_CubesOverlapDiagonally()
    {
        Vector3[] a = UnitCube();

        // Shifted by (+0.5, +0.5, 0): minimum penetration is 0.5 along x or y
        Vector3[] b = [.. a.Select(p => p + new Vector3(0.5f, 0.5f, 0))];

        double distance = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 contactNormal);

        Assert.True(distance < 0);
        Assert.Equal(0.5, -distance, 6);
        Assert.Equal(1.0, Math.Max(Math.Abs(contactNormal.X), Math.Abs(contactNormal.Y)), 6);
        Assert.Equal(0.0, contactNormal.Z, 6);
    }

    [Fact]
    public void Should_ComputeMinimumDistance_When_GivenVector3Arrays()
    {
        Vector3[] a = UnitCube();
        Vector3[] b = [.. a.Select(p => p + new Vector3(3f, 0, 0))];

        double distance = OpenGjk.ComputeMinimumDistance(a, b);

        Assert.Equal(2.0, distance, 6);
    }

    [Fact]
    public void Should_ComputeWitnessPoints_When_GivenVector3Arrays()
    {
        Vector3[] a = UnitCube();
        Vector3[] b = [.. a.Select(p => p + new Vector3(3f, 0, 0))];

        double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

        Assert.Equal(2.0, distance, 6);
        Assert.Equal(1.0, closestA.X, 5);
        Assert.Equal(3.0, closestB.X, 5);
        Assert.Equal(distance, (closestB - closestA).Length(), 5);
    }

    [Fact]
    public void Should_ReturnConsistentWitnessPoints_For_RandomSeparatedClouds()
    {
        var random = new Random(12345);

        for (int iteration = 0; iteration < 200; iteration++)
        {
            Vector3[] a = RandomCloud(random, new Vector3(0, 0, 0));
            Vector3[] b = RandomCloud(random, new Vector3(5, 0, 0));

            double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

            Assert.True(IsFinite(distance), $"iteration {iteration}: distance is not finite");
            Assert.True(IsFinite(closestA.X) && IsFinite(closestA.Y) && IsFinite(closestA.Z),
                $"iteration {iteration}: closestA is not finite");
            Assert.True(IsFinite(closestB.X) && IsFinite(closestB.Y) && IsFinite(closestB.Z),
                $"iteration {iteration}: closestB is not finite");
            Assert.Equal(distance, (closestB - closestA).Length(), 3);
        }

        static Vector3[] RandomCloud(Random random, Vector3 center) =>
            [.. Enumerable.Range(0, 20).Select(_ => center + new Vector3(
                (float)random.NextDouble(),
                (float)random.NextDouble(),
                (float)random.NextDouble()))];
    }

    [Fact]
    public void Should_ReturnConsistentWitnessPoints_When_PolytopeHasDuplicatedVertices()
    {
        // Duplicated vertices can degenerate the simplex (identical simplex corners),
        // which exercises the degenerate branches of the witness computation.
        Vector3[] a = [new(0, 0, 0), new(0, 0, 0), new(1, 0, 0), new(1, 0, 0), new(0.5f, 1, 0), new(0.5f, 1, 0)];
        Vector3[] b = [new(3, 0, 0), new(3, 0, 0), new(4, 0, 0), new(4, 0, 0), new(3.5f, 1, 0), new(3.5f, 1, 0)];

        double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 closestA, out Vector3 closestB);

        Assert.Equal(2.0, distance, 6);
        Assert.Equal(distance, (closestB - closestA).Length(), 5);
    }

    [Fact]
    public void Should_ThrowArgumentOutOfRangeException_When_PrecisionIsNegative()
    {
        Vector3[] a = UnitCube();
        Vector3[] b = [.. a.Select(p => p + new Vector3(3f, 0, 0))];
        Vector2[] a2 = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b2 = [.. a2.Select(p => p + new Vector2(3f, 0))];

        Assert.Throws<ArgumentOutOfRangeException>(() => OpenGjk.HasCollision(a, b, precision: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => OpenGjk.HasCollision(a2, b2, precision: -1));
    }

    [Fact]
    public void Should_PointContactNormalFromFirstBodyToSecond_When_Separated()
    {
        Vector3[] a = UnitCube();
        Vector3[] b = [.. a.Select(p => p + new Vector3(3f, 0, 0))];

        OpenGjk.ComputeCollisionInformation(a, b, out Vector3 nAb);
        OpenGjk.ComputeCollisionInformation(b, a, out Vector3 nBa);

        Assert.Equal(1.0, nAb.X, 5);
        Assert.Equal(-1.0, nBa.X, 5);
    }

    [Fact]
    public void Should_ComputeMinimumDistance_When_GivenVector2Arrays()
    {
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [.. a.Select(p => p + new Vector2(3f, 0))];

        double distance = OpenGjk.ComputeMinimumDistance(a, b);

        Assert.Equal(2.0, distance, 6);
    }

    [Fact]
    public void Should_ComputeWitnessPoints_When_GivenVector2Arrays()
    {
        Vector2[] a = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] b = [.. a.Select(p => p + new Vector2(3f, 0))];

        double distance = OpenGjk.ComputeMinimumDistance(a, b, out Vector2 closestA, out Vector2 closestB);

        Assert.Equal(2.0, distance, 6);
        Assert.Equal(1.0, closestA.X, 5);
        Assert.Equal(3.0, closestB.X, 5);
        Assert.Equal(distance, (closestB - closestA).Length(), 5);
    }

    [Fact]
    public void Should_ThrowArgumentException_When_CoordRowIsShorterThanThree()
    {
        var invalid = new GkPolytope { NumPoints = 1, Coord = [[1.0, 2.0]] };
        var valid = new GkPolytope { NumPoints = 1, Coord = [[0.0, 0.0, 0.0]] };

        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(invalid, valid));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(valid, invalid));
    }

    [Fact]
    public void Should_ComputeAccuratePenetrationDepth_For_RandomOverlappingBoxes()
    {
        var random = new Random(888);
        int tested = 0;

        for (int i = 0; i < 500; i++)
        {
            Vector3 min1 = NextVector(random) - new Vector3(0.5f);
            Vector3 max1 = min1 + NextVector(random) + new Vector3(0.5f);
            Vector3 min2 = NextVector(random) - new Vector3(0.5f);
            Vector3 max2 = min2 + NextVector(random) + new Vector3(0.5f);

            double ox = Math.Min(max1.X, max2.X) - Math.Max(min1.X, min2.X);
            double oy = Math.Min(max1.Y, max2.Y) - Math.Max(min1.Y, min2.Y);
            double oz = Math.Min(max1.Z, max2.Z) - Math.Max(min1.Z, min2.Z);
            if (ox <= 0.01 || oy <= 0.01 || oz <= 0.01)
            {
                continue; // Skip separated or near-touching pairs
            }

            tested++;

            // Minimum translation to separate along an axis is the distance needed to
            // push one box fully past the other, not the projected overlap length
            // (they differ when one box's projection contains the other's).
            double sx = Math.Min(max1.X - min2.X, max2.X - min1.X);
            double sy = Math.Min(max1.Y - min2.Y, max2.Y - min1.Y);
            double sz = Math.Min(max1.Z - min2.Z, max2.Z - min1.Z);
            double expected = Math.Min(sx, Math.Min(sy, sz));

            double actual = OpenGjk.ComputeCollisionInformation(Box(min1, max1), Box(min2, max2), out _);

            Assert.True(Math.Abs(-actual - expected) < 1e-3,
                $"pair {i}: depth = {-actual}, expected = {expected}");
        }

        Assert.True(tested > 100, $"only {tested} overlapping pairs generated");

        static Vector3 NextVector(Random random) => new(
            (float)random.NextDouble(),
            (float)random.NextDouble(),
            (float)random.NextDouble());

        static Vector3[] Box(Vector3 min, Vector3 max) =>
        [
            new(min.X, min.Y, min.Z), new(max.X, min.Y, min.Z), new(min.X, max.Y, min.Z), new(max.X, max.Y, min.Z),
            new(min.X, min.Y, max.Z), new(max.X, min.Y, max.Z), new(min.X, max.Y, max.Z), new(max.X, max.Y, max.Z),
        ];
    }

    [Fact]
    public void Should_ThrowArgumentException_When_CoordinatesAreNotFinite()
    {
        Vector3[] cube = UnitCube();
        Vector3[] withNaN = [new(float.NaN, 0, 0), new(1, 0, 0), new(0, 1, 0), new(0, 0, 1)];
        Vector3[] withInfinity = [new(float.PositiveInfinity, 0, 0), new(1, 0, 0), new(0, 1, 0), new(0, 0, 1)];
        Vector2[] withNaN2D = [new(float.NaN, 0), new(1, 0), new(0, 1)];

        Assert.Throws<ArgumentException>(() => OpenGjk.HasCollision(withNaN, cube));
        Assert.Throws<ArgumentException>(() => OpenGjk.HasCollision(cube, withInfinity));
        Assert.Throws<ArgumentException>(() => OpenGjk.HasCollision(withNaN2D, [new(0, 0), new(1, 1)]));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(withNaN, cube));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeCollisionInformation(withInfinity, cube, out _));
    }

    [Fact]
    public void Should_ThrowArgumentNullException_When_InputIsNull()
    {
        Vector3[] cube = UnitCube();

        Assert.Throws<ArgumentNullException>(() => OpenGjk.HasCollision(null!, cube));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.HasCollision(cube, null!));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.HasCollision((Vector2[])null!, [new(0, 0)]));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeCollisionInformation(null!, cube, out _));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeCollisionInformation(cube, null!, out _));
    }

    [Fact]
    public void Should_ThrowArgumentException_When_InputIsEmpty()
    {
        Vector3[] cube = UnitCube();

        Assert.Throws<ArgumentException>(() => OpenGjk.HasCollision([], cube));
        Assert.Throws<ArgumentException>(() => OpenGjk.HasCollision(cube, []));
        Assert.Throws<ArgumentException>(() => OpenGjk.HasCollision([new(0, 0)], Array.Empty<Vector2>()));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeCollisionInformation([], cube, out _));
    }

    [Theory]
    [InlineData(0.5f)]
    [InlineData(1.0f)]
    [InlineData(1.5f)]
    [InlineData(1.9f)]
    public void Should_ReturnPenetrationDepth_When_SphereHullsOverlap(float centerDistance)
    {
        // Regression: GJK finishes with a degenerate simplex whose plane passes exactly
        // through the origin; EPA used to mis-orient that face's normal (inward) and
        // report ~0 penetration depth regardless of the actual overlap.
        Vector3[] sphere1 = SphereHull(0, 0, 0, 1);
        Vector3[] sphere2 = SphereHull(centerDistance, 0, 0, 1);

        double distance = OpenGjk.ComputeCollisionInformation(sphere1, sphere2, out Vector3 normal);

        // The discrete hull is slightly smaller than the ideal sphere, so allow ~5%.
        double expectedPenetration = 2 - centerDistance;
        Assert.True(distance < 0, $"expected collision, got distance {distance}");
        Assert.Equal(expectedPenetration, -distance, expectedPenetration * 0.05);
        Assert.Equal(1, normal.Length(), 1e-3);
        Assert.True(normal.X > 0.9f, $"normal should point from body1 toward body2, got {normal}");
    }

    private static Vector3[] SphereHull(float cx, float cy, float cz, float r, int slices = 24, int stacks = 12)
    {
        var points = new List<Vector3>();
        for (int i = 0; i <= stacks; i++)
        {
            double phi = Math.PI * i / stacks;
            for (int j = 0; j < slices; j++)
            {
                double theta = 2 * Math.PI * j / slices;
                points.Add(new Vector3(
                    cx + (float)(r * Math.Sin(phi) * Math.Cos(theta)),
                    cy + (float)(r * Math.Sin(phi) * Math.Sin(theta)),
                    cz + (float)(r * Math.Cos(phi))));
            }
        }

        return [.. points];
    }

    [Fact]
    public void Should_ThrowArgumentException_When_NumPointsExceedsCoordLength()
    {
        var invalid = new GkPolytope { NumPoints = 10, Coord = [[0, 0, 0]] };
        var valid = new GkPolytope { NumPoints = 1, Coord = [[5, 0, 0]] };

        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(invalid, valid));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(valid, invalid));
    }

    private static Vector3[] UnitCube() =>
    [
        new(0, 0, 0),
        new(1, 0, 0),
        new(0, 1, 0),
        new(1, 1, 0),
        new(0, 0, 1),
        new(1, 0, 1),
        new(0, 1, 1),
        new(1, 1, 1),
    ];

    private static GkPolytope ToGkPolytope(Vector3[] points) => new()
    {
        NumPoints = points.Length,
        Coord = [.. points.Select(p => (double[])[p.X, p.Y, p.Z])],
    };

    // double.IsFinite/float.IsFinite are unavailable on net48
    private static bool IsFinite(double value) => !double.IsNaN(value) && !double.IsInfinity(value);

    private static bool IsFinite(float value) => !float.IsNaN(value) && !float.IsInfinity(value);
}
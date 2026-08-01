using System.Numerics;

namespace OpenGJKSharp.Tests;

public class DoublePrecisionQueryTests
{
    [Fact]
    public void Should_ResolveGap_When_CoordinatesAreTooLargeForSinglePrecision()
    {
        // At 1e6 the spacing of float is about 6.25e-2, so a 1e-3 gap vanishes in the
        // Vector3 overloads. This was the one thing only the obsolete CsFunction could do.
        double[] a = [1e6, 0, 0];
        double[] b = [1e6 + 1e-3, 0, 0];

        double distance = OpenGjk.ComputeMinimumDistance(a, b);

        Assert.Equal(1e-3, distance, 10);
        Assert.Equal(0.0, OpenGjk.ComputeMinimumDistance([new Vector3(1e6f, 0, 0)], [new Vector3((float)(1e6 + 1e-3), 0, 0)]));
    }

    [Fact]
    public void Should_KeepDoublePrecision_When_DistanceIsNotRepresentableAsFloat()
    {
        double[] a = [0, 0, 0];
        double[] b = [1.0 / 3.0, 0, 0];

        double distance = OpenGjk.ComputeMinimumDistance(a, b);

        // Bit-identical to the double 1/3, which differs from the float 0.33333334.
        Assert.Equal(1.0 / 3.0, distance);
    }

    [Fact]
    public void Should_MatchObsoleteCsFunction_When_SameBodiesAreQueried()
    {
        var random = new Random(20260804);

        for (int i = 0; i < 100; i++)
        {
            int countA = 1 + random.Next(8);
            int countB = 1 + random.Next(8);
            double[] a = RandomCoordinates(random, countA);
            double[] b = RandomCoordinates(random, countB);

#pragma warning disable CS0618 // the replaced API is exactly what this test compares against
            double expected = OpenGjk.CsFunction(countA, ToLegacyLayout(a), countB, ToLegacyLayout(b));
#pragma warning restore CS0618
            double actual = OpenGjk.ComputeMinimumDistance(a, b);

            Assert.Equal(expected, actual);
        }
    }

    [Fact]
    public void Should_MatchVectorOverload_When_CoordinatesAreExactlyRepresentable()
    {
        // Halves and quarters survive the float conversion, so both paths must agree exactly.
        double[] a = [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1];
        double[] b = [3.5, 0, 0, 4.5, 0, 0, 3.5, 1, 0, 3.5, 0, 1];
        Vector3[] va = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(0, 0, 1)];
        Vector3[] vb = [new(3.5f, 0, 0), new(4.5f, 0, 0), new(3.5f, 1, 0), new(3.5f, 0, 1)];

        Assert.Equal(OpenGjk.ComputeMinimumDistance(va, vb), OpenGjk.ComputeMinimumDistance(a, b));
        Assert.Equal(OpenGjk.HasCollision(va, vb), OpenGjk.HasCollision(a, b));

        double expectedInfo = OpenGjk.ComputeCollisionInformation(va, vb, out Vector3 expectedNormal);
        Span<double> normal = stackalloc double[3];
        double actualInfo = OpenGjk.ComputeCollisionInformation(a, b, normal);

        Assert.Equal(expectedInfo, actualInfo);
        Assert.Equal(expectedNormal.X, (float)normal[0]);
        Assert.Equal(expectedNormal.Y, (float)normal[1]);
        Assert.Equal(expectedNormal.Z, (float)normal[2]);
    }

    [Fact]
    public void Should_ReportWitnessPoints_When_BodiesAreSeparated()
    {
        double[] a = [0, 0, 0];
        double[] b = [2.5, 0, 0, 3.5, 0, 0];

        Span<double> closestA = stackalloc double[3];
        Span<double> closestB = stackalloc double[3];
        double distance = OpenGjk.ComputeMinimumDistance(a, b, closestA, closestB);

        Assert.Equal(2.5, distance, 10);
        Assert.Equal(0.0, closestA[0], 10);
        Assert.Equal(2.5, closestB[0], 10);
    }

    [Fact]
    public void Should_ReportPenetration_When_BodiesOverlap()
    {
        // Two unit cubes sharing their y/z extent and overlapping by 0.25 along x,
        // so the minimum translation is unambiguously along x.
        double[] a = Cube(0, 1);
        double[] b = Cube(0.75, 1.75);


        Span<double> normal = stackalloc double[3];
        double info = OpenGjk.ComputeCollisionInformation(a, b, normal);

        Assert.True(info < 0, $"expected collision, got {info}");
        Assert.Equal(0.25, -info, 6);
        Assert.Equal(1.0, Math.Abs(normal[0]), 6);
    }

    [Fact]
    public void Should_MatchAllocatingOverload_When_WorkspaceIsReused()
    {
        var workspace = new GjkWorkspace();
        var random = new Random(20260805);

        for (int i = 0; i < 100; i++)
        {
            double[] a = RandomCoordinates(random, 1 + random.Next(8));
            double[] b = RandomCoordinates(random, 1 + random.Next(8));

            Assert.Equal(OpenGjk.ComputeMinimumDistance(a, b), OpenGjk.ComputeMinimumDistance(a, b, workspace));
            Assert.Equal(OpenGjk.HasCollision(a, b), OpenGjk.HasCollision(a, b, workspace));
        }
    }

    [Fact]
    public void Should_Throw_When_CoordinatesAreInvalid()
    {
        double[] valid = [0, 0, 0];

        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance([], valid));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(valid, [1, 2]));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(valid, [double.NaN, 0, 0]));
        Assert.Throws<ArgumentException>(() => OpenGjk.ComputeMinimumDistance(valid, [double.PositiveInfinity, 0, 0]));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeMinimumDistance(valid, valid, null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => OpenGjk.HasCollision(valid, valid, -1));

        Assert.Throws<ArgumentException>(() =>
        {
            Span<double> tooShort = stackalloc double[2];
            OpenGjk.ComputeCollisionInformation(valid, valid, tooShort);
        });
    }

    // A unit box spanning [minX, maxX] along x and [0, 1] along y and z.
    private static double[] Cube(double minX, double maxX)
    {
        return
        [
            minX, 0, 0, maxX, 0, 0, minX, 1, 0, maxX, 1, 0,
            minX, 0, 1, maxX, 0, 1, minX, 1, 1, maxX, 1, 1,
        ];
    }

    private static double[] RandomCoordinates(Random random, int count)
    {
        var coordinates = new double[count * 3];
        for (int i = 0; i < coordinates.Length; i++)
        {
            coordinates[i] = random.NextDouble() * 6 - 3;
        }

        return coordinates;
    }

    private static double[,] ToLegacyLayout(double[] coordinates)
    {
        int count = coordinates.Length / 3;
        var legacy = new double[3, count];
        for (int i = 0; i < count; i++)
        {
            legacy[0, i] = coordinates[i * 3];
            legacy[1, i] = coordinates[(i * 3) + 1];
            legacy[2, i] = coordinates[(i * 3) + 2];
        }

        return legacy;
    }
}

using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

/// <summary>
/// Contact points reported by the <c>TryGetPenetration</c> overloads with <c>pointA</c>/<c>pointB</c>:
/// each on its own body's surface, separated along the normal by the depth.
/// </summary>
public class ContactPointTests
{
    private const double _tolerance = 1e-5;

    private static readonly Quaternion _rotation =
        Quaternion.CreateFromAxisAngle(Vector3.Normalize(new Vector3(1, 2, 3)), 0.4f);

    private static readonly Dictionary<string, (ISupportMappable A, ISupportMappable B)> _shapePairs = new()
    {
        ["box-box"] = (new Box(Vector3.Zero, Vector3.One), new Box(new Vector3(1.5f, 0.3f, 0.2f), Vector3.One)),
        ["rotated-box-box"] = (
            new TransformedShape(new Box(Vector3.Zero, new Vector3(1, 0.5f, 0.5f)), _rotation, Vector3.Zero),
            new Box(new Vector3(1.2f, 0.1f, 0), new Vector3(0.5f, 0.5f, 0.5f))),
        ["sphere-sphere"] = (new Sphere(Vector3.Zero, 1), new Sphere(new Vector3(1.5f, 0.2f, -0.1f), 1)),
        ["sphere-box"] = (new Sphere(new Vector3(0.3f, 1.8f, 0), 1), new Box(Vector3.Zero, new Vector3(2, 1, 2))),
        ["sphere-inside-box"] = (new Sphere(new Vector3(0.2f, 0, 0), 0.5f), new Box(Vector3.Zero, new Vector3(2, 2, 2))),
        ["capsule-box"] = (new Capsule(new Vector3(-1, 1.3f, 0), new Vector3(1, 1.3f, 0), 0.5f),
            new Box(Vector3.Zero, new Vector3(2, 1, 2))),
        ["capsule-sphere"] = (new Capsule(new Vector3(-1, 0, 0), new Vector3(1, 0, 0), 0.5f),
            new Sphere(new Vector3(0.3f, 1.2f, 0), 1)),
        ["hull-hull"] = (
            new ConvexHull([new(1, 0, 0), new(-1, 0, 0), new(0, 1, 0), new(0, -1, 0), new(0, 0, 1), new(0, 0, -1)]),
            new ConvexHull([new(0.6f, 0.1f, 0.1f), new(2, 0, 0), new(1, 1.5f, 0), new(1, 0, 1.5f)])),
        ["transformed-sphere-box"] = (new TransformedShape(new Sphere(Vector3.Zero, 1), _rotation, new Vector3(0, 1.7f, 0)),
            new Box(Vector3.Zero, new Vector3(2, 1, 2))),
        ["concentric-spheres"] = (new Sphere(Vector3.Zero, 1), new Sphere(Vector3.Zero, 0.5f)),
        ["coplanar-hulls"] = (
            new ConvexHull([new(0, 0, 0), new(1, 0, 0), new(1, 1, 0), new(0, 1, 0)]),
            new ConvexHull([new(0.75f, 0, 0), new(1.75f, 0, 0), new(1.75f, 1, 0), new(0.75f, 1, 0)])),
        ["collinear-hulls"] = (new ConvexHull([new(0, 0, 0), new(2, 0, 0)]), new ConvexHull([new(1, 0, 0), new(3, 0, 0)])),
    };

    public static TheoryData<string> ShapePairs => [.. _shapePairs.Keys];

    private static readonly Dictionary<string, (Vector3[] A, Vector3[] B)> _vertexPairs = new()
    {
        ["cube-cube"] = (Cube(Vector3.Zero, 1), Cube(new Vector3(1.5f, 0.3f, 0.2f), 1)),
        ["cube-tetrahedron"] = (Cube(Vector3.Zero, 1), [new(0.7f, 0.2f, 0.1f), new(2, 0, 0), new(1, 1.5f, 0), new(1, 0, 1.5f)]),
        ["coplanar-quads"] = (
            [new(0, 0, 0), new(1, 0, 0), new(1, 1, 0), new(0, 1, 0)],
            [new(0.75f, 0, 0), new(1.75f, 0, 0), new(1.75f, 1, 0), new(0.75f, 1, 0)]),
        ["collinear-segments"] = ([new(0, 0, 0), new(2, 0, 0)], [new(1, 0, 0), new(3, 0, 0)]),
        ["nested-segments"] = ([new(0, 0, 0), new(10, 0, 0)], [new(4, 0, 0), new(6, 0, 0)]),
    };

    public static TheoryData<string> VertexPairs => [.. _vertexPairs.Keys];

    [Theory]
    [MemberData(nameof(ShapePairs))]
    public void Should_PlaceContactPointsOnSurfaces_When_ShapesOverlap(string pair)
    {
        (ISupportMappable a, ISupportMappable b) = _shapePairs[pair];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal,
            out Vector3 pointA, out Vector3 pointB));

        Assert.True(depth > 0, $"expected a positive depth, got {depth}");
        AssertDepthRelation(depth, normal, pointA, pointB);
        AssertOnSurface(a, pointA, normal);
        AssertOnSurface(b, pointB, -normal);
    }

    [Theory]
    [MemberData(nameof(VertexPairs))]
    public void Should_PlaceContactPointsOnSurfaces_When_VertexArraysOverlap(string pair)
    {
        (Vector3[] a, Vector3[] b) = _vertexPairs[pair];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal,
            out Vector3 pointA, out Vector3 pointB));

        Assert.True(depth > 0, $"expected a positive depth, got {depth}");
        AssertDepthRelation(depth, normal, pointA, pointB);
        AssertOnSurface(a, pointA, normal);
        AssertOnSurface(b, pointB, -normal);
    }

    [Fact]
    public void Should_FaceEachOtherInsideContactRegion_When_BoxFacesOverlap()
    {
        // The overlap is 0.5 deep along x and spans y in [-0.7, 1], z in [-0.8, 1].
        (ISupportMappable a, ISupportMappable b) = _shapePairs["box-box"];

        Assert.True(OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal,
            out Vector3 pointA, out Vector3 pointB));

        Assert.Equal(0.5, depth, 6);
        Assert.Equal(1.0, normal.X, 6);
        Assert.Equal(1.0, pointA.X, 6);
        Assert.Equal(0.5, pointB.X, 6);
        Assert.Equal(pointA.Y, pointB.Y, 6);
        Assert.Equal(pointA.Z, pointB.Z, 6);
        Assert.InRange(pointA.Y, -0.7f - 1e-6f, 1 + 1e-6f);
        Assert.InRange(pointA.Z, -0.8f - 1e-6f, 1 + 1e-6f);
    }

    [Fact]
    public void Should_ReturnPointsOnSpheres_When_SpheresOverlap()
    {
        var a = new Sphere(Vector3.Zero, 1);
        var b = new Sphere(new Vector3(1.5f, 0.2f, -0.1f), 1);

        Assert.True(OpenGjk.TryGetPenetration(a, b, out _, out _, out Vector3 pointA, out Vector3 pointB));

        Assert.Equal(1.0, (pointA - a.Center).Length(), 5);
        Assert.Equal(1.0, (pointB - b.Center).Length(), 5);
    }

    [Fact]
    public void Should_ReportSameDepthAndNormal_As_OverloadsWithoutPoints()
    {
        foreach ((ISupportMappable a, ISupportMappable b) in _shapePairs.Values)
        {
            bool hit = OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal);
            Assert.Equal(hit, OpenGjk.TryGetPenetration(a, b, out double depthWithPoints,
                out Vector3 normalWithPoints, out _, out _));
            Assert.Equal(depth, depthWithPoints);
            Assert.Equal(normal, normalWithPoints);
        }

        foreach ((Vector3[] a, Vector3[] b) in _vertexPairs.Values)
        {
            bool hit = OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal);
            Assert.Equal(hit, OpenGjk.TryGetPenetration(a, b, out double depthWithPoints,
                out Vector3 normalWithPoints, out _, out _));
            Assert.Equal(depth, depthWithPoints);
            Assert.Equal(normal, normalWithPoints);
        }
    }

    [Fact]
    public void Should_PlaceContactPointsOnSurfaces_When_PolygonsOverlap()
    {
        var random = new Random(20261101);

        for (int i = 0; i < 200; i++)
        {
            Vector2[] a = RandomConvexPolygon(random, Vector2.Zero);
            Vector2[] b = RandomConvexPolygon(random, new Vector2(
                (float)(random.NextDouble() * 2 - 1),
                (float)(random.NextDouble() * 2 - 1)));

            if (!OpenGjk.TryGetPenetration(a, b, out double depth, out Vector2 normal,
                out Vector2 pointA, out Vector2 pointB) || depth == 0)
            {
                continue;
            }

            AssertDepthRelation(depth, new Vector3(normal, 0), new Vector3(pointA, 0), new Vector3(pointB, 0));
            AssertOnSurface(a, pointA, normal);
            AssertOnSurface(b, pointB, -normal);
        }
    }

    [Fact]
    public void Should_PlaceContactPointsOnSurfaces_When_DoubleCoordinatesOverlap()
    {
        double[] a = Flatten(Cube(Vector3.Zero, 1));
        double[] b = Flatten(Cube(new Vector3(1.5f, 0.3f, 0.2f), 1));
        double[] normal = new double[3];
        double[] pointA = new double[3];
        double[] pointB = new double[3];

        Assert.True(OpenGjk.TryGetPenetration(a, b, normal, pointA, pointB, out double depth));

        Assert.Equal(0.5, depth, 12);
        double separation = 0;
        for (int c = 0; c < 3; c++)
        {
            separation += (pointA[c] - pointB[c]) * normal[c];
        }

        Assert.Equal(depth, separation, 12);
        Assert.Equal(1.0, pointA[0], 12);
        Assert.Equal(0.5, pointB[0], 12);
    }

    [Fact]
    public void Should_ReturnClosestPoints_When_BodiesAreApart()
    {
        var sphere = new Sphere(Vector3.Zero, 1);
        var box = new Box(new Vector3(4, 0.5f, 0), Vector3.One);
        OpenGjk.ComputeMinimumDistance(sphere, box, out Vector3 closestA, out Vector3 closestB);
        Assert.False(OpenGjk.TryGetPenetration(sphere, box, out double depth, out _, out Vector3 pointA, out Vector3 pointB));
        Assert.Equal(0.0, depth);
        Assert.Equal(closestA, pointA);
        Assert.Equal(closestB, pointB);

        Vector3[] cubeA = Cube(Vector3.Zero, 1);
        Vector3[] cubeB = Cube(new Vector3(3, 1, 0.5f), 1);
        OpenGjk.ComputeMinimumDistance(cubeA, cubeB, out closestA, out closestB);
        Assert.False(OpenGjk.TryGetPenetration(cubeA, cubeB, out _, out _, out pointA, out pointB));
        Assert.Equal(closestA, pointA);
        Assert.Equal(closestB, pointB);

        Vector2[] squareA = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] squareB = [new(2, 0.5f), new(3, 0.5f), new(3, 1.5f), new(2, 1.5f)];
        OpenGjk.ComputeMinimumDistance(squareA, squareB, out Vector2 closest2A, out Vector2 closest2B);
        Assert.False(OpenGjk.TryGetPenetration(squareA, squareB, out _, out _, out Vector2 point2A, out Vector2 point2B));
        Assert.Equal(closest2A, point2A);
        Assert.Equal(closest2B, point2B);

        double[] coordsA = Flatten(cubeA);
        double[] coordsB = Flatten(cubeB);
        double[] expectedA = new double[3];
        double[] expectedB = new double[3];
        double[] actualA = new double[3];
        double[] actualB = new double[3];
        OpenGjk.ComputeMinimumDistance(coordsA, coordsB, expectedA, expectedB);
        Assert.False(OpenGjk.TryGetPenetration(coordsA, coordsB, new double[3], actualA, actualB, out _));
        Assert.Equal(expectedA, actualA);
        Assert.Equal(expectedB, actualB);
    }

    [Fact]
    public void Should_MatchAllocatingOverload_When_WorkspaceIsReused()
    {
        var workspace = new GjkWorkspace();

        foreach ((ISupportMappable a, ISupportMappable b) in _shapePairs.Values)
        {
            OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal, out Vector3 pointA, out Vector3 pointB);
            OpenGjk.TryGetPenetration(a, b, workspace, out double depthW, out Vector3 normalW,
                out Vector3 pointAW, out Vector3 pointBW);
            Assert.Equal((depth, normal, pointA, pointB), (depthW, normalW, pointAW, pointBW));
        }

        foreach ((Vector3[] a, Vector3[] b) in _vertexPairs.Values)
        {
            OpenGjk.TryGetPenetration(a, b, out double depth, out Vector3 normal, out Vector3 pointA, out Vector3 pointB);
            OpenGjk.TryGetPenetration(a, b, workspace, out double depthW, out Vector3 normalW,
                out Vector3 pointAW, out Vector3 pointBW);
            Assert.Equal((depth, normal, pointA, pointB), (depthW, normalW, pointAW, pointBW));

            double[] coordsA = Flatten(a);
            double[] coordsB = Flatten(b);
            double[] n = new double[3], pa = new double[3], pb = new double[3];
            double[] nW = new double[3], paW = new double[3], pbW = new double[3];
            OpenGjk.TryGetPenetration(coordsA, coordsB, n, pa, pb, out double depthD);
            OpenGjk.TryGetPenetration(coordsA, coordsB, workspace, nW, paW, pbW, out double depthDW);
            Assert.Equal(depthD, depthDW);
            Assert.Equal(n, nW);
            Assert.Equal(pa, paW);
            Assert.Equal(pb, pbW);
        }

        var random = new Random(7);
        for (int i = 0; i < 20; i++)
        {
            Vector2[] a = RandomConvexPolygon(random, Vector2.Zero);
            Vector2[] b = RandomConvexPolygon(random, new Vector2(0.5f, 0.2f));
            OpenGjk.TryGetPenetration(a, b, out double depth, out Vector2 normal, out Vector2 pointA, out Vector2 pointB);
            OpenGjk.TryGetPenetration(a, b, workspace, out double depthW, out Vector2 normalW,
                out Vector2 pointAW, out Vector2 pointBW);
            Assert.Equal((depth, normal, pointA, pointB), (depthW, normalW, pointAW, pointBW));
        }
    }

    [Fact]
    public void Should_Throw_When_ArgumentsAreInvalid()
    {
        double[] coords = [0, 0, 0];
        double[] output = new double[3];

        Assert.Throws<ArgumentNullException>(
            () => OpenGjk.TryGetPenetration(coords, coords, null!, output, output, output, out _));
        Assert.Throws<ArgumentException>(
            () => OpenGjk.TryGetPenetration(coords, coords, output, new double[2], output, out _));
        Assert.Throws<ArgumentException>(
            () => OpenGjk.TryGetPenetration(coords, coords, output, output, new double[2], out _));
        Assert.Throws<ArgumentNullException>(
            () => OpenGjk.TryGetPenetration(null!, new Sphere(Vector3.Zero, 1), out _, out _, out _, out _));
    }

#if NET8_0_OR_GREATER
    [Fact]
    public void Should_StopAllocating_When_ContactQueriesReuseAWorkspace()
    {
        var workspace = new GjkWorkspace();
        (ISupportMappable sphere, ISupportMappable box) = _shapePairs["sphere-box"];
        Vector3[] cubeA = Cube(Vector3.Zero, 1);
        Vector3[] cubeB = Cube(new Vector3(1.5f, 0.3f, 0.2f), 1);
        Vector2[] squareA = [new(0, 0), new(1, 0), new(1, 1), new(0, 1)];
        Vector2[] squareB = [new(0.75f, 0.2f), new(1.75f, 0.2f), new(1.75f, 1.2f), new(0.75f, 1.2f)];
        double[] coordsA = Flatten(cubeA);
        double[] coordsB = Flatten(cubeB);
        double[] normal = new double[3], pointA = new double[3], pointB = new double[3];

        void Run()
        {
            OpenGjk.TryGetPenetration(sphere, box, workspace, out _, out _, out _, out _);
            OpenGjk.TryGetPenetration(cubeA, cubeB, workspace, out _, out _, out _, out _);
            OpenGjk.TryGetPenetration(squareA, squareB, workspace, out _, out _, out _, out _);
            OpenGjk.TryGetPenetration(coordsA, coordsB, workspace, normal, pointA, pointB, out _);
        }

        for (int i = 0; i < 20; i++)
        {
            Run();
        }

        long before = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < 100; i++)
        {
            Run();
        }
        long bytes = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.True(bytes == 0, $"expected no allocation with a reused workspace, got {bytes} bytes");
    }
#endif

    private static void AssertDepthRelation(double depth, Vector3 normal, Vector3 pointA, Vector3 pointB)
    {
        double separation = Vector3.Dot(pointA - pointB, normal);
        Assert.True(Math.Abs(separation - depth) <= _tolerance * (1 + depth),
            $"dot(pointA - pointB, normal) = {separation}, depth = {depth}");
    }

    // On the support plane facing outward and inside the body, which together mean on its surface.
    private static void AssertOnSurface(ISupportMappable shape, Vector3 point, Vector3 outward)
    {
        float plane = Vector3.Dot(shape.GetSupportPoint(outward), outward);
        Assert.Equal(plane, Vector3.Dot(point, outward), _tolerance);
        Assert.True(OpenGjk.ComputeMinimumDistance(new ConvexHull([point]), shape) <= _tolerance,
            $"{point} is outside the shape");
    }

    private static void AssertOnSurface(Vector3[] body, Vector3 point, Vector3 outward)
    {
        float plane = body.Max(v => Vector3.Dot(v, outward));
        Assert.Equal(plane, Vector3.Dot(point, outward), _tolerance);
        Assert.True(OpenGjk.ComputeMinimumDistance([point], body) <= _tolerance, $"{point} is outside the body");
    }

    private static void AssertOnSurface(Vector2[] polygon, Vector2 point, Vector2 outward)
    {
        float line = polygon.Max(v => Vector2.Dot(v, outward));
        Assert.Equal(line, Vector2.Dot(point, outward), _tolerance);
        Assert.True(OpenGjk.ComputeMinimumDistance([point], polygon) <= _tolerance, $"{point} is outside the polygon");
    }

    private static Vector3[] Cube(Vector3 center, float half)
    {
        var vertices = new Vector3[8];
        for (int i = 0; i < 8; i++)
        {
            vertices[i] = center + new Vector3(
                (i & 1) == 0 ? -half : half,
                (i & 2) == 0 ? -half : half,
                (i & 4) == 0 ? -half : half);
        }

        return vertices;
    }

    private static double[] Flatten(Vector3[] vertices)
    {
        return [.. vertices.SelectMany(v => new double[] { v.X, v.Y, v.Z })];
    }

    private static Vector2[] RandomConvexPolygon(Random random, Vector2 center)
    {
        // Points on a circle in angle order are convex by construction.
        int count = random.Next(3, 9);
        double[] angles = [.. Enumerable.Range(0, count).Select(_ => random.NextDouble() * 2 * Math.PI).OrderBy(x => x)];
        double radius = 0.5 + random.NextDouble();

        return [.. angles.Select(t => center + new Vector2((float)(radius * Math.Cos(t)), (float)(radius * Math.Sin(t))))];
    }
}

using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class GjkWorkspaceTests
{
    [Fact]
    public void Should_MatchAllocatingOverload_When_ShapeQueriesReuseWorkspace()
    {
        var workspace = new GjkWorkspace();
        var random = new Random(20260801);

        for (int i = 0; i < 200; i++)
        {
            var a = new Sphere(RandomPoint(random), 0.5f + (float)random.NextDouble());
            var b = new Box(RandomPoint(random), new Vector3(0.4f, 0.7f, 1.1f));

            double expected = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 expectedA, out Vector3 expectedB);
            double actual = OpenGjk.ComputeMinimumDistance(a, b, workspace, out Vector3 actualA, out Vector3 actualB);

            Assert.Equal(expected, actual);
            Assert.Equal(expectedA, actualA);
            Assert.Equal(expectedB, actualB);
            Assert.Equal(OpenGjk.HasCollision(a, b), OpenGjk.HasCollision(a, b, workspace));

            double expectedInfo = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 expectedNormal);
            double actualInfo = OpenGjk.ComputeCollisionInformation(a, b, workspace, out Vector3 actualNormal);

            Assert.Equal(expectedInfo, actualInfo);
            Assert.Equal(expectedNormal, actualNormal);
        }
    }

    [Fact]
    public void Should_MatchAllocatingOverload_When_VertexQueriesReuseWorkspace()
    {
        var workspace = new GjkWorkspace();
        var random = new Random(20260802);

        for (int i = 0; i < 200; i++)
        {
            // Alternating vertex counts also exercise growing and shrinking the reused buffers.
            Vector3[] a = RandomCloud(random, 1 + (i % 8));
            Vector3[] b = RandomCloud(random, 1 + ((i * 3) % 8));

            double expected = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 expectedA, out Vector3 expectedB);
            double actual = OpenGjk.ComputeMinimumDistance(a, b, workspace, out Vector3 actualA, out Vector3 actualB);

            Assert.Equal(expected, actual);
            Assert.Equal(expectedA, actualA);
            Assert.Equal(expectedB, actualB);
            Assert.Equal(OpenGjk.HasCollision(a, b), OpenGjk.HasCollision(a, b, workspace));

            double expectedInfo = OpenGjk.ComputeCollisionInformation(a, b, out Vector3 expectedNormal);
            double actualInfo = OpenGjk.ComputeCollisionInformation(a, b, workspace, out Vector3 actualNormal);

            Assert.Equal(expectedInfo, actualInfo);
            Assert.Equal(expectedNormal, actualNormal);
        }
    }

    [Fact]
    public void Should_MatchAllocatingOverload_When_PolygonQueriesReuseWorkspace()
    {
        var workspace = new GjkWorkspace();
        Vector2[] a = [new(0, 0), new(1, 0), new(0, 1)];
        Vector2[] b = [new(3, 0), new(4, 0), new(3, 1)];

        double expected = OpenGjk.ComputeMinimumDistance(a, b, out Vector2 expectedA, out Vector2 expectedB);
        double actual = OpenGjk.ComputeMinimumDistance(a, b, workspace, out Vector2 actualA, out Vector2 actualB);

        Assert.Equal(expected, actual);
        Assert.Equal(expectedA, actualA);
        Assert.Equal(expectedB, actualB);
        Assert.Equal(OpenGjk.HasCollision(a, b), OpenGjk.HasCollision(a, b, workspace));
    }

    [Fact]
    public void Should_StayConsistent_When_WorkspaceIsSharedByDifferentQueryKinds()
    {
        var workspace = new GjkWorkspace();
        var sphere = new Sphere(new Vector3(0, 0, 0), 1f);
        var box = new Box(new Vector3(4, 0, 0), new Vector3(1, 1, 1));
        Vector3[] cloudA = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(0, 0, 1)];
        Vector3[] cloudB = [new(5, 0, 0), new(6, 0, 0), new(5, 1, 0), new(5, 0, 1)];

        // Interleaving the two body kinds must not leak state between queries.
        for (int i = 0; i < 5; i++)
        {
            Assert.Equal(2.0, OpenGjk.ComputeMinimumDistance(sphere, box, workspace), 4);
            Assert.Equal(4.0, OpenGjk.ComputeMinimumDistance(cloudA, cloudB, workspace), 4);
        }
    }

    [Fact]
    public void Should_Throw_When_WorkspaceIsNull()
    {
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var b = new Sphere(new Vector3(3, 0, 0), 1f);
        Vector3[] cloud = [new(0, 0, 0)];

        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeMinimumDistance(a, b, null!));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.HasCollision(a, b, null!));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeCollisionInformation(a, b, null!, out _));
        Assert.Throws<ArgumentNullException>(() => OpenGjk.ComputeMinimumDistance(cloud, cloud, null!));
    }

#if NET8_0_OR_GREATER
    [Fact]
    public void Should_StopAllocating_When_WorkspaceIsReusedAcrossQueries()
    {
        var workspace = new GjkWorkspace();
        var a = new Sphere(new Vector3(0, 0, 0), 1f);
        var separated = new Box(new Vector3(4, 0, 0), new Vector3(1, 1, 1));
        var overlapping = new Box(new Vector3(1.2f, 0, 0), new Vector3(1, 1, 1));
        Vector3[] cloudA = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(0, 0, 1)];
        Vector3[] cloudB = [new(5, 0, 0), new(6, 0, 0), new(5, 1, 0), new(5, 0, 1)];

        // Warm up: the first queries grow the buffers the workspace then keeps.
        for (int i = 0; i < 20; i++)
        {
            OpenGjk.ComputeMinimumDistance(a, separated, workspace);
            OpenGjk.ComputeCollisionInformation(a, overlapping, workspace, out _);
            OpenGjk.ComputeMinimumDistance(cloudA, cloudB, workspace);
            OpenGjk.HasCollision(a, separated, workspace);
            OpenGjk.HasCollision(cloudA, cloudB, workspace);
        }

        const int queries = 100;
        long before = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < queries; i++)
        {
            OpenGjk.ComputeMinimumDistance(a, separated, workspace);
            OpenGjk.ComputeCollisionInformation(a, overlapping, workspace, out _);
            OpenGjk.ComputeMinimumDistance(cloudA, cloudB, workspace);
            OpenGjk.HasCollision(a, separated, workspace);
            OpenGjk.HasCollision(cloudA, cloudB, workspace);
        }
        long reusedBytes = GC.GetAllocatedBytesForCurrentThread() - before;

        before = GC.GetAllocatedBytesForCurrentThread();
        for (int i = 0; i < queries; i++)
        {
            OpenGjk.ComputeMinimumDistance(a, separated);
            OpenGjk.ComputeCollisionInformation(a, overlapping, out _);
            OpenGjk.ComputeMinimumDistance(cloudA, cloudB);
        }
        long allocatingBytes = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.True(reusedBytes == 0, $"expected no allocation with a reused workspace, got {reusedBytes} bytes");
        Assert.True(allocatingBytes > 0, "the allocating overloads are expected to allocate");
    }
#endif

    private static Vector3 RandomPoint(Random random)
    {
        return new Vector3(
            (float)(random.NextDouble() * 6 - 3),
            (float)(random.NextDouble() * 6 - 3),
            (float)(random.NextDouble() * 6 - 3));
    }

    private static Vector3[] RandomCloud(Random random, int count)
    {
        var points = new Vector3[count];
        for (int i = 0; i < count; i++)
        {
            points[i] = RandomPoint(random);
        }

        return points;
    }
}

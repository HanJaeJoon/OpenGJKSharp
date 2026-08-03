using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

public class GjkWorkspaceGuardTests
{
    private static readonly Vector3[] _cloudA = [new(0, 0, 0), new(1, 0, 0), new(0, 1, 0), new(0, 0, 1)];
    private static readonly Vector3[] _cloudB = [new(5, 0, 0), new(6, 0, 0), new(5, 1, 0), new(5, 0, 1)];
    private static readonly Vector2[] _polygonA = [new(0, 0), new(1, 0), new(0, 1)];
    private static readonly Vector2[] _polygonB = [new(5, 0), new(6, 0), new(5, 1)];
    private static readonly double[] _coordsA = [0, 0, 0, 1, 0, 0, 0, 1, 0];
    private static readonly double[] _coordsB = [5, 0, 0, 6, 0, 0, 5, 1, 0];

    /// <summary>
    /// Every public overload that takes a workspace, so that no entry point can be left unguarded.
    /// </summary>
    public static TheoryData<string, Action<GjkWorkspace>> WorkspaceOverloads()
    {
        var shapeA = new Sphere(new Vector3(0, 0, 0), 1f);
        var shapeB = new Box(new Vector3(5, 0, 0), new Vector3(1, 1, 1));

        return new TheoryData<string, Action<GjkWorkspace>>
        {
            { "HasCollision(Vector3[])", w => OpenGjk.HasCollision(_cloudA, _cloudB, w) },
            { "HasCollision(Vector2[])", w => OpenGjk.HasCollision(_polygonA, _polygonB, w) },
            { "HasCollision(ReadOnlySpan<double>)", w => OpenGjk.HasCollision(_coordsA, _coordsB, w) },
            { "HasCollision(ISupportMappable)", w => OpenGjk.HasCollision(shapeA, shapeB, w) },
            { "ComputeMinimumDistance(Vector3[])", w => OpenGjk.ComputeMinimumDistance(_cloudA, _cloudB, w) },
            {
                "ComputeMinimumDistance(Vector3[], witnesses)",
                w => OpenGjk.ComputeMinimumDistance(_cloudA, _cloudB, w, out Vector3 _, out Vector3 _)
            },
            { "ComputeMinimumDistance(Vector2[])", w => OpenGjk.ComputeMinimumDistance(_polygonA, _polygonB, w) },
            {
                "ComputeMinimumDistance(Vector2[], witnesses)",
                w => OpenGjk.ComputeMinimumDistance(_polygonA, _polygonB, w, out Vector2 _, out Vector2 _)
            },
            { "ComputeMinimumDistance(ReadOnlySpan<double>)", w => OpenGjk.ComputeMinimumDistance(_coordsA, _coordsB, w) },
            {
                "ComputeMinimumDistance(ReadOnlySpan<double>, witnesses)",
                w => OpenGjk.ComputeMinimumDistance(_coordsA, _coordsB, w, new double[3], new double[3])
            },
            { "ComputeMinimumDistance(ISupportMappable)", w => OpenGjk.ComputeMinimumDistance(shapeA, shapeB, w) },
            {
                "ComputeMinimumDistance(ISupportMappable, witnesses)",
                w => OpenGjk.ComputeMinimumDistance(shapeA, shapeB, w, out Vector3 _, out Vector3 _)
            },
            {
                "ComputeCollisionInformation(Vector3[])",
                w => OpenGjk.ComputeCollisionInformation(_cloudA, _cloudB, w, out Vector3 _)
            },
            {
                "ComputeCollisionInformation(Vector2[])",
                w => OpenGjk.ComputeCollisionInformation(_polygonA, _polygonB, w, out Vector2 _)
            },
            {
                "ComputeCollisionInformation(ReadOnlySpan<double>)",
                w => OpenGjk.ComputeCollisionInformation(_coordsA, _coordsB, w, new double[3])
            },
            {
                "ComputeCollisionInformation(ISupportMappable)",
                w => OpenGjk.ComputeCollisionInformation(shapeA, shapeB, w, out Vector3 _)
            },
            {
                "TryGetPenetration(Vector3[])",
                w => OpenGjk.TryGetPenetration(_cloudA, _cloudB, w, out double _, out Vector3 _)
            },
            {
                "TryGetPenetration(Vector2[])",
                w => OpenGjk.TryGetPenetration(_polygonA, _polygonB, w, out double _, out Vector2 _)
            },
            {
                "TryGetPenetration(ISupportMappable)",
                w => OpenGjk.TryGetPenetration(shapeA, shapeB, w, out double _, out Vector3 _)
            },
            {
                "TryGetPenetration(ReadOnlySpan<double>)",
                w => OpenGjk.TryGetPenetration(_coordsA, _coordsB, w, new double[3], out double _)
            },
        };
    }

    [Theory]
    [MemberData(nameof(WorkspaceOverloads))]
    public void Should_Throw_When_WorkspaceIsAlreadyInUse(string entryPoint, Action<GjkWorkspace> query)
    {
        var workspace = new GjkWorkspace();

        using var holder = new WorkspaceHolder(workspace);

        var exception = Assert.Throws<InvalidOperationException>(() => query(workspace));
        Assert.Contains("GjkWorkspace", exception.Message);
        Assert.Contains("thread", exception.Message, StringComparison.OrdinalIgnoreCase);

        // Releasing the holding query must leave the workspace usable again.
        holder.Release();
        query(workspace);

        Assert.NotEmpty(entryPoint);
    }

    [Fact]
    public void Should_Throw_When_QueryReentersWithTheSameWorkspace()
    {
        var workspace = new GjkWorkspace();
        var other = new Sphere(new Vector3(5, 0, 0), 1f);
        var reentrant = new CallbackShape(_ => OpenGjk.ComputeMinimumDistance(_cloudA, _cloudB, workspace));

        // A custom GetSupportPoint that queries with the workspace it was called from
        // is the single-threaded shape of the same misuse.
        Assert.Throws<InvalidOperationException>(
            () => OpenGjk.ComputeMinimumDistance(reentrant, other, workspace));

        // The outer query released the workspace on the way out.
        Assert.Equal(4.0, OpenGjk.ComputeMinimumDistance(_cloudA, _cloudB, workspace), 4);
    }

    [Fact]
    public void Should_ReleaseWorkspace_When_QueryThrows()
    {
        var workspace = new GjkWorkspace();
        var other = new Sphere(new Vector3(5, 0, 0), 1f);
        var faulty = new CallbackShape(_ => throw new NotSupportedException("boom"));

        Assert.Throws<NotSupportedException>(() => OpenGjk.ComputeMinimumDistance(faulty, other, workspace));

        // A query that failed inside the library must not leave the workspace locked.
        Assert.Equal(4.0, OpenGjk.ComputeMinimumDistance(_cloudA, _cloudB, workspace), 4);
    }

    [Fact]
    public void Should_StayCorrect_When_EachThreadHasItsOwnWorkspace()
    {
        const int threads = 8;
        const int queries = 500;

        double[] expected =
        [
            OpenGjk.ComputeMinimumDistance(_cloudA, _cloudB),
            OpenGjk.ComputeCollisionInformation(_cloudA, _cloudA, out _),
        ];

        Parallel.For(0, threads, thread =>
        {
            _ = thread;
            var workspace = new GjkWorkspace();
            for (int i = 0; i < queries; i++)
            {
                Assert.Equal(expected[0], OpenGjk.ComputeMinimumDistance(_cloudA, _cloudB, workspace));
                Assert.Equal(expected[1], OpenGjk.ComputeCollisionInformation(_cloudA, _cloudA, workspace, out _));
            }
        });
    }

    /// <summary>Occupies a workspace on a background thread until released.</summary>
    private sealed class WorkspaceHolder : IDisposable
    {
        private readonly ManualResetEventSlim _entered = new();
        private readonly ManualResetEventSlim _release = new();
        private readonly Task _query;

        public WorkspaceHolder(GjkWorkspace workspace)
        {
            var blocking = new CallbackShape(_ =>
            {
                _entered.Set();
                _release.Wait();
            });
            var other = new Sphere(new Vector3(5, 0, 0), 1f);

            _query = Task.Run(() => OpenGjk.ComputeMinimumDistance(blocking, other, workspace));

            Assert.True(_entered.Wait(TimeSpan.FromSeconds(30)), "the holding query never started");
        }

        public void Release()
        {
            _release.Set();
            Assert.True(_query.Wait(TimeSpan.FromSeconds(30)), "the holding query never finished");
        }

        public void Dispose()
        {
            _release.Set();
            _entered.Dispose();
            _release.Dispose();
        }
    }

    /// <summary>A point shape that runs an arbitrary action inside <c>GetSupportPoint</c>.</summary>
    private sealed class CallbackShape : ISupportMappable
    {
        private readonly Action<Vector3> _onSupport;

        public CallbackShape(Action<Vector3> onSupport) => _onSupport = onSupport;

        public Vector3 GetSupportPoint(Vector3 direction)
        {
            _onSupport(direction);
            return Vector3.Zero;
        }

        public Vector3 GetCenter() => Vector3.Zero;
    }
}

using OpenGJKSharp.Shapes;
using System.Numerics;

namespace OpenGJKSharp.Tests;

/// <summary>
/// HasCollision stops GJK early; its answer must still be the converged distance compared
/// against the precision.
/// </summary>
public class HasCollisionEquivalenceTests
{
    private static readonly double[] _precisions = [0, 1e-6, 0.25];

    [Fact]
    public void Should_MatchMinimumDistance_When_VertexArraysCollideTouchOrSeparate()
    {
        var random = new Random(20261101);
        var workspace = new GjkWorkspace();

        for (int i = 0; i < 300; i++)
        {
            Vector3[] a = RandomCloud(random, 4 + (i % 12), Origin(i));
            Vector3[] b = RandomCloud(random, 4 + ((i * 5) % 12), Origin(i) + RandomPoint(random, 4));
            double d = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 ca, out Vector3 cb);
            Vector3 offset = TouchingOffset(random, d, ca, cb, i);
            b = Array.ConvertAll(b, q => q + offset);

            double distance = OpenGjk.ComputeMinimumDistance(a, b);
            foreach (double p in _precisions)
            {
                Assert.Equal(distance <= p, OpenGjk.HasCollision(a, b, p));
                Assert.Equal(distance <= p, OpenGjk.HasCollision(a, b, workspace, p));
            }
        }
    }

    [Fact]
    public void Should_MatchMinimumDistance_When_PolygonsCollideTouchOrSeparate()
    {
        var random = new Random(20261102);
        var workspace = new GjkWorkspace();

        for (int i = 0; i < 300; i++)
        {
            Vector2[] a = Flatten(RandomCloud(random, 3 + (i % 10), Origin(i)));
            Vector2[] b = Flatten(RandomCloud(random, 3 + ((i * 7) % 10), Origin(i) + RandomPoint(random, 4)));
            double d = OpenGjk.ComputeMinimumDistance(a, b, out Vector2 ca, out Vector2 cb);
            Vector3 offset = TouchingOffset(random, d, new Vector3(ca, 0), new Vector3(cb, 0), i);
            b = Array.ConvertAll(b, q => q + new Vector2(offset.X, offset.Y));

            double distance = OpenGjk.ComputeMinimumDistance(a, b);
            foreach (double p in _precisions)
            {
                Assert.Equal(distance <= p, OpenGjk.HasCollision(a, b, p));
                Assert.Equal(distance <= p, OpenGjk.HasCollision(a, b, workspace, p));
            }
        }
    }

    [Fact]
    public void Should_MatchMinimumDistance_When_ShapesCollideTouchOrSeparate()
    {
        var random = new Random(20261103);
        var workspace = new GjkWorkspace();

        for (int i = 0; i < 300; i++)
        {
            ISupportMappable a = new TransformedShape(RandomShape(random, i), Quaternion.Identity, Origin(i));
            ISupportMappable b = new TransformedShape(RandomShape(random, i + 1), Quaternion.Identity, Origin(i) + RandomPoint(random, 4));
            double d = OpenGjk.ComputeMinimumDistance(a, b, out Vector3 ca, out Vector3 cb);
            b = new TransformedShape(b, Quaternion.Identity, TouchingOffset(random, d, ca, cb, i));

            double distance = OpenGjk.ComputeMinimumDistance(a, b);
            foreach (double p in _precisions)
            {
                Assert.Equal(distance <= p, OpenGjk.HasCollision(a, b, p));
                Assert.Equal(distance <= p, OpenGjk.HasCollision(a, b, workspace, p));
            }
        }
    }

    // Rounding grows with the coordinates, so some pairs sit far from the origin.
    private static Vector3 Origin(int i)
    {
        return i % 5 == 4 ? new Vector3(1000, -500, 250) : Vector3.Zero;
    }

    // Every third pair is moved until its gap is close to one of the precisions, so the
    // early-exit bounds are exercised right at the threshold as well as far from it.
    private static Vector3 TouchingOffset(Random random, double distance, Vector3 closestA, Vector3 closestB, int i)
    {
        if (i % 3 != 0 || distance <= 0)
        {
            return Vector3.Zero;
        }

        double[] gaps = [0, 5e-7, 1e-6, 2e-6, 0.25];
        double gap = gaps[random.Next(gaps.Length)];
        return (closestA - closestB) * (float)(1 - gap / distance);
    }

    private static ISupportMappable RandomShape(Random random, int kind)
    {
        Vector3 center = RandomPoint(random, 0.5);
        float size = 0.3f + (float)random.NextDouble();
        return (kind % 4) switch
        {
            0 => new Sphere(center, size),
            1 => new Box(center, new Vector3(size, 0.5f * size, 0.8f), Quaternion.Normalize(new Quaternion(RandomPoint(random, 1), 1))),
            2 => new Capsule(center, center + RandomPoint(random, 1), 0.5f * size),
            _ => new ConvexHull(RandomCloud(random, 8, center)),
        };
    }

    private static Vector3 RandomPoint(Random random, double scale)
    {
        return new Vector3(
            (float)((random.NextDouble() * 2 - 1) * scale),
            (float)((random.NextDouble() * 2 - 1) * scale),
            (float)((random.NextDouble() * 2 - 1) * scale));
    }

    private static Vector3[] RandomCloud(Random random, int count, Vector3 center)
    {
        var points = new Vector3[count];
        for (int i = 0; i < count; i++)
        {
            points[i] = center + RandomPoint(random, 1);
        }

        return points;
    }

    private static Vector2[] Flatten(Vector3[] points)
    {
        return Array.ConvertAll(points, q => new Vector2(q.X, q.Y));
    }
}

using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// The convex hull of a set of world-space points.
/// The hull is never computed explicitly; the support mapping scans the vertices.
/// </summary>
public sealed class ConvexHull : ISupportMappable
{
    private readonly Vector3[] _points;
    private readonly Vector3 _centroid;

    /// <summary>
    /// Creates a convex hull from a set of points. The array is copied.
    /// </summary>
    /// <param name="points">The points spanning the hull. Must contain at least one finite point.</param>
    /// <exception cref="ArgumentNullException"><paramref name="points"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="points"/> is empty or contains a non-finite point.</exception>
    public ConvexHull(Vector3[] points)
    {
        if (points is null)
        {
            throw new ArgumentNullException(nameof(points));
        }

        if (points.Length == 0)
        {
            throw new ArgumentException("At least one point is required.", nameof(points));
        }

        var copy = new Vector3[points.Length];
        var sum = Vector3.Zero;
        for (int i = 0; i < points.Length; i++)
        {
            ShapeArgument.EnsureFinite(points[i], nameof(points), $"Point at index {i}");
            copy[i] = points[i];
            sum += points[i];
        }

        _points = copy;
        _centroid = sum / points.Length;
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        Vector3 best = _points[0];
        float bestDot = Vector3.Dot(best, direction);

        for (int i = 1; i < _points.Length; i++)
        {
            float dot = Vector3.Dot(_points[i], direction);
            if (dot > bestDot)
            {
                bestDot = dot;
                best = _points[i];
            }
        }

        return best;
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return _centroid;
    }
}

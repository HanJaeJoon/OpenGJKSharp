using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// A capsule: the set of points within <see cref="Radius"/> of the world-space
/// segment from <see cref="P0"/> to <see cref="P1"/>.
/// </summary>
public sealed class Capsule : ISupportMappable
{
    /// <summary>The first endpoint of the capsule axis.</summary>
    public Vector3 P0 { get; }

    /// <summary>The second endpoint of the capsule axis.</summary>
    public Vector3 P1 { get; }

    /// <summary>The radius of the capsule.</summary>
    public float Radius { get; }

    /// <summary>
    /// Creates a capsule.
    /// </summary>
    /// <param name="p0">The first endpoint of the axis. Must be finite.</param>
    /// <param name="p1">The second endpoint of the axis. Must be finite.
    /// May equal <paramref name="p0"/>, in which case the capsule degenerates to a sphere.</param>
    /// <param name="radius">The radius of the capsule. Must be positive and finite.</param>
    /// <exception cref="ArgumentException"><paramref name="p0"/> or <paramref name="p1"/> is not finite.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="radius"/> is not positive and finite.</exception>
    public Capsule(Vector3 p0, Vector3 p1, float radius)
    {
        ShapeArgument.EnsureFinite(p0, nameof(p0), "P0");
        ShapeArgument.EnsureFinite(p1, nameof(p1), "P1");
        ShapeArgument.EnsurePositiveFinite(radius, nameof(radius), "Radius");

        P0 = p0;
        P1 = p1;
        Radius = radius;
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        Vector3 endpoint = Vector3.Dot(P0, direction) >= Vector3.Dot(P1, direction) ? P0 : P1;

        return endpoint + Radius * Vector3.Normalize(direction);
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return (P0 + P1) * 0.5f;
    }
}

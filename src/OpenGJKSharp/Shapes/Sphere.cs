using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// A sphere defined by a world-space center and radius.
/// </summary>
public sealed class Sphere : ISupportMappable
{
    /// <summary>The center of the sphere.</summary>
    public Vector3 Center { get; }

    /// <summary>The radius of the sphere.</summary>
    public float Radius { get; }

    /// <summary>
    /// Creates a sphere.
    /// </summary>
    /// <param name="center">The center of the sphere. Must be finite.</param>
    /// <param name="radius">The radius of the sphere. Must be positive and finite.</param>
    /// <exception cref="ArgumentException"><paramref name="center"/> is not finite.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="radius"/> is not positive and finite.</exception>
    public Sphere(Vector3 center, float radius)
    {
        ShapeArgument.EnsureFinite(center, nameof(center), "Center");
        ShapeArgument.EnsurePositiveFinite(radius, nameof(radius), "Radius");

        Center = center;
        Radius = radius;
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        return Center + Radius * Vector3.Normalize(direction);
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return Center;
    }
}

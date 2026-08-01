using System.Numerics;

namespace OpenGJKSharp;

/// <summary>
/// A convex shape described by its support mapping.
/// Implementations must describe a convex shape, and <see cref="GetSupportPoint"/>
/// must be exact for the shape - GJK/EPA accuracy depends on it.
/// </summary>
public interface ISupportMappable
{
    /// <summary>
    /// Returns the furthest point of the shape in the given direction (world space).
    /// The direction is not required to be normalized and is never zero.
    /// </summary>
    /// <param name="direction">The (non-zero) search direction.</param>
    /// <returns>The furthest point of the shape along <paramref name="direction"/>.</returns>
    Vector3 GetSupportPoint(Vector3 direction);

    /// <summary>
    /// Returns an interior point of the shape (typically the center).
    /// Used as the initial guess for GJK and as the seed for EPA.
    /// </summary>
    /// <returns>An interior point of the shape.</returns>
    Vector3 GetCenter();
}

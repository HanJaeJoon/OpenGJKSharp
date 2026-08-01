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
    /// </summary>
    /// <param name="direction">
    /// The search direction. Guaranteed to be non-zero, so implementations that need a unit
    /// vector may normalize it without a zero check. Its length is unspecified: do not assume
    /// it is normalized, and do not use the length for anything.
    /// </param>
    /// <returns>The furthest point of the shape along <paramref name="direction"/>.</returns>
    /// <remarks>
    /// Ties may be broken arbitrarily - any point that maximizes the dot product with
    /// <paramref name="direction"/> is a valid answer.
    /// </remarks>
    Vector3 GetSupportPoint(Vector3 direction);

    /// <summary>
    /// Returns an interior point of the shape (typically the center).
    /// Used as the initial guess for GJK and as the seed for EPA.
    /// </summary>
    /// <returns>An interior point of the shape.</returns>
    /// <remarks>
    /// The point must lie inside the shape; returning a point outside it makes query results
    /// undefined. Which interior point you return does not affect correctness - it only seeds
    /// the first search direction, so a poor choice costs iterations rather than accuracy
    /// (results may still differ within the convergence tolerance, as they do for any change
    /// in iteration path). When the exact centroid is expensive, any cheap interior point
    /// (for example the midpoint of the bounding box) is a good answer.
    /// </remarks>
    Vector3 GetCenter();
}

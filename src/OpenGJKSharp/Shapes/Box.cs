using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// An oriented box defined by a world-space center, half extents, and orientation.
/// </summary>
public sealed class Box : ISupportMappable
{
    /// <summary>The center of the box.</summary>
    public Vector3 Center { get; }

    /// <summary>The half extents of the box along its local axes.</summary>
    public Vector3 HalfExtents { get; }

    /// <summary>The orientation of the box (normalized).</summary>
    public Quaternion Orientation { get; }

    /// <summary>
    /// Creates an axis-aligned box (identity orientation).
    /// </summary>
    /// <param name="center">The center of the box. Must be finite.</param>
    /// <param name="halfExtents">The half extents of the box. Every component must be positive and finite.</param>
    /// <exception cref="ArgumentException">An argument is not finite or a half extent is not positive.</exception>
    public Box(Vector3 center, Vector3 halfExtents)
        : this(center, halfExtents, Quaternion.Identity)
    {
    }

    /// <summary>
    /// Creates an oriented box.
    /// </summary>
    /// <param name="center">The center of the box. Must be finite.</param>
    /// <param name="halfExtents">The half extents of the box. Every component must be positive and finite.</param>
    /// <param name="orientation">The orientation of the box. Must be finite with non-zero length; normalized internally.</param>
    /// <exception cref="ArgumentException">An argument is not finite, a half extent is not positive,
    /// or <paramref name="orientation"/> has zero length.</exception>
    public Box(Vector3 center, Vector3 halfExtents, Quaternion orientation)
    {
        ShapeArgument.EnsureFinite(center, nameof(center), "Center");

        if (!ShapeArgument.IsFinite(halfExtents.X) || !ShapeArgument.IsFinite(halfExtents.Y) || !ShapeArgument.IsFinite(halfExtents.Z)
            || halfExtents.X <= 0 || halfExtents.Y <= 0 || halfExtents.Z <= 0)
        {
            throw new ArgumentException($"Half extents must be positive and finite: {halfExtents}.", nameof(halfExtents));
        }

        ShapeArgument.EnsureRotation(orientation, nameof(orientation), "Orientation");

        Center = center;
        HalfExtents = halfExtents;
        Orientation = Quaternion.Normalize(orientation);
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        // Rotate the direction into local space, pick the corner per axis, rotate back.
        Vector3 localDirection = Vector3.Transform(direction, Quaternion.Conjugate(Orientation));
        var localSupport = new Vector3(
            localDirection.X >= 0 ? HalfExtents.X : -HalfExtents.X,
            localDirection.Y >= 0 ? HalfExtents.Y : -HalfExtents.Y,
            localDirection.Z >= 0 ? HalfExtents.Z : -HalfExtents.Z);

        return Center + Vector3.Transform(localSupport, Orientation);
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return Center;
    }
}

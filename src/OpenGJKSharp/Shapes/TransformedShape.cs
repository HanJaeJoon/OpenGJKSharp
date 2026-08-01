using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// Places another shape at a pose (rotation and translation), so a shape can be defined once
/// in local space and moved without rebuilding it.
/// </summary>
/// <remarks>
/// <para>
/// The wrapped shape is queried in its own space: the search direction is rotated into local
/// space, and the resulting support point is rotated and translated back into world space.
/// Wrapping therefore costs no accuracy - it is the same transform the shape would apply itself.
/// </para>
/// <para>
/// Unlike the other built-in shapes this one is mutable: <see cref="SetPose"/> updates the pose
/// in place, so a physics loop can reuse the same instance (and, with a
/// <see cref="GjkWorkspace"/>, run queries without allocating at all). Do not change the pose
/// while a query is running.
/// </para>
/// <example>
/// <code>
/// var hull = new ConvexHull(localPoints);   // built once, in local space
/// var body = new TransformedShape(hull);
///
/// foreach (var frame in frames)
/// {
///     body.SetPose(frame.Rotation, frame.Position);
///     double distance = OpenGjk.ComputeMinimumDistance(body, ground, workspace);
/// }
/// </code>
/// </example>
/// </remarks>
public sealed class TransformedShape : ISupportMappable
{
    private Quaternion _inverseRotation;

    /// <summary>
    /// Wraps a shape with an identity pose.
    /// </summary>
    /// <param name="shape">The shape to place. Its own coordinates are treated as local space.</param>
    /// <exception cref="ArgumentNullException"><paramref name="shape"/> is null.</exception>
    public TransformedShape(ISupportMappable shape)
        : this(shape, Quaternion.Identity, Vector3.Zero)
    {
    }

    /// <summary>
    /// Wraps a shape and places it at the given pose.
    /// </summary>
    /// <param name="shape">The shape to place. Its own coordinates are treated as local space.</param>
    /// <param name="rotation">The rotation applied to the shape. Must be finite with non-zero length; normalized internally.</param>
    /// <param name="translation">The translation applied after the rotation. Must be finite.</param>
    /// <exception cref="ArgumentNullException"><paramref name="shape"/> is null.</exception>
    /// <exception cref="ArgumentException">An argument is not finite, or <paramref name="rotation"/> has zero length.</exception>
    public TransformedShape(ISupportMappable shape, Quaternion rotation, Vector3 translation)
    {
        Shape = shape ?? throw new ArgumentNullException(nameof(shape));
        SetPose(rotation, translation);
    }

    /// <summary>The wrapped shape, defined in local space.</summary>
    public ISupportMappable Shape { get; }

    /// <summary>The current rotation (normalized).</summary>
    public Quaternion Rotation { get; private set; }

    /// <summary>The current translation, applied after <see cref="Rotation"/>.</summary>
    public Vector3 Translation { get; private set; }

    /// <summary>
    /// Moves the shape to a new pose in place.
    /// </summary>
    /// <param name="rotation">The rotation applied to the shape. Must be finite with non-zero length; normalized internally.</param>
    /// <param name="translation">The translation applied after the rotation. Must be finite.</param>
    /// <exception cref="ArgumentException">An argument is not finite, or <paramref name="rotation"/> has zero length.</exception>
    public void SetPose(Quaternion rotation, Vector3 translation)
    {
        ShapeArgument.EnsureFinite(translation, nameof(translation), "Translation");
        ShapeArgument.EnsureRotation(rotation, nameof(rotation));

        Rotation = Quaternion.Normalize(rotation);
        _inverseRotation = Quaternion.Conjugate(Rotation);
        Translation = translation;
    }

    /// <inheritdoc />
    public Vector3 GetSupportPoint(Vector3 direction)
    {
        // world support(d) = R * localSupport(R^-1 * d) + t
        Vector3 localDirection = Vector3.Transform(direction, _inverseRotation);
        Vector3 localSupport = Shape.GetSupportPoint(localDirection);

        return Vector3.Transform(localSupport, Rotation) + Translation;
    }

    /// <inheritdoc />
    public Vector3 GetCenter()
    {
        return Vector3.Transform(Shape.GetCenter(), Rotation) + Translation;
    }
}

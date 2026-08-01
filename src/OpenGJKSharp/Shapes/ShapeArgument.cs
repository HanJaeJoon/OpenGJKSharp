using System.Numerics;

namespace OpenGJKSharp.Shapes;

/// <summary>
/// Shared constructor-argument validation for the built-in shapes.
/// </summary>
internal static class ShapeArgument
{
    // float.IsFinite is unavailable on netstandard2.0
    public static bool IsFinite(float value)
    {
        return !float.IsNaN(value) && !float.IsInfinity(value);
    }

    public static void EnsureFinite(Vector3 value, string paramName, string description)
    {
        if (!IsFinite(value.X) || !IsFinite(value.Y) || !IsFinite(value.Z))
        {
            throw new ArgumentException($"{description} contains a non-finite coordinate: {value}.", paramName);
        }
    }

    public static void EnsurePositiveFinite(float value, string paramName, string description)
    {
        if (!IsFinite(value) || value <= 0)
        {
            throw new ArgumentOutOfRangeException(paramName, value, $"{description} must be positive and finite.");
        }
    }
}

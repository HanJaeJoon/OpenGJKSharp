namespace OpenGJKSharp.Tests;

/// <summary>
/// Numeric helpers the tests need on every target framework.
/// </summary>
internal static class Numbers
{
    /// <summary>
    /// Whether a value carries a negative sign, telling -0 apart from 0.
    /// <c>double.IsNegative</c> would do, but it does not exist on .NET Framework, where the
    /// tests run against the netstandard2.0 assembly.
    /// </summary>
    public static bool IsNegative(double value)
    {
        return BitConverter.DoubleToInt64Bits(value) < 0;
    }
}

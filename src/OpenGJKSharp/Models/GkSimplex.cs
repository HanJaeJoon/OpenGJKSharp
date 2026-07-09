namespace OpenGJKSharp.Models;

/// <summary>
/// Data structure for simplex.
/// The simplex is updated at each GJK-iteration.
/// For the first iteration this value is a guess - and this guess not irrelevant.
/// </summary>
public class GkSimplex
{
    /// <summary>
    /// Number of points defining the simplex.
    /// </summary>
    public int NVrtx { get; set; }

    /// <summary>
    /// Coordinates of the points of the simplex.
    /// </summary>
    public double[,] Vrtx { get; } = new double[4, 3];

    /// <summary>
    /// Indices of the points of the simplex.
    /// </summary>
    public int[,] VrtxIdx { get; } = new int[4, 2];

    /// <summary>
    /// Witness points (closest points on each body).
    /// After calling <see cref="OpenGJKSharp.ComputeMinimumDistance(GkPolytope, GkPolytope, GkSimplex)"/>:
    /// Witnesses[0, *] contains the closest point on the first body,
    /// Witnesses[1, *] contains the closest point on the second body.
    /// These are computed using barycentric coordinates from the final simplex vertices.
    /// </summary>
    public double[,] Witnesses { get; } = new double[2, 3];
}

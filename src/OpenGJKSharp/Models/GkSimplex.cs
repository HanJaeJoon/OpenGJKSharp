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
    public int NVrtx;

    /// <summary>
    /// Coordinates of the points of the simplex.
    /// </summary>
    public double[,] Vrtx = new double[4, 3];

    /// <summary>
    /// Indices of the points of the simplex.
    /// </summary>
    public int[,] VrtxIdx = new int[4, 2];

    /// <summary>
    /// Coordinates of the witness points.
    /// </summary>
    public double[,] Witnesses = new double[2, 3];
}
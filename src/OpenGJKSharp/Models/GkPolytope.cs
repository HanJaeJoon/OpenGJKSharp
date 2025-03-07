namespace OpenGJKSharp.Models;

/// <summary>
/// Data structure for convex polytopes.
/// Polytopes are three-dimensional shapes and the GJK algorithm works directly on their convex-hull.
/// However the convex-hull is never computed explicitly,
/// instead each GJK-iteration employs a support function that has a cost linearly dependent on the number of points defining the polytope.
/// </summary>
public class GkPolytope
{
    /// <summary>
    /// Number of points defining the polytope.
    /// </summary>
    public int NumPoints;

    /// <summary>
    /// Furthest point returned by the support function and updated at each GJK-iteration.
    /// For the first iteration this value is a guess - and this guess not irrelevant.
    /// </summary>
    public double[] S = new double[3];

    /// <summary>
    /// Index of the furthest point returned by the support function.
    /// </summary>
    public int SIdx;

    /// <summary>
    /// Coordinates of the points of the polytope.
    /// This is owned by user who manages and garbage-collects the memory for these coordinates.
    /// </summary>
    public required double[][] Coord;
}

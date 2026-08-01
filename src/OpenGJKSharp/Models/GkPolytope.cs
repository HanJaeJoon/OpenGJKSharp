using System.Numerics;

namespace OpenGJKSharp.Models;

/// <summary>
/// Data structure for convex polytopes.
/// Polytopes are three-dimensional shapes and the GJK algorithm works directly on their convex-hull.
/// However the convex-hull is never computed explicitly,
/// instead each GJK-iteration employs a support function that has a cost linearly dependent on the number of points defining the polytope.
/// </summary>
internal sealed class GkPolytope
{
    /// <summary>
    /// Number of points defining the polytope.
    /// </summary>
    public int NumPoints { get; set; }

    /// <summary>
    /// Furthest point returned by the support function and updated at each GJK-iteration.
    /// For the first iteration this value is a guess - and this guess not irrelevant.
    /// </summary>
    public double[] S { get; } = new double[3];

    /// <summary>
    /// Index of the furthest point returned by the support function.
    /// </summary>
    public int SIdx { get; set; }

    /// <summary>
    /// Coordinates of the points of the polytope, laid out as [point][x, y, z].
    /// </summary>
    public double[][] Coord { get; set; } = [];

    /// <summary>
    /// Fills this polytope from a vertex array, growing <see cref="Coord"/> only when
    /// the existing buffer is too small so repeated queries stay allocation-free.
    /// </summary>
    public void Load(Vector3[] points)
    {
        EnsureCapacity(points.Length);

        for (int i = 0; i < points.Length; i++)
        {
            double[] coordinates = Coord[i];
            coordinates[0] = points[i].X;
            coordinates[1] = points[i].Y;
            coordinates[2] = points[i].Z;
        }

        NumPoints = points.Length;
    }

    /// <summary>
    /// Fills this polytope from a 2D vertex array, placing every point on the z = 0 plane.
    /// </summary>
    public void Load(Vector2[] points)
    {
        EnsureCapacity(points.Length);

        for (int i = 0; i < points.Length; i++)
        {
            double[] coordinates = Coord[i];
            coordinates[0] = points[i].X;
            coordinates[1] = points[i].Y;
            coordinates[2] = 0;
        }

        NumPoints = points.Length;
    }

    /// <summary>
    /// Fills this polytope from consecutive (x, y, z) coordinate triples, keeping full
    /// double precision (the <see cref="Vector3"/> overloads quantize to single precision).
    /// </summary>
    public void Load(ReadOnlySpan<double> coordinates)
    {
        int count = coordinates.Length / 3;
        EnsureCapacity(count);

        for (int i = 0; i < count; i++)
        {
            double[] point = Coord[i];
            point[0] = coordinates[i * 3];
            point[1] = coordinates[(i * 3) + 1];
            point[2] = coordinates[(i * 3) + 2];
        }

        NumPoints = count;
    }

    private void EnsureCapacity(int count)
    {
        if (Coord.Length < count)
        {
            var grown = new double[count][];
            Array.Copy(Coord, grown, Coord.Length);
            Coord = grown;
        }

        for (int i = 0; i < count; i++)
        {
            if (Coord[i] is not { Length: >= 3 })
            {
                Coord[i] = new double[3];
            }
        }
    }
}

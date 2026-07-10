using System.Diagnostics;
using OpenGJKSharp.Models;
using System.Numerics;

namespace OpenGJKSharp;

/// <summary>
/// GJK (Gilbert-Johnson-Keerthi) collision detection between convex polytopes,
/// with EPA (Expanding Polytope Algorithm) support for penetration depth.
/// </summary>
public static class OpenGjk
{
    private const double _gkEpsilon = 2.2204460492503131e-16; // DBL_EPSILON
    private const double _precision = 1e-6;

    /// <summary>
    /// Determines whether two convex polytopes collide.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="precision">Distance below which the bodies are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(Vector3[] a, Vector3[] b, double precision = _precision)
    {
        ValidateVertices(a, b);
        if (precision < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(precision), precision, "Precision must not be negative.");
        }

        double distance = ComputeMinimumDistance(ToGkPolytope(a), ToGkPolytope(b));

        return distance <= precision;
    }

    /// <summary>
    /// Determines whether two convex 2D polygons collide (treated as flat shapes on the z = 0 plane).
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="precision">Distance below which the polygons are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(Vector2[] a, Vector2[] b, double precision = _precision)
    {
        ValidateVertices(a, b);
        if (precision < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(precision), precision, "Precision must not be negative.");
        }

        double distance = ComputeMinimumDistance(ToGkPolytope(a), ToGkPolytope(b));

        return distance <= precision;
    }

    /// <summary>
    /// Computes collision information between two convex polytopes.
    /// Returns the minimum distance when the bodies are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first body toward the second.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    public static double ComputeCollisionInformation(Vector3[] a, Vector3[] b, out Vector3 contactNormal)
    {
        ValidateVertices(a, b);

        var bd1 = ToGkPolytope(a);
        var bd2 = ToGkPolytope(b);
        var simplex = new GkSimplex();

        double distance = ComputeMinimumDistance(bd1, bd2, simplex);

        double[] normal = new double[3];
        ComputeCollisionInformation(bd1, bd2, simplex, ref distance, normal);

        contactNormal = new Vector3((float)normal[0], (float)normal[1], (float)normal[2]);
        return distance;
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes.
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    public static double ComputeMinimumDistance(Vector3[] a, Vector3[] b)
    {
        return ComputeMinimumDistance(a, b, out _, out _);
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes,
    /// along with the closest point on each body (witness points).
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="closestA">Closest point on the first body.</param>
    /// <param name="closestB">Closest point on the second body.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    public static double ComputeMinimumDistance(Vector3[] a, Vector3[] b, out Vector3 closestA, out Vector3 closestB)
    {
        ValidateVertices(a, b);

        var simplex = new GkSimplex();
        double distance = ComputeMinimumDistance(ToGkPolytope(a), ToGkPolytope(b), simplex);

        closestA = new Vector3(
            (float)simplex.Witnesses[0, 0],
            (float)simplex.Witnesses[0, 1],
            (float)simplex.Witnesses[0, 2]);
        closestB = new Vector3(
            (float)simplex.Witnesses[1, 0],
            (float)simplex.Witnesses[1, 1],
            (float)simplex.Witnesses[1, 2]);

        return distance;
    }

    /// <summary>
    /// Computes the minimum distance between two convex 2D polygons
    /// (treated as flat shapes on the z = 0 plane).
    /// Returns 0 when the polygons collide.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <returns>The minimum distance between the two polygons.</returns>
    public static double ComputeMinimumDistance(Vector2[] a, Vector2[] b)
    {
        return ComputeMinimumDistance(a, b, out _, out _);
    }

    /// <summary>
    /// Computes the minimum distance between two convex 2D polygons
    /// (treated as flat shapes on the z = 0 plane),
    /// along with the closest point on each polygon (witness points).
    /// Returns 0 when the polygons collide.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="closestA">Closest point on the first polygon.</param>
    /// <param name="closestB">Closest point on the second polygon.</param>
    /// <returns>The minimum distance between the two polygons.</returns>
    public static double ComputeMinimumDistance(Vector2[] a, Vector2[] b, out Vector2 closestA, out Vector2 closestB)
    {
        ValidateVertices(a, b);

        var simplex = new GkSimplex();
        double distance = ComputeMinimumDistance(ToGkPolytope(a), ToGkPolytope(b), simplex);

        closestA = new Vector2(
            (float)simplex.Witnesses[0, 0],
            (float)simplex.Witnesses[0, 1]);
        closestB = new Vector2(
            (float)simplex.Witnesses[1, 0],
            (float)simplex.Witnesses[1, 1]);

        return distance;
    }

    private static void ValidateVertices(Vector3[] a, Vector3[] b)
    {
        ValidateNotNullOrEmpty(a, b);

        for (int i = 0; i < a.Length; i++)
        {
            if (!float.IsFinite(a[i].X) || !float.IsFinite(a[i].Y) || !float.IsFinite(a[i].Z))
            {
                throw new ArgumentException($"Vertex at index {i} contains a non-finite coordinate: {a[i]}.", nameof(a));
            }
        }

        for (int i = 0; i < b.Length; i++)
        {
            if (!float.IsFinite(b[i].X) || !float.IsFinite(b[i].Y) || !float.IsFinite(b[i].Z))
            {
                throw new ArgumentException($"Vertex at index {i} contains a non-finite coordinate: {b[i]}.", nameof(b));
            }
        }
    }

    private static void ValidateVertices(Vector2[] a, Vector2[] b)
    {
        ValidateNotNullOrEmpty(a, b);

        for (int i = 0; i < a.Length; i++)
        {
            if (!float.IsFinite(a[i].X) || !float.IsFinite(a[i].Y))
            {
                throw new ArgumentException($"Vertex at index {i} contains a non-finite coordinate: {a[i]}.", nameof(a));
            }
        }

        for (int i = 0; i < b.Length; i++)
        {
            if (!float.IsFinite(b[i].X) || !float.IsFinite(b[i].Y))
            {
                throw new ArgumentException($"Vertex at index {i} contains a non-finite coordinate: {b[i]}.", nameof(b));
            }
        }
    }

    private static void ValidateNotNullOrEmpty<T>(T[] a, T[] b)
    {
        if (a is null)
        {
            throw new ArgumentNullException(nameof(a));
        }

        if (b is null)
        {
            throw new ArgumentNullException(nameof(b));
        }

        if (a.Length == 0)
        {
            throw new ArgumentException("At least one vertex is required.", nameof(a));
        }

        if (b.Length == 0)
        {
            throw new ArgumentException("At least one vertex is required.", nameof(b));
        }
    }

    private static void ValidatePolytope(GkPolytope polytope, string paramName)
    {
        if (polytope is null)
        {
            throw new ArgumentNullException(paramName);
        }

        if (polytope.NumPoints < 1 || polytope.NumPoints > polytope.Coord.Length)
        {
            throw new ArgumentException(
                $"NumPoints ({polytope.NumPoints}) must be between 1 and Coord.Length ({polytope.Coord.Length}).",
                paramName);
        }

        for (int i = 0; i < polytope.NumPoints; i++)
        {
            if (polytope.Coord[i] is not { Length: >= 3 })
            {
                throw new ArgumentException(
                    $"Coord[{i}] must contain at least 3 elements (x, y, z).",
                    paramName);
            }
        }
    }

    private static GkPolytope ToGkPolytope(Vector3[] points)
    {
        var coordinates = new double[points.Length][];

        for (int i = 0; i < points.Length; i++)
        {
            var point = points[i];
            coordinates[i] = [point.X, point.Y, point.Z];
        }

        return new GkPolytope()
        {
            NumPoints = points.Length,
            Coord = coordinates,
        };
    }

    private static GkPolytope ToGkPolytope(Vector2[] points)
    {
        var coordinates = new double[points.Length][];

        for (int i = 0; i < points.Length; i++)
        {
            var point = points[i];
            coordinates[i] = [point.X, point.Y, 0];
        }

        return new GkPolytope()
        {
            NumPoints = points.Length,
            Coord = coordinates,
        };
    }

    #region OpenGJK ported from https://github.com/MattiaMontanari/openGJK

    /// <summary>
    /// Computes the minimum distance between two bodies given as [3, n] coordinate arrays.
    /// Kept for compatibility with the original openGJK C interface.
    /// </summary>
    /// <param name="nCoordsA">Number of points of the first body.</param>
    /// <param name="inCoordsA">Coordinates of the first body as a [3, n] array (x/y/z rows).</param>
    /// <param name="nCoordsB">Number of points of the second body.</param>
    /// <param name="inCoordsB">Coordinates of the second body as a [3, n] array (x/y/z rows).</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    [Obsolete("Use ComputeMinimumDistance(Vector3[], Vector3[]) or ComputeMinimumDistance(GkPolytope, GkPolytope) instead.")]
    public static double CsFunction(int nCoordsA, double[,] inCoordsA, int nCoordsB, double[,] inCoordsB)
    {
        double[][] pinCoordsA = new double[nCoordsA][];
        for (int j = 0; j < nCoordsA; j++)
        {
            pinCoordsA[j] = [inCoordsA[0, j], inCoordsA[1, j], inCoordsA[2, j]];
        }

        double[][] pinCoordsB = new double[nCoordsB][];
        for (int j = 0; j < nCoordsB; j++)
        {
            pinCoordsB[j] = [inCoordsB[0, j], inCoordsB[1, j], inCoordsB[2, j]];
        }

        var bd1 = new GkPolytope()
        {
            Coord = pinCoordsA,
            NumPoints = nCoordsA,
        };
        var bd2 = new GkPolytope()
        {
            Coord = pinCoordsB,
            NumPoints = nCoordsB
        };

        // Compute minimum distance using GJK algorithm
        double distance = ComputeMinimumDistance(bd1, bd2);

        return distance;
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes using the GJK algorithm.
    /// </summary>
    /// <param name="bd1">The first body.</param>
    /// <param name="bd2">The second body.</param>
    /// <returns>The minimum distance between the two bodies (0 when colliding).</returns>
    public static double ComputeMinimumDistance(GkPolytope bd1, GkPolytope bd2)
    {
        return ComputeMinimumDistance(bd1, bd2, new GkSimplex());
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes using the GJK algorithm.
    /// After the call, <paramref name="s"/> holds the final simplex and the witness points
    /// (closest point on each body) in <see cref="GkSimplex.Witnesses"/>.
    /// </summary>
    /// <param name="bd1">The first body.</param>
    /// <param name="bd2">The second body.</param>
    /// <param name="s">Simplex updated during the GJK iterations.</param>
    /// <returns>The minimum distance between the two bodies (0 when colliding).</returns>
    public static double ComputeMinimumDistance(GkPolytope bd1, GkPolytope bd2, GkSimplex s)
    {
        ValidatePolytope(bd1, nameof(bd1));
        ValidatePolytope(bd2, nameof(bd2));
        if (s is null)
        {
            throw new ArgumentNullException(nameof(s));
        }

        uint k = 0; // Iteration counter
        const int mk = 25; // Maximum number of GJK iterations
        const double eps_rel = _gkEpsilon * 1e4; // Tolerance on relative
        const double eps_tot = _gkEpsilon * 1e2; // Tolerance on absolute distance
        const double eps_rel2 = eps_rel * eps_rel;
        double[] w = new double[3];
        int[] w_idx = new int[2];
        double[] v = new double[3];
        double[] vMinus = new double[3];
        double norm2WMax = 0;

        // Initialize search direction
        v[0] = bd1.Coord[0][0] - bd2.Coord[0][0];
        v[1] = bd1.Coord[0][1] - bd2.Coord[0][1];
        v[2] = bd1.Coord[0][2] - bd2.Coord[0][2];

        // Initialize simplex
        s.NVrtx = 1;
        for (int t = 0; t < 3; ++t)
        {
            s.Vrtx[0, t] = v[t];
        }
        s.VrtxIdx[0, 0] = 0;
        s.VrtxIdx[0, 1] = 0;

        for (int t = 0; t < 3; ++t)
        {
            bd1.S[t] = bd1.Coord[0][t];
        }
        bd1.SIdx = 0;

        for (int t = 0; t < 3; ++t)
        {
            bd2.S[t] = bd2.Coord[0][t];
        }
        bd2.SIdx = 0;

        // Begin GJK iteration
        do
        {
            k++;

            // Update negative search direction
            for (int t = 0; t < 3; ++t)
            {
                vMinus[t] = -v[t];
            }

            // Support function
            Support(bd1, vMinus);
            Support(bd2, v);
            for (int t = 0; t < 3; ++t)
            {
                w[t] = bd1.S[t] - bd2.S[t];
            }
            w_idx[0] = bd1.SIdx;
            w_idx[1] = bd2.SIdx;

            // Test first exit condition
            double exceedTol_rel = Norm2(v) - DotProduct(v, w);
            if (exceedTol_rel <= (eps_rel * Norm2(v)) || exceedTol_rel < eps_tot)
            {
                break;
            }

            if (Norm2(v) < eps_rel2)
            {
                break;
            }

            // Add new vertex to simplex
            int i = s.NVrtx;
            for (int t = 0; t < 3; ++t)
            {
                s.Vrtx[i, t] = w[t];
            }
            s.VrtxIdx[i, 0] = w_idx[0];
            s.VrtxIdx[i, 1] = w_idx[1];
            s.NVrtx++;

            // Invoke distance sub-algorithm
            SubAlgorithm(s, v);

            // Test
            for (int jj = 0; jj < s.NVrtx; jj++)
            {
                double testNorm = Norm2([s.Vrtx[jj, 0], s.Vrtx[jj, 1], s.Vrtx[jj, 2]]);
                if (testNorm > norm2WMax)
                {
                    norm2WMax = testNorm;
                }
            }

            if (Norm2(v) <= (eps_tot * eps_tot * norm2WMax))
            {
                break;
            }

        } while (s.NVrtx != 4 && k != mk);

        if (k == mk)
        {
            Debug.WriteLine("\n * * * * * * * * * * * * MAXIMUM ITERATION NUMBER REACHED!!!  * * * * * * * * * * * * * * \n");
        }

        ComputeWitnesses(bd1, bd2, s);
        return Math.Sqrt(Norm2(v));
    }

    private static void Support(GkPolytope body, double[] v)
    {
        double s, maxS;
        double[] vrt;
        int better = -1;

        maxS = DotProduct(body.S, v);

        for (int i = 0; i < body.NumPoints; ++i)
        {
            vrt = body.Coord[i];
            s = DotProduct(vrt, v);
            if (s > maxS)
            {
                maxS = s;
                better = i;
            }
        }

        if (better != -1)
        {
            body.S[0] = body.Coord[better][0];
            body.S[1] = body.Coord[better][1];
            body.S[2] = body.Coord[better][2];
            body.SIdx = better;
        }
    }

    private static double DotProduct(double[] a, double[] b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    private static double Norm2(double[] a)
    {
        return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    }

    private static void SubAlgorithm(GkSimplex s, double[] v)
    {
        switch (s.NVrtx)
        {
            case 4:
                S3D(s, v);
                break;
            case 3:
                S2D(s, v);
                break;
            case 2:
                S1D(s, v);
                break;
            default:
                Debug.WriteLine("\nERROR:\t invalid simplex\n");
                break;
        }
    }

    private static void S1D(GkSimplex s, double[] v)
    {
        double[] s1p = [s.Vrtx[1, 0], s.Vrtx[1, 1], s.Vrtx[1, 2]];
        double[] s2p = [s.Vrtx[0, 0], s.Vrtx[0, 1], s.Vrtx[0, 2]];

        if (Hff1(s1p, s2p))
        {
            ProjectOnLine(s1p, s2p, v);
            return;
        }
        else
        {
            S1DRegion1(s, v);
            return;
        }
    }

    private static bool Hff1(double[] p, double[] q)
    {
        double tmp = 0;

        for (int i = 0; i < 3; i++)
        {
            tmp += p[i] * p[i] - p[i] * q[i];
        }

        if (tmp > 0)
        {
            return true;
        }
        return false;
    }

    private static void ProjectOnLine(double[] p, double[] q, double[] v)
    {
        double[] pq = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
        double tmp = DotProduct(p, pq) / DotProduct(pq, pq);

        for (int i = 0; i < 3; i++)
        {
            v[i] = p[i] - pq[i] * tmp;
        }
    }

    private static void S1DRegion1(GkSimplex s, double[] v)
    {
        v[0] = s.Vrtx[1, 0];
        v[1] = s.Vrtx[1, 1];
        v[2] = s.Vrtx[1, 2];

        s.NVrtx = 1;

        s.Vrtx[0, 0] = s.Vrtx[1, 0];
        s.Vrtx[0, 1] = s.Vrtx[1, 1];
        s.Vrtx[0, 2] = s.Vrtx[1, 2];

        s.VrtxIdx[0, 0] = s.VrtxIdx[1, 0];
        s.VrtxIdx[0, 1] = s.VrtxIdx[1, 1];
    }

    private static void S2D(GkSimplex s, double[] v)
    {
        double[] s1p = [s.Vrtx[2, 0], s.Vrtx[2, 1], s.Vrtx[2, 2]];
        double[] s2p = [s.Vrtx[1, 0], s.Vrtx[1, 1], s.Vrtx[1, 2]];
        double[] s3p = [s.Vrtx[0, 0], s.Vrtx[0, 1], s.Vrtx[0, 2]];

        bool hff1f_s12 = Hff1(s1p, s2p);
        bool hff1f_s13 = Hff1(s1p, s3p);

        if (hff1f_s12)
        {
            bool hff2f_23 = !Hff2(s1p, s2p, s3p);
            if (hff2f_23)
            {
                if (hff1f_s13)
                {
                    bool hff2f_32 = !Hff2(s1p, s3p, s2p);
                    if (hff2f_32)
                    {
                        ProjectOnPlane(s1p, s2p, s3p, v);
                        return;
                    }
                    else
                    {
                        ProjectOnLine(s1p, s3p, v);
                        S2DRegion13(s);
                        return;
                    }
                }
                else
                {
                    ProjectOnPlane(s1p, s2p, s3p, v);
                    return;
                }
            }
            else
            {
                ProjectOnLine(s1p, s2p, v);
                S2DRegion12(s);
                return;
            }
        }
        else if (hff1f_s13)
        {
            bool hff2f_32 = !Hff2(s1p, s3p, s2p);
            if (hff2f_32)
            {
                ProjectOnPlane(s1p, s2p, s3p, v);
                return;
            }
            else
            {
                ProjectOnLine(s1p, s3p, v);
                S2DRegion13(s);
                return;
            }
        }
        else
        {
            S2DRegion1(s, v);
            return;
        }
    }

    private static bool Hff2(double[] p, double[] q, double[] r)
    {
        double[] nTmp = new double[3];
        double[] n = new double[3];
        double[] pq = new double[3];
        double[] pr = new double[3];

        for (int i = 0; i < 3; i++)
        {
            pq[i] = q[i] - p[i];
        }
        for (int i = 0; i < 3; i++)
        {
            pr[i] = r[i] - p[i];
        }

        CrossProduct(pq, pr, nTmp);
        CrossProduct(pq, nTmp, n);

        return DotProduct(p, n) < 0;
    }

    private static void CrossProduct(double[] a, double[] b, double[] c)
    {
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
    }

    private static void ProjectOnPlane(double[] p, double[] q, double[] r, double[] v)
    {
        double[] n = new double[3];
        double[] pq = new double[3];
        double[] pr = new double[3];

        for (int i = 0; i < 3; i++)
        {
            pq[i] = p[i] - q[i];
        }
        for (int i = 0; i < 3; i++)
        {
            pr[i] = p[i] - r[i];
        }

        CrossProduct(pq, pr, n);
        double tmp = DotProduct(n, p) / DotProduct(n, n);

        for (int i = 0; i < 3; i++)
        {
            v[i] = n[i] * tmp;
        }
    }

    private static void S2DRegion13(GkSimplex s)
    {
        s.NVrtx = 2;

        s.Vrtx[1, 0] = s.Vrtx[2, 0];
        s.Vrtx[1, 1] = s.Vrtx[2, 1];
        s.Vrtx[1, 2] = s.Vrtx[2, 2];

        s.VrtxIdx[1, 0] = s.VrtxIdx[2, 0];
        s.VrtxIdx[1, 1] = s.VrtxIdx[2, 1];
    }

    private static void S2DRegion12(GkSimplex s)
    {
        s.NVrtx = 2;

        s.Vrtx[0, 0] = s.Vrtx[2, 0];
        s.Vrtx[0, 1] = s.Vrtx[2, 1];
        s.Vrtx[0, 2] = s.Vrtx[2, 2];

        s.VrtxIdx[0, 0] = s.VrtxIdx[2, 0];
        s.VrtxIdx[0, 1] = s.VrtxIdx[2, 1];
    }

    private static void S2DRegion1(GkSimplex s, double[] v)
    {
        v[0] = s.Vrtx[2, 0];
        v[1] = s.Vrtx[2, 1];
        v[2] = s.Vrtx[2, 2];

        s.NVrtx = 1;

        s.Vrtx[0, 0] = s.Vrtx[2, 0];
        s.Vrtx[0, 1] = s.Vrtx[2, 1];
        s.Vrtx[0, 2] = s.Vrtx[2, 2];

        s.VrtxIdx[0, 0] = s.VrtxIdx[2, 0];
        s.VrtxIdx[0, 1] = s.VrtxIdx[2, 1];
    }

    private static void S3D(GkSimplex s, double[] v)
    {
        double[] s1 = new double[3];
        double[] s2 = new double[3];
        double[] s3 = new double[3];
        double[] s4 = new double[3];
        double[] s1s2 = new double[3];
        double[] s1s3 = new double[3];
        double[] s1s4 = new double[3];
        double[] si = new double[3];
        double[] sj = new double[3];
        double[] sk = new double[3];
        int[] s1Idx = new int[2];
        int[] s2Idx = new int[2];
        int[] s3Idx = new int[2];
        int[] siIdx = new int[2];
        int[] sjIdx = new int[2];
        int[] skIdx = new int[2];
        bool testLineThree, testLineFour;
        int testPlaneTwo, testPlaneThree, testPlaneFour, dotTotal;
        int i, j, k;

        GetVrtxIdx(s, s1, s1Idx, 3);
        GetVrtxIdx(s, s2, s2Idx, 2);
        GetVrtxIdx(s, s3, s3Idx, 1);
        GetVrtx(s, s4, 0);
        CalculateEdgeVector(s, s1s2, s2);
        CalculateEdgeVector(s, s1s3, s3);
        CalculateEdgeVector(s, s1s4, s4);

        bool[] hff1Tests = new bool[3];
        hff1Tests[2] = Hff1(s1, s2);
        hff1Tests[1] = Hff1(s1, s3);
        hff1Tests[0] = Hff1(s1, s4);
        testLineThree = Hff1(s1, s3);
        testLineFour = Hff1(s1, s4);

        dotTotal = (Hff1(s1, s2) ? 1 : 0) + (testLineThree ? 1 : 0) + (testLineFour ? 1 : 0);
        if (dotTotal == 0)
        {
            S3DRegion1(s, v, s1, s1Idx);
            return;
        }

        double det134 = Determinant(s1s3, s1s4, s1s2);
        int sss = det134 <= 0 ? 1 : 0;

        testPlaneTwo = (Hff3(s1, s3, s4) ? 1 : 0) - sss;
        testPlaneTwo *= testPlaneTwo;
        testPlaneThree = (Hff3(s1, s4, s2) ? 1 : 0) - sss;
        testPlaneThree *= testPlaneThree;
        testPlaneFour = (Hff3(s1, s2, s3) ? 1 : 0) - sss;
        testPlaneFour *= testPlaneFour;

        switch (testPlaneTwo + testPlaneThree + testPlaneFour)
        {
            case 3:
                S3DRegion1234(s, v);
                break;

            case 2:
                s.NVrtx = 3;
                if (testPlaneTwo == 0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        s.Vrtx[2, i] = s.Vrtx[3, i];
                    }
                    for (i = 0; i < 2; i++)
                    {
                        s.VrtxIdx[2, i] = s.VrtxIdx[3, i];
                    }
                }
                else if (testPlaneThree == 0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        s.Vrtx[1, i] = s2[i];
                        s.Vrtx[2, i] = s.Vrtx[3, i];
                    }
                    for (i = 0; i < 2; i++)
                    {
                        s.VrtxIdx[1, i] = s2Idx[i];
                        s.VrtxIdx[2, i] = s.VrtxIdx[3, i];
                    }
                }
                else if (testPlaneFour == 0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        s.Vrtx[0, i] = s3[i];
                        s.Vrtx[1, i] = s2[i];
                        s.Vrtx[2, i] = s.Vrtx[3, i];
                    }
                    for (i = 0; i < 2; i++)
                    {
                        s.VrtxIdx[0, i] = s3Idx[i];
                        s.VrtxIdx[1, i] = s2Idx[i];
                        s.VrtxIdx[2, i] = s.VrtxIdx[3, i];
                    }
                }
                S2D(s, v);
                break;

            case 1:
                s.NVrtx = 3;
                if (testPlaneTwo != 0)
                {
                    k = 2;
                    i = 1;
                    j = 0;
                }
                else if (testPlaneThree != 0)
                {
                    k = 1;
                    i = 0;
                    j = 2;
                }
                else
                {
                    k = 0;
                    i = 2;
                    j = 1;
                }

                GetVrtxIdx(s, si, siIdx, i);
                GetVrtxIdx(s, sj, sjIdx, j);
                GetVrtxIdx(s, sk, skIdx, k);

                if (dotTotal == 1)
                {
                    if (hff1Tests[k])
                    {
                        if (!Hff2(s1, sk, si))
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, si, sk, v);
                        }
                        else if (!Hff2(s1, sk, sj))
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sj, sk, v);
                        }
                        else
                        {
                            Select1(s, sk, skIdx);
                            ProjectOnLine(s1, sk, v);
                        }
                    }
                    else if (hff1Tests[i])
                    {
                        if (!Hff2(s1, si, sk))
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, si, sk, v);
                        }
                        else
                        {
                            Select1(s, si, siIdx);
                            ProjectOnLine(s1, si, v);
                        }
                    }
                    else
                    {
                        if (!Hff2(s1, sj, sk))
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sj, sk, v);
                        }
                        else
                        {
                            Select1(s, sj, sjIdx);
                            ProjectOnLine(s1, sj, v);
                        }
                    }
                }
                else if (dotTotal == 2)
                {
                    if (hff1Tests[i])
                    {
                        if (!Hff2(s1, sk, si))
                        {
                            if (!Hff2(s1, si, sk))
                            {
                                Select1(s, si, siIdx, sk, skIdx);
                                ProjectOnPlane(s1, si, sk, v);
                            }
                            else
                            {
                                Select1(s, sk, skIdx);
                                ProjectOnLine(s1, sk, v);
                            }
                        }
                        else
                        {
                            if (!Hff2(s1, sk, sj))
                            {
                                Select1(s, sj, sjIdx, sk, skIdx);
                                ProjectOnPlane(s1, sj, sk, v);
                            }
                            else
                            {
                                Select1(s, sk, skIdx);
                                ProjectOnLine(s1, sk, v);
                            }
                        }
                    }
                    else if (hff1Tests[j])
                    {
                        if (!Hff2(s1, sk, sj))
                        {
                            if (!Hff2(s1, sj, sk))
                            {
                                Select1(s, sj, sjIdx, sk, skIdx);
                                ProjectOnPlane(s1, sj, sk, v);
                            }
                            else
                            {
                                Select1(s, sj, sjIdx);
                                ProjectOnLine(s1, sj, v);
                            }
                        }
                        else
                        {
                            if (!Hff2(s1, sk, si))
                            {
                                Select1(s, si, siIdx, sk, skIdx);
                                ProjectOnPlane(s1, si, sk, v);
                            }
                            else
                            {
                                Select1(s, sk, skIdx);
                                ProjectOnLine(s1, sk, v);
                            }
                        }
                    }
                }
                else if (dotTotal == 3)
                {
                    bool hff2Ik = Hff2(s1, si, sk);
                    bool hff2Jk = Hff2(s1, sj, sk);
                    bool hff2Ki = Hff2(s1, sk, si);
                    bool hff2Kj = Hff2(s1, sk, sj);

                    if (!hff2Ki && !hff2Kj)
                    {
                        Debug.WriteLine("\n\n UNEXPECTED VALUES!!! \n\n");
                    }
                    if (hff2Ki && hff2Kj)
                    {
                        Select1(s, sk, skIdx);
                        ProjectOnLine(s1, sk, v);
                    }
                    else if (hff2Ki)
                    {
                        if (hff2Jk)
                        {
                            Select1(s, sj, sjIdx);
                            ProjectOnLine(s1, sj, v);
                        }
                        else
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sk, sj, v);
                        }
                    }
                    else
                    {
                        if (hff2Ik)
                        {
                            Select1(s, si, siIdx);
                            ProjectOnLine(s1, si, v);
                        }
                        else
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, sk, si, v);
                        }
                    }
                }
                break;

            case 0:
                if (dotTotal == 1)
                {
                    if (testLineThree)
                    {
                        k = 2;
                        i = 1;
                        j = 0;
                    }
                    else if (testLineFour)
                    {
                        k = 1;
                        i = 0;
                        j = 2;
                    }
                    else
                    {
                        k = 0;
                        i = 2;
                        j = 1;
                    }
                    GetVrtxIdx(s, si, siIdx, i);
                    GetVrtxIdx(s, sj, sjIdx, j);
                    GetVrtxIdx(s, sk, skIdx, k);

                    if (!Hff2(s1, si, sj))
                    {
                        Select1(s, si, siIdx, sj, sjIdx);
                        ProjectOnPlane(s1, si, sj, v);
                    }
                    else if (!Hff2(s1, si, sk))
                    {
                        Select1(s, si, siIdx, sk, skIdx);
                        ProjectOnPlane(s1, si, sk, v);
                    }
                    else
                    {
                        Select1(s, si, siIdx);
                        ProjectOnLine(s1, si, v);
                    }
                }
                else if (dotTotal == 2)
                {
                    s.NVrtx = 3;
                    if (!testLineThree)
                    {
                        k = 2;
                        i = 1;
                        j = 0;
                    }
                    else if (!testLineFour)
                    {
                        k = 1;
                        i = 0;
                        j = 2;
                    }
                    else
                    {
                        k = 0;
                        i = 2;
                        j = 1;
                    }
                    GetVrtxIdx(s, si, siIdx, i);
                    GetVrtxIdx(s, sj, sjIdx, j);
                    GetVrtxIdx(s, sk, skIdx, k);

                    if (!Hff2(s1, sj, sk))
                    {
                        if (!Hff2(s1, sk, sj))
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sj, sk, v);
                        }
                        else if (!Hff2(s1, sk, si))
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, sk, si, v);
                        }
                        else
                        {
                            Select1(s, sk, skIdx);
                            ProjectOnLine(s1, sk, v);
                        }
                    }
                    else if (!Hff2(s1, sj, si))
                    {
                        Select1(s, si, siIdx, sj, sjIdx);
                        ProjectOnPlane(s1, si, sj, v);
                    }
                    else
                    {
                        Select1(s, sj, sjIdx);
                        ProjectOnLine(s1, sj, v);
                    }
                }
                break;

            default:
                Debug.WriteLine("\nERROR:\tunhandled");
                break;
        }
    }

    private static void GetVrtxIdx(GkSimplex s, double[] point, int[] index, int location)
    {
        point[0] = s.Vrtx[location, 0];
        point[1] = s.Vrtx[location, 1];
        point[2] = s.Vrtx[location, 2];

        index[0] = s.VrtxIdx[location, 0];
        index[1] = s.VrtxIdx[location, 1];
    }

    private static void GetVrtx(GkSimplex s, double[] point, int location)
    {
        point[0] = s.Vrtx[location, 0];
        point[1] = s.Vrtx[location, 1];
        point[2] = s.Vrtx[location, 2];
    }

    private static void CalculateEdgeVector(GkSimplex s, double[] p1p2, double[] p2)
    {
        p1p2[0] = p2[0] - s.Vrtx[3, 0];
        p1p2[1] = p2[1] - s.Vrtx[3, 1];
        p1p2[2] = p2[2] - s.Vrtx[3, 2];
    }

    private static void S3DRegion1(GkSimplex s, double[] v, double[] s1, int[] s1Idx)
    {
        v[0] = s1[0];
        v[1] = s1[1];
        v[2] = s1[2];

        s.NVrtx = 1;

        s.Vrtx[0, 0] = s1[0];
        s.Vrtx[0, 1] = s1[1];
        s.Vrtx[0, 2] = s1[2];

        s.VrtxIdx[0, 0] = s1Idx[0];
        s.VrtxIdx[0, 1] = s1Idx[1];
    }

    private static double Determinant(double[] p, double[] q, double[] r)
    {
        return p[0] * (q[1] * r[2] - r[1] * q[2])
            - p[1] * (q[0] * r[2] - r[0] * q[2])
            + p[2] * (q[0] * r[1] - r[0] * q[1]);
    }

    private static bool Hff3(double[] p, double[] q, double[] r)
    {
        double[] n = new double[3];
        CrossProduct(q, r, n);
        return DotProduct(p, n) <= 0; // discard s if true
    }

    private static void S3DRegion1234(GkSimplex s, double[] v)
    {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
        s.NVrtx = 4;
    }

    private static void Select1(GkSimplex s, double[] si, int[] siIdx, double[] some, int[] someIdx)
    {
        s.NVrtx = 3;

        for (int t = 0; t < 3; t++)
            s.Vrtx[2, t] = s.Vrtx[3, t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[2, t] = s.VrtxIdx[3, t];
        for (int t = 0; t < 3; t++)
            s.Vrtx[1, t] = si[t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[1, t] = siIdx[t];
        for (int t = 0; t < 3; t++)
            s.Vrtx[0, t] = some[t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[0, t] = someIdx[t];
    }

    private static void Select1(GkSimplex s, double[] some, int[] someIdx)
    {
        s.NVrtx = 2;

        for (int t = 0; t < 3; t++)
            s.Vrtx[1, t] = s.Vrtx[3, t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[1, t] = s.VrtxIdx[3, t];
        for (int t = 0; t < 3; t++)
            s.Vrtx[0, t] = some[t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[0, t] = someIdx[t];
    }

    private static void ComputeWitnesses(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        switch (smp.NVrtx)
        {
            case 4:
                W3D(bd1, bd2, smp);
                break;
            case 3:
                W2D(bd1, bd2, smp);
                break;
            case 2:
                W1D(bd1, bd2, smp);
                break;
            case 1:
                W0D(bd1, bd2, smp);
                break;
            default:
                Debug.WriteLine("\nERROR:\t invalid simplex\n");
                break;
        }
    }

    private static void W0D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        var w01 = bd2.Coord[smp.VrtxIdx[0, 1]];

        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t];
            smp.Witnesses[1, t] = w01[t];
        }
    }

    private static void W1D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var pq = new double[3];
        var po = new double[3];

        double[] p = [smp.Vrtx[0, 0], smp.Vrtx[0, 1], smp.Vrtx[0, 2]];
        double[] q = [smp.Vrtx[1, 0], smp.Vrtx[1, 1], smp.Vrtx[1, 2]];

        for (int t = 0; t < 3; t++)
        {
            pq[t] = q[t] - p[t];
            po[t] = -p[t];
        }

        // Compute barycentric coordinates via matrix inversion
        // (in the linear case the matrix is 1x1 thus simplified)
        double det = DotProduct(pq, pq);
        if (det == 0.0)
        {
            // Degenerate case
            W0D(bd1, bd2, smp);
            return;
        }

        double a1 = DotProduct(pq, po) / det;
        double a0 = 1.0 - a1;

        // Compute witness points
        double[] w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        double[] w01 = bd2.Coord[smp.VrtxIdx[0, 1]];
        double[] w10 = bd1.Coord[smp.VrtxIdx[1, 0]];
        double[] w11 = bd2.Coord[smp.VrtxIdx[1, 1]];
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1;
        }
    }

    private static void W2D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var pq = new double[3];
        var pr = new double[3];
        var po = new double[3];

        double[] p = [smp.Vrtx[0, 0], smp.Vrtx[0, 1], smp.Vrtx[0, 2]];
        double[] q = [smp.Vrtx[1, 0], smp.Vrtx[1, 1], smp.Vrtx[1, 2]];
        double[] r = [smp.Vrtx[2, 0], smp.Vrtx[2, 1], smp.Vrtx[2, 2]];

        for (int t = 0; t < 3; t++)
        {
            pq[t] = q[t] - p[t];
            pr[t] = r[t] - p[t];
            po[t] = -p[t];
        }

        double t00 = DotProduct(pq, pq);
        double t01 = DotProduct(pq, pr);
        double t11 = DotProduct(pr, pr);
        double det = t00 * t11 - t01 * t01;
        if (det == 0.0)
        {
            // Degenerate case
            W1D(bd1, bd2, smp);
            return;
        }

        double b0 = DotProduct(pq, po);
        double b1 = DotProduct(pr, po);
        double i00 = t11 / det;
        double i01 = -t01 / det;
        double i11 = t00 / det;
        double a1 = i00 * b0 + i01 * b1;
        double a2 = i01 * b0 + i11 * b1;
        double a0 = 1.0 - a1 - a2;

        // check if the origin is very close to one of the edges of the
        // simplex. In this case, a 1D projection will be more accurate.
        if (a0 < _gkEpsilon)
        {
            smp.NVrtx = 2;
            smp.Vrtx[0, 0] = smp.Vrtx[2, 0];
            smp.Vrtx[0, 1] = smp.Vrtx[2, 1];
            smp.Vrtx[0, 2] = smp.Vrtx[2, 2];
            smp.VrtxIdx[0, 0] = smp.VrtxIdx[2, 0];
            smp.VrtxIdx[0, 1] = smp.VrtxIdx[2, 1];
            W1D(bd1, bd2, smp);
            return;
        }
        else if (a1 < _gkEpsilon)
        {
            smp.NVrtx = 2;
            smp.Vrtx[1, 0] = smp.Vrtx[2, 0];
            smp.Vrtx[1, 1] = smp.Vrtx[2, 1];
            smp.Vrtx[1, 2] = smp.Vrtx[2, 2];
            smp.VrtxIdx[1, 0] = smp.VrtxIdx[2, 0];
            smp.VrtxIdx[1, 1] = smp.VrtxIdx[2, 1];
            W1D(bd1, bd2, smp);
            return;
        }
        else if (a2 < _gkEpsilon)
        {
            smp.NVrtx = 2;
            W1D(bd1, bd2, smp);
            return;
        }

        // Compute witness points
        // This is done by blending the source points using
        // the barycentric coordinates
        double[] w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        double[] w01 = bd2.Coord[smp.VrtxIdx[0, 1]];
        double[] w10 = bd1.Coord[smp.VrtxIdx[1, 0]];
        double[] w11 = bd2.Coord[smp.VrtxIdx[1, 1]];
        double[] w20 = bd1.Coord[smp.VrtxIdx[2, 0]];
        double[] w21 = bd2.Coord[smp.VrtxIdx[2, 1]];
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1 + w20[t] * a2;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1 + w21[t] * a2;
        }
    }

    private static void W3D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var pq = new double[3];
        var pr = new double[3];
        var ps = new double[3];
        var po = new double[3];

        double[] p = [smp.Vrtx[0, 0], smp.Vrtx[0, 1], smp.Vrtx[0, 2]];
        double[] q = [smp.Vrtx[1, 0], smp.Vrtx[1, 1], smp.Vrtx[1, 2]];
        double[] r = [smp.Vrtx[2, 0], smp.Vrtx[2, 1], smp.Vrtx[2, 2]];
        double[] s = [smp.Vrtx[3, 0], smp.Vrtx[3, 1], smp.Vrtx[3, 2]];

        for (int t = 0; t < 3; t++)
        {
            pq[t] = q[t] - p[t];
            pr[t] = r[t] - p[t];
            ps[t] = s[t] - p[t];
            po[t] = -p[t];
        }

        double t00 = DotProduct(pq, pq);
        double t01 = DotProduct(pq, pr);
        double t02 = DotProduct(pq, ps);
        double t11 = DotProduct(pr, pr);
        double t12 = DotProduct(pr, ps);
        double t22 = DotProduct(ps, ps);
        double det00 = t11 * t22 - t12 * t12;
        double det01 = t01 * t22 - t02 * t12;
        double det02 = t01 * t12 - t02 * t11;
        double det = t00 * det00 - t01 * det01 + t02 * det02;
        if (det == 0.0)
        {
            // Degenerate case
            W2D(bd1, bd2, smp);
            return;
        }

        double b0 = DotProduct(pq, po);
        double b1 = DotProduct(pr, po);
        double b2 = DotProduct(ps, po);

        // inverse matrix
        // (the matrix is symmetric, so we can use the cofactor matrix)
        double det11 = t00 * t22 - t02 * t02;
        double det12 = t00 * t12 - t01 * t02;
        double det22 = t00 * t11 - t01 * t01;
        double i00 = det00 / det;
        double i01 = -det01 / det;
        double i02 = det02 / det;
        double i11 = det11 / det;
        double i12 = -det12 / det;
        double i22 = det22 / det;

        double a1 = i00 * b0 + i01 * b1 + i02 * b2;
        double a2 = i01 * b0 + i11 * b1 + i12 * b2;
        double a3 = i02 * b0 + i12 * b1 + i22 * b2;
        double a0 = 1.0 - a1 - a2 - a3;

        // check if the origin is very close to one of the faces of the
        // simplex. In this case, a 2D projection will be more accurate.
        if (a0 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            smp.Vrtx[0, 0] = smp.Vrtx[3, 0];
            smp.Vrtx[0, 1] = smp.Vrtx[3, 1];
            smp.Vrtx[0, 2] = smp.Vrtx[3, 2];
            smp.VrtxIdx[0, 0] = smp.VrtxIdx[3, 0];
            smp.VrtxIdx[0, 1] = smp.VrtxIdx[3, 1];
            W2D(bd1, bd2, smp);
            return;
        }
        else if (a1 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            smp.Vrtx[1, 0] = smp.Vrtx[3, 0];
            smp.Vrtx[1, 1] = smp.Vrtx[3, 1];
            smp.Vrtx[1, 2] = smp.Vrtx[3, 2];
            smp.VrtxIdx[1, 0] = smp.VrtxIdx[3, 0];
            smp.VrtxIdx[1, 1] = smp.VrtxIdx[3, 1];
            W2D(bd1, bd2, smp);
            return;
        }
        else if (a2 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            smp.Vrtx[2, 0] = smp.Vrtx[3, 0];
            smp.Vrtx[2, 1] = smp.Vrtx[3, 1];
            smp.Vrtx[2, 2] = smp.Vrtx[3, 2];
            smp.VrtxIdx[2, 0] = smp.VrtxIdx[3, 0];
            smp.VrtxIdx[2, 1] = smp.VrtxIdx[3, 1];
            W2D(bd1, bd2, smp);
            return;
        }
        else if (a3 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            W2D(bd1, bd2, smp);
            return;
        }

        // Compute witness points
        // This is done by blending the original points using
        // the barycentric coordinates
        double[] w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        double[] w01 = bd2.Coord[smp.VrtxIdx[0, 1]];
        double[] w10 = bd1.Coord[smp.VrtxIdx[1, 0]];
        double[] w11 = bd2.Coord[smp.VrtxIdx[1, 1]];
        double[] w20 = bd1.Coord[smp.VrtxIdx[2, 0]];
        double[] w21 = bd2.Coord[smp.VrtxIdx[2, 1]];
        double[] w30 = bd1.Coord[smp.VrtxIdx[3, 0]];
        double[] w31 = bd2.Coord[smp.VrtxIdx[3, 1]];
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1 + w20[t] * a2 + w30[t] * a3;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1 + w21[t] * a2 + w31[t] * a3;
        }
    }

    //*******************************************************************************************
    // EPA Algorithm
    //*******************************************************************************************

    // Compute face normal and distance of face from origin.
    // Winding is already fixed at face creation time, so the cross product
    // direction is trusted directly - no centroid-based orientation check needed.
    private static void ComputeFaceNormalDistance(EpaPolytope poly, int faceIdx, double[] centroid)
    {
        EpaFace face = poly.Faces[faceIdx];

        double[] v0 = poly.Vertices[face.V[0]];
        double[] v1 = poly.Vertices[face.V[1]];
        double[] v2 = poly.Vertices[face.V[2]];

        double[] e0 = new double[3];
        double[] e1 = new double[3];
        for (int i = 0; i < 3; i++)
        {
            e0[i] = v1[i] - v0[i];
            e1[i] = v2[i] - v0[i];
        }

        CrossProduct(e0, e1, face.Normal);
        double normSq = Norm2(face.Normal);

        if (normSq > _gkEpsilon * _gkEpsilon)
        {
            double norm = Math.Sqrt(normSq);
            for (int i = 0; i < 3; i++)
            {
                face.Normal[i] /= norm;
            }

            face.Distance = DotProduct(face.Normal, v0);

            // Orient the normal away from the polytope interior. Deciding by the sign
            // of Distance alone mis-orients faces whose plane passes through the origin
            // (Distance ~ +-epsilon): a tiny negative value flipped the normal inward,
            // the support point in that direction duplicated an existing vertex, and
            // EPA stalled reporting zero penetration depth.
            double[] toCentroid =
            [
                centroid[0] - v0[0],
                centroid[1] - v0[1],
                centroid[2] - v0[2],
            ];
            double side = DotProduct(face.Normal, toCentroid);
            bool flip = Math.Abs(side) > _gkEpsilon
                ? side > 0
                : face.Distance < 0;
            if (flip)
            {
                for (int i = 0; i < 3; i++)
                {
                    face.Normal[i] = -face.Normal[i];
                }
                face.Distance = -face.Distance;
            }

            // The origin sits inside the polytope, so an outward-oriented face can only
            // report a negative distance through numerical noise - clamp it so the face
            // stays eligible in the closest-face search.
            if (face.Distance < 0 && face.Distance > -_precision)
            {
                face.Distance = 0;
            }
        }
        else
        {
            face.Valid = false;
            face.Distance = 1e10;
        }
    }

    // Check if a face is visible from a point (point is on positive side of face).
    // Needed to determine which faces to restructure when vertex is added.
    private static bool IsFaceVisible(EpaPolytope poly, int faceIdx, double[] point)
    {
        EpaFace face = poly.Faces[faceIdx];
        if (!face.Valid)
        {
            return false;
        }

        double[] v0 = poly.Vertices[face.V[0]];
        double[] diff = new double[3];
        for (int i = 0; i < 3; i++)
        {
            diff[i] = point[i] - v0[i];
        }

        return DotProduct(face.Normal, diff) > _gkEpsilon;
    }

    // Initialize EPA polytope from GJK simplex (should be a tetrahedron)
    private static void InitEpaPolytope(EpaPolytope poly, GkSimplex simplex, double[] centroid)
    {
        // Copy vertices from simplex
        poly.NumVertices = 4;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                poly.Vertices[i][j] = simplex.Vrtx[i, j];
            }
            poly.VertexIndices[i][0] = simplex.VrtxIdx[i, 0];
            poly.VertexIndices[i][1] = simplex.VrtxIdx[i, 1];
        }

        // Compute centroid of the tetrahedron
        centroid[0] = centroid[1] = centroid[2] = 0;

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                centroid[j] += poly.Vertices[i][j] * 0.25;
            }
        }

        // Create 4 faces of tetrahedron
        // set up the faces and then fix winding based on normal direction
        // Face 0: vertices 0, 1, 2
        poly.Faces[0].V[0] = 0;
        poly.Faces[0].V[1] = 1;
        poly.Faces[0].V[2] = 2;
        poly.Faces[0].Valid = true;

        // Face 1: vertices 0, 3, 1
        poly.Faces[1].V[0] = 0;
        poly.Faces[1].V[1] = 3;
        poly.Faces[1].V[2] = 1;
        poly.Faces[1].Valid = true;

        // Face 2: vertices 0, 2, 3
        poly.Faces[2].V[0] = 0;
        poly.Faces[2].V[1] = 2;
        poly.Faces[2].V[2] = 3;
        poly.Faces[2].Valid = true;

        // Face 3: vertices 1, 3, 2
        poly.Faces[3].V[0] = 1;
        poly.Faces[3].V[1] = 3;
        poly.Faces[3].V[2] = 2;
        poly.Faces[3].Valid = true;

        // Copy vertex indices for witness point computation
        for (int f = 0; f < 4; f++)
        {
            for (int v = 0; v < 3; v++)
            {
                int vi = poly.Faces[f].V[v];
                poly.Faces[f].VIdx[v][0] = simplex.VrtxIdx[vi, 0];
                poly.Faces[f].VIdx[v][1] = simplex.VrtxIdx[vi, 1];
            }
        }

        // Compute normals and fix winding
        for (int f = 0; f < 4; f++)
        {
            double[] v0 = poly.Vertices[poly.Faces[f].V[0]];
            double[] v1 = poly.Vertices[poly.Faces[f].V[1]];
            double[] v2 = poly.Vertices[poly.Faces[f].V[2]];

            double[] e0 = new double[3];
            double[] e1 = new double[3];
            double[] normal = new double[3];

            for (int i = 0; i < 3; i++)
            {
                e0[i] = v1[i] - v0[i];
                e1[i] = v2[i] - v0[i];
            }
            CrossProduct(e0, e1, normal);

            // If normal points toward centroid need to flip the winding
            double[] toCentroid = new double[3];
            for (int i = 0; i < 3; i++)
            {
                toCentroid[i] = centroid[i] - v0[i];
            }
            if (DotProduct(normal, toCentroid) > 0)
            {
                (poly.Faces[f].V[1], poly.Faces[f].V[2]) = (poly.Faces[f].V[2], poly.Faces[f].V[1]);
                (poly.Faces[f].VIdx[1], poly.Faces[f].VIdx[2]) = (poly.Faces[f].VIdx[2], poly.Faces[f].VIdx[1]);
            }
        }

        poly.MaxFaceIndex = 4;
    }

    // barycentric coordinate compute closest point on triangle to origin
    private static void ComputeBarycentricOrigin(double[] v0, double[] v1, double[] v2,
        out double a0, out double a1, out double a2)
    {
        // Compute vectors
        double[] e0 = new double[3];
        double[] e1 = new double[3];

        for (int i = 0; i < 3; i++)
        {
            e0[i] = v1[i] - v0[i];
            e1[i] = v2[i] - v0[i];
        }

        // Compute dot products for barycentric coords
        double d00 = DotProduct(e0, e0);
        double d01 = DotProduct(e0, e1);
        double d11 = DotProduct(e1, e1);
        double d20 = -DotProduct(v0, e0);
        double d21 = -DotProduct(v0, e1);

        double denom = d00 * d11 - d01 * d01;

        if (Math.Abs(denom) < _gkEpsilon)
        {
            // Degenerate
            a0 = a1 = a2 = 1.0 / 3.0;
            return;
        }

        double invDenom = 1.0 / denom;
        double u = (d11 * d20 - d01 * d21) * invDenom;
        double v = (d00 * d21 - d01 * d20) * invDenom;
        double w = 1.0 - u - v;

        // Clamp to triangle
        if (w < 0)
        {
            // Origin projects outside edge v1-v2
            // Project onto edge v1-v2
            double[] e12 = new double[3];
            for (int i = 0; i < 3; i++)
            {
                e12[i] = v2[i] - v1[i];
            }
            double t = -DotProduct(v1, e12) / DotProduct(e12, e12);
            t = Math.Clamp(t, 0.0, 1.0);
            a0 = 0;
            a1 = 1.0 - t;
            a2 = t;
        }
        else if (u < 0)
        {
            // Origin projects outside edge v0-v2
            double t = -DotProduct(v0, e1) / DotProduct(e1, e1);
            t = Math.Clamp(t, 0.0, 1.0);
            a0 = 1.0 - t;
            a1 = 0;
            a2 = t;
        }
        else if (v < 0)
        {
            // Origin projects outside edge v0-v1
            double t = -DotProduct(v0, e0) / DotProduct(e0, e0);
            t = Math.Clamp(t, 0.0, 1.0);
            a0 = 1.0 - t;
            a1 = t;
            a2 = 0;
        }
        else
        {
            // Inside triangle
            a0 = w;
            a1 = u;
            a2 = v;
        }
    }

    // Support function for EPA, basically the GJK one but only cares about the Minkowski difference point
    private static void SupportEpa(GkPolytope body1, GkPolytope body2,
        double[] direction, double[] result, int[] resultIdx)
    {
        double localMax1 = -1e10;
        double localMax2 = -1e10;
        int localBest1 = -1;
        int localBest2 = -1;

        // Search body1
        for (int i = 0; i < body1.NumPoints; i++)
        {
            double s = DotProduct(body1.Coord[i], direction);
            if (s > localMax1)
            {
                localMax1 = s;
                localBest1 = i;
            }
        }

        // Search body2 in opposite direction
        for (int i = 0; i < body2.NumPoints; i++)
        {
            double s = DotProduct(body2.Coord[i], direction);
            if (-s > localMax2)
            {
                localMax2 = -s;
                localBest2 = i;
            }
        }

        // Compute Minkowski difference point
        if (localBest1 >= 0 && localBest2 >= 0)
        {
            for (int i = 0; i < 3; i++)
            {
                result[i] = body1.Coord[localBest1][i] - body2.Coord[localBest2][i];
            }
            resultIdx[0] = localBest1;
            resultIdx[1] = localBest2;
        }
    }

    private static void SetContactNormal(GkSimplex smp, double[] contactNormal)
    {
        double[] d =
        [
            smp.Witnesses[1, 0] - smp.Witnesses[0, 0],
            smp.Witnesses[1, 1] - smp.Witnesses[0, 1],
            smp.Witnesses[1, 2] - smp.Witnesses[0, 2],
        ];
        double n = Math.Sqrt(Norm2(d));
        if (n > _gkEpsilon)
        {
            contactNormal[0] = d[0] / n;
            contactNormal[1] = d[1] / n;
            contactNormal[2] = d[2] / n;
        }
        else
        {
            contactNormal[0] = 1.0;
            contactNormal[1] = 0.0;
            contactNormal[2] = 0.0;
        }
    }

    //*******************************************************************************************
    // Entry point to the EPA Implementation
    //*******************************************************************************************

    /// <summary>
    /// Runs the EPA algorithm on the final GJK simplex to compute the penetration depth
    /// and contact normal when the bodies collide.
    /// </summary>
    /// <param name="bd1">The first body.</param>
    /// <param name="bd2">The second body.</param>
    /// <param name="simplex">The simplex produced by <see cref="ComputeMinimumDistance(GkPolytope, GkPolytope, GkSimplex)"/>.</param>
    /// <param name="distance">On input the GJK distance; on output negative penetration depth when colliding.</param>
    /// <param name="contactNormal">Receives the contact normal (length 3).</param>
    public static void ComputeCollisionInformation(GkPolytope bd1, GkPolytope bd2,
        GkSimplex simplex, ref double distance, double[] contactNormal)
    {
        const double epsSq = _gkEpsilon * _gkEpsilon;

        // if distance isn't 0 didn't detect collision - skip EPA
        if (distance > _gkEpsilon)
        {
            SetContactNormal(simplex, contactNormal);
            return;
        }

        // If GJK returned a degenerate simplex, rebuild it properly for EPA
        if (simplex.NVrtx != 4)
        {
            if (simplex.NVrtx == 1)
            {
                // Grow simplex from a single point: fire a support in some direction.
                // We use current simplex point for new direction for the
                // support; if this does not produce a new point, treat penetration as 0.
                double[] newVertex = new double[3];
                int[] newVertexIdx = new int[2];

                double[] dir = [simplex.Vrtx[0, 0], simplex.Vrtx[0, 1], simplex.Vrtx[0, 2]];
                if (Norm2(dir) < epsSq)
                {
                    // Zero initial direction (e.g. both bodies share their first vertex):
                    // fall back to an arbitrary axis so the simplex can still grow.
                    dir = [1.0, 0.0, 0.0];
                }

                SupportEpa(bd1, bd2, dir, newVertex, newVertexIdx);

                // Check if this is a new point relative to the existing simplex vertex.
                double dx = newVertex[0] - simplex.Vrtx[0, 0];
                double dy = newVertex[1] - simplex.Vrtx[0, 1];
                double dz = newVertex[2] - simplex.Vrtx[0, 2];
                bool isNew = dx * dx + dy * dy + dz * dz >= epsSq;

                if (isNew)
                {
                    int idx = simplex.NVrtx;
                    for (int c = 0; c < 3; ++c)
                    {
                        simplex.Vrtx[idx, c] = newVertex[c];
                    }
                    simplex.VrtxIdx[idx, 0] = newVertexIdx[0];
                    simplex.VrtxIdx[idx, 1] = newVertexIdx[1];
                    simplex.NVrtx = 2;
                }
                else
                {
                    // No new support point means penetration depth effectively zero.
                    distance = 0;
                    for (int c = 0; c < 3; ++c)
                    {
                        simplex.Witnesses[0, c] = bd1.Coord[newVertexIdx[0]][c];
                        simplex.Witnesses[1, c] = bd2.Coord[newVertexIdx[1]][c];
                    }
                    SetContactNormal(simplex, contactNormal);
                    return;
                }
            }
            if (simplex.NVrtx == 2)
            {
                // Grow simplex from an edge: fire a support in a direction perpendicular
                // to the edge. If this does not produce a new point, treat penetration as 0.
                double[] dir = new double[3];
                double[] newVertex = new double[3];
                int[] newVertexIdx = new int[2];

                double[] edge = new double[3];
                for (int c = 0; c < 3; ++c)
                {
                    edge[c] = simplex.Vrtx[1, c] - simplex.Vrtx[0, c];
                }

                // Build a perpendicular
                double[] axis = [1.0, 0.0, 0.0];
                double edgeNorm = Math.Sqrt(Norm2(edge));
                if (edgeNorm > _gkEpsilon && Math.Abs(edge[0]) > 0.9 * edgeNorm)
                {
                    axis[0] = 0.0;
                    axis[1] = 1.0;
                    axis[2] = 0.0;
                }

                // dir = edge x axis
                CrossProduct(edge, axis, dir);
                if (Norm2(dir) < _gkEpsilon)
                {
                    // Fallback axis
                    axis[0] = 0.0;
                    axis[1] = 0.0;
                    axis[2] = 1.0;
                    CrossProduct(edge, axis, dir);
                }

                SupportEpa(bd1, bd2, dir, newVertex, newVertexIdx);

                // Check if this is a new point relative to both existing simplex vertices.
                bool isNew = true;
                for (int vtx = 0; vtx < simplex.NVrtx; ++vtx)
                {
                    double dx = newVertex[0] - simplex.Vrtx[vtx, 0];
                    double dy = newVertex[1] - simplex.Vrtx[vtx, 1];
                    double dz = newVertex[2] - simplex.Vrtx[vtx, 2];
                    if (dx * dx + dy * dy + dz * dz < epsSq)
                    {
                        isNew = false;
                        break;
                    }
                }

                if (isNew)
                {
                    int idx = simplex.NVrtx;
                    for (int c = 0; c < 3; ++c)
                    {
                        simplex.Vrtx[idx, c] = newVertex[c];
                    }
                    simplex.VrtxIdx[idx, 0] = newVertexIdx[0];
                    simplex.VrtxIdx[idx, 1] = newVertexIdx[1];
                    simplex.NVrtx = 3;
                }
                else
                {
                    // No new support point means penetration depth effectively zero.
                    distance = 0;
                    for (int c = 0; c < 3; ++c)
                    {
                        simplex.Witnesses[0, c] = bd1.Coord[newVertexIdx[0]][c];
                        simplex.Witnesses[1, c] = bd2.Coord[newVertexIdx[1]][c];
                    }
                    SetContactNormal(simplex, contactNormal);
                    return;
                }
            }
            if (simplex.NVrtx == 3)
            {
                // Grow simplex from a triangle: fire a support in the direction of the
                // triangle normal. If this does not produce a new point, treat penetration as 0.
                double[] dir = new double[3];
                double[] newVertex = new double[3];
                int[] newVertexIdx = new int[2];

                double[] e0 = new double[3];
                double[] e1 = new double[3];
                for (int c = 0; c < 3; ++c)
                {
                    e0[c] = simplex.Vrtx[1, c] - simplex.Vrtx[0, c];
                    e1[c] = simplex.Vrtx[2, c] - simplex.Vrtx[0, c];
                }
                // dir = e0 x e1 (normal to the triangle)
                CrossProduct(e0, e1, dir);

                SupportEpa(bd1, bd2, dir, newVertex, newVertexIdx);

                // Check if this is a new point relative to all three existing simplex vertices.
                bool isNew = true;
                for (int vtx = 0; vtx < simplex.NVrtx; ++vtx)
                {
                    double dx = newVertex[0] - simplex.Vrtx[vtx, 0];
                    double dy = newVertex[1] - simplex.Vrtx[vtx, 1];
                    double dz = newVertex[2] - simplex.Vrtx[vtx, 2];
                    if (dx * dx + dy * dy + dz * dz < epsSq)
                    {
                        isNew = false;
                        break;
                    }
                }

                if (isNew)
                {
                    int idx = simplex.NVrtx;
                    for (int c = 0; c < 3; ++c)
                    {
                        simplex.Vrtx[idx, c] = newVertex[c];
                    }
                    simplex.VrtxIdx[idx, 0] = newVertexIdx[0];
                    simplex.VrtxIdx[idx, 1] = newVertexIdx[1];
                    simplex.NVrtx = 4;
                }
                else
                {
                    // Try opposite direction
                    dir[0] = -dir[0];
                    dir[1] = -dir[1];
                    dir[2] = -dir[2];
                }

                // If first direction didn't work, try opposite
                if (simplex.NVrtx == 3)
                {
                    SupportEpa(bd1, bd2, dir, newVertex, newVertexIdx);
                    isNew = true;
                    for (int vtx = 0; vtx < simplex.NVrtx; ++vtx)
                    {
                        double dx = newVertex[0] - simplex.Vrtx[vtx, 0];
                        double dy = newVertex[1] - simplex.Vrtx[vtx, 1];
                        double dz = newVertex[2] - simplex.Vrtx[vtx, 2];
                        if (dx * dx + dy * dy + dz * dz < epsSq)
                        {
                            isNew = false;
                            break;
                        }
                    }

                    if (isNew)
                    {
                        int idx = simplex.NVrtx;
                        for (int c = 0; c < 3; ++c)
                        {
                            simplex.Vrtx[idx, c] = newVertex[c];
                        }
                        simplex.VrtxIdx[idx, 0] = newVertexIdx[0];
                        simplex.VrtxIdx[idx, 1] = newVertexIdx[1];
                        simplex.NVrtx = 4;
                    }
                    else
                    {
                        distance = 0;
                        for (int c = 0; c < 3; ++c)
                        {
                            simplex.Witnesses[0, c] = bd1.Coord[newVertexIdx[0]][c];
                            simplex.Witnesses[1, c] = bd2.Coord[newVertexIdx[1]][c];
                        }
                        SetContactNormal(simplex, contactNormal);
                        return;
                    }
                }
            }

            // If we still don't have 4 vertices, abort
            if (simplex.NVrtx != 4)
            {
                distance = 0;
                // Set witness points from best available simplex vertex
                int best = simplex.NVrtx > 0 ? simplex.NVrtx - 1 : 0;
                for (int c = 0; c < 3; ++c)
                {
                    simplex.Witnesses[0, c] = bd1.Coord[simplex.VrtxIdx[best, 0]][c];
                    simplex.Witnesses[1, c] = bd2.Coord[simplex.VrtxIdx[best, 1]][c];
                }
                SetContactNormal(simplex, contactNormal);
                return;
            }
        }

        // On to actual EPA alg with a valid tetrahedron simplex
        // Initialize EPA polytope from simplex
        var poly = new EpaPolytope();
        double[] centroid = new double[3];
        InitEpaPolytope(poly, simplex, centroid);

        // EPA iteration parameters
        const int maxIterations = 64;
        const double tolerance = _gkEpsilon * 1e2;
        int iteration = 0;

        // Main EPA loop
        while (iteration < maxIterations && poly.NumVertices < EpaPolytope.MaxEpaVertices - 1)
        {
            iteration++;

            // Recompute normals & distances for assigned faces
            for (int i = 0; i < poly.MaxFaceIndex; ++i)
            {
                if (poly.Faces[i].Valid)
                {
                    ComputeFaceNormalDistance(poly, i, centroid);
                }
            }

            // Find the closest face
            int closestFace = -1;
            double closestDistance = 1e10;

            for (int i = 0; i < poly.MaxFaceIndex; ++i)
            {
                if (!poly.Faces[i].Valid)
                {
                    continue;
                }
                if (poly.Faces[i].Distance >= 0 && poly.Faces[i].Distance < closestDistance)
                {
                    closestDistance = poly.Faces[i].Distance;
                    closestFace = i;
                }
            }

            if (closestFace < 0)
            {
                break;
            }

            EpaFace closest = poly.Faces[closestFace];

            // Get support point in direction of closest face normal
            double[] newVertex = new double[3];
            int[] newVertexIdx = new int[2];
            SupportEpa(bd1, bd2, closest.Normal, newVertex, newVertexIdx);

            // Check termination condition: if distance to new vertex along normal
            // is not more than tolerance further than closest face
            double distToNew = DotProduct(closest.Normal, newVertex);
            double improvement = distToNew - closestDistance;

            if (improvement < tolerance)
            {
                // Converged, compute witness points with bary coords
                ComputeBarycentricOrigin(poly.Vertices[closest.V[0]],
                    poly.Vertices[closest.V[1]],
                    poly.Vertices[closest.V[2]], out double a0, out double a1, out double a2);
                for (int i = 0; i < 3; i++)
                {
                    simplex.Witnesses[0, i] = bd1.Coord[closest.VIdx[0][0]][i] * a0
                        + bd1.Coord[closest.VIdx[1][0]][i] * a1
                        + bd1.Coord[closest.VIdx[2][0]][i] * a2;
                    simplex.Witnesses[1, i] = bd2.Coord[closest.VIdx[0][1]][i] * a0
                        + bd2.Coord[closest.VIdx[1][1]][i] * a1
                        + bd2.Coord[closest.VIdx[2][1]][i] * a2;
                    contactNormal[i] = closest.Normal[i];
                }
                distance = -closestDistance;
                break;
            }

            // Check if new vertex is duplicate
            bool isDuplicate = false;
            for (int i = 0; i < poly.NumVertices; i++)
            {
                double dx = newVertex[0] - poly.Vertices[i][0];
                double dy = newVertex[1] - poly.Vertices[i][1];
                double dz = newVertex[2] - poly.Vertices[i][2];
                if (dx * dx + dy * dy + dz * dz < epsSq)
                {
                    isDuplicate = true;
                    break;
                }
            }

            if (isDuplicate)
            {
                // Can't make progress, use current best
                ComputeBarycentricOrigin(poly.Vertices[closest.V[0]],
                    poly.Vertices[closest.V[1]],
                    poly.Vertices[closest.V[2]], out double a0, out double a1, out double a2);
                for (int i = 0; i < 3; i++)
                {
                    simplex.Witnesses[0, i] = bd1.Coord[closest.VIdx[0][0]][i] * a0
                        + bd1.Coord[closest.VIdx[1][0]][i] * a1
                        + bd1.Coord[closest.VIdx[2][0]][i] * a2;
                    simplex.Witnesses[1, i] = bd2.Coord[closest.VIdx[0][1]][i] * a0
                        + bd2.Coord[closest.VIdx[1][1]][i] * a1
                        + bd2.Coord[closest.VIdx[2][1]][i] * a2;
                    contactNormal[i] = closest.Normal[i];
                }
                distance = -closestDistance;
                break;
            }

            // Add new vertex to polytope
            int newVertexId = poly.NumVertices;
            for (int i = 0; i < 3; i++)
            {
                poly.Vertices[newVertexId][i] = newVertex[i];
            }
            poly.VertexIndices[newVertexId][0] = newVertexIdx[0];
            poly.VertexIndices[newVertexId][1] = newVertexIdx[1];
            poly.NumVertices++;

            // Update centroid incrementally (running mean)
            double invN = 1.0 / poly.NumVertices;
            for (int i = 0; i < 3; i++)
            {
                centroid[i] += (newVertex[i] - centroid[i]) * invN;
            }

            // Find horizon edges: collect edges from faces being removed this iteration
            // only, then mark them invalid. Collecting from ALL invalid faces (including
            // ones from previous iterations) would pull in stale interior edges.
            var edges = new EpaEdge[EpaPolytope.MaxEpaFaces * 3];
            int numEdges = 0;

            for (int f = 0; f < poly.MaxFaceIndex; f++)
            {
                if (!poly.Faces[f].Valid)
                {
                    continue;
                }
                if (!IsFaceVisible(poly, f, newVertex))
                {
                    continue;
                }

                // Collect edges before invalidating the face
                for (int e = 0; e < 3; e++)
                {
                    if (numEdges >= EpaPolytope.MaxEpaFaces * 3)
                    {
                        break;
                    }

                    int va = e;
                    int vb = (e + 1) % 3;
                    edges[numEdges] = new EpaEdge
                    {
                        V1 = poly.Faces[f].V[va],
                        V2 = poly.Faces[f].V[vb],
                        VIdx1 = [poly.Faces[f].VIdx[va][0], poly.Faces[f].VIdx[va][1]],
                        VIdx2 = [poly.Faces[f].VIdx[vb][0], poly.Faces[f].VIdx[vb][1]],
                        Valid = true,
                    };
                    numEdges++;
                }

                poly.Faces[f].Valid = false;
            }

            // Remove duplicate edges (edges shared by two removed faces)
            for (int i = 0; i < numEdges; i++)
            {
                if (!edges[i].Valid)
                {
                    continue;
                }

                for (int j = i + 1; j < numEdges; j++)
                {
                    if (!edges[j].Valid)
                    {
                        continue;
                    }

                    // Check if same edge (either direction)
                    if ((edges[i].V1 == edges[j].V1 && edges[i].V2 == edges[j].V2) ||
                        (edges[i].V1 == edges[j].V2 && edges[i].V2 == edges[j].V1))
                    {
                        edges[i].Valid = false;
                        edges[j].Valid = false;
                    }
                }
            }

            // Create new faces from horizon edges
            for (int i = 0; i < numEdges; i++)
            {
                if (!edges[i].Valid)
                {
                    continue;
                }

                // Find next available face slot
                int newFaceIdx = -1;
                for (int j = 0; j < EpaPolytope.MaxEpaFaces; j++)
                {
                    if (!poly.Faces[j].Valid)
                    {
                        newFaceIdx = j;
                        break;
                    }
                }

                if (newFaceIdx < 0 || newFaceIdx >= EpaPolytope.MaxEpaFaces)
                {
                    break;
                }

                // Create new face: edge horizon vertices + new vertex
                EpaFace newFace = poly.Faces[newFaceIdx];
                newFace.V[0] = edges[i].V1;
                newFace.V[1] = edges[i].V2;
                newFace.V[2] = newVertexId;

                newFace.VIdx[0][0] = edges[i].VIdx1[0];
                newFace.VIdx[0][1] = edges[i].VIdx1[1];
                newFace.VIdx[1][0] = edges[i].VIdx2[0];
                newFace.VIdx[1][1] = edges[i].VIdx2[1];
                newFace.VIdx[2][0] = newVertexIdx[0];
                newFace.VIdx[2][1] = newVertexIdx[1];

                newFace.Valid = true;

                // Check winding and fix if necessary
                double[] fv0 = poly.Vertices[newFace.V[0]];
                double[] fv1 = poly.Vertices[newFace.V[1]];
                double[] fv2 = poly.Vertices[newFace.V[2]];

                double[] fe0 = new double[3];
                double[] fe1 = new double[3];
                double[] fnormal = new double[3];

                for (int c = 0; c < 3; c++)
                {
                    fe0[c] = fv1[c] - fv0[c];
                    fe1[c] = fv2[c] - fv0[c];
                }
                CrossProduct(fe0, fe1, fnormal);

                // If normal points toward centroid flip winding
                double[] toCent = new double[3];
                for (int c = 0; c < 3; c++)
                {
                    toCent[c] = centroid[c] - fv0[c];
                }
                if (DotProduct(fnormal, toCent) > 0)
                {
                    (newFace.V[1], newFace.V[2]) = (newFace.V[2], newFace.V[1]);
                    (newFace.VIdx[1], newFace.VIdx[2]) = (newFace.VIdx[2], newFace.VIdx[1]);
                }

                // Update max face index
                if (newFaceIdx >= poly.MaxFaceIndex)
                {
                    poly.MaxFaceIndex = newFaceIdx + 1;
                }
            }
        }

        // If we exited due to max iterations, recompute closest face and use it
        if (iteration >= maxIterations)
        {
            // Find closest face and compute result
            for (int i = 0; i < poly.MaxFaceIndex; ++i)
            {
                if (!poly.Faces[i].Valid)
                {
                    continue;
                }
                ComputeFaceNormalDistance(poly, i, centroid);
            }

            int closestFace = -1;
            double closestDistance = 1e10;
            for (int i = 0; i < poly.MaxFaceIndex; ++i)
            {
                if (!poly.Faces[i].Valid)
                {
                    continue;
                }
                if (poly.Faces[i].Distance >= 0 && poly.Faces[i].Distance < closestDistance)
                {
                    closestDistance = poly.Faces[i].Distance;
                    closestFace = i;
                }
            }

            if (closestFace >= 0)
            {
                EpaFace closest = poly.Faces[closestFace];
                ComputeBarycentricOrigin(poly.Vertices[closest.V[0]],
                    poly.Vertices[closest.V[1]],
                    poly.Vertices[closest.V[2]], out double a0, out double a1, out double a2);
                for (int i = 0; i < 3; i++)
                {
                    simplex.Witnesses[0, i] = bd1.Coord[closest.VIdx[0][0]][i] * a0
                        + bd1.Coord[closest.VIdx[1][0]][i] * a1
                        + bd1.Coord[closest.VIdx[2][0]][i] * a2;
                    simplex.Witnesses[1, i] = bd2.Coord[closest.VIdx[0][1]][i] * a0
                        + bd2.Coord[closest.VIdx[1][1]][i] * a1
                        + bd2.Coord[closest.VIdx[2][1]][i] * a2;
                    contactNormal[i] = closest.Normal[i];
                }
                distance = -closestDistance;
            }
        }
    }

    #endregion
}

using OpenGJKSharp.Models;
using System.Numerics;

namespace OpenGJKSharp;

public static partial class OpenGjk
{
    // Relative to the largest coordinate seen, so it sits above the single-precision
    // quantization of the support points (about 6e-8 relative) at any scale.
    private const double _castTolerance = 1e-6;

    // Twice the worst case measured over random rays against a capsule (31 iterations; rays
    // against polytopes need far fewer). Each advancement along the ray adds iterations to what
    // a plain distance query needs, and curved shapes hit at a grazing angle use the most.
    private const int _castMaxIterations = 64;

    /// <summary>
    /// Casts a ray against a support-mapped convex shape.
    /// </summary>
    /// <param name="shape">The shape to cast against.</param>
    /// <param name="origin">Start of the ray.</param>
    /// <param name="translation">
    /// The ray as a displacement: it covers the segment from <paramref name="origin"/> to
    /// <paramref name="origin"/> + <paramref name="translation"/>. A zero translation tests
    /// whether <paramref name="origin"/> lies in the shape.
    /// </param>
    /// <param name="fraction">
    /// On a hit, the fraction of <paramref name="translation"/> in [0, 1] at which the ray first
    /// touches the shape: the hit point is <paramref name="origin"/> + fraction * <paramref name="translation"/>.
    /// 0 when the ray misses.
    /// </param>
    /// <param name="normal">
    /// On a hit, the unit outward surface normal at the hit point. The zero vector when the ray
    /// misses, or when <paramref name="origin"/> already lies in the shape (or on its surface).
    /// </param>
    /// <returns><c>true</c> when the ray touches the shape.</returns>
    /// <remarks>
    /// <para>
    /// Uses the GJK ray cast (Gino van den Bergen, "Ray Casting against General Convex Objects
    /// with Application to Continuous Collision Detection", 2004): GJK runs on the distance
    /// between the shape and a point that advances along the ray to each separating plane it
    /// finds, so the reported hit point never lies past the true one.
    /// </para>
    /// <para>
    /// The normal is that of the last separating plane. It is exact on a flat face; on a curved
    /// surface it lags slightly behind the true normal, by about 1e-3 for unit-sized shapes.
    /// </para>
    /// <para>
    /// A ray that starts inside the shape (or on its surface) reports a hit with
    /// <paramref name="fraction"/> 0 and a zero <paramref name="normal"/>: there is no entry point
    /// to report a normal for.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="shape"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="origin"/> or <paramref name="translation"/> has a non-finite coordinate.</exception>
    public static bool TryCastRay(ISupportMappable shape, Vector3 origin, Vector3 translation,
        out double fraction, out Vector3 normal)
    {
        return TryCastRay(shape, origin, translation, new GjkWorkspace(), out fraction, out normal);
    }

    /// <summary>
    /// Casts a ray against a support-mapped convex shape, reusing the buffers of
    /// <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="shape">The shape to cast against.</param>
    /// <param name="origin">Start of the ray.</param>
    /// <param name="translation">
    /// The ray as a displacement: it covers the segment from <paramref name="origin"/> to
    /// <paramref name="origin"/> + <paramref name="translation"/>. A zero translation tests
    /// whether <paramref name="origin"/> lies in the shape.
    /// </param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="fraction">
    /// On a hit, the fraction of <paramref name="translation"/> in [0, 1] at which the ray first
    /// touches the shape: the hit point is <paramref name="origin"/> + fraction * <paramref name="translation"/>.
    /// 0 when the ray misses.
    /// </param>
    /// <param name="normal">
    /// On a hit, the unit outward surface normal at the hit point. The zero vector when the ray
    /// misses, or when <paramref name="origin"/> already lies in the shape (or on its surface).
    /// </param>
    /// <returns><c>true</c> when the ray touches the shape.</returns>
    /// <inheritdoc cref="TryCastRay(ISupportMappable, Vector3, Vector3, out double, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="shape"/> or <paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="origin"/> or <paramref name="translation"/> has a non-finite coordinate.</exception>
    public static bool TryCastRay(ISupportMappable shape, Vector3 origin, Vector3 translation, GjkWorkspace workspace,
        out double fraction, out Vector3 normal)
    {
        if (shape is null)
        {
            throw new ArgumentNullException(nameof(shape));
        }

        ValidateFinite(origin, nameof(origin));
        ValidateFinite(translation, nameof(translation));
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            // The ray becomes a cast from the origin against shape - {origin}.
            ReadOnlySpan<double> start = stackalloc double[] { origin.X, origin.Y, origin.Z };
            workspace.PolytopeA.Load(start);
            workspace.PolytopeBodyA.Reset(workspace.PolytopeA);
            workspace.MappedBodyB.Reset(shape);

            double[] n = workspace.ContactNormal;
            bool hit = Cast(workspace.MappedBodyB, workspace.PolytopeBodyA, translation, workspace.Simplex, n,
                out fraction);

            // The cast normal points from the ray into the shape; the surface normal is its opposite.
            normal = -new Vector3((float)n[0], (float)n[1], (float)n[2]);
            return hit;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Sweeps a support-mapped convex shape along a straight line and finds where it first
    /// touches another one.
    /// </summary>
    /// <param name="a">The moving shape, at its start pose.</param>
    /// <param name="b">The static shape.</param>
    /// <param name="translation">
    /// The displacement of <paramref name="a"/> over the sweep. A zero translation tests whether
    /// the shapes overlap at the start pose.
    /// </param>
    /// <param name="fraction">
    /// On a hit, the fraction of <paramref name="translation"/> in [0, 1] at which the shapes
    /// first touch: <paramref name="a"/> moved by fraction * <paramref name="translation"/> touches
    /// <paramref name="b"/>. 0 when they never touch.
    /// </param>
    /// <param name="normal">
    /// On a hit, the unit contact normal at first contact, pointing from <paramref name="a"/>
    /// toward <paramref name="b"/>. The zero vector when the shapes never touch, or when they
    /// already overlap (or touch) at the start pose.
    /// </param>
    /// <param name="point">
    /// On a hit, the contact point on the surface of <paramref name="b"/>. When the shapes already
    /// overlap at the start pose, a point inside both of them. The zero vector when they never touch.
    /// </param>
    /// <returns><c>true</c> when <paramref name="a"/> touches <paramref name="b"/> anywhere along the sweep.</returns>
    /// <remarks>
    /// <para>
    /// The sweep reduces to a ray cast from the origin against the Minkowski difference
    /// <paramref name="b"/> - <paramref name="a"/>, solved with the GJK ray cast (Gino van den
    /// Bergen, 2004). The reported fraction never lies past the true time of impact, so moving
    /// <paramref name="a"/> by it never leaves the shapes overlapping.
    /// </para>
    /// <para>
    /// The normal is that of the last separating plane. It is exact where a flat face is hit;
    /// where a curved surface is involved it lags slightly behind the true normal, by about 1e-3
    /// for unit-sized shapes, and the contact point by a similar amount.
    /// </para>
    /// <para>
    /// Shapes that already overlap (or touch) at the start pose report a hit with
    /// <paramref name="fraction"/> 0 and a zero <paramref name="normal"/>; use
    /// <see cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)"/>
    /// to separate them.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="translation"/> has a non-finite coordinate.</exception>
    public static bool TryCastShape(ISupportMappable a, ISupportMappable b, Vector3 translation,
        out double fraction, out Vector3 normal, out Vector3 point)
    {
        return TryCastShape(a, b, translation, new GjkWorkspace(), out fraction, out normal, out point);
    }

    /// <summary>
    /// Sweeps a support-mapped convex shape along a straight line and finds where it first
    /// touches another one, reusing the buffers of <paramref name="workspace"/> instead of
    /// allocating them.
    /// </summary>
    /// <param name="a">The moving shape, at its start pose.</param>
    /// <param name="b">The static shape.</param>
    /// <param name="translation">
    /// The displacement of <paramref name="a"/> over the sweep. A zero translation tests whether
    /// the shapes overlap at the start pose.
    /// </param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="fraction">
    /// On a hit, the fraction of <paramref name="translation"/> in [0, 1] at which the shapes
    /// first touch: <paramref name="a"/> moved by fraction * <paramref name="translation"/> touches
    /// <paramref name="b"/>. 0 when they never touch.
    /// </param>
    /// <param name="normal">
    /// On a hit, the unit contact normal at first contact, pointing from <paramref name="a"/>
    /// toward <paramref name="b"/>. The zero vector when the shapes never touch, or when they
    /// already overlap (or touch) at the start pose.
    /// </param>
    /// <param name="point">
    /// On a hit, the contact point on the surface of <paramref name="b"/>. When the shapes already
    /// overlap at the start pose, a point inside both of them. The zero vector when they never touch.
    /// </param>
    /// <returns><c>true</c> when <paramref name="a"/> touches <paramref name="b"/> anywhere along the sweep.</returns>
    /// <inheritdoc cref="TryCastShape(ISupportMappable, ISupportMappable, Vector3, out double, out Vector3, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="a"/>, <paramref name="b"/> or <paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="translation"/> has a non-finite coordinate.</exception>
    public static bool TryCastShape(ISupportMappable a, ISupportMappable b, Vector3 translation, GjkWorkspace workspace,
        out double fraction, out Vector3 normal, out Vector3 point)
    {
        ValidateShapes(a, b);
        ValidateFinite(translation, nameof(translation));
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            // The cast runs on b - a, so b takes the first body slot.
            workspace.LoadShapes(b, a, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;
            double[] n = workspace.ContactNormal;

            if (!Cast(bd1, bd2, translation, simplex, n, out fraction))
            {
                normal = Vector3.Zero;
                point = Vector3.Zero;
                return false;
            }

            normal = new Vector3((float)n[0], (float)n[1], (float)n[2]);
            point = new Vector3(
                (float)simplex.Witnesses[0, 0],
                (float)simplex.Witnesses[0, 1],
                (float)simplex.Witnesses[0, 2]);
            return true;
        }
        finally
        {
            workspace.Release();
        }
    }

    // GJK ray cast (van den Bergen 2004) from the origin along translation against the
    // Minkowski difference bd1 - bd2. Simplex vertices are stored relative to the current ray
    // point x, so GJK measures the distance from x to the difference. Whenever the search
    // direction yields a plane separating x from the difference, x jumps to that plane - no
    // point before it can touch - and the simplex vertices move with it. On a hit, normal holds
    // the unit normal of the last such plane, pointing from the ray into the difference (zero
    // when x never moved), and Witnesses[0] the matching point of bd1; on a miss normal is zero.
    private static bool Cast(GkBody bd1, GkBody bd2, Vector3 translation, GkSimplex s, double[] normal,
        out double fraction)
    {
        double rX = translation.X, rY = translation.Y, rZ = translation.Z;
        GkScratch scratch = s.Scratch;
        double[] v = scratch.V;
        double[] vMinus = scratch.VMinus;
        double[] w = scratch.W;
        Span<double> weights = stackalloc double[4];
        double lambda = 0;
        double xX = 0, xY = 0, xZ = 0;
        double nX = 0, nY = 0, nZ = 0;

        normal[0] = 0;
        normal[1] = 0;
        normal[2] = 0;

        bd1.InitSupport();
        bd2.InitSupport();

        // The interior points give an initial direction from x toward the difference.
        for (int t = 0; t < 3; t++)
        {
            v[t] = bd1.S[t] - bd2.S[t];
        }

        double scale2 = Math.Max(Norm2(bd1.S), Norm2(bd2.S));
        s.NVrtx = 0;

        for (int k = 0; k < _castMaxIterations; k++)
        {
            for (int t = 0; t < 3; t++)
            {
                vMinus[t] = -v[t];
            }

            bd1.Support(vMinus);
            bd2.Support(v);
            scale2 = Math.Max(scale2, Math.Max(Norm2(bd1.S), Norm2(bd2.S)));

            w[0] = bd1.S[0] - bd2.S[0] - xX;
            w[1] = bd1.S[1] - bd2.S[1] - xY;
            w[2] = bd1.S[2] - bd2.S[2] - xZ;

            double vw = DotProduct(v, w);
            if (vw > 0)
            {
                // The plane through the support point with normal v separates x from the difference.
                double vr = (v[0] * rX) + (v[1] * rY) + (v[2] * rZ);
                if (vr <= 0)
                {
                    fraction = 0;
                    return false;
                }

                lambda += vw / vr;
                if (lambda > 1)
                {
                    fraction = 0;
                    return false;
                }

                double dX = (lambda * rX) - xX, dY = (lambda * rY) - xY, dZ = (lambda * rZ) - xZ;
                xX += dX;
                xY += dY;
                xZ += dZ;
                w[0] -= dX;
                w[1] -= dY;
                w[2] -= dZ;
                for (int i = 0; i < s.NVrtx; i++)
                {
                    s.Vrtx[i, 0] -= dX;
                    s.Vrtx[i, 1] -= dY;
                    s.Vrtx[i, 2] -= dZ;
                }

                nX = v[0];
                nY = v[1];
                nZ = v[2];
            }

            int last = s.NVrtx;
            for (int t = 0; t < 3; t++)
            {
                s.Vrtx[last, t] = w[t];
            }

            s.VrtxIdx[last, 0] = bd1.SIdx;
            s.VrtxIdx[last, 1] = bd2.SIdx;
            s.NVrtx++;

            ClosestOnSimplex(s, v, weights);

            if (s.NVrtx == 4 || Norm2(v) <= _castTolerance * _castTolerance * scale2)
            {
                break;
            }
        }

        // Reaching the iteration cap still reports a hit at the last lambda: it never lies
        // past the true hit, which is the safe answer for continuous collision detection.
        double length = Math.Sqrt(Norm2(nX, nY, nZ));
        if (length > 0)
        {
            normal[0] = nX / length;
            normal[1] = nY / length;
            normal[2] = nZ / length;
        }

        for (int t = 0; t < 3; t++)
        {
            double sum = 0;
            for (int i = 0; i < s.NVrtx; i++)
            {
                sum += weights[i] * bd1.GetPoint(s.VrtxIdx[i, 0])[t];
            }

            s.Witnesses[0, t] = sum;
        }

        fraction = lambda;
        return true;
    }

    // Closest point to the origin of the simplex, found by trying every subset of its vertices;
    // reduces the simplex to the subset it lies on, with its barycentric weights in weights.
    // SubAlgorithm cannot be used here: it assumes the newest vertex belongs to the closest
    // feature, which stops holding once the ray point has moved under the older vertices.
    private static void ClosestOnSimplex(GkSimplex s, double[] v, Span<double> weights)
    {
        Span<double> candidate = stackalloc double[4];
        double bestNorm = double.PositiveInfinity;
        int bestMask = 1;

        for (int mask = 1; mask < (1 << s.NVrtx); mask++)
        {
            if (!TryProjectOnSubset(s, mask, candidate, out double pX, out double pY, out double pZ))
            {
                continue;
            }

            double norm = Norm2(pX, pY, pZ);
            if (norm < bestNorm)
            {
                bestNorm = norm;
                bestMask = mask;
                v[0] = pX;
                v[1] = pY;
                v[2] = pZ;
                candidate.CopyTo(weights);
            }
        }

        int count = 0;
        for (int i = 0; i < s.NVrtx; i++)
        {
            if ((bestMask & (1 << i)) == 0)
            {
                continue;
            }

            for (int t = 0; t < 3; t++)
            {
                s.Vrtx[count, t] = s.Vrtx[i, t];
            }

            s.VrtxIdx[count, 0] = s.VrtxIdx[i, 0];
            s.VrtxIdx[count, 1] = s.VrtxIdx[i, 1];
            weights[count] = weights[i];
            count++;
        }

        s.NVrtx = count;
    }

    // Projects the origin onto the affine hull of the vertices selected by mask. Succeeds when
    // the subset is not degenerate and the projection lies strictly inside it; weights then
    // holds its barycentric coordinates, indexed by vertex.
    private static bool TryProjectOnSubset(GkSimplex s, int mask, Span<double> weights,
        out double pX, out double pY, out double pZ)
    {
        Span<int> index = stackalloc int[4];
        int count = 0;
        for (int i = 0; i < 4; i++)
        {
            weights[i] = 0;
            if ((mask & (1 << i)) != 0)
            {
                index[count++] = i;
            }
        }

        int o = index[0];
        double oX = s.Vrtx[o, 0], oY = s.Vrtx[o, 1], oZ = s.Vrtx[o, 2];
        pX = oX;
        pY = oY;
        pZ = oZ;
        weights[o] = 1;

        // Least squares in the edges e_j = y_j - y_0: solve G mu = c with G_jk = e_j . e_k
        // and c_j = -e_j . y_0 by Cramer's rule.
        int d = count - 1;
        Span<double> e = stackalloc double[9];
        Span<double> g = stackalloc double[9];
        Span<double> c = stackalloc double[3];
        for (int j = 0; j < d; j++)
        {
            int vertex = index[j + 1];
            e[(j * 3) + 0] = s.Vrtx[vertex, 0] - oX;
            e[(j * 3) + 1] = s.Vrtx[vertex, 1] - oY;
            e[(j * 3) + 2] = s.Vrtx[vertex, 2] - oZ;
            c[j] = -((e[j * 3] * oX) + (e[(j * 3) + 1] * oY) + (e[(j * 3) + 2] * oZ));
        }

        for (int j = 0; j < d; j++)
        {
            for (int k = 0; k < d; k++)
            {
                g[(j * 3) + k] = (e[j * 3] * e[k * 3]) + (e[(j * 3) + 1] * e[(k * 3) + 1]) + (e[(j * 3) + 2] * e[(k * 3) + 2]);
            }
        }

        if (d == 0)
        {
            return true;
        }

        double det = DeterminantOf(g, d, c, -1);
        double diagonal = d == 1 ? g[0] : d == 2 ? g[0] * g[4] : g[0] * g[4] * g[8];
        if (!(det > 1e-12 * diagonal))
        {
            return false;
        }

        double first = 1;
        for (int j = 0; j < d; j++)
        {
            double mu = DeterminantOf(g, d, c, j) / det;
            if (!(mu > 0))
            {
                return false;
            }

            first -= mu;
            weights[index[j + 1]] = mu;
            pX += mu * e[j * 3];
            pY += mu * e[(j * 3) + 1];
            pZ += mu * e[(j * 3) + 2];
        }

        weights[o] = first;
        return first > 0;
    }

    // Determinant of the leading d x d block of the row-major 3x3 matrix m, with column
    // replaced by c (Cramer's rule) unless replaced is negative.
    private static double DeterminantOf(Span<double> m, int d, Span<double> c, int replaced)
    {
        Span<double> a = stackalloc double[9];
        m.CopyTo(a);
        for (int row = 0; replaced >= 0 && row < d; row++)
        {
            a[(row * 3) + replaced] = c[row];
        }

        return d switch
        {
            1 => a[0],
            2 => (a[0] * a[4]) - (a[1] * a[3]),
            _ => (a[0] * ((a[4] * a[8]) - (a[5] * a[7])))
                - (a[1] * ((a[3] * a[8]) - (a[5] * a[6])))
                + (a[2] * ((a[3] * a[7]) - (a[4] * a[6]))),
        };
    }

    private static void ValidateFinite(Vector3 value, string paramName)
    {
        if (!IsFinite(value.X) || !IsFinite(value.Y) || !IsFinite(value.Z))
        {
            throw new ArgumentException($"Value contains a non-finite coordinate: {value}.", paramName);
        }
    }
}

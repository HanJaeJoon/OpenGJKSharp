using System.Diagnostics;
using OpenGJKSharp.Models;
using System.Numerics;

namespace OpenGJKSharp;

/// <summary>
/// GJK (Gilbert-Johnson-Keerthi) collision detection between convex polytopes,
/// with EPA (Expanding Polytope Algorithm) support for penetration depth.
/// </summary>
/// <remarks>
/// <para>
/// The following holds for every query, whatever the input kind (vertex arrays, 2D polygons,
/// double-precision coordinates, or <see cref="ISupportMappable"/> shapes):
/// </para>
/// <list type="bullet">
/// <item><description>
/// <b>Distance on collision.</b> Colliding bodies report a distance of exactly 0: the
/// residual left by the iteration is measured against the size of the bodies and flushed to
/// zero. What remains is input quantization - with the <see cref="Vector3"/> overloads at
/// large coordinates, a pair meant to touch exactly may end up slightly separated or slightly
/// overlapping once the coordinates are rounded to single precision.
/// </description></item>
/// <item><description>
/// <b>Witness points on collision.</b> The closest points are only meaningful while the
/// bodies are separated, where they lie on the respective surfaces. Once the bodies overlap
/// there is no closest pair: both witness points collapse onto a common point inside the
/// intersection, which is not on either surface. Use <c>ComputeCollisionInformation</c> for
/// penetration depth and contact normal instead, or the <c>TryGetPenetration</c> overloads with
/// <c>pointA</c>/<c>pointB</c> for contact points on the surfaces.
/// </description></item>
/// <item><description>
/// <b>Contact normal.</b> It is a unit vector pointing from the first body toward the second,
/// so pushing the second body along it (or the first against it) resolves the penetration.
/// For exactly concentric bodies the direction is mathematically arbitrary - every direction
/// is a valid minimum translation - and an unspecified unit vector is returned.
/// </description></item>
/// <item><description>
/// <b>Penetration depth.</b> Exact for polyhedra. Curved shapes cannot be resolved exactly by
/// a finite polytope, so their depth is reported as a slight overestimate (within about 0.3%
/// for spheres) rather than an underestimate - pushing the bodies apart by it always resolves
/// the penetration.
/// </description></item>
/// <item><description>
/// <b>Precision.</b> Internal math runs in double precision. The <see cref="Vector3"/> and
/// <see cref="Vector2"/> overloads quantize their input to single precision; the
/// <see cref="ReadOnlySpan{T}"/> of <see cref="double"/> overloads do not.
/// </description></item>
/// </list>
/// </remarks>
public static partial class OpenGjk
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
        return HasCollision(a, b, new GjkWorkspace(), precision);
    }

    /// <summary>
    /// Determines whether two convex polytopes collide, reusing the buffers of
    /// <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="precision">Distance below which the bodies are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(Vector3[] a, Vector3[] b, GjkWorkspace workspace, double precision = _precision)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);
        ValidatePrecision(precision);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            double distance = ComputeMinimumDistance(bd1, bd2, workspace.Simplex, precision);

            return distance <= precision;
        }
        finally
        {
            workspace.Release();
        }
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
        return HasCollision(a, b, new GjkWorkspace(), precision);
    }

    /// <summary>
    /// Determines whether two convex 2D polygons collide (treated as flat shapes on the z = 0 plane),
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="precision">Distance below which the polygons are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(Vector2[] a, Vector2[] b, GjkWorkspace workspace, double precision = _precision)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);
        ValidatePrecision(precision);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            double distance = ComputeMinimumDistance(bd1, bd2, workspace.Simplex, precision);

            return distance <= precision;
        }
        finally
        {
            workspace.Release();
        }
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
        return ComputeCollisionInformation(a, b, new GjkWorkspace(), out contactNormal);
    }

    /// <summary>
    /// Computes collision information between two convex polytopes,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns the minimum distance when the bodies are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first body toward the second.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static double ComputeCollisionInformation(Vector3[] a, Vector3[] b, GjkWorkspace workspace, out Vector3 contactNormal)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;

            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

            double[] normal = workspace.ContactNormal;
            ComputeCollisionInformation(bd1, bd2, simplex, workspace.Epa, ref distance, normal);

            contactNormal = new Vector3((float)normal[0], (float)normal[1], (float)normal[2]);
            return distance;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Computes collision information between two convex 2D polygons
    /// (treated as flat shapes on the z = 0 plane).
    /// Returns the minimum distance when the polygons are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first polygon toward the second.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    /// <remarks>
    /// The depth and the normal are measured in the plane of the polygons, which is the whole
    /// problem for a 2D pair.
    /// </remarks>
    public static double ComputeCollisionInformation(Vector2[] a, Vector2[] b, out Vector2 contactNormal)
    {
        return ComputeCollisionInformation(a, b, new GjkWorkspace(), out contactNormal);
    }

    /// <summary>
    /// Computes collision information between two convex 2D polygons
    /// (treated as flat shapes on the z = 0 plane),
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns the minimum distance when the polygons are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first polygon toward the second.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    /// <inheritdoc cref="ComputeCollisionInformation(Vector2[], Vector2[], out Vector2)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static double ComputeCollisionInformation(Vector2[] a, Vector2[] b, GjkWorkspace workspace,
        out Vector2 contactNormal)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;

            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

            double[] normal = workspace.ContactNormal;
            ComputeCollisionInformation(bd1, bd2, simplex, workspace.Epa, ref distance, normal);

            contactNormal = new Vector2((float)normal[0], (float)normal[1]);
            return distance;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Determines whether two convex 2D polygons (treated as flat shapes on the z = 0 plane)
    /// collide and, if so, by how much.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the polygons are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first polygon toward the second.</param>
    /// <returns><c>true</c> when the polygons overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)" path="/remarks"/>
    public static bool TryGetPenetration(Vector2[] a, Vector2[] b, out double depth, out Vector2 contactNormal)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), out depth, out contactNormal);
    }

    /// <summary>
    /// Determines whether two convex 2D polygons (treated as flat shapes on the z = 0 plane)
    /// collide and, if so, by how much, reusing the buffers of <paramref name="workspace"/>
    /// instead of allocating them.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the polygons are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first polygon toward the second.</param>
    /// <returns><c>true</c> when the polygons overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static bool TryGetPenetration(Vector2[] a, Vector2[] b, GjkWorkspace workspace,
        out double depth, out Vector2 contactNormal)
    {
        double info = ComputeCollisionInformation(a, b, workspace, out contactNormal);

        return ToPenetration(info, out depth);
    }

    /// <summary>
    /// Determines whether two convex 2D polygons (treated as flat shapes on the z = 0 plane)
    /// collide and, if so, by how much and where they touch.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the polygons are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first polygon toward the second.</param>
    /// <param name="pointA">Contact point on the first polygon, or its closest point when the polygons are apart.</param>
    /// <param name="pointB">Contact point on the second polygon, or its closest point when the polygons are apart.</param>
    /// <returns><c>true</c> when the polygons overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3, out Vector3, out Vector3)" path="/remarks"/>
    public static bool TryGetPenetration(Vector2[] a, Vector2[] b, out double depth, out Vector2 contactNormal,
        out Vector2 pointA, out Vector2 pointB)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), out depth, out contactNormal, out pointA, out pointB);
    }

    /// <summary>
    /// Determines whether two convex 2D polygons (treated as flat shapes on the z = 0 plane)
    /// collide and, if so, by how much and where they touch, reusing the buffers of
    /// <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the polygons are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first polygon toward the second.</param>
    /// <param name="pointA">Contact point on the first polygon, or its closest point when the polygons are apart.</param>
    /// <param name="pointB">Contact point on the second polygon, or its closest point when the polygons are apart.</param>
    /// <returns><c>true</c> when the polygons overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3, out Vector3, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static bool TryGetPenetration(Vector2[] a, Vector2[] b, GjkWorkspace workspace,
        out double depth, out Vector2 contactNormal, out Vector2 pointA, out Vector2 pointB)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            double info = ComputeContact(bd1, bd2, workspace);

            double[] normal = workspace.ContactNormal;
            double[,] points = workspace.Simplex.Witnesses;
            contactNormal = new Vector2((float)normal[0], (float)normal[1]);
            pointA = new Vector2((float)points[0, 0], (float)points[0, 1]);
            pointB = new Vector2((float)points[1, 0], (float)points[1, 1]);

            return ToPenetration(info, out depth);
        }
        finally
        {
            workspace.Release();
        }
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
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static double ComputeMinimumDistance(Vector3[] a, Vector3[] b, GjkWorkspace workspace)
    {
        return ComputeMinimumDistance(a, b, workspace, out _, out _);
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
        return ComputeMinimumDistance(a, b, new GjkWorkspace(), out closestA, out closestB);
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes along with the closest point
    /// on each body (witness points), reusing the buffers of <paramref name="workspace"/>
    /// instead of allocating them.
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="closestA">Closest point on the first body.</param>
    /// <param name="closestB">Closest point on the second body.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static double ComputeMinimumDistance(Vector3[] a, Vector3[] b, GjkWorkspace workspace, out Vector3 closestA, out Vector3 closestB)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;
            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

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
        finally
        {
            workspace.Release();
        }
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
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns 0 when the polygons collide.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <returns>The minimum distance between the two polygons.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static double ComputeMinimumDistance(Vector2[] a, Vector2[] b, GjkWorkspace workspace)
    {
        return ComputeMinimumDistance(a, b, workspace, out _, out _);
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
        return ComputeMinimumDistance(a, b, new GjkWorkspace(), out closestA, out closestB);
    }

    /// <summary>
    /// Computes the minimum distance between two convex 2D polygons
    /// (treated as flat shapes on the z = 0 plane) along with the closest point on each polygon
    /// (witness points), reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns 0 when the polygons collide.
    /// </summary>
    /// <param name="a">Vertices of the first polygon.</param>
    /// <param name="b">Vertices of the second polygon.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="closestA">Closest point on the first polygon.</param>
    /// <param name="closestB">Closest point on the second polygon.</param>
    /// <returns>The minimum distance between the two polygons.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static double ComputeMinimumDistance(Vector2[] a, Vector2[] b, GjkWorkspace workspace, out Vector2 closestA, out Vector2 closestB)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;
            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

            closestA = new Vector2(
                (float)simplex.Witnesses[0, 0],
                (float)simplex.Witnesses[0, 1]);
            closestB = new Vector2(
                (float)simplex.Witnesses[1, 0],
                (float)simplex.Witnesses[1, 1]);

            return distance;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Determines whether two convex polytopes given in double precision collide.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="precision">Distance below which the bodies are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, or a coordinate is not finite.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(ReadOnlySpan<double> a, ReadOnlySpan<double> b, double precision = _precision)
    {
        return HasCollision(a, b, new GjkWorkspace(), precision);
    }

    /// <summary>
    /// Determines whether two convex polytopes given in double precision collide,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="precision">Distance below which the bodies are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, or a coordinate is not finite.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(ReadOnlySpan<double> a, ReadOnlySpan<double> b, GjkWorkspace workspace, double precision = _precision)
    {
        ValidateCoordinates(a, b);
        ValidateWorkspace(workspace);
        ValidatePrecision(precision);

        workspace.Occupy();
        try
        {
            LoadCoordinates(workspace, a, b, out GkBody bd1, out GkBody bd2);
            double distance = ComputeMinimumDistance(bd1, bd2, workspace.Simplex, precision);

            return distance <= precision;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes given in double precision.
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, or a coordinate is not finite.</exception>
    public static double ComputeMinimumDistance(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        return ComputeMinimumDistance(a, b, new GjkWorkspace());
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes given in double precision,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, or a coordinate is not finite.</exception>
    public static double ComputeMinimumDistance(ReadOnlySpan<double> a, ReadOnlySpan<double> b, GjkWorkspace workspace)
    {
        ValidateCoordinates(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            LoadCoordinates(workspace, a, b, out GkBody bd1, out GkBody bd2);

            return ComputeMinimumDistance(bd1, bd2, workspace.Simplex);
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes given in double precision,
    /// along with the closest point on each body (witness points).
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="closestA">Receives the closest point on the first body. Must hold at least 3 elements.</param>
    /// <param name="closestB">Receives the closest point on the second body. Must hold at least 3 elements.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or an output span is shorter than 3 elements.</exception>
    public static double ComputeMinimumDistance(ReadOnlySpan<double> a, ReadOnlySpan<double> b,
        Span<double> closestA, Span<double> closestB)
    {
        return ComputeMinimumDistance(a, b, new GjkWorkspace(), closestA, closestB);
    }

    /// <summary>
    /// Computes the minimum distance between two convex polytopes given in double precision,
    /// along with the closest point on each body (witness points),
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns 0 when the bodies collide.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="closestA">Receives the closest point on the first body. Must hold at least 3 elements.</param>
    /// <param name="closestB">Receives the closest point on the second body. Must hold at least 3 elements.</param>
    /// <returns>The minimum distance between the two bodies.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or an output span is shorter than 3 elements.</exception>
    public static double ComputeMinimumDistance(ReadOnlySpan<double> a, ReadOnlySpan<double> b,
        GjkWorkspace workspace, Span<double> closestA, Span<double> closestB)
    {
        ValidateCoordinates(a, b);
        ValidateWorkspace(workspace);
        ValidateOutput(closestA, nameof(closestA));
        ValidateOutput(closestB, nameof(closestB));

        workspace.Occupy();
        try
        {
            LoadCoordinates(workspace, a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;
            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

            for (int t = 0; t < 3; t++)
            {
                closestA[t] = simplex.Witnesses[0, t];
                closestB[t] = simplex.Witnesses[1, t];
            }

            return distance;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Computes collision information between two convex polytopes given in double precision.
    /// Returns the minimum distance when the bodies are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="contactNormal">Receives the contact normal, pointing from the first body toward the second. Must hold at least 3 elements.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or <paramref name="contactNormal"/> is shorter than 3 elements.</exception>
    public static double ComputeCollisionInformation(ReadOnlySpan<double> a, ReadOnlySpan<double> b, Span<double> contactNormal)
    {
        return ComputeCollisionInformation(a, b, new GjkWorkspace(), contactNormal);
    }

    /// <summary>
    /// Computes collision information between two convex polytopes given in double precision,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns the minimum distance when the bodies are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="contactNormal">Receives the contact normal, pointing from the first body toward the second. Must hold at least 3 elements.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or <paramref name="contactNormal"/> is shorter than 3 elements.</exception>
    public static double ComputeCollisionInformation(ReadOnlySpan<double> a, ReadOnlySpan<double> b,
        GjkWorkspace workspace, Span<double> contactNormal)
    {
        ValidateCoordinates(a, b);
        ValidateWorkspace(workspace);
        ValidateOutput(contactNormal, nameof(contactNormal));

        workspace.Occupy();
        try
        {
            LoadCoordinates(workspace, a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;

            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

            double[] normal = workspace.ContactNormal;
            ComputeCollisionInformation(bd1, bd2, simplex, workspace.Epa, ref distance, normal);

            for (int t = 0; t < 3; t++)
            {
                contactNormal[t] = normal[t];
            }

            return distance;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Determines whether two support-mapped convex shapes collide.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="precision">Distance below which the shapes are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(ISupportMappable a, ISupportMappable b, double precision = _precision)
    {
        return HasCollision(a, b, new GjkWorkspace(), precision);
    }

    /// <summary>
    /// Determines whether two support-mapped convex shapes collide,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="precision">Distance below which the shapes are considered colliding. Must not be negative.</param>
    /// <returns><c>true</c> when the minimum distance is within <paramref name="precision"/>.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/>, <paramref name="b"/> or <paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="precision"/> is negative.</exception>
    public static bool HasCollision(ISupportMappable a, ISupportMappable b, GjkWorkspace workspace, double precision = _precision)
    {
        ValidateShapes(a, b);
        ValidateWorkspace(workspace);
        ValidatePrecision(precision);

        workspace.Occupy();
        try
        {
            workspace.LoadShapes(a, b, out GkBody bd1, out GkBody bd2);
            double distance = ComputeMinimumDistance(bd1, bd2, workspace.Simplex, precision);

            return distance <= precision;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Computes the minimum distance between two support-mapped convex shapes.
    /// Returns 0 when the shapes collide.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <returns>The minimum distance between the two shapes.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
    public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b)
    {
        return ComputeMinimumDistance(a, b, out _, out _);
    }

    /// <summary>
    /// Computes the minimum distance between two support-mapped convex shapes,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns 0 when the shapes collide.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <returns>The minimum distance between the two shapes.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/>, <paramref name="b"/> or <paramref name="workspace"/> is null.</exception>
    public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b, GjkWorkspace workspace)
    {
        return ComputeMinimumDistance(a, b, workspace, out _, out _);
    }

    /// <summary>
    /// Computes the minimum distance between two support-mapped convex shapes,
    /// along with the closest point on each shape (witness points).
    /// Returns 0 when the shapes collide.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="closestA">Closest point on the first shape.</param>
    /// <param name="closestB">Closest point on the second shape.</param>
    /// <returns>The minimum distance between the two shapes.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
    public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b, out Vector3 closestA, out Vector3 closestB)
    {
        return ComputeMinimumDistance(a, b, new GjkWorkspace(), out closestA, out closestB);
    }

    /// <summary>
    /// Computes the minimum distance between two support-mapped convex shapes along with the
    /// closest point on each shape (witness points), reusing the buffers of
    /// <paramref name="workspace"/> instead of allocating them.
    /// Returns 0 when the shapes collide.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="closestA">Closest point on the first shape.</param>
    /// <param name="closestB">Closest point on the second shape.</param>
    /// <returns>The minimum distance between the two shapes.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/>, <paramref name="b"/> or <paramref name="workspace"/> is null.</exception>
    public static double ComputeMinimumDistance(ISupportMappable a, ISupportMappable b, GjkWorkspace workspace, out Vector3 closestA, out Vector3 closestB)
    {
        ValidateShapes(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            workspace.LoadShapes(a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;
            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

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
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Computes collision information between two support-mapped convex shapes.
    /// Returns the minimum distance when the shapes are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first shape toward the second.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
    public static double ComputeCollisionInformation(ISupportMappable a, ISupportMappable b, out Vector3 contactNormal)
    {
        return ComputeCollisionInformation(a, b, new GjkWorkspace(), out contactNormal);
    }

    /// <summary>
    /// Computes collision information between two support-mapped convex shapes,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// Returns the minimum distance when the shapes are separated (positive value),
    /// or the penetration depth computed by the EPA algorithm when they collide (negative value).
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first shape toward the second.</param>
    /// <returns>Positive separation distance, or negative penetration depth on collision.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="a"/>, <paramref name="b"/> or <paramref name="workspace"/> is null.</exception>
    public static double ComputeCollisionInformation(ISupportMappable a, ISupportMappable b, GjkWorkspace workspace, out Vector3 contactNormal)
    {
        ValidateShapes(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            // The same body instances must serve both the GJK and EPA passes:
            // the simplex records ids into each body's point store.
            workspace.LoadShapes(a, b, out GkBody bd1, out GkBody bd2);
            GkSimplex simplex = workspace.Simplex;

            double distance = ComputeMinimumDistance(bd1, bd2, simplex);

            double[] normal = workspace.ContactNormal;
            ComputeCollisionInformation(bd1, bd2, simplex, workspace.Epa, ref distance, normal);

            contactNormal = new Vector3((float)normal[0], (float)normal[1], (float)normal[2]);
            return distance;
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Determines whether two support-mapped convex shapes collide and, if so, by how much.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the shapes are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first shape toward the second.</param>
    /// <returns><c>true</c> when the shapes overlap or touch exactly.</returns>
    /// <remarks>
    /// The same query as <see cref="ComputeCollisionInformation(ISupportMappable, ISupportMappable, out Vector3)"/>,
    /// reported without a sign convention to interpret: move the second shape by
    /// <paramref name="depth"/> along <paramref name="contactNormal"/> (or the first against it)
    /// to separate the pair. When the shapes are apart the normal still points from the first
    /// toward the second, so it can be used as a separation direction.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
    public static bool TryGetPenetration(ISupportMappable a, ISupportMappable b, out double depth, out Vector3 contactNormal)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), out depth, out contactNormal);
    }

    /// <summary>
    /// Determines whether two support-mapped convex shapes collide and, if so, by how much,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the shapes are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first shape toward the second.</param>
    /// <returns><c>true</c> when the shapes overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="a"/>, <paramref name="b"/> or <paramref name="workspace"/> is null.</exception>
    public static bool TryGetPenetration(ISupportMappable a, ISupportMappable b, GjkWorkspace workspace,
        out double depth, out Vector3 contactNormal)
    {
        double info = ComputeCollisionInformation(a, b, workspace, out contactNormal);

        return ToPenetration(info, out depth);
    }

    /// <summary>
    /// Determines whether two support-mapped convex shapes collide and, if so, by how much and
    /// where they touch.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the shapes are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first shape toward the second.</param>
    /// <param name="pointA">Contact point on the first shape, or its closest point when the shapes are apart.</param>
    /// <param name="pointB">Contact point on the second shape, or its closest point when the shapes are apart.</param>
    /// <returns><c>true</c> when the shapes overlap or touch exactly.</returns>
    /// <remarks>
    /// <para>
    /// The same query as <see cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)"/>
    /// (same depth and normal), also reporting the contact points.
    /// </para>
    /// <para>
    /// While the bodies overlap, <paramref name="pointA"/> is the deepest point of the first body
    /// along <paramref name="contactNormal"/> and <paramref name="pointB"/> the deepest point of the
    /// second body against it. Both lie on the surface of their own body, and they are separated
    /// along the normal by exactly the depth: <c>dot(pointA - pointB, contactNormal) == depth</c>,
    /// up to rounding. Where flat faces meet, the two points face each other inside the contact
    /// region; on a curved surface the point is the body's support point along the normal. When
    /// the normal is arbitrary (concentric bodies) the points follow whatever normal is returned.
    /// </para>
    /// <para>
    /// When the bodies are apart, the points are their closest points - the witness points
    /// <c>ComputeMinimumDistance</c> reports - and when they touch exactly, both are the touching
    /// point. Unlike the witness points of <c>ComputeMinimumDistance</c>, which collapse inside
    /// the overlap, they are meaningful in every case.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="a"/> or <paramref name="b"/> is null.</exception>
    public static bool TryGetPenetration(ISupportMappable a, ISupportMappable b, out double depth,
        out Vector3 contactNormal, out Vector3 pointA, out Vector3 pointB)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), out depth, out contactNormal, out pointA, out pointB);
    }

    /// <summary>
    /// Determines whether two support-mapped convex shapes collide and, if so, by how much and
    /// where they touch, reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">The first shape.</param>
    /// <param name="b">The second shape.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the shapes are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first shape toward the second.</param>
    /// <param name="pointA">Contact point on the first shape, or its closest point when the shapes are apart.</param>
    /// <param name="pointB">Contact point on the second shape, or its closest point when the shapes are apart.</param>
    /// <returns><c>true</c> when the shapes overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3, out Vector3, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="a"/>, <paramref name="b"/> or <paramref name="workspace"/> is null.</exception>
    public static bool TryGetPenetration(ISupportMappable a, ISupportMappable b, GjkWorkspace workspace,
        out double depth, out Vector3 contactNormal, out Vector3 pointA, out Vector3 pointB)
    {
        ValidateShapes(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            workspace.LoadShapes(a, b, out GkBody bd1, out GkBody bd2);
            double info = ComputeContact(bd1, bd2, workspace);

            ReadContact(workspace, out contactNormal, out pointA, out pointB);
            return ToPenetration(info, out depth);
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Determines whether two convex polytopes collide and, if so, by how much.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first body toward the second.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)" path="/remarks"/>
    public static bool TryGetPenetration(Vector3[] a, Vector3[] b, out double depth, out Vector3 contactNormal)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), out depth, out contactNormal);
    }

    /// <summary>
    /// Determines whether two convex polytopes collide and, if so, by how much,
    /// reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first body toward the second.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static bool TryGetPenetration(Vector3[] a, Vector3[] b, GjkWorkspace workspace,
        out double depth, out Vector3 contactNormal)
    {
        double info = ComputeCollisionInformation(a, b, workspace, out contactNormal);

        return ToPenetration(info, out depth);
    }

    /// <summary>
    /// Determines whether two convex polytopes collide and, if so, by how much and where they touch.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first body toward the second.</param>
    /// <param name="pointA">Contact point on the first body, or its closest point when the bodies are apart.</param>
    /// <param name="pointB">Contact point on the second body, or its closest point when the bodies are apart.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3, out Vector3, out Vector3)" path="/remarks"/>
    public static bool TryGetPenetration(Vector3[] a, Vector3[] b, out double depth, out Vector3 contactNormal,
        out Vector3 pointA, out Vector3 pointB)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), out depth, out contactNormal, out pointA, out pointB);
    }

    /// <summary>
    /// Determines whether two convex polytopes collide and, if so, by how much and where they
    /// touch, reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">Vertices of the first body.</param>
    /// <param name="b">Vertices of the second body.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <param name="contactNormal">Contact normal, pointing from the first body toward the second.</param>
    /// <param name="pointA">Contact point on the first body, or its closest point when the bodies are apart.</param>
    /// <param name="pointB">Contact point on the second body, or its closest point when the bodies are apart.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3, out Vector3, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    public static bool TryGetPenetration(Vector3[] a, Vector3[] b, GjkWorkspace workspace,
        out double depth, out Vector3 contactNormal, out Vector3 pointA, out Vector3 pointB)
    {
        ValidateVertices(a, b);
        ValidateWorkspace(workspace);

        workspace.Occupy();
        try
        {
            LoadVertices(workspace, a, b, out GkBody bd1, out GkBody bd2);
            double info = ComputeContact(bd1, bd2, workspace);

            ReadContact(workspace, out contactNormal, out pointA, out pointB);
            return ToPenetration(info, out depth);
        }
        finally
        {
            workspace.Release();
        }
    }

    /// <summary>
    /// Determines whether two convex polytopes given in double precision collide and,
    /// if so, by how much.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="contactNormal">Receives the contact normal, pointing from the first body toward the second. Must hold at least 3 elements.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or <paramref name="contactNormal"/> is shorter than 3 elements.</exception>
    public static bool TryGetPenetration(ReadOnlySpan<double> a, ReadOnlySpan<double> b,
        Span<double> contactNormal, out double depth)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), contactNormal, out depth);
    }

    /// <summary>
    /// Determines whether two convex polytopes given in double precision collide and, if so,
    /// by how much, reusing the buffers of <paramref name="workspace"/> instead of allocating them.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="contactNormal">Receives the contact normal, pointing from the first body toward the second. Must hold at least 3 elements.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or <paramref name="contactNormal"/> is shorter than 3 elements.</exception>
    public static bool TryGetPenetration(ReadOnlySpan<double> a, ReadOnlySpan<double> b,
        GjkWorkspace workspace, Span<double> contactNormal, out double depth)
    {
        double info = ComputeCollisionInformation(a, b, workspace, contactNormal);

        return ToPenetration(info, out depth);
    }

    /// <summary>
    /// Determines whether two convex polytopes given in double precision collide and,
    /// if so, by how much and where they touch.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="contactNormal">Receives the contact normal, pointing from the first body toward the second. Must hold at least 3 elements.</param>
    /// <param name="pointA">Receives the contact point on the first body, or its closest point when the bodies are apart. Must hold at least 3 elements.</param>
    /// <param name="pointB">Receives the contact point on the second body, or its closest point when the bodies are apart. Must hold at least 3 elements.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3, out Vector3, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or an output span is shorter than 3 elements.</exception>
    public static bool TryGetPenetration(ReadOnlySpan<double> a, ReadOnlySpan<double> b,
        Span<double> contactNormal, Span<double> pointA, Span<double> pointB, out double depth)
    {
        return TryGetPenetration(a, b, new GjkWorkspace(), contactNormal, pointA, pointB, out depth);
    }

    /// <summary>
    /// Determines whether two convex polytopes given in double precision collide and, if so,
    /// by how much and where they touch, reusing the buffers of <paramref name="workspace"/>
    /// instead of allocating them.
    /// </summary>
    /// <param name="a">Coordinates of the first body as consecutive (x, y, z) triples.</param>
    /// <param name="b">Coordinates of the second body as consecutive (x, y, z) triples.</param>
    /// <param name="workspace">Scratch memory reused across queries. Must not be shared by concurrent queries.</param>
    /// <param name="contactNormal">Receives the contact normal, pointing from the first body toward the second. Must hold at least 3 elements.</param>
    /// <param name="pointA">Receives the contact point on the first body, or its closest point when the bodies are apart. Must hold at least 3 elements.</param>
    /// <param name="pointB">Receives the contact point on the second body, or its closest point when the bodies are apart. Must hold at least 3 elements.</param>
    /// <param name="depth">The penetration depth as a positive value, or 0 when the bodies are apart or touch exactly.</param>
    /// <returns><c>true</c> when the bodies overlap or touch exactly.</returns>
    /// <inheritdoc cref="TryGetPenetration(ISupportMappable, ISupportMappable, out double, out Vector3, out Vector3, out Vector3)" path="/remarks"/>
    /// <exception cref="ArgumentNullException"><paramref name="workspace"/> is null.</exception>
    /// <exception cref="ArgumentException">A span is empty, its length is not a multiple of 3, a coordinate is not finite,
    /// or an output span is shorter than 3 elements.</exception>
    public static bool TryGetPenetration(ReadOnlySpan<double> a, ReadOnlySpan<double> b, GjkWorkspace workspace,
        Span<double> contactNormal, Span<double> pointA, Span<double> pointB, out double depth)
    {
        ValidateCoordinates(a, b);
        ValidateWorkspace(workspace);
        ValidateOutput(contactNormal, nameof(contactNormal));
        ValidateOutput(pointA, nameof(pointA));
        ValidateOutput(pointB, nameof(pointB));

        workspace.Occupy();
        try
        {
            LoadCoordinates(workspace, a, b, out GkBody bd1, out GkBody bd2);
            double info = ComputeContact(bd1, bd2, workspace);

            double[] normal = workspace.ContactNormal;
            double[,] points = workspace.Simplex.Witnesses;
            for (int t = 0; t < 3; t++)
            {
                contactNormal[t] = normal[t];
                pointA[t] = points[0, t];
                pointB[t] = points[1, t];
            }

            return ToPenetration(info, out depth);
        }
        finally
        {
            workspace.Release();
        }
    }

    // Collision information with the witnesses turned into contact points, for the
    // TryGetPenetration overloads that report them.
    private static double ComputeContact(GkBody bd1, GkBody bd2, GjkWorkspace workspace)
    {
        GkSimplex simplex = workspace.Simplex;
        double distance = ComputeMinimumDistance(bd1, bd2, simplex);

        double[] normal = workspace.ContactNormal;
        EpaPolytope poly = workspace.Epa;
        ComputeCollisionInformation(bd1, bd2, simplex, poly, ref distance, normal);
        SetContactPoints(bd1, bd2, simplex, poly, distance, normal);

        return distance;
    }

    private static void ReadContact(GjkWorkspace workspace, out Vector3 contactNormal,
        out Vector3 pointA, out Vector3 pointB)
    {
        double[] normal = workspace.ContactNormal;
        double[,] points = workspace.Simplex.Witnesses;
        contactNormal = new Vector3((float)normal[0], (float)normal[1], (float)normal[2]);
        pointA = new Vector3((float)points[0, 0], (float)points[0, 1], (float)points[0, 2]);
        pointB = new Vector3((float)points[1, 0], (float)points[1, 1], (float)points[1, 2]);
    }

    // Collision information is signed - negative depth on collision, positive separation
    // distance otherwise - which is exactly the convention that is easy to get backwards.
    private static bool ToPenetration(double collisionInformation, out double depth)
    {
        if (collisionInformation > 0)
        {
            depth = 0;
            return false;
        }

        // Exactly touching reports 0, and negating it would hand back -0.
        depth = collisionInformation < 0 ? -collisionInformation : 0;
        return true;
    }

    private static void ValidateShapes(ISupportMappable a, ISupportMappable b)
    {
        if (a is null)
        {
            throw new ArgumentNullException(nameof(a));
        }

        if (b is null)
        {
            throw new ArgumentNullException(nameof(b));
        }
    }

    private static void ValidateVertices(Vector3[] a, Vector3[] b)
    {
        ValidateNotNullOrEmpty(a, b);

        for (int i = 0; i < a.Length; i++)
        {
            if (!IsFinite(a[i].X) || !IsFinite(a[i].Y) || !IsFinite(a[i].Z))
            {
                throw new ArgumentException($"Vertex at index {i} contains a non-finite coordinate: {a[i]}.", nameof(a));
            }
        }

        for (int i = 0; i < b.Length; i++)
        {
            if (!IsFinite(b[i].X) || !IsFinite(b[i].Y) || !IsFinite(b[i].Z))
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
            if (!IsFinite(a[i].X) || !IsFinite(a[i].Y))
            {
                throw new ArgumentException($"Vertex at index {i} contains a non-finite coordinate: {a[i]}.", nameof(a));
            }
        }

        for (int i = 0; i < b.Length; i++)
        {
            if (!IsFinite(b[i].X) || !IsFinite(b[i].Y))
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

    private static void ValidateWorkspace(GjkWorkspace workspace)
    {
        if (workspace is null)
        {
            throw new ArgumentNullException(nameof(workspace));
        }
    }

    private static void ValidatePrecision(double precision)
    {
        if (precision < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(precision), precision, "Precision must not be negative.");
        }
    }

    private static void ValidateCoordinates(ReadOnlySpan<double> a, ReadOnlySpan<double> b)
    {
        ValidateCoordinates(a, nameof(a));
        ValidateCoordinates(b, nameof(b));
    }

    private static void ValidateCoordinates(ReadOnlySpan<double> coordinates, string paramName)
    {
        if (coordinates.Length == 0)
        {
            throw new ArgumentException("At least one vertex is required.", paramName);
        }

        if (coordinates.Length % 3 != 0)
        {
            throw new ArgumentException(
                $"Length must be a multiple of 3 (consecutive x, y, z triples), but was {coordinates.Length}.",
                paramName);
        }

        for (int i = 0; i < coordinates.Length; i++)
        {
            if (double.IsNaN(coordinates[i]) || double.IsInfinity(coordinates[i]))
            {
                throw new ArgumentException(
                    $"Coordinate at index {i} is not finite: {coordinates[i]}.",
                    paramName);
            }
        }
    }

    private static void ValidateOutput(Span<double> destination, string paramName)
    {
        if (destination.Length < 3)
        {
            throw new ArgumentException(
                $"Must hold at least 3 elements (x, y, z), but holds {destination.Length}.",
                paramName);
        }
    }

    private static void LoadCoordinates(GjkWorkspace workspace, ReadOnlySpan<double> a, ReadOnlySpan<double> b,
        out GkBody bd1, out GkBody bd2)
    {
        workspace.PolytopeA.Load(a);
        workspace.PolytopeB.Load(b);
        workspace.LoadPolytopes(out bd1, out bd2);
    }

    private static void LoadVertices(GjkWorkspace workspace, Vector3[] a, Vector3[] b, out GkBody bd1, out GkBody bd2)
    {
        workspace.PolytopeA.Load(a);
        workspace.PolytopeB.Load(b);
        workspace.LoadPolytopes(out bd1, out bd2);
    }

    private static void LoadVertices(GjkWorkspace workspace, Vector2[] a, Vector2[] b, out GkBody bd1, out GkBody bd2)
    {
        workspace.PolytopeA.Load(a);
        workspace.PolytopeB.Load(b);
        workspace.LoadPolytopes(out bd1, out bd2);
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
    [Obsolete("Use ComputeMinimumDistance(ReadOnlySpan<double>, ReadOnlySpan<double>) for double precision input (coordinates as consecutive x, y, z triples instead of a [3, n] array), or ComputeMinimumDistance(Vector3[], Vector3[]). Add the GjkWorkspace overload to avoid per-query allocation.")]
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

        var bd1 = new GkPolytope
        {
            Coord = pinCoordsA,
            NumPoints = nCoordsA,
        };
        var bd2 = new GkPolytope
        {
            Coord = pinCoordsB,
            NumPoints = nCoordsB,
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
    internal static double ComputeMinimumDistance(GkPolytope bd1, GkPolytope bd2)
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
    internal static double ComputeMinimumDistance(GkPolytope bd1, GkPolytope bd2, GkSimplex s)
    {
        ValidatePolytope(bd1, nameof(bd1));
        ValidatePolytope(bd2, nameof(bd2));
        if (s is null)
        {
            throw new ArgumentNullException(nameof(s));
        }

        return ComputeMinimumDistance(new GkPolytopeBody(bd1), new GkPolytopeBody(bd2), s);
    }

    /// <summary>
    /// Runs GJK on two bodies. With a non-negative <paramref name="collisionPrecision"/> the loop
    /// stops as soon as <c>distance &lt;= collisionPrecision</c> is settled: it then returns 0 or
    /// <see cref="double.PositiveInfinity"/>, leaves the witnesses unset, and the result is only
    /// good for that comparison.
    /// </summary>
    internal static double ComputeMinimumDistance(GkBody bd1, GkBody bd2, GkSimplex s, double collisionPrecision = -1)
    {
        uint k = 0; // Iteration counter
        const int mk = 25; // Maximum number of GJK iterations
        const double eps_rel = _gkEpsilon * 1e4; // Tolerance on relative
        const double eps_tot = _gkEpsilon * 1e2; // Tolerance on absolute distance
        const double eps_rel2 = eps_rel * eps_rel;
        GkScratch scratch = s.Scratch;
        double[] w = scratch.W;
        int[] w_idx = scratch.WIdx;
        double[] v = scratch.V;
        double[] vMinus = scratch.VMinus;
        double norm2WMax = 0;
        double norm2Scale = 0;

        // Mapped shapes take their direction and return their support point as float, so
        // the early-exit bounds of HasCollision carry float-sized error for them.
        double boundTolerance = bd1 is GkMappedBody || bd2 is GkMappedBody ? 1e-5 : eps_rel;

        bd1.InitSupport();
        bd2.InitSupport();

        // Initialize search direction from the two interior points
        v[0] = bd1.S[0] - bd2.S[0];
        v[1] = bd1.S[1] - bd2.S[1];
        v[2] = bd1.S[2] - bd2.S[2];

        if (Norm2(v) < eps_rel2)
        {
            // The interior points coincide (concentric bodies), so they say nothing about
            // where to look. Any direction works as a start.
            v[0] = 1;
            v[1] = 0;
            v[2] = 0;
        }

        // Initialize the simplex with the first support point along that direction. The
        // interior points must not become simplex vertices themselves: EPA later builds its
        // polytope out of these vertices and needs boundary points, or it inherits a face
        // through the origin - which is exactly what an interior vertex produces - and its
        // closest-face search stalls at distance zero.
        vMinus[0] = -v[0];
        vMinus[1] = -v[1];
        vMinus[2] = -v[2];
        bd1.Support(vMinus);
        bd2.Support(v);

        s.NVrtx = 1;
        for (int t = 0; t < 3; ++t)
        {
            v[t] = bd1.S[t] - bd2.S[t];
            s.Vrtx[0, t] = v[t];
        }
        s.VrtxIdx[0, 0] = bd1.SIdx;
        s.VrtxIdx[0, 1] = bd2.SIdx;

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
            bd1.Support(vMinus);
            bd2.Support(v);
            for (int t = 0; t < 3; ++t)
            {
                w[t] = bd1.S[t] - bd2.S[t];
            }
            w_idx[0] = bd1.SIdx;
            w_idx[1] = bd2.SIdx;

            if (collisionPrecision >= 0)
            {
                // dot(v, w) / |v| is a lower bound on the distance and |v| an upper bound. The
                // margin keeps either verdict clear of rounding and of the zero clamp below, so
                // the answer matches the converged distance compared against the precision.
                // It scales with the coordinates, since rounding does, not with the distance.
                norm2Scale = Math.Max(Math.Max(norm2Scale, Norm2(v)), Math.Max(Norm2(w), Math.Max(Norm2(bd1.S), Norm2(bd2.S))));
                double normV = Math.Sqrt(Norm2(v));
                double margin = boundTolerance * Math.Sqrt(norm2Scale);
                if (DotProduct(v, w) > (collisionPrecision + margin) * normV)
                {
                    return double.PositiveInfinity;
                }

                if (normV + margin <= collisionPrecision)
                {
                    return 0;
                }
            }

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
                double testNorm = Norm2(s.Vrtx[jj, 0], s.Vrtx[jj, 1], s.Vrtx[jj, 2]);
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

        if (collisionPrecision < 0)
        {
            ComputeWitnesses(bd1, bd2, s);
        }

        // Overlapping bodies leave the search direction at the origin up to rounding, so the
        // raw result is noise (about 1e-16 at unit scale) rather than the documented 0. The
        // noise floor grows with the coordinates, so it is measured against the size of the
        // simplex - the same relative criterion the loop above terminates on.
        double norm2 = Norm2(v);
        for (int jj = 0; jj < s.NVrtx; jj++)
        {
            double vertexNorm = Norm2(s.Vrtx[jj, 0], s.Vrtx[jj, 1], s.Vrtx[jj, 2]);
            if (vertexNorm > norm2WMax)
            {
                norm2WMax = vertexNorm;
            }
        }

        if (norm2 <= eps_tot * eps_tot * norm2WMax || norm2 <= _gkEpsilon * _gkEpsilon)
        {
            return 0;
        }

        return Math.Sqrt(norm2);
    }

    private static double DotProduct(double[] a, double[] b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    private static double Norm2(double[] a)
    {
        return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    }

    private static double Norm2(double x, double y, double z)
    {
        return x * x + y * y + z * z;
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
        double[] s1p = CopyVrtx(s, 1, s.Scratch.L1);
        double[] s2p = CopyVrtx(s, 0, s.Scratch.L2);

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
        double pqX = p[0] - q[0], pqY = p[1] - q[1], pqZ = p[2] - q[2];
        double tmp = (p[0] * pqX + p[1] * pqY + p[2] * pqZ) / Norm2(pqX, pqY, pqZ);

        v[0] = p[0] - pqX * tmp;
        v[1] = p[1] - pqY * tmp;
        v[2] = p[2] - pqZ * tmp;
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
        double[] s1p = CopyVrtx(s, 2, s.Scratch.T1);
        double[] s2p = CopyVrtx(s, 1, s.Scratch.T2);
        double[] s3p = CopyVrtx(s, 0, s.Scratch.T3);

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
        double pqX = q[0] - p[0], pqY = q[1] - p[1], pqZ = q[2] - p[2];
        double prX = r[0] - p[0], prY = r[1] - p[1], prZ = r[2] - p[2];

        // nTmp = pq x pr, n = pq x nTmp
        double tX = pqY * prZ - pqZ * prY;
        double tY = pqZ * prX - pqX * prZ;
        double tZ = pqX * prY - pqY * prX;

        double nX = pqY * tZ - pqZ * tY;
        double nY = pqZ * tX - pqX * tZ;
        double nZ = pqX * tY - pqY * tX;

        return p[0] * nX + p[1] * nY + p[2] * nZ < 0;
    }

    private static void CrossProduct(double[] a, double[] b, double[] c)
    {
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
    }

    private static void ProjectOnPlane(double[] p, double[] q, double[] r, double[] v)
    {
        double pqX = p[0] - q[0], pqY = p[1] - q[1], pqZ = p[2] - q[2];
        double prX = p[0] - r[0], prY = p[1] - r[1], prZ = p[2] - r[2];

        double nX = pqY * prZ - pqZ * prY;
        double nY = pqZ * prX - pqX * prZ;
        double nZ = pqX * prY - pqY * prX;

        double tmp = (nX * p[0] + nY * p[1] + nZ * p[2]) / Norm2(nX, nY, nZ);

        v[0] = nX * tmp;
        v[1] = nY * tmp;
        v[2] = nZ * tmp;
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
        GkScratch scratch = s.Scratch;
        double[] s1 = scratch.S1;
        double[] s2 = scratch.S2;
        double[] s3 = scratch.S3;
        double[] s4 = scratch.S4;
        double[] s1s2 = scratch.E12;
        double[] s1s3 = scratch.E13;
        double[] s1s4 = scratch.E14;
        double[] si = scratch.Si;
        double[] sj = scratch.Sj;
        double[] sk = scratch.Sk;
        int[] s1Idx = scratch.S1Idx;
        int[] s2Idx = scratch.S2Idx;
        int[] s3Idx = scratch.S3Idx;
        int[] siIdx = scratch.SiIdx;
        int[] sjIdx = scratch.SjIdx;
        int[] skIdx = scratch.SkIdx;
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

        bool[] hff1Tests = scratch.Hff1Tests;
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

    // Copies a simplex vertex into a scratch buffer; the sub-algorithms need the values
    // to survive the simplex rewrites they perform.
    private static double[] CopyVrtx(GkSimplex s, int location, double[] destination)
    {
        destination[0] = s.Vrtx[location, 0];
        destination[1] = s.Vrtx[location, 1];
        destination[2] = s.Vrtx[location, 2];

        return destination;
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
        double nX = q[1] * r[2] - q[2] * r[1];
        double nY = q[2] * r[0] - q[0] * r[2];
        double nZ = q[0] * r[1] - q[1] * r[0];

        return p[0] * nX + p[1] * nY + p[2] * nZ <= 0; // discard s if true
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

    private static void ComputeWitnesses(GkBody bd1, GkBody bd2, GkSimplex smp)
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

    private static void W0D(GkBody bd1, GkBody bd2, GkSimplex smp)
    {
        var w00 = bd1.GetPoint(smp.VrtxIdx[0, 0]);
        var w01 = bd2.GetPoint(smp.VrtxIdx[0, 1]);

        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t];
            smp.Witnesses[1, t] = w01[t];
        }
    }

    private static void W1D(GkBody bd1, GkBody bd2, GkSimplex smp)
    {
        double pX = smp.Vrtx[0, 0], pY = smp.Vrtx[0, 1], pZ = smp.Vrtx[0, 2];
        double qX = smp.Vrtx[1, 0], qY = smp.Vrtx[1, 1], qZ = smp.Vrtx[1, 2];

        double pqX = qX - pX, pqY = qY - pY, pqZ = qZ - pZ;

        // Compute barycentric coordinates via matrix inversion
        // (in the linear case the matrix is 1x1 thus simplified)
        double det = Norm2(pqX, pqY, pqZ);
        if (det == 0.0)
        {
            // Degenerate case
            W0D(bd1, bd2, smp);
            return;
        }

        double a1 = (pqX * -pX + pqY * -pY + pqZ * -pZ) / det;
        double a0 = 1.0 - a1;

        // Compute witness points
        double[] w00 = bd1.GetPoint(smp.VrtxIdx[0, 0]);
        double[] w01 = bd2.GetPoint(smp.VrtxIdx[0, 1]);
        double[] w10 = bd1.GetPoint(smp.VrtxIdx[1, 0]);
        double[] w11 = bd2.GetPoint(smp.VrtxIdx[1, 1]);
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1;
        }
    }

    private static void W2D(GkBody bd1, GkBody bd2, GkSimplex smp)
    {
        double pX = smp.Vrtx[0, 0], pY = smp.Vrtx[0, 1], pZ = smp.Vrtx[0, 2];
        double qX = smp.Vrtx[1, 0], qY = smp.Vrtx[1, 1], qZ = smp.Vrtx[1, 2];
        double rX = smp.Vrtx[2, 0], rY = smp.Vrtx[2, 1], rZ = smp.Vrtx[2, 2];

        double pqX = qX - pX, pqY = qY - pY, pqZ = qZ - pZ;
        double prX = rX - pX, prY = rY - pY, prZ = rZ - pZ;

        double t00 = Norm2(pqX, pqY, pqZ);
        double t01 = pqX * prX + pqY * prY + pqZ * prZ;
        double t11 = Norm2(prX, prY, prZ);
        double det = t00 * t11 - t01 * t01;
        if (det == 0.0)
        {
            // Degenerate case
            W1D(bd1, bd2, smp);
            return;
        }

        double b0 = pqX * -pX + pqY * -pY + pqZ * -pZ;
        double b1 = prX * -pX + prY * -pY + prZ * -pZ;
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
        double[] w00 = bd1.GetPoint(smp.VrtxIdx[0, 0]);
        double[] w01 = bd2.GetPoint(smp.VrtxIdx[0, 1]);
        double[] w10 = bd1.GetPoint(smp.VrtxIdx[1, 0]);
        double[] w11 = bd2.GetPoint(smp.VrtxIdx[1, 1]);
        double[] w20 = bd1.GetPoint(smp.VrtxIdx[2, 0]);
        double[] w21 = bd2.GetPoint(smp.VrtxIdx[2, 1]);
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1 + w20[t] * a2;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1 + w21[t] * a2;
        }
    }

    private static void W3D(GkBody bd1, GkBody bd2, GkSimplex smp)
    {
        double pX = smp.Vrtx[0, 0], pY = smp.Vrtx[0, 1], pZ = smp.Vrtx[0, 2];
        double qX = smp.Vrtx[1, 0], qY = smp.Vrtx[1, 1], qZ = smp.Vrtx[1, 2];
        double rX = smp.Vrtx[2, 0], rY = smp.Vrtx[2, 1], rZ = smp.Vrtx[2, 2];
        double sX = smp.Vrtx[3, 0], sY = smp.Vrtx[3, 1], sZ = smp.Vrtx[3, 2];

        double pqX = qX - pX, pqY = qY - pY, pqZ = qZ - pZ;
        double prX = rX - pX, prY = rY - pY, prZ = rZ - pZ;
        double psX = sX - pX, psY = sY - pY, psZ = sZ - pZ;

        double t00 = Norm2(pqX, pqY, pqZ);
        double t01 = pqX * prX + pqY * prY + pqZ * prZ;
        double t02 = pqX * psX + pqY * psY + pqZ * psZ;
        double t11 = Norm2(prX, prY, prZ);
        double t12 = prX * psX + prY * psY + prZ * psZ;
        double t22 = Norm2(psX, psY, psZ);
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

        double b0 = pqX * -pX + pqY * -pY + pqZ * -pZ;
        double b1 = prX * -pX + prY * -pY + prZ * -pZ;
        double b2 = psX * -pX + psY * -pY + psZ * -pZ;

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
        double[] w00 = bd1.GetPoint(smp.VrtxIdx[0, 0]);
        double[] w01 = bd2.GetPoint(smp.VrtxIdx[0, 1]);
        double[] w10 = bd1.GetPoint(smp.VrtxIdx[1, 0]);
        double[] w11 = bd2.GetPoint(smp.VrtxIdx[1, 1]);
        double[] w20 = bd1.GetPoint(smp.VrtxIdx[2, 0]);
        double[] w21 = bd2.GetPoint(smp.VrtxIdx[2, 1]);
        double[] w30 = bd1.GetPoint(smp.VrtxIdx[3, 0]);
        double[] w31 = bd2.GetPoint(smp.VrtxIdx[3, 1]);
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

        double e0X = v1[0] - v0[0], e0Y = v1[1] - v0[1], e0Z = v1[2] - v0[2];
        double e1X = v2[0] - v0[0], e1Y = v2[1] - v0[1], e1Z = v2[2] - v0[2];

        face.Normal[0] = e0Y * e1Z - e0Z * e1Y;
        face.Normal[1] = e0Z * e1X - e0X * e1Z;
        face.Normal[2] = e0X * e1Y - e0Y * e1X;
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
            double side = face.Normal[0] * (centroid[0] - v0[0])
                + face.Normal[1] * (centroid[1] - v0[1])
                + face.Normal[2] * (centroid[2] - v0[2]);
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
        double distance = face.Normal[0] * (point[0] - v0[0])
            + face.Normal[1] * (point[1] - v0[1])
            + face.Normal[2] * (point[2] - v0[2]);

        return distance > _gkEpsilon;
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

            double e0X = v1[0] - v0[0], e0Y = v1[1] - v0[1], e0Z = v1[2] - v0[2];
            double e1X = v2[0] - v0[0], e1Y = v2[1] - v0[1], e1Z = v2[2] - v0[2];

            double nX = e0Y * e1Z - e0Z * e1Y;
            double nY = e0Z * e1X - e0X * e1Z;
            double nZ = e0X * e1Y - e0Y * e1X;

            // If normal points toward centroid need to flip the winding
            double towardCentroid = nX * (centroid[0] - v0[0])
                + nY * (centroid[1] - v0[1])
                + nZ * (centroid[2] - v0[2]);
            if (towardCentroid > 0)
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
        double e0X = v1[0] - v0[0], e0Y = v1[1] - v0[1], e0Z = v1[2] - v0[2];
        double e1X = v2[0] - v0[0], e1Y = v2[1] - v0[1], e1Z = v2[2] - v0[2];

        // Compute dot products for barycentric coords
        double d00 = e0X * e0X + e0Y * e0Y + e0Z * e0Z;
        double d01 = e0X * e1X + e0Y * e1Y + e0Z * e1Z;
        double d11 = e1X * e1X + e1Y * e1Y + e1Z * e1Z;
        double d20 = -(v0[0] * e0X + v0[1] * e0Y + v0[2] * e0Z);
        double d21 = -(v0[0] * e1X + v0[1] * e1Y + v0[2] * e1Z);

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
            double e12X = v2[0] - v1[0], e12Y = v2[1] - v1[1], e12Z = v2[2] - v1[2];
            double t = -(v1[0] * e12X + v1[1] * e12Y + v1[2] * e12Z)
                / (e12X * e12X + e12Y * e12Y + e12Z * e12Z);
            t = Clamp01(t);
            a0 = 0;
            a1 = 1.0 - t;
            a2 = t;
        }
        else if (u < 0)
        {
            // Origin projects outside edge v0-v2
            double t = -(v0[0] * e1X + v0[1] * e1Y + v0[2] * e1Z) / d11;
            t = Clamp01(t);
            a0 = 1.0 - t;
            a1 = 0;
            a2 = t;
        }
        else if (v < 0)
        {
            // Origin projects outside edge v0-v1
            double t = -(v0[0] * e0X + v0[1] * e0Y + v0[2] * e0Z) / d00;
            t = Clamp01(t);
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

    // float.IsFinite is unavailable on netstandard2.0
    private static bool IsFinite(float value)
    {
        return !float.IsNaN(value) && !float.IsInfinity(value);
    }

    // Math.Clamp is unavailable on netstandard2.0
    private static double Clamp01(double value)
    {
        if (value < 0.0)
        {
            return 0.0;
        }

        return value > 1.0 ? 1.0 : value;
    }

    // Support function for EPA: fresh (stateless) support of each body,
    // body2 searched in the opposite direction, combined into the Minkowski difference point.
    private static void SupportEpa(GkBody body1, GkBody body2, EpaPolytope poly,
        double[] direction, double[] result, int[] resultIdx)
    {
        double[] negDirection = poly.NegDirection;
        negDirection[0] = -direction[0];
        negDirection[1] = -direction[1];
        negDirection[2] = -direction[2];

        double[] p1 = poly.SupportA;
        double[] p2 = poly.SupportB;

        body1.SupportInto(direction, p1, out int idx1);
        body2.SupportInto(negDirection, p2, out int idx2);

        for (int i = 0; i < 3; i++)
        {
            result[i] = p1[i] - p2[i];
        }
        resultIdx[0] = idx1;
        resultIdx[1] = idx2;
    }

    private static void SetContactNormal(GkSimplex smp, double[] contactNormal)
    {
        double dX = smp.Witnesses[1, 0] - smp.Witnesses[0, 0];
        double dY = smp.Witnesses[1, 1] - smp.Witnesses[0, 1];
        double dZ = smp.Witnesses[1, 2] - smp.Witnesses[0, 2];

        double n = Math.Sqrt(Norm2(dX, dY, dZ));
        if (n > _gkEpsilon)
        {
            contactNormal[0] = dX / n;
            contactNormal[1] = dY / n;
            contactNormal[2] = dZ / n;
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
    internal static void ComputeCollisionInformation(GkPolytope bd1, GkPolytope bd2,
        GkSimplex simplex, ref double distance, double[] contactNormal)
    {
        ComputeCollisionInformation(new GkPolytopeBody(bd1), new GkPolytopeBody(bd2), simplex,
            new EpaPolytope(), ref distance, contactNormal);
    }

    internal static void ComputeCollisionInformation(GkBody bd1, GkBody bd2,
        GkSimplex simplex, EpaPolytope poly, ref double distance, double[] contactNormal)
    {
        const double epsSq = _gkEpsilon * _gkEpsilon;

        // Support point of the current step and the ids behind it. Both buffers live on the
        // polytope so a reused workspace runs the whole EPA pass without allocating.
        double[] newVertex = poly.NewVertex;
        int[] newVertexIdx = poly.NewVertexIdx;

        // if distance isn't 0 didn't detect collision - skip EPA
        if (distance > _gkEpsilon)
        {
            SetContactNormal(simplex, contactNormal);
            return;
        }

        // GJK can stop with fewer than four simplex vertices, or with four that enclose no
        // volume, and EPA needs a real tetrahedron. Growing the simplex with support queries
        // also establishes how many dimensions the Minkowski difference actually has: when it
        // cannot be grown any further, that is proof the difference is flat or thinner rather
        // than a failure to find a point.
        int dimension = PrepareEpaSimplex(bd1, bd2, simplex, poly);
        if (dimension < 4)
        {
            ComputeDegenerateCollisionInformation(bd1, bd2, simplex, poly, dimension,
                ref distance, contactNormal);
            return;
        }

        // On to actual EPA alg with a valid tetrahedron simplex
        // Initialize EPA polytope from simplex
        poly.Reset();
        double[] centroid = poly.Centroid;
        InitEpaPolytope(poly, simplex, centroid);

        // EPA iteration parameters
        const int maxIterations = 64;
        const double tolerance = _gkEpsilon * 1e2;
        int iteration = 0;

        // Horizon edges are collected into a pre-allocated buffer; this was allocated per
        // iteration before, which dominated the allocation cost of a colliding query.
        EpaEdge[] edges = poly.Edges;

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
            SupportEpa(bd1, bd2, poly, closest.Normal, newVertex, newVertexIdx);

            // Check termination condition: if distance to new vertex along normal
            // is not more than tolerance further than closest face
            double distToNew = DotProduct(closest.Normal, newVertex);
            double improvement = distToNew - closestDistance;

            // Every face normal tried yields a support plane touching both shapes, so its
            // distance is an upper bound on the penetration depth. Keep the tightest one: on
            // a curved shape the closest face of the final polytope is a poor direction
            // estimate, while the best plane seen is the closest the search ever came to the
            // true minimum translation.
            if (distToNew < poly.BestDepth)
            {
                RecordBestFace(poly, closest, distToNew);
            }

            if (improvement < tolerance)
            {
                // Converged
                Finish(bd1, bd2, poly, simplex, closestDistance, ref distance, contactNormal);
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
                Finish(bd1, bd2, poly, simplex, closestDistance, ref distance, contactNormal);
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
                    if (numEdges >= EpaPolytope.MaxEpaEdges)
                    {
                        break;
                    }

                    int va = e;
                    int vb = (e + 1) % 3;
                    edges[numEdges].Set(
                        poly.Faces[f].V[va],
                        poly.Faces[f].V[vb],
                        poly.Faces[f].VIdx[va],
                        poly.Faces[f].VIdx[vb]);
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

                double fe0X = fv1[0] - fv0[0], fe0Y = fv1[1] - fv0[1], fe0Z = fv1[2] - fv0[2];
                double fe1X = fv2[0] - fv0[0], fe1Y = fv2[1] - fv0[1], fe1Z = fv2[2] - fv0[2];

                double fnX = fe0Y * fe1Z - fe0Z * fe1Y;
                double fnY = fe0Z * fe1X - fe0X * fe1Z;
                double fnZ = fe0X * fe1Y - fe0Y * fe1X;

                // If normal points toward centroid flip winding
                double fTowardCentroid = fnX * (centroid[0] - fv0[0])
                    + fnY * (centroid[1] - fv0[1])
                    + fnZ * (centroid[2] - fv0[2]);
                if (fTowardCentroid > 0)
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
                // The polytope ran out of budget - this is where curved shapes end up, since
                // no finite polytope resolves them. The final face was never measured against
                // the shapes, so do that before picking the best plane.
                EpaFace closest = poly.Faces[closestFace];
                SupportEpa(bd1, bd2, poly, closest.Normal, newVertex, newVertexIdx);
                double distToNew = DotProduct(closest.Normal, newVertex);
                if (distToNew < poly.BestDepth)
                {
                    RecordBestFace(poly, closest, distToNew);
                }

                Finish(bd1, bd2, poly, simplex, closestDistance, ref distance, contactNormal);
            }
        }
    }

    // Remembers a face and the support-plane distance along its normal, so the result can be
    // built from it after later iterations have recycled its face slot.
    // Relative tolerances for the affine independence tests below. A new simplex vertex must
    // add a dimension by a wide enough margin that rounding cannot be mistaken for extent, and
    // by a small enough one that genuinely thin geometry still counts as three-dimensional:
    // a pair 1e-6 thick at unit scale is six orders of magnitude above these.
    private const double _sinTolerance = 1e-12;
    private const double _planeTolerance = 1e-12;

    // Makes the simplex an affinely independent set and grows it toward a tetrahedron with
    // support queries, returning the number of vertices reached. That number is the dimension
    // of the affine hull of the Minkowski difference: a stage that cannot be grown proves the
    // difference has no extent in the missing direction, since support queries in both
    // directions along it return points the simplex already has.
    //
    // GJK's own simplex cannot be trusted for this. Two unit squares overlapping in the z = 0
    // plane leave three *collinear* vertices on the x axis even though their difference is a
    // two-dimensional polygon, and a normal taken from those three points is meaningless.
    private static int PrepareEpaSimplex(GkBody bd1, GkBody bd2, GkSimplex simplex, EpaPolytope poly)
    {
        // Drop vertices that add no dimension (duplicates, collinear triples, coplanar quads).
        int count = 0;
        for (int vtx = 0; vtx < simplex.NVrtx; ++vtx)
        {
            if (!AddsDimension(simplex, count, simplex.Vrtx[vtx, 0], simplex.Vrtx[vtx, 1], simplex.Vrtx[vtx, 2]))
            {
                continue;
            }

            if (count != vtx)
            {
                for (int c = 0; c < 3; ++c)
                {
                    simplex.Vrtx[count, c] = simplex.Vrtx[vtx, c];
                }

                simplex.VrtxIdx[count, 0] = simplex.VrtxIdx[vtx, 0];
                simplex.VrtxIdx[count, 1] = simplex.VrtxIdx[vtx, 1];
            }

            count++;
        }

        simplex.NVrtx = count;

        while (count is > 0 and < 4 && TryGrowSimplex(bd1, bd2, simplex, poly, count))
        {
            count++;
            simplex.NVrtx = count;
        }

        return count;
    }

    // Fires support queries along directions that would leave the current simplex's affine hull
    // and appends the first point that does. Returns false when no direction produces one.
    private static bool TryGrowSimplex(GkBody bd1, GkBody bd2, GkSimplex simplex,
        EpaPolytope poly, int count)
    {
        double[] direction = poly.Direction;
        double[] newVertex = poly.NewVertex;
        int[] newVertexIdx = poly.NewVertexIdx;

        // Up to four candidate directions, all perpendicular to the current hull so that a
        // point beyond it is the one thing they can find.
        int candidates = SetupGrowthDirections(simplex, poly, count);

        for (int candidate = 0; candidate < candidates; ++candidate)
        {
            SelectGrowthDirection(poly, candidate, direction);
            if (Norm2(direction) <= _gkEpsilon * _gkEpsilon)
            {
                continue;
            }

            SupportEpa(bd1, bd2, poly, direction, newVertex, newVertexIdx);

            if (!AddsDimension(simplex, count, newVertex[0], newVertex[1], newVertex[2]))
            {
                continue;
            }

            for (int c = 0; c < 3; ++c)
            {
                simplex.Vrtx[count, c] = newVertex[c];
            }

            simplex.VrtxIdx[count, 0] = newVertexIdx[0];
            simplex.VrtxIdx[count, 1] = newVertexIdx[1];

            return true;
        }

        return false;
    }

    // Fills the growth direction pair (and its negatives) for the current simplex size and
    // returns how many candidates that makes.
    private static int SetupGrowthDirections(GkSimplex simplex, EpaPolytope poly, int count)
    {
        double[] first = poly.TempA;
        double[] second = poly.TempB;

        if (count == 1)
        {
            // Any direction can find a second point; the vertex itself points away from the
            // origin, which is where the difference has extent if it has any.
            first[0] = simplex.Vrtx[0, 0];
            first[1] = simplex.Vrtx[0, 1];
            first[2] = simplex.Vrtx[0, 2];

            if (Norm2(first) <= _gkEpsilon * _gkEpsilon)
            {
                first[0] = 1;
                first[1] = 0;
                first[2] = 0;
            }

            second[0] = 0;
            second[1] = 1;
            second[2] = 0;

            return 4;
        }

        double[] edge = poly.TempC;
        for (int c = 0; c < 3; ++c)
        {
            edge[c] = simplex.Vrtx[1, c] - simplex.Vrtx[0, c];
        }

        if (count == 2)
        {
            // Two independent perpendiculars of the edge: a point off the line has to have a
            // component along one of them.
            double edgeNorm = Math.Sqrt(Norm2(edge));
            double[] axis = poly.PlaneNormal;
            axis[0] = 1;
            axis[1] = 0;
            axis[2] = 0;
            if (edgeNorm > _gkEpsilon && Math.Abs(edge[0]) > 0.9 * edgeNorm)
            {
                axis[0] = 0;
                axis[1] = 1;
                axis[2] = 0;
            }

            CrossProduct(edge, axis, first);
            CrossProduct(edge, first, second);

            return 4;
        }

        // Three vertices: the only direction left is the triangle normal.
        double[] other = poly.PlaneNormal;
        for (int c = 0; c < 3; ++c)
        {
            other[c] = simplex.Vrtx[2, c] - simplex.Vrtx[0, c];
        }

        CrossProduct(edge, other, first);
        second[0] = 0;
        second[1] = 0;
        second[2] = 0;

        return 2;
    }

    // Candidates alternate between the two directions and their negatives.
    private static void SelectGrowthDirection(EpaPolytope poly, int candidate, double[] direction)
    {
        double[] source = candidate < 2 ? poly.TempA : poly.TempB;
        double sign = candidate % 2 == 0 ? 1 : -1;

        for (int c = 0; c < 3; ++c)
        {
            direction[c] = source[c] * sign;
        }
    }

    // Whether a point extends the affine hull of the first <paramref name="count"/> simplex
    // vertices, and so may be appended as vertex number <paramref name="count"/>.
    private static bool AddsDimension(GkSimplex simplex, int count, double x, double y, double z)
    {
        if (count == 0)
        {
            return true;
        }

        double e0X = x - simplex.Vrtx[0, 0];
        double e0Y = y - simplex.Vrtx[0, 1];
        double e0Z = z - simplex.Vrtx[0, 2];
        double e0Sq = Norm2(e0X, e0Y, e0Z);

        if (count == 1)
        {
            return e0Sq > _gkEpsilon * _gkEpsilon;
        }

        double aX = simplex.Vrtx[1, 0] - simplex.Vrtx[0, 0];
        double aY = simplex.Vrtx[1, 1] - simplex.Vrtx[0, 1];
        double aZ = simplex.Vrtx[1, 2] - simplex.Vrtx[0, 2];
        double aSq = Norm2(aX, aY, aZ);

        Cross(aX, aY, aZ, e0X, e0Y, e0Z, out double nX, out double nY, out double nZ);

        if (count == 2)
        {
            // |a x e| = |a| |e| sin(angle), so testing against the product of the lengths
            // tests the angle rather than an area, which would depend on the size of the bodies.
            return Norm2(nX, nY, nZ) > aSq * e0Sq * (_sinTolerance * _sinTolerance);
        }

        double bX = simplex.Vrtx[2, 0] - simplex.Vrtx[0, 0];
        double bY = simplex.Vrtx[2, 1] - simplex.Vrtx[0, 1];
        double bZ = simplex.Vrtx[2, 2] - simplex.Vrtx[0, 2];

        Cross(aX, aY, aZ, bX, bY, bZ, out double pX, out double pY, out double pZ);
        double planeNorm = Math.Sqrt(Norm2(pX, pY, pZ));
        if (planeNorm <= 0)
        {
            return false;
        }

        // Distance from the triangle's plane, relative to the size of the triangle.
        double height = Math.Abs((pX * e0X) + (pY * e0Y) + (pZ * e0Z)) / planeNorm;
        double scaleSq = aSq > Norm2(bX, bY, bZ) ? aSq : Norm2(bX, bY, bZ);

        return height > Math.Sqrt(scaleSq) * _planeTolerance;
    }

    private static void Cross(double aX, double aY, double aZ, double bX, double bY, double bZ,
        out double x, out double y, out double z)
    {
        x = (aY * bZ) - (aZ * bY);
        y = (aZ * bX) - (aX * bZ);
        z = (aX * bY) - (aY * bX);
    }

    // Reached when the Minkowski difference has no volume: EPA proved that support queries in
    // the missing direction return points the simplex already has, so the difference is a flat
    // polygon (dimension 3 vertices), a segment (2), or a single point (1).
    //
    // Such a difference has the origin on its boundary, so the penetration depth measured in 3D
    // is 0 - moving the bodies apart perpendicular to the plane separates them by any amount.
    // That is useless to a caller working with flat geometry, who wants the minimum translation
    // *within* the plane, so that is what is reported. Note the discontinuity this creates: a
    // pair with a small but non-zero thickness is not degenerate, keeps the 3D path, and gets
    // its (smaller) out-of-plane depth.
    private static void ComputeDegenerateCollisionInformation(GkBody bd1, GkBody bd2,
        GkSimplex simplex, EpaPolytope poly, int dimension, ref double distance, double[] contactNormal)
    {
        poly.BestDepth = double.MaxValue;

        if (dimension >= 3)
        {
            double[] e0 = poly.TempA;
            double[] e1 = poly.TempB;
            for (int c = 0; c < 3; ++c)
            {
                e0[c] = simplex.Vrtx[1, c] - simplex.Vrtx[0, c];
                e1[c] = simplex.Vrtx[2, c] - simplex.Vrtx[0, c];
            }

            CrossProduct(e0, e1, poly.PlaneNormal);
            ExpandPlanarPolygon(bd1, bd2, simplex, poly, ref distance, contactNormal);
            return;
        }

        if (dimension == 2)
        {
            ExpandAlongLine(bd1, bd2, simplex, poly, ref distance, contactNormal);
            return;
        }

        // A single point: the bodies touch at one point and 0 is the honest depth.
        for (int c = 0; c < 3; ++c)
        {
            simplex.Witnesses[0, c] = bd1.GetPoint(simplex.VrtxIdx[0, 0])[c];
            simplex.Witnesses[1, c] = bd2.GetPoint(simplex.VrtxIdx[0, 1])[c];
        }

        distance = 0;
        SetContactNormal(simplex, contactNormal);
    }

    // Depth of a pair whose Minkowski difference is a segment through the origin: the origin
    // leaves it at one end or the other, so the depth is the distance to the nearer end.
    private static void ExpandAlongLine(GkBody bd1, GkBody bd2, GkSimplex simplex,
        EpaPolytope poly, ref double distance, double[] contactNormal)
    {
        double[] direction = poly.Direction;

        // Any two distinct simplex vertices span the line.
        double lengthSq = 0;
        for (int vtx = 1; vtx < simplex.NVrtx; ++vtx)
        {
            for (int c = 0; c < 3; ++c)
            {
                direction[c] = simplex.Vrtx[vtx, c] - simplex.Vrtx[0, c];
            }

            lengthSq = Norm2(direction);
            if (lengthSq > _gkEpsilon * _gkEpsilon)
            {
                break;
            }
        }

        if (lengthSq <= _gkEpsilon * _gkEpsilon)
        {
            // A single point: the bodies touch and 0 is the honest answer.
            SetDegenerateResult(bd1, bd2, simplex, poly, 0, ref distance, contactNormal);
            return;
        }

        double length = Math.Sqrt(lengthSq);
        for (int c = 0; c < 3; ++c)
        {
            direction[c] /= length;
        }

        double[] support = poly.NewVertex;
        int[] supportIdx = poly.NewVertexIdx;
        double[] negative = poly.TempC;

        SupportEpa(bd1, bd2, poly, direction, support, supportIdx);
        double forward = DotProduct(direction, support);
        int forwardA = supportIdx[0];
        int forwardB = supportIdx[1];

        for (int c = 0; c < 3; ++c)
        {
            negative[c] = -direction[c];
        }

        SupportEpa(bd1, bd2, poly, negative, support, supportIdx);
        double backward = DotProduct(negative, support);

        // The origin lies inside the segment, so both extents are non-negative up to rounding.
        double depth = forward < backward ? forward : backward;
        if (depth < 0)
        {
            depth = 0;
        }

        if (forward >= backward)
        {
            for (int c = 0; c < 3; ++c)
            {
                direction[c] = negative[c];
            }
        }
        else
        {
            supportIdx[0] = forwardA;
            supportIdx[1] = forwardB;
        }

        for (int c = 0; c < 3; ++c)
        {
            contactNormal[c] = direction[c];
        }

        SetDegenerateResult(bd1, bd2, simplex, poly, depth, ref distance, contactNormal);
    }

    // Depth of a pair whose Minkowski difference is a flat polygon: the two-dimensional EPA.
    // The polygon starts as the triangle EPA could not grow and is expanded toward the origin's
    // nearest edge until no support point beyond that edge exists.
    private static void ExpandPlanarPolygon(GkBody bd1, GkBody bd2, GkSimplex simplex,
        EpaPolytope poly, ref double distance, double[] contactNormal)
    {
        const double tolerance = _gkEpsilon * 1e2;
        const double epsSq = _gkEpsilon * _gkEpsilon;
        const int maxIterations = 64;

        double[] planeNormal = poly.PlaneNormal;
        double planeLength = Math.Sqrt(Norm2(planeNormal));
        for (int c = 0; c < 3; ++c)
        {
            planeNormal[c] /= planeLength;
        }

        // The triangle is wound counter-clockwise about the plane normal by construction
        // (the normal came from its own edge cross product), so edge x normal points outward.
        double[][] vertices = poly.PlanarVertices;
        int[][] indices = poly.PlanarVertexIndices;
        for (int vtx = 0; vtx < 3; ++vtx)
        {
            for (int c = 0; c < 3; ++c)
            {
                vertices[vtx][c] = simplex.Vrtx[vtx, c];
            }

            indices[vtx][0] = simplex.VrtxIdx[vtx, 0];
            indices[vtx][1] = simplex.VrtxIdx[vtx, 1];
        }

        poly.NumPlanarVertices = 3;

        double[] support = poly.NewVertex;
        int[] supportIdx = poly.NewVertexIdx;
        double[] edge = poly.TempA;
        double[] outward = poly.TempB;
        double edgeDistance = 0;

        for (int iteration = 0; iteration < maxIterations; ++iteration)
        {
            int count = poly.NumPlanarVertices;

            // Find the edge whose supporting line is closest to the origin.
            int closest = -1;
            double closestDistance = double.MaxValue;
            for (int i = 0; i < count; ++i)
            {
                int j = i + 1 == count ? 0 : i + 1;
                if (!EdgeOutwardNormal(vertices[i], vertices[j], planeNormal, edge, outward))
                {
                    continue;
                }

                // The origin is inside the polygon, so this is non-negative up to rounding.
                double signedDistance = DotProduct(outward, vertices[i]);
                if (signedDistance < 0)
                {
                    signedDistance = 0;
                }

                if (signedDistance < closestDistance)
                {
                    closestDistance = signedDistance;
                    closest = i;
                }
            }

            if (closest < 0)
            {
                break;
            }

            int next = closest + 1 == count ? 0 : closest + 1;
            EdgeOutwardNormal(vertices[closest], vertices[next], planeNormal, edge, outward);
            edgeDistance = closestDistance;

            SupportEpa(bd1, bd2, poly, outward, support, supportIdx);
            double supportDistance = DotProduct(outward, support);

            // Every outward normal tried gives a support line touching both bodies, so its
            // distance bounds the depth from above; keep the tightest one, exactly as the
            // three-dimensional pass does for support planes.
            if (supportDistance < poly.BestDepth)
            {
                RecordBestEdge(poly, vertices[closest], vertices[next],
                    indices[closest], indices[next], outward, supportDistance);
            }

            if (supportDistance - edgeDistance < tolerance)
            {
                break;
            }

            bool isDuplicate = false;
            for (int i = 0; i < count; ++i)
            {
                double dx = support[0] - vertices[i][0];
                double dy = support[1] - vertices[i][1];
                double dz = support[2] - vertices[i][2];
                if (dx * dx + dy * dy + dz * dz < epsSq)
                {
                    isDuplicate = true;
                    break;
                }
            }

            if (isDuplicate || count == EpaPolytope.MaxPlanarVertices)
            {
                break;
            }

            // Split the closest edge with the new support point, keeping the winding.
            for (int i = count; i > next; --i)
            {
                double[] source = vertices[i - 1];
                double[] target = vertices[i];
                target[0] = source[0];
                target[1] = source[1];
                target[2] = source[2];
                indices[i][0] = indices[i - 1][0];
                indices[i][1] = indices[i - 1][1];
            }

            for (int c = 0; c < 3; ++c)
            {
                vertices[next][c] = support[c];
            }

            indices[next][0] = supportIdx[0];
            indices[next][1] = supportIdx[1];
            poly.NumPlanarVertices = count + 1;
        }

        // Same bound arithmetic as the three-dimensional pass: the edge distance is a lower
        // bound on the depth, the tightest support line an upper one.
        double depth = poly.BestDepth > edgeDistance ? poly.BestDepth : edgeDistance;
        if (depth < 0 || double.IsInfinity(depth) || double.IsNaN(depth))
        {
            depth = 0;
        }

        for (int c = 0; c < 3; ++c)
        {
            contactNormal[c] = poly.PlanarBestNormal[c];
        }

        SetPlanarWitnesses(bd1, bd2, simplex, poly);
        distance = depth > 0 ? -depth : 0;
    }

    // In-plane outward normal of the edge from <paramref name="from"/> to <paramref name="to"/>,
    // for a polygon wound counter-clockwise about <paramref name="planeNormal"/>.
    // Returns false for an edge too short to have a direction.
    private static bool EdgeOutwardNormal(double[] from, double[] to, double[] planeNormal,
        double[] edge, double[] outward)
    {
        for (int c = 0; c < 3; ++c)
        {
            edge[c] = to[c] - from[c];
        }

        CrossProduct(edge, planeNormal, outward);
        double length = Math.Sqrt(Norm2(outward));
        if (length <= _gkEpsilon)
        {
            return false;
        }

        for (int c = 0; c < 3; ++c)
        {
            outward[c] /= length;
        }

        return true;
    }

    private static void RecordBestEdge(EpaPolytope poly, double[] from, double[] to,
        int[] fromIndices, int[] toIndices, double[] outward, double supportDistance)
    {
        poly.BestDepth = supportDistance;

        for (int c = 0; c < 3; ++c)
        {
            poly.PlanarBestEdge[0][c] = from[c];
            poly.PlanarBestEdge[1][c] = to[c];
            poly.PlanarBestNormal[c] = outward[c];
        }

        poly.PlanarBestEdgeIndices[0][0] = fromIndices[0];
        poly.PlanarBestEdgeIndices[0][1] = fromIndices[1];
        poly.PlanarBestEdgeIndices[1][0] = toIndices[0];
        poly.PlanarBestEdgeIndices[1][1] = toIndices[1];
    }

    // Witness points from the closest edge, weighted by where the origin projects onto it.
    // Both collapse onto a common point inside the overlap, as they do for the 3D pass.
    private static void SetPlanarWitnesses(GkBody bd1, GkBody bd2, GkSimplex simplex, EpaPolytope poly)
    {
        double[] from = poly.PlanarBestEdge[0];
        double[] to = poly.PlanarBestEdge[1];

        double edgeX = to[0] - from[0];
        double edgeY = to[1] - from[1];
        double edgeZ = to[2] - from[2];
        double lengthSq = Norm2(edgeX, edgeY, edgeZ);

        double t = 0;
        if (lengthSq > _gkEpsilon * _gkEpsilon)
        {
            t = -(from[0] * edgeX + from[1] * edgeY + from[2] * edgeZ) / lengthSq;
            t = t < 0 ? 0 : t > 1 ? 1 : t;
        }

        double w0 = 1 - t;
        int[] fromIndices = poly.PlanarBestEdgeIndices[0];
        int[] toIndices = poly.PlanarBestEdgeIndices[1];

        for (int c = 0; c < 3; ++c)
        {
            simplex.Witnesses[0, c] = bd1.GetPoint(fromIndices[0])[c] * w0
                + bd1.GetPoint(toIndices[0])[c] * t;
            simplex.Witnesses[1, c] = bd2.GetPoint(fromIndices[1])[c] * w0
                + bd2.GetPoint(toIndices[1])[c] * t;
        }
    }

    private static void SetDegenerateResult(GkBody bd1, GkBody bd2, GkSimplex simplex,
        EpaPolytope poly, double depth, ref double distance, double[] contactNormal)
    {
        int[] supportIdx = poly.NewVertexIdx;
        for (int c = 0; c < 3; ++c)
        {
            simplex.Witnesses[0, c] = bd1.GetPoint(supportIdx[0])[c];
            simplex.Witnesses[1, c] = bd2.GetPoint(supportIdx[1])[c];
        }

        if (depth <= 0)
        {
            distance = 0;
            SetContactNormal(simplex, contactNormal);
            return;
        }

        distance = -depth;
    }

    private static void RecordBestFace(EpaPolytope poly, EpaFace face, double supportDistance)
    {
        poly.BestDepth = supportDistance;

        for (int i = 0; i < 3; i++)
        {
            poly.BestNormal[i] = face.Normal[i];
            poly.BestVertices[i] = face.V[i];
            poly.BestVertexIndices[i][0] = face.VIdx[i][0];
            poly.BestVertexIndices[i][1] = face.VIdx[i][1];
        }
    }

    // Builds the EPA result from the best face seen.
    //
    // On depth: a face of the polytope is inscribed in the Minkowski difference, so its
    // distance from the origin is a lower bound on the penetration depth - exact for a
    // polyhedron, but for a curved shape short by the chord of the face (up to 5.6% for a
    // sphere, which no iteration budget removes). The support plane along a face normal is
    // the matching upper bound, equal to the depth when the normal is the true minimum
    // translation direction and larger otherwise, so the smallest one seen is the best
    // estimate available. Taking the larger of the two bounds keeps the result at least as
    // deep as the polytope proves it to be.
    private static void Finish(GkBody bd1, GkBody bd2, EpaPolytope poly, GkSimplex simplex,
        double faceDistance, ref double distance, double[] contactNormal)
    {
        ComputeBarycentricOrigin(poly.Vertices[poly.BestVertices[0]],
            poly.Vertices[poly.BestVertices[1]],
            poly.Vertices[poly.BestVertices[2]], out double a0, out double a1, out double a2);

        for (int i = 0; i < 3; i++)
        {
            simplex.Witnesses[0, i] = bd1.GetPoint(poly.BestVertexIndices[0][0])[i] * a0
                + bd1.GetPoint(poly.BestVertexIndices[1][0])[i] * a1
                + bd1.GetPoint(poly.BestVertexIndices[2][0])[i] * a2;
            simplex.Witnesses[1, i] = bd2.GetPoint(poly.BestVertexIndices[0][1])[i] * a0
                + bd2.GetPoint(poly.BestVertexIndices[1][1])[i] * a1
                + bd2.GetPoint(poly.BestVertexIndices[2][1])[i] * a2;
            contactNormal[i] = poly.BestNormal[i];
        }

        distance = -(poly.BestDepth > faceDistance ? poly.BestDepth : faceDistance);
    }

    // Turns the witnesses of a penetrating pair into contact points: the deepest point of each
    // body along the contact normal. The passes above leave the barycentric image of the
    // origin's projection onto the best face (or edge), which lies on a body's support plane
    // only where that face is flush with the body - true for polyhedra, but on a curved shape
    // the face is a chord and its image sits inside. A witness off the plane by more than
    // rounding is therefore replaced with the support point itself, which is on the surface by
    // definition; one within rounding is moved onto the plane, keeping its place in the contact
    // region. Either way dot(pointA - pointB, normal) becomes the support distance along the
    // normal, which is the depth every pass reports.
    private static void SetContactPoints(GkBody bd1, GkBody bd2, GkSimplex simplex, EpaPolytope poly,
        double distance, double[] contactNormal)
    {
        // Apart or exactly touching: the witnesses already lie on both surfaces.
        if (!(distance < 0))
        {
            return;
        }

        SupportEpa(bd1, bd2, poly, contactNormal, poly.NewVertex, poly.NewVertexIdx);
        MoveOntoSupportPlane(simplex, 0, poly.SupportA, contactNormal);
        MoveOntoSupportPlane(simplex, 1, poly.SupportB, contactNormal);
    }

    private static void MoveOntoSupportPlane(GkSimplex simplex, int body, double[] support, double[] normal)
    {
        double offset = 0;
        double scale = 1;
        for (int c = 0; c < 3; ++c)
        {
            offset += (support[c] - simplex.Witnesses[body, c]) * normal[c];
            scale += Math.Abs(support[c]);
        }

        // Single-precision input leaves the vertices of a rotated face off its plane by about
        // 1e-7 relative, so this tolerance has to sit above that to recognize a flush face.
        bool flush = Math.Abs(offset) <= scale * _precision;
        for (int c = 0; c < 3; ++c)
        {
            simplex.Witnesses[body, c] = flush
                ? simplex.Witnesses[body, c] + offset * normal[c]
                : support[c];
        }
    }

    #endregion
}

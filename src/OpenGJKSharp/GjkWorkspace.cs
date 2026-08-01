using OpenGJKSharp.Models;

namespace OpenGJKSharp;

/// <summary>
/// Reusable scratch memory for GJK and EPA queries.
/// </summary>
/// <remarks>
/// <para>
/// Every <see cref="OpenGjk"/> query needs working buffers (simplex, support point store,
/// EPA polytope). The overloads that do not take a workspace allocate those buffers per call,
/// which is fine for occasional queries but shows up as steady garbage in a physics loop.
/// Passing the same workspace to repeated queries reuses the buffers instead: after the first
/// few calls have grown them, further queries of the same kind allocate nothing.
/// </para>
/// <para>
/// A workspace holds no result state - everything a query produces is returned or written to
/// its <c>out</c> parameters - so it can be reused immediately after a query returns.
/// The EPA buffers are created lazily on the first colliding query.
/// </para>
/// <para>
/// Instances are <b>not thread-safe</b> and must not be shared by concurrent queries.
/// Use one workspace per thread (for example a <c>[ThreadStatic]</c> field or one per worker).
/// </para>
/// <example>
/// <code>
/// var workspace = new GjkWorkspace();
/// foreach (var (a, b) in pairs)
/// {
///     double distance = OpenGjk.ComputeMinimumDistance(a, b, workspace);
/// }
/// </code>
/// </example>
/// </remarks>
public sealed class GjkWorkspace
{
    private EpaPolytope? _epa;

    /// <summary>
    /// Creates an empty workspace. Buffers are allocated on first use and grow as needed.
    /// </summary>
    public GjkWorkspace()
    {
    }

    internal GkSimplex Simplex { get; } = new();

    internal GkMappedBody MappedBodyA { get; } = new();

    internal GkMappedBody MappedBodyB { get; } = new();

    internal GkPolytope PolytopeA { get; } = new();

    internal GkPolytope PolytopeB { get; } = new();

    internal GkPolytopeBody PolytopeBodyA { get; } = new();

    internal GkPolytopeBody PolytopeBodyB { get; } = new();

    /// <summary>Buffer receiving the contact normal before it is converted to a vector.</summary>
    internal double[] ContactNormal { get; } = new double[3];

    /// <summary>EPA buffers, created on the first colliding query.</summary>
    internal EpaPolytope Epa => _epa ??= new EpaPolytope();

    /// <summary>Binds both vertex-cloud bodies to the workspace polytopes and returns them.</summary>
    internal void LoadPolytopes(out GkBody bodyA, out GkBody bodyB)
    {
        PolytopeBodyA.Reset(PolytopeA);
        PolytopeBodyB.Reset(PolytopeB);
        bodyA = PolytopeBodyA;
        bodyB = PolytopeBodyB;
    }

    /// <summary>Binds both support-mapped bodies to the given shapes and returns them.</summary>
    internal void LoadShapes(ISupportMappable a, ISupportMappable b, out GkBody bodyA, out GkBody bodyB)
    {
        MappedBodyA.Reset(a);
        MappedBodyB.Reset(b);
        bodyA = MappedBodyA;
        bodyB = MappedBodyB;
    }
}

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
/// A query that finds the workspace already occupied throws <see cref="InvalidOperationException"/>
/// rather than reading and writing the buffers underneath the query that owns them, so sharing one
/// by mistake fails loudly instead of returning wrong results. The same check catches re-entering a
/// query with the workspace it was called from, for example from a custom
/// <see cref="ISupportMappable.GetSupportPoint"/>.
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
    private int _inUse;

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

    /// <summary>
    /// Takes exclusive ownership of the buffers for the duration of one query.
    /// Every query that touches the workspace must pair this with <see cref="Release"/>
    /// in a <c>finally</c> block, so a failing query does not leave it occupied.
    /// </summary>
    /// <exception cref="InvalidOperationException">Another query already owns the workspace.</exception>
    internal void Occupy()
    {
        if (Interlocked.CompareExchange(ref _inUse, 1, 0) != 0)
        {
            throw new InvalidOperationException(
                "This GjkWorkspace is already in use by another query. A GjkWorkspace is not " +
                "thread-safe: use one per thread instead of sharing one across threads, and do not " +
                "start a query with the workspace of a query that is still running (for example " +
                "from inside a custom GetSupportPoint).");
        }
    }

    /// <summary>Releases the ownership taken by <see cref="Occupy"/>.</summary>
    internal void Release()
    {
        Volatile.Write(ref _inUse, 0);
    }

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

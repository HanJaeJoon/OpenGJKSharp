namespace OpenGJKSharp.Models;

/// <summary>
/// Internal abstraction over the two body kinds the GJK/EPA loops operate on:
/// vertex clouds (<see cref="GkPolytope"/>) and support-mapped shapes (<see cref="ISupportMappable"/>).
/// Support points are identified by an integer id so the index-based simplex/EPA
/// bookkeeping works for both kinds; <see cref="GetPoint"/> resolves ids back to
/// coordinates for witness reconstruction.
/// </summary>
internal abstract class GkBody
{
    /// <summary>Current support point, seeded by <see cref="InitSupport"/> and updated by <see cref="Support"/>.</summary>
    public abstract double[] S { get; }

    /// <summary>Id of the current support point, resolvable via <see cref="GetPoint"/>.</summary>
    public abstract int SIdx { get; }

    /// <summary>Seeds <see cref="S"/>/<see cref="SIdx"/> with the initial point (first vertex or shape center).</summary>
    public abstract void InitSupport();

    /// <summary>Stateful GJK support: moves <see cref="S"/>/<see cref="SIdx"/> to the furthest point along <paramref name="v"/>.</summary>
    public abstract void Support(double[] v);

    /// <summary>Stateless support used by EPA: writes the furthest point along <paramref name="direction"/> into <paramref name="result"/>.</summary>
    public abstract void SupportInto(double[] direction, double[] result, out int idx);

    /// <summary>Resolves a support point id to its coordinates.</summary>
    public abstract double[] GetPoint(int idx);
}

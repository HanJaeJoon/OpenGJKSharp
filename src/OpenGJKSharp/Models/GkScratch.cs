namespace OpenGJKSharp.Models;

/// <summary>
/// Pre-allocated working buffers for the GJK main loop and its distance sub-algorithm.
/// The sub-algorithms copy simplex vertices into local arrays before rewriting the simplex;
/// those copies live here so a reused <see cref="GkSimplex"/> makes repeated queries
/// allocation-free. Each set of buffers belongs to exactly one routine, so routines that
/// call each other (S3D -> S2D) never alias.
/// </summary>
internal sealed class GkScratch
{
    /// <summary>Search direction of the GJK main loop.</summary>
    public double[] V { get; } = new double[3];

    /// <summary>Negated search direction of the GJK main loop.</summary>
    public double[] VMinus { get; } = new double[3];

    /// <summary>Latest Minkowski difference support point of the GJK main loop.</summary>
    public double[] W { get; } = new double[3];

    /// <summary>Body point ids behind <see cref="W"/>.</summary>
    public int[] WIdx { get; } = new int[2];

    /// <summary>Vertex copies used by S1D.</summary>
    public double[] L1 { get; } = new double[3];

    /// <inheritdoc cref="L1"/>
    public double[] L2 { get; } = new double[3];

    /// <summary>Vertex copies used by S2D.</summary>
    public double[] T1 { get; } = new double[3];

    /// <inheritdoc cref="T1"/>
    public double[] T2 { get; } = new double[3];

    /// <inheritdoc cref="T1"/>
    public double[] T3 { get; } = new double[3];

    /// <summary>Vertex copies used by S3D.</summary>
    public double[] S1 { get; } = new double[3];

    /// <inheritdoc cref="S1"/>
    public double[] S2 { get; } = new double[3];

    /// <inheritdoc cref="S1"/>
    public double[] S3 { get; } = new double[3];

    /// <inheritdoc cref="S1"/>
    public double[] S4 { get; } = new double[3];

    /// <inheritdoc cref="S1"/>
    public double[] Si { get; } = new double[3];

    /// <inheritdoc cref="S1"/>
    public double[] Sj { get; } = new double[3];

    /// <inheritdoc cref="S1"/>
    public double[] Sk { get; } = new double[3];

    /// <summary>Edge vectors used by S3D.</summary>
    public double[] E12 { get; } = new double[3];

    /// <inheritdoc cref="E12"/>
    public double[] E13 { get; } = new double[3];

    /// <inheritdoc cref="E12"/>
    public double[] E14 { get; } = new double[3];

    /// <summary>Body point ids used by S3D.</summary>
    public int[] S1Idx { get; } = new int[2];

    /// <inheritdoc cref="S1Idx"/>
    public int[] S2Idx { get; } = new int[2];

    /// <inheritdoc cref="S1Idx"/>
    public int[] S3Idx { get; } = new int[2];

    /// <inheritdoc cref="S1Idx"/>
    public int[] SiIdx { get; } = new int[2];

    /// <inheritdoc cref="S1Idx"/>
    public int[] SjIdx { get; } = new int[2];

    /// <inheritdoc cref="S1Idx"/>
    public int[] SkIdx { get; } = new int[2];

    /// <summary>Results of <see cref="OpenGjk"/>'s Hff1 tests, indexed by simplex vertex.</summary>
    public bool[] Hff1Tests { get; } = new bool[3];
}

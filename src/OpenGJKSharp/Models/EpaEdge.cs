namespace OpenGJKSharp.Models;

/// <summary>
/// Structure for horizon edge collection.
/// Instances are pre-allocated by <see cref="EpaPolytope"/> and refilled per iteration.
/// </summary>
internal sealed class EpaEdge
{
    /// <summary>
    /// First vertex index in polytope.
    /// </summary>
    public int V1 { get; set; }

    /// <summary>
    /// Second vertex index in polytope.
    /// </summary>
    public int V2 { get; set; }

    /// <summary>
    /// Original vertex indices of <see cref="V1"/> for witness computation.
    /// </summary>
    public int[] VIdx1 { get; } = new int[2];

    /// <summary>
    /// Original vertex indices of <see cref="V2"/> for witness computation.
    /// </summary>
    public int[] VIdx2 { get; } = new int[2];

    public bool Valid { get; set; }

    /// <summary>Refills this edge in place, keeping the index arrays.</summary>
    public void Set(int v1, int v2, int[] vIdx1, int[] vIdx2)
    {
        V1 = v1;
        V2 = v2;
        VIdx1[0] = vIdx1[0];
        VIdx1[1] = vIdx1[1];
        VIdx2[0] = vIdx2[0];
        VIdx2[1] = vIdx2[1];
        Valid = true;
    }
}

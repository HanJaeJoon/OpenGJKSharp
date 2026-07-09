namespace OpenGJKSharp.Models;

/// <summary>
/// Structure for horizon edge collection.
/// </summary>
internal class EpaEdge
{
    /// <summary>
    /// First vertex index in polytope.
    /// </summary>
    public int V1 { get; init; }

    /// <summary>
    /// Second vertex index in polytope.
    /// </summary>
    public int V2 { get; init; }

    /// <summary>
    /// Original vertex indices of <see cref="V1"/> for witness computation.
    /// </summary>
    public int[] VIdx1 { get; init; } = new int[2];

    /// <summary>
    /// Original vertex indices of <see cref="V2"/> for witness computation.
    /// </summary>
    public int[] VIdx2 { get; init; } = new int[2];

    public bool Valid { get; set; }
}

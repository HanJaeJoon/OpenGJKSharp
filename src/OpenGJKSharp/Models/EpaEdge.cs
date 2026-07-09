namespace OpenGJKSharp.Models;

/// <summary>
/// Structure for horizon edge collection.
/// </summary>
internal class EpaEdge
{
    /// <summary>
    /// Vertex indices in polytope.
    /// </summary>
    public int V1, V2;

    /// <summary>
    /// Original vertex indices for witness computation.
    /// </summary>
    public int[] VIdx1 = new int[2];

    /// <summary>
    /// Original vertex indices for witness computation.
    /// </summary>
    public int[] VIdx2 = new int[2];

    public bool Valid;
}

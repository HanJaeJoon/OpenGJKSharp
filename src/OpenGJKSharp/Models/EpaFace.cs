namespace OpenGJKSharp.Models;

/// <summary>
/// Face structure for EPA polytope.
/// Each face is a triangle with 3 vertex indices.
/// </summary>
internal class EpaFace
{
    /// <summary>
    /// Vertex indices in the polytope.
    /// </summary>
    public int[] V { get; } = new int[3];

    /// <summary>
    /// Original vertex indices from original polytopes for witness computation [vertex][body].
    /// </summary>
    public int[][] VIdx { get; } = [new int[2], new int[2], new int[2]];

    /// <summary>
    /// Face normal (pointing outward from origin).
    /// </summary>
    public double[] Normal { get; } = new double[3];

    /// <summary>
    /// Distance from origin to face plane.
    /// </summary>
    public double Distance { get; set; }

    /// <summary>
    /// Whether this face is still valid (not removed).
    /// </summary>
    public bool Valid { get; set; }
}

namespace OpenGJKSharp.Models;

/// <summary>
/// Polytope structure for EPA.
/// </summary>
internal class EpaPolytope
{
    public const int MaxEpaFaces = 128;

    public const int MaxEpaVertices = MaxEpaFaces + 4;

    /// <summary>
    /// Vertex coordinates in the Minkowski difference.
    /// </summary>
    public double[][] Vertices { get; } = new double[MaxEpaVertices][];

    /// <summary>
    /// Original vertex indices [vertex][body].
    /// </summary>
    public int[][] VertexIndices { get; } = new int[MaxEpaVertices][];

    public int NumVertices { get; set; }

    public EpaFace[] Faces { get; } = new EpaFace[MaxEpaFaces];

    /// <summary>
    /// Highest face index in use (for iteration bounds).
    /// </summary>
    public int MaxFaceIndex { get; set; }

    public EpaPolytope()
    {
        for (int i = 0; i < MaxEpaVertices; i++)
        {
            Vertices[i] = new double[3];
            VertexIndices[i] = new int[2];
        }

        for (int i = 0; i < MaxEpaFaces; i++)
        {
            Faces[i] = new EpaFace();
        }
    }
}

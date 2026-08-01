namespace OpenGJKSharp.Models;

/// <summary>
/// Polytope structure for EPA.
/// Everything the EPA loop needs (vertices, faces, horizon edges and support scratch)
/// is allocated once here, so an instance reused across queries
/// (see <see cref="GjkWorkspace"/>) runs without further allocation.
/// </summary>
internal sealed class EpaPolytope
{
    public const int MaxEpaFaces = 128;

    public const int MaxEpaVertices = MaxEpaFaces + 4;

    public const int MaxEpaEdges = MaxEpaFaces * 3;

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
    /// Horizon edges collected while a new vertex is inserted.
    /// </summary>
    public EpaEdge[] Edges { get; } = new EpaEdge[MaxEpaEdges];

    /// <summary>
    /// Highest face index in use (for iteration bounds).
    /// </summary>
    public int MaxFaceIndex { get; set; }

    /// <summary>Centroid of the polytope, kept for face orientation tests.</summary>
    public double[] Centroid { get; } = new double[3];

    /// <summary>Support point of the current EPA iteration.</summary>
    public double[] NewVertex { get; } = new double[3];

    /// <summary>Body point ids behind <see cref="NewVertex"/>.</summary>
    public int[] NewVertexIdx { get; } = new int[2];

    /// <summary>Per-body support results, filled by the EPA support helper.</summary>
    public double[] SupportA { get; } = new double[3];

    /// <inheritdoc cref="SupportA"/>
    public double[] SupportB { get; } = new double[3];

    /// <summary>Negated search direction, filled by the EPA support helper.</summary>
    public double[] NegDirection { get; } = new double[3];

    /// <summary>Search direction used while a degenerate GJK simplex is grown to a tetrahedron.</summary>
    public double[] Direction { get; } = new double[3];

    /// <summary>General purpose vectors used while a degenerate GJK simplex is grown.</summary>
    public double[] TempA { get; } = new double[3];

    /// <inheritdoc cref="TempA"/>
    public double[] TempB { get; } = new double[3];

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

        for (int i = 0; i < MaxEpaEdges; i++)
        {
            Edges[i] = new EpaEdge();
        }
    }

    /// <summary>Rewinds the polytope so it can be rebuilt for another query.</summary>
    public void Reset()
    {
        NumVertices = 0;
        MaxFaceIndex = 0;

        for (int i = 0; i < MaxEpaFaces; i++)
        {
            Faces[i].Valid = false;
        }
    }
}

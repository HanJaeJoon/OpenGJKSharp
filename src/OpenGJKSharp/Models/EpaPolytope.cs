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

    /// <summary>
    /// Tightest penetration depth seen so far, that is the smallest support-plane distance
    /// over all face normals tried. <see cref="double.MaxValue"/> while nothing has been tried.
    /// </summary>
    public double BestDepth { get; set; }

    /// <summary>Face normal that produced <see cref="BestDepth"/>.</summary>
    public double[] BestNormal { get; } = new double[3];

    /// <summary>Polytope vertices of that face, kept because face slots get recycled.</summary>
    public int[] BestVertices { get; } = new int[3];

    /// <summary>Body point ids of that face's vertices [vertex][body].</summary>
    public int[][] BestVertexIndices { get; } = [new int[2], new int[2], new int[2]];

    /// <summary>
    /// Maximum number of vertices of the expanding polygon used when the Minkowski difference
    /// is flat and no tetrahedron exists.
    /// </summary>
    public const int MaxPlanarVertices = 64;

    /// <summary>
    /// Vertices of the expanding polygon, in cyclic order around the shared plane.
    /// </summary>
    public double[][] PlanarVertices { get; } = new double[MaxPlanarVertices][];

    /// <summary>Body point ids behind <see cref="PlanarVertices"/> [vertex][body].</summary>
    public int[][] PlanarVertexIndices { get; } = new int[MaxPlanarVertices][];

    public int NumPlanarVertices { get; set; }

    /// <summary>Normal of the plane the flat Minkowski difference lies in (not normalized).</summary>
    public double[] PlaneNormal { get; } = new double[3];

    /// <summary>Endpoints of the polygon edge that produced <see cref="BestDepth"/>.</summary>
    public double[][] PlanarBestEdge { get; } = [new double[3], new double[3]];

    /// <summary>Body point ids of those endpoints [endpoint][body].</summary>
    public int[][] PlanarBestEdgeIndices { get; } = [new int[2], new int[2]];

    /// <summary>In-plane outward normal that produced <see cref="BestDepth"/>.</summary>
    public double[] PlanarBestNormal { get; } = new double[3];

    /// <summary>Search direction used while a degenerate GJK simplex is grown to a tetrahedron.</summary>
    public double[] Direction { get; } = new double[3];

    /// <summary>General purpose vectors used while a degenerate GJK simplex is grown.</summary>
    public double[] TempA { get; } = new double[3];

    /// <inheritdoc cref="TempA"/>
    public double[] TempB { get; } = new double[3];

    /// <inheritdoc cref="TempA"/>
    public double[] TempC { get; } = new double[3];

    public EpaPolytope()
    {
        for (int i = 0; i < MaxEpaVertices; i++)
        {
            Vertices[i] = new double[3];
            VertexIndices[i] = new int[2];
        }

        for (int i = 0; i < MaxPlanarVertices; i++)
        {
            PlanarVertices[i] = new double[3];
            PlanarVertexIndices[i] = new int[2];
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
        BestDepth = double.MaxValue;

        for (int i = 0; i < MaxEpaFaces; i++)
        {
            Faces[i].Valid = false;
        }
    }
}

using System.Numerics;

namespace OpenGJKSharp.Models;

/// <summary>
/// Adapts a user-provided <see cref="ISupportMappable"/> shape to the internal body contract.
/// Every point the shape returns is appended to a growing store, and its store position acts
/// as the point id tracked by the simplex/EPA bookkeeping; <see cref="GetPoint"/> resolves ids
/// back to coordinates for witness reconstruction. The same instance must therefore be used
/// for the GJK and EPA passes of a single query.
/// <para>
/// <see cref="Reset"/> rewinds the store without releasing it, so an instance reused across
/// queries (see <see cref="GjkWorkspace"/>) keeps recycling the same coordinate arrays.
/// </para>
/// </summary>
internal sealed class GkMappedBody : GkBody
{
    private readonly List<double[]> _points = [];
    private ISupportMappable _shape = null!;
    private int _count;
    private int _sIdx;

    public GkMappedBody()
    {
    }

    public GkMappedBody(ISupportMappable shape)
    {
        _shape = shape;
    }

    public override double[] S { get; } = new double[3];

    public override int SIdx => _sIdx;

    /// <summary>Points this body at another shape and rewinds the point store for a new query.</summary>
    public void Reset(ISupportMappable shape)
    {
        _shape = shape;
        _count = 0;
        _sIdx = 0;
    }

    public override void InitSupport()
    {
        _sIdx = AddPoint(_shape.GetCenter(), S);
    }

    public override void Support(double[] v)
    {
        SupportInto(v, S, out _sIdx);
    }

    public override void SupportInto(double[] direction, double[] result, out int idx)
    {
        Vector3 point = _shape.GetSupportPoint(ToUnitDirection(direction));
        idx = AddPoint(point, result);
    }

    public override double[] GetPoint(int idx)
    {
        return _points[idx];
    }

    private int AddPoint(Vector3 point, double[] destination)
    {
        double[] coordinates;
        if (_count < _points.Count)
        {
            coordinates = _points[_count];
        }
        else
        {
            coordinates = new double[3];
            _points.Add(coordinates);
        }

        coordinates[0] = point.X;
        coordinates[1] = point.Y;
        coordinates[2] = point.Z;

        destination[0] = coordinates[0];
        destination[1] = coordinates[1];
        destination[2] = coordinates[2];

        return _count++;
    }

    // Normalizes in double before converting to float so tiny directions survive the
    // conversion. The interface contract promises implementations a non-zero direction,
    // so the zero vector (e.g. concentric shapes make the initial direction
    // centerA - centerB vanish) is replaced with an arbitrary axis here.
    private static Vector3 ToUnitDirection(double[] direction)
    {
        double normSquared = direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2];
        if (normSquared <= 0)
        {
            return Vector3.UnitX;
        }

        double inverseNorm = 1.0 / Math.Sqrt(normSquared);
        return new Vector3(
            (float)(direction[0] * inverseNorm),
            (float)(direction[1] * inverseNorm),
            (float)(direction[2] * inverseNorm));
    }
}

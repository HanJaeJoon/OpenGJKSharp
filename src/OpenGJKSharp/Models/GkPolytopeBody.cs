namespace OpenGJKSharp.Models;

/// <summary>
/// Vertex-cloud body backed by a <see cref="GkPolytope"/>.
/// Delegates support state to the wrapped polytope's <see cref="GkPolytope.S"/>/<see cref="GkPolytope.SIdx"/>
/// so their publicly observable behavior is unchanged.
/// </summary>
internal sealed class GkPolytopeBody : GkBody
{
    private readonly GkPolytope _polytope;

    public GkPolytopeBody(GkPolytope polytope)
    {
        _polytope = polytope;
    }

    public override double[] S => _polytope.S;

    public override int SIdx => _polytope.SIdx;

    public override void InitSupport()
    {
        for (int t = 0; t < 3; ++t)
        {
            _polytope.S[t] = _polytope.Coord[0][t];
        }
        _polytope.SIdx = 0;
    }

    public override void Support(double[] v)
    {
        // Stateful scan: starts from the current support point and only a strictly
        // better vertex replaces it (ties keep the previous point/index).
        double maxS = Dot(_polytope.S, v);
        int better = -1;

        for (int i = 0; i < _polytope.NumPoints; ++i)
        {
            double s = Dot(_polytope.Coord[i], v);
            if (s > maxS)
            {
                maxS = s;
                better = i;
            }
        }

        if (better != -1)
        {
            _polytope.S[0] = _polytope.Coord[better][0];
            _polytope.S[1] = _polytope.Coord[better][1];
            _polytope.S[2] = _polytope.Coord[better][2];
            _polytope.SIdx = better;
        }
    }

    public override void SupportInto(double[] direction, double[] result, out int idx)
    {
        // Fresh scan, first strict max wins (matches the previous EPA scan).
        int best = 0;
        double bestDot = Dot(_polytope.Coord[0], direction);

        for (int i = 1; i < _polytope.NumPoints; ++i)
        {
            double s = Dot(_polytope.Coord[i], direction);
            if (s > bestDot)
            {
                bestDot = s;
                best = i;
            }
        }

        for (int t = 0; t < 3; ++t)
        {
            result[t] = _polytope.Coord[best][t];
        }
        idx = best;
    }

    public override double[] GetPoint(int idx)
    {
        return _polytope.Coord[idx];
    }

    private static double Dot(double[] a, double[] b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
}

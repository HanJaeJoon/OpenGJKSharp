using OpenGJKSharp.Models;

namespace OpenGJKSharp;

public static class OpenGJKSharp
{
    private const double _gkEpsilon = 2.2204460492503131e-16; // DBL_EPSILON

    /// <summary>
    /// Invoke this function from C# applications
    /// </summary>
    /// <param name="nCoordsA"></param>
    /// <param name="inCoordsA"></param>
    /// <param name="nCoordsB"></param>
    /// <param name="inCoordsB"></param>
    /// <returns></returns>
    public static double CsFunction(int nCoordsA, double[,] inCoordsA, int nCoordsB, double[,] inCoordsB)
    {
        double[][] pinCoordsA = new double[nCoordsA][];
        for (int j = 0; j < nCoordsA; j++)
        {
            pinCoordsA[j] = [inCoordsA[0, j], inCoordsA[1, j], inCoordsA[2, j]];
        }

        double[][] pinCoordsB = new double[nCoordsB][];
        for (int j = 0; j < nCoordsB; j++)
        {
            pinCoordsB[j] = [inCoordsB[0, j], inCoordsB[1, j], inCoordsB[2, j]];
        }

        var bd1 = new GkPolytope()
        {
            Coord = pinCoordsA,
            NumPoints = nCoordsA,
        };
        var bd2 = new GkPolytope()
        {
            Coord = pinCoordsB,
            NumPoints = nCoordsB
        };

        var s = new GkSimplex();

        // Compute squared distance using GJK algorithm
        double distance = ComputeMinimumDistance(bd1, bd2, s);

        return distance;
    }

    public static double ComputeMinimumDistance(GkPolytope bd1, GkPolytope bd2, GkSimplex s)
    {
        uint k = 0; // Iteration counter
        const int mk = 25; // Maximum number of GJK iterations
        const double eps_rel = _gkEpsilon * 1e4; // Tolerance on relative
        const double eps_tot = _gkEpsilon * 1e2; // Tolerance on absolute distance
        const double eps_rel2 = eps_rel * eps_rel;
        double[] w = new double[3];
        int[] w_idx = new int[2];
        double[] v = new double[3];
        double[] vMinus = new double[3];
        double norm2WMax = 0;

        // Initialize search direction
        v[0] = bd1.Coord[0][0] - bd2.Coord[0][0];
        v[1] = bd1.Coord[0][1] - bd2.Coord[0][1];
        v[2] = bd1.Coord[0][2] - bd2.Coord[0][2];

        // Initialize simplex
        s.NVrtx = 1;
        for (int t = 0; t < 3; ++t)
        {
            s.Vrtx[0, t] = v[t];
        }
        s.VrtxIdx[0, 0] = 0;
        s.VrtxIdx[0, 1] = 0;

        for (int t = 0; t < 3; ++t)
        {
            bd1.S[t] = bd1.Coord[0][t];
        }
        bd1.SIdx = 0;

        for (int t = 0; t < 3; ++t)
        {
            bd2.S[t] = bd2.Coord[0][t];
        }
        bd2.SIdx = 0;

        // Begin GJK iteration
        do
        {
            k++;

            // Update negative search direction
            for (int t = 0; t < 3; ++t)
            {
                vMinus[t] = -v[t];
            }

            // Support function
            Support(bd1, vMinus);
            Support(bd2, v);
            for (int t = 0; t < 3; ++t)
            {
                w[t] = bd1.S[t] - bd2.S[t];
            }
            w_idx[0] = bd1.SIdx;
            w_idx[1] = bd2.SIdx;

            // Test first exit condition
            double exceedTol_rel = Norm2(v) - DotProduct(v, w);
            if (exceedTol_rel <= (eps_rel * Norm2(v)) || exceedTol_rel < eps_tot)
            {
                break;
            }

            if (Norm2(v) < eps_rel2)
            {
                break;
            }

            // Add new vertex to simplex
            int i = s.NVrtx;
            for (int t = 0; t < 3; ++t)
            {
                s.Vrtx[i, t] = w[t];
            }
            s.VrtxIdx[i, 0] = w_idx[0];
            s.VrtxIdx[i, 1] = w_idx[1];
            s.NVrtx++;

            // Invoke distance sub-algorithm
            SubAlgorithm(s, v);

            // Test
            for (int jj = 0; jj < s.NVrtx; jj++)
            {
                double testNorm = Norm2([s.Vrtx[jj, 0], s.Vrtx[jj, 1], s.Vrtx[jj, 2]]);
                if (testNorm > norm2WMax)
                {
                    norm2WMax = testNorm;
                }
            }

            if (Norm2(v) <= (eps_tot * eps_tot * norm2WMax))
            {
                break;
            }

        } while (s.NVrtx != 4 && k != mk);

        if (k == mk)
        {
            Console.WriteLine("\n * * * * * * * * * * * * MAXIMUM ITERATION NUMBER REACHED!!!  * * * * * * * * * * * * * * \n");
        }

        ComputeWitnesses(bd1, bd2, s);
        return Math.Sqrt(Norm2(v));
    }

    private static void Support(GkPolytope body, double[] v)
    {
        double s, maxS;
        double[] vrt;
        int better = -1;

        maxS = DotProduct(body.S, v);

        for (int i = 0; i < body.NumPoints; ++i)
        {
            vrt = body.Coord[i];
            s = DotProduct(vrt, v);
            if (s > maxS)
            {
                maxS = s;
                better = i;
            }
        }

        if (better != -1)
        {
            body.S[0] = body.Coord[better][0];
            body.S[1] = body.Coord[better][1];
            body.S[2] = body.Coord[better][2];
            body.SIdx = better;
        }
    }

    private static double DotProduct(double[] a, double[] b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    private static double Norm2(double[] a)
    {
        return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
    }

    private static void SubAlgorithm(GkSimplex s, double[] v)
    {
        switch (s.NVrtx)
        {
            case 4:
                S3D(s, v);
                break;
            case 3:
                S2D(s, v);
                break;
            case 2:
                S1D(s, v);
                break;
            default:
                Console.WriteLine("\nERROR:\t invalid simplex\n");
                break;
        }
    }

    private static void S1D(GkSimplex s, double[] v)
    {
        double[] s1p = [s.Vrtx[1, 0], s.Vrtx[1, 1], s.Vrtx[1, 2]];
        double[] s2p = [s.Vrtx[0, 0], s.Vrtx[0, 1], s.Vrtx[0, 2]];

        if (Hff1(s1p, s2p))
        {
            ProjectOnLine(s1p, s2p, v);
            return;
        }
        else
        {
            S1DRegion1(s, v);
            return;
        }
    }

    private static bool Hff1(double[] p, double[] q)
    {
        double tmp = 0;

        for (int i = 0; i < 3; i++)
        {
            tmp += p[i] * p[i] - p[i] * q[i];
        }

        if (tmp > 0)
        {
            return true;
        }
        return false;
    }

    private static void ProjectOnLine(double[] p, double[] q, double[] v)
    {
        double[] pq = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
        double tmp = DotProduct(p, pq) / DotProduct(pq, pq);

        for (int i = 0; i < 3; i++)
        {
            v[i] = p[i] - pq[i] * tmp;
        }
    }

    private static void S1DRegion1(GkSimplex s, double[] v)
    {
        v[0] = s.Vrtx[1, 0];
        v[1] = s.Vrtx[1, 1];
        v[2] = s.Vrtx[1, 2];

        s.NVrtx = 1;

        s.Vrtx[0, 0] = s.Vrtx[1, 0];
        s.Vrtx[0, 1] = s.Vrtx[1, 1];
        s.Vrtx[0, 2] = s.Vrtx[1, 2];

        s.VrtxIdx[0, 0] = s.VrtxIdx[1, 0];
        s.VrtxIdx[0, 1] = s.VrtxIdx[1, 1];
    }

    private static void S2D(GkSimplex s, double[] v)
    {
        double[] s1p = [s.Vrtx[2, 0], s.Vrtx[2, 1], s.Vrtx[2, 2]];
        double[] s2p = [s.Vrtx[1, 0], s.Vrtx[1, 1], s.Vrtx[1, 2]];
        double[] s3p = [s.Vrtx[0, 0], s.Vrtx[0, 1], s.Vrtx[0, 2]];

        bool hff1f_s12 = Hff1(s1p, s2p);
        bool hff1f_s13 = Hff1(s1p, s3p);

        if (hff1f_s12)
        {
            bool hff2f_23 = !Hff2(s1p, s2p, s3p);
            if (hff2f_23)
            {
                if (hff1f_s13)
                {
                    bool hff2f_32 = !Hff2(s1p, s3p, s2p);
                    if (hff2f_32)
                    {
                        ProjectOnPlane(s1p, s2p, s3p, v);
                        return;
                    }
                    else
                    {
                        ProjectOnLine(s1p, s3p, v);
                        S2DRegion13(s);
                        return;
                    }
                }
                else
                {
                    ProjectOnPlane(s1p, s2p, s3p, v);
                    return;
                }
            }
            else
            {
                ProjectOnLine(s1p, s2p, v);
                S2DRegion12(s);
                return;
            }
        }
        else if (hff1f_s13)
        {
            bool hff2f_32 = !Hff2(s1p, s3p, s2p);
            if (hff2f_32)
            {
                ProjectOnPlane(s1p, s2p, s3p, v);
                return;
            }
            else
            {
                ProjectOnLine(s1p, s3p, v);
                S2DRegion13(s);
                return;
            }
        }
        else
        {
            S2DRegion1(s, v);
            return;
        }
    }

    private static bool Hff2(double[] p, double[] q, double[] r)
    {
        double[] nTmp = new double[3];
        double[] n = new double[3];
        double[] pq = new double[3];
        double[] pr = new double[3];

        for (int i = 0; i < 3; i++)
        {
            pq[i] = q[i] - p[i];
        }
        for (int i = 0; i < 3; i++)
        {
            pr[i] = r[i] - p[i];
        }

        CrossProduct(pq, pr, nTmp);
        CrossProduct(pq, nTmp, n);

        return DotProduct(p, n) < 0;
    }

    private static void CrossProduct(double[] a, double[] b, double[] c)
    {
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
    }

    private static void ProjectOnPlane(double[] p, double[] q, double[] r, double[] v)
    {
        double[] n = new double[3];
        double[] pq = new double[3];
        double[] pr = new double[3];

        for (int i = 0; i < 3; i++)
        {
            pq[i] = p[i] - q[i];
        }
        for (int i = 0; i < 3; i++)
        {
            pr[i] = p[i] - r[i];
        }

        CrossProduct(pq, pr, n);
        double tmp = DotProduct(n, p) / DotProduct(n, n);

        for (int i = 0; i < 3; i++)
        {
            v[i] = n[i] * tmp;
        }
    }

    private static void S2DRegion13(GkSimplex s)
    {
        s.NVrtx = 2;

        s.Vrtx[1, 0] = s.Vrtx[2, 0];
        s.Vrtx[1, 1] = s.Vrtx[2, 1];
        s.Vrtx[1, 2] = s.Vrtx[2, 2];

        s.VrtxIdx[1, 0] = s.VrtxIdx[2, 0];
        s.VrtxIdx[1, 1] = s.VrtxIdx[2, 1];
    }

    private static void S2DRegion12(GkSimplex s)
    {
        s.NVrtx = 2;

        s.Vrtx[0, 0] = s.Vrtx[2, 0];
        s.Vrtx[0, 1] = s.Vrtx[2, 1];
        s.Vrtx[0, 2] = s.Vrtx[2, 2];

        s.VrtxIdx[0, 0] = s.VrtxIdx[2, 0];
        s.VrtxIdx[0, 1] = s.VrtxIdx[2, 1];
    }

    private static void S2DRegion1(GkSimplex s, double[] v)
    {
        v[0] = s.Vrtx[2, 0];
        v[1] = s.Vrtx[2, 1];
        v[2] = s.Vrtx[2, 2];

        s.NVrtx = 1;

        s.Vrtx[0, 0] = s.Vrtx[2, 0];
        s.Vrtx[0, 1] = s.Vrtx[2, 1];
        s.Vrtx[0, 2] = s.Vrtx[2, 2];

        s.VrtxIdx[0, 0] = s.VrtxIdx[2, 0];
        s.VrtxIdx[0, 1] = s.VrtxIdx[2, 1];
    }

    private static void S3D(GkSimplex s, double[] v)
    {
        double[] s1 = new double[3];
        double[] s2 = new double[3];
        double[] s3 = new double[3];
        double[] s4 = new double[3];
        double[] s1s2 = new double[3];
        double[] s1s3 = new double[3];
        double[] s1s4 = new double[3];
        double[] si = new double[3];
        double[] sj = new double[3];
        double[] sk = new double[3];
        int[] s1Idx = new int[2];
        int[] s2Idx = new int[2];
        int[] s3Idx = new int[2];
        int[] siIdx = new int[2];
        int[] sjIdx = new int[2];
        int[] skIdx = new int[2];
        bool testLineThree, testLineFour;
        int testPlaneTwo, testPlaneThree, testPlaneFour, dotTotal;
        int i, j, k;

        GetVrtxIdx(s, s1, s1Idx, 3);
        GetVrtxIdx(s, s2, s2Idx, 2);
        GetVrtxIdx(s, s3, s3Idx, 1);
        GetVrtx(s, s4, 0);
        CalculateEdgeVector(s, s1s2, s2);
        CalculateEdgeVector(s, s1s3, s3);
        CalculateEdgeVector(s, s1s4, s4);

        bool[] hff1Tests = new bool[3];
        hff1Tests[2] = Hff1(s1, s2);
        hff1Tests[1] = Hff1(s1, s3);
        hff1Tests[0] = Hff1(s1, s4);
        testLineThree = Hff1(s1, s3);
        testLineFour = Hff1(s1, s4);

        dotTotal = (Hff1(s1, s2) ? 1 : 0) + (testLineThree ? 1 : 0) + (testLineFour ? 1 : 0);
        if (dotTotal == 0)
        {
            S3DRegion1(s, v, s1, s1Idx);
            return;
        }

        double det134 = Determinant(s1s3, s1s4, s1s2);
        int sss = det134 <= 0 ? 1 : 0;

        testPlaneTwo = (Hff3(s1, s3, s4) ? 1 : 0) - sss;
        testPlaneTwo *= testPlaneTwo;
        testPlaneThree = (Hff3(s1, s4, s2) ? 1 : 0) - sss;
        testPlaneThree *= testPlaneThree;
        testPlaneFour = (Hff3(s1, s2, s3) ? 1 : 0) - sss;
        testPlaneFour *= testPlaneFour;

        switch (testPlaneTwo + testPlaneThree + testPlaneFour)
        {
            case 3:
                S3DRegion1234(s, v);
                break;

            case 2:
                s.NVrtx = 3;
                if (testPlaneTwo == 0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        s.Vrtx[2, i] = s.Vrtx[3, i];
                    }
                    for (i = 0; i < 2; i++)
                    {
                        s.VrtxIdx[2, i] = s.VrtxIdx[3, i];
                    }
                }
                else if (testPlaneThree == 0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        s.Vrtx[1, i] = s2[i];
                        s.Vrtx[2, i] = s.Vrtx[3, i];
                    }
                    for (i = 0; i < 2; i++)
                    {
                        s.VrtxIdx[1, i] = s2Idx[i];
                        s.VrtxIdx[2, i] = s.VrtxIdx[3, i];
                    }
                }
                else if (testPlaneFour == 0)
                {
                    for (i = 0; i < 3; i++)
                    {
                        s.Vrtx[0, i] = s3[i];
                        s.Vrtx[1, i] = s2[i];
                        s.Vrtx[2, i] = s.Vrtx[3, i];
                    }
                    for (i = 0; i < 2; i++)
                    {
                        s.VrtxIdx[0, i] = s3Idx[i];
                        s.VrtxIdx[1, i] = s2Idx[i];
                        s.VrtxIdx[2, i] = s.VrtxIdx[3, i];
                    }
                }
                S2D(s, v);
                break;

            case 1:
                s.NVrtx = 3;
                if (testPlaneTwo != 0)
                {
                    k = 2;
                    i = 1;
                    j = 0;
                }
                else if (testPlaneThree != 0)
                {
                    k = 1;
                    i = 0;
                    j = 2;
                }
                else
                {
                    k = 0;
                    i = 2;
                    j = 1;
                }

                GetVrtxIdx(s, si, siIdx, i);
                GetVrtxIdx(s, sj, sjIdx, j);
                GetVrtxIdx(s, sk, skIdx, k);

                if (dotTotal == 1)
                {
                    if (hff1Tests[k])
                    {
                        if (!Hff2(s1, sk, si))
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, si, sk, v);
                        }
                        else if (!Hff2(s1, sk, sj))
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sj, sk, v);
                        }
                        else
                        {
                            Select1(s, sk, skIdx);
                            ProjectOnLine(s1, sk, v);
                        }
                    }
                    else if (hff1Tests[i])
                    {
                        if (!Hff2(s1, si, sk))
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, si, sk, v);
                        }
                        else
                        {
                            Select1(s, si, siIdx);
                            ProjectOnLine(s1, si, v);
                        }
                    }
                    else
                    {
                        if (!Hff2(s1, sj, sk))
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sj, sk, v);
                        }
                        else
                        {
                            Select1(s, sj, sjIdx);
                            ProjectOnLine(s1, sj, v);
                        }
                    }
                }
                else if (dotTotal == 2)
                {
                    if (hff1Tests[i])
                    {
                        if (!Hff2(s1, sk, si))
                        {
                            if (!Hff2(s1, si, sk))
                            {
                                Select1(s, si, siIdx, sk, skIdx);
                                ProjectOnPlane(s1, si, sk, v);
                            }
                            else
                            {
                                Select1(s, sk, skIdx);
                                ProjectOnLine(s1, sk, v);
                            }
                        }
                        else
                        {
                            if (!Hff2(s1, sk, sj))
                            {
                                Select1(s, sj, sjIdx, sk, skIdx);
                                ProjectOnPlane(s1, sj, sk, v);
                            }
                            else
                            {
                                Select1(s, sk, skIdx);
                                ProjectOnLine(s1, sk, v);
                            }
                        }
                    }
                    else if (hff1Tests[j])
                    {
                        if (!Hff2(s1, sk, sj))
                        {
                            if (!Hff2(s1, sj, sk))
                            {
                                Select1(s, sj, sjIdx, sk, skIdx);
                                ProjectOnPlane(s1, sj, sk, v);
                            }
                            else
                            {
                                Select1(s, sj, sjIdx);
                                ProjectOnLine(s1, sj, v);
                            }
                        }
                        else
                        {
                            if (!Hff2(s1, sk, si))
                            {
                                Select1(s, si, siIdx, sk, skIdx);
                                ProjectOnPlane(s1, si, sk, v);
                            }
                            else
                            {
                                Select1(s, sk, skIdx);
                                ProjectOnLine(s1, sk, v);
                            }
                        }
                    }
                }
                else if (dotTotal == 3)
                {
                    bool hff2Ik = Hff2(s1, si, sk);
                    bool hff2Jk = Hff2(s1, sj, sk);
                    bool hff2Ki = Hff2(s1, sk, si);
                    bool hff2Kj = Hff2(s1, sk, sj);

                    if (!hff2Ki && !hff2Kj)
                    {
                        Console.WriteLine("\n\n UNEXPECTED VALUES!!! \n\n");
                    }
                    if (hff2Ki && hff2Kj)
                    {
                        Select1(s, sk, skIdx);
                        ProjectOnLine(s1, sk, v);
                    }
                    else if (hff2Ki)
                    {
                        if (hff2Jk)
                        {
                            Select1(s, sj, sjIdx);
                            ProjectOnLine(s1, sj, v);
                        }
                        else
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sk, sj, v);
                        }
                    }
                    else
                    {
                        if (hff2Ik)
                        {
                            Select1(s, si, siIdx);
                            ProjectOnLine(s1, si, v);
                        }
                        else
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, sk, si, v);
                        }
                    }
                }
                break;

            case 0:
                if (dotTotal == 1)
                {
                    if (testLineThree)
                    {
                        k = 2;
                        i = 1;
                        j = 0;
                    }
                    else if (testLineFour)
                    {
                        k = 1;
                        i = 0;
                        j = 2;
                    }
                    else
                    {
                        k = 0;
                        i = 2;
                        j = 1;
                    }
                    GetVrtxIdx(s, si, siIdx, i);
                    GetVrtxIdx(s, sj, sjIdx, j);
                    GetVrtxIdx(s, sk, skIdx, k);

                    if (!Hff2(s1, si, sj))
                    {
                        Select1(s, si, siIdx, sj, sjIdx);
                        ProjectOnPlane(s1, si, sj, v);
                    }
                    else if (!Hff2(s1, si, sk))
                    {
                        Select1(s, si, siIdx, sk, skIdx);
                        ProjectOnPlane(s1, si, sk, v);
                    }
                    else
                    {
                        Select1(s, si, siIdx);
                        ProjectOnLine(s1, si, v);
                    }
                }
                else if (dotTotal == 2)
                {
                    s.NVrtx = 3;
                    if (!testLineThree)
                    {
                        k = 2;
                        i = 1;
                        j = 0;
                    }
                    else if (!testLineFour)
                    {
                        k = 1;
                        i = 0;
                        j = 2;
                    }
                    else
                    {
                        k = 0;
                        i = 2;
                        j = 1;
                    }
                    GetVrtxIdx(s, si, siIdx, i);
                    GetVrtxIdx(s, sj, sjIdx, j);
                    GetVrtxIdx(s, sk, skIdx, k);

                    if (!Hff2(s1, sj, sk))
                    {
                        if (!Hff2(s1, sk, sj))
                        {
                            Select1(s, sj, sjIdx, sk, skIdx);
                            ProjectOnPlane(s1, sj, sk, v);
                        }
                        else if (!Hff2(s1, sk, si))
                        {
                            Select1(s, si, siIdx, sk, skIdx);
                            ProjectOnPlane(s1, sk, si, v);
                        }
                        else
                        {
                            Select1(s, sk, skIdx);
                            ProjectOnLine(s1, sk, v);
                        }
                    }
                    else if (!Hff2(s1, sj, si))
                    {
                        Select1(s, si, siIdx, sj, sjIdx);
                        ProjectOnPlane(s1, si, sj, v);
                    }
                    else
                    {
                        Select1(s, sj, sjIdx);
                        ProjectOnLine(s1, sj, v);
                    }
                }
                break;

            default:
                Console.WriteLine("\nERROR:\tunhandled");
                break;
        }
    }

    private static void GetVrtxIdx(GkSimplex s, double[] point, int[] index, int location)
    {
        point[0] = s.Vrtx[location, 0];
        point[1] = s.Vrtx[location, 1];
        point[2] = s.Vrtx[location, 2];

        index[0] = s.VrtxIdx[location, 0];
        index[1] = s.VrtxIdx[location, 1];
    }

    private static void GetVrtx(GkSimplex s, double[] point, int location)
    {
        point[0] = s.Vrtx[location, 0];
        point[1] = s.Vrtx[location, 1];
        point[2] = s.Vrtx[location, 2];
    }

    private static void CalculateEdgeVector(GkSimplex s, double[] p1p2, double[] p2)
    {
        p1p2[0] = p2[0] - s.Vrtx[3, 0];
        p1p2[1] = p2[1] - s.Vrtx[3, 1];
        p1p2[2] = p2[2] - s.Vrtx[3, 2];
    }

    private static void S3DRegion1(GkSimplex s, double[] v, double[] s1, int[] s1Idx)
    {
        v[0] = s1[0];
        v[1] = s1[1];
        v[2] = s1[2];

        s.NVrtx = 1;

        s.Vrtx[0, 0] = s1[0];
        s.Vrtx[0, 1] = s1[1];
        s.Vrtx[0, 2] = s1[2];

        s.VrtxIdx[0, 0] = s1Idx[0];
        s.VrtxIdx[0, 1] = s1Idx[1];
    }

    private static double Determinant(double[] p, double[] q, double[] r)
    {
        return p[0] * (q[1] * r[2] - r[1] * q[2])
            - p[1] * (q[0] * r[2] - r[0] * q[2])
            + p[2] * (q[0] * r[1] - r[0] * q[1]);
    }

    private static bool Hff3(double[] p, double[] q, double[] r)
    {
        double[] n = new double[3];
        double[] pq = new double[3];
        double[] pr = new double[3];

        for (int i = 0; i < 3; i++)
        {
            pq[i] = q[i] - p[i];
        }
        for (int i = 0; i < 3; i++)
        {
            pr[i] = r[i] - p[i];
        }

        CrossProduct(pq, pr, n);
        return DotProduct(p, n) <= 0;
    }

    private static void S3DRegion1234(GkSimplex s, double[] v)
    {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
        s.NVrtx = 4;
    }

    private static void Select1(GkSimplex s, double[] si, int[] siIdx, double[] some, int[] someIdx)
    {
        s.NVrtx = 3;

        for (int t = 0; t < 3; t++)
            s.Vrtx[2, t] = s.Vrtx[3, t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[2, t] = s.VrtxIdx[3, t];
        for (int t = 0; t < 3; t++)
            s.Vrtx[1, t] = si[t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[1, t] = siIdx[t];
        for (int t = 0; t < 3; t++)
            s.Vrtx[0, t] = some[t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[0, t] = someIdx[t];
    }

    private static void Select1(GkSimplex s, double[] some, int[] someIdx)
    {
        s.NVrtx = 2;

        for (int t = 0; t < 3; t++)
            s.Vrtx[1, t] = s.Vrtx[3, t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[1, t] = s.VrtxIdx[3, t];
        for (int t = 0; t < 3; t++)
            s.Vrtx[0, t] = some[t];
        for (int t = 0; t < 2; t++)
            s.VrtxIdx[0, t] = someIdx[t];
    }

    private static void ComputeWitnesses(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        switch (smp.NVrtx)
        {
            case 4:
                W3D(bd1, bd2, smp);
                break;
            case 3:
                W2D(bd1, bd2, smp);
                break;
            case 2:
                W1D(bd1, bd2, smp);
                break;
            case 1:
                W0D(bd1, bd2, smp);
                break;
            default:
                Console.WriteLine("\nERROR:\t invalid simplex\n");
                break;
        }
    }

    private static void W0D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        var w01 = bd2.Coord[smp.VrtxIdx[0, 1]];

        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t];
            smp.Witnesses[1, t] = w01[t];
        }
    }

    private static void W1D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var pq = new double[3];
        var po = new double[3];

        double[] p = [smp.Vrtx[0, 0], smp.Vrtx[0, 1], smp.Vrtx[0, 2]];
        double[] q = [smp.Vrtx[1, 0], smp.Vrtx[1, 1], smp.Vrtx[1, 2]];

        for (int t = 0; t < 3; t++)
        {
            pq[t] = q[t] - p[t];
            po[t] = -p[t];
        }

        // Compute barycentric coordinates via matrix inversion
        // (in the linear case the matrix is 1x1 thus simplified)
        double det = DotProduct(pq, pq);
        if (det == 0.0)
        {
            // Degenerate case
            W0D(bd1, bd2, smp);
        }

        double a1 = DotProduct(pq, po) / det;
        double a0 = 1.0 - a1;

        // Compute witness points
        double[] w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        double[] w01 = bd2.Coord[smp.VrtxIdx[0, 1]];
        double[] w10 = bd1.Coord[smp.VrtxIdx[1, 0]];
        double[] w11 = bd2.Coord[smp.VrtxIdx[1, 1]];
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1;
        }
    }

    private static void W2D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var pq = new double[3];
        var pr = new double[3];
        var po = new double[3];

        double[] p = [smp.Vrtx[0, 0], smp.Vrtx[0, 1], smp.Vrtx[0, 2]];
        double[] q = [smp.Vrtx[1, 0], smp.Vrtx[1, 1], smp.Vrtx[1, 2]];
        double[] r = [smp.Vrtx[2, 0], smp.Vrtx[2, 1], smp.Vrtx[2, 2]];

        for (int t = 0; t < 3; t++)
        {
            pq[t] = q[t] - p[t];
            pr[t] = r[t] - p[t];
            po[t] = -p[t];
        }

        double t00 = DotProduct(pq, pq);
        double t01 = DotProduct(pq, pr);
        double t11 = DotProduct(pr, pr);
        double det = t00 * t11 - t01 * t01;
        if (det == 0.0)
        {
            // Degenerate case
            W1D(bd1, bd2, smp);
        }

        double b0 = DotProduct(pq, po);
        double b1 = DotProduct(pr, po);
        double i00 = t11 / det;
        double i01 = -t01 / det;
        double i11 = t00 / det;
        double a1 = i00 * b0 + i01 * b1;
        double a2 = i01 * b0 + i11 * b1;
        double a0 = 1.0 - a1 - a2;

        // check if the origin is very close to one of the edges of the
        // simplex. In this case, a 1D projection will be more accurate.
        if (a0 < _gkEpsilon)
        {
            smp.NVrtx = 2;
            smp.Vrtx[0, 0] = smp.Vrtx[2, 0];
            smp.Vrtx[0, 1] = smp.Vrtx[2, 1];
            smp.Vrtx[0, 2] = smp.Vrtx[2, 2];
            smp.VrtxIdx[0, 0] = smp.VrtxIdx[2, 0];
            smp.VrtxIdx[0, 1] = smp.VrtxIdx[2, 1];
            W1D(bd1, bd2, smp);
        }
        else if (a1 < _gkEpsilon)
        {
            smp.NVrtx = 2;
            smp.Vrtx[1, 0] = smp.Vrtx[2, 0];
            smp.Vrtx[1, 1] = smp.Vrtx[2, 1];
            smp.Vrtx[1, 2] = smp.Vrtx[2, 2];
            smp.VrtxIdx[1, 0] = smp.VrtxIdx[2, 0];
            smp.VrtxIdx[1, 1] = smp.VrtxIdx[2, 1];
            W1D(bd1, bd2, smp);
        }
        else if (a2 < _gkEpsilon)
        {
            smp.NVrtx = 2;
            W1D(bd1, bd2, smp);
        }

        // Compute witness points
        // This is done by blending the source points using
        // the barycentric coordinates
        double[] w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        double[] w01 = bd2.Coord[smp.VrtxIdx[0, 1]];
        double[] w10 = bd1.Coord[smp.VrtxIdx[1, 0]];
        double[] w11 = bd2.Coord[smp.VrtxIdx[1, 1]];
        double[] w20 = bd1.Coord[smp.VrtxIdx[2, 0]];
        double[] w21 = bd2.Coord[smp.VrtxIdx[2, 1]];
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1 + w20[t] * a2;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1 + w21[t] * a2;
        }
    }

    private static void W3D(GkPolytope bd1, GkPolytope bd2, GkSimplex smp)
    {
        var pq = new double[3];
        var pr = new double[3];
        var ps = new double[3];
        var po = new double[3];

        double[] p = [smp.Vrtx[0, 0], smp.Vrtx[0, 1], smp.Vrtx[0, 2]];
        double[] q = [smp.Vrtx[1, 0], smp.Vrtx[1, 1], smp.Vrtx[1, 2]];
        double[] r = [smp.Vrtx[2, 0], smp.Vrtx[2, 1], smp.Vrtx[2, 2]];
        double[] s = [smp.Vrtx[3, 0], smp.Vrtx[3, 1], smp.Vrtx[3, 2]];

        for (int t = 0; t < 3; t++)
        {
            pq[t] = q[t] - p[t];
            pr[t] = r[t] - p[t];
            ps[t] = s[t] - p[t];
            po[t] = -p[t];
        }

        double t00 = DotProduct(pq, pq);
        double t01 = DotProduct(pq, pr);
        double t02 = DotProduct(pq, ps);
        double t11 = DotProduct(pr, pr);
        double t12 = DotProduct(pr, ps);
        double t22 = DotProduct(ps, ps);
        double det00 = t11 * t22 - t12 * t12;
        double det01 = t01 * t22 - t02 * t12;
        double det02 = t01 * t12 - t02 * t11;
        double det = t00 * det00 - t01 * det01 + t02 * det02;
        if (det == 0.0)
        {
            // Degenerate case
            W2D(bd1, bd2, smp);
        }

        double b0 = DotProduct(pq, po);
        double b1 = DotProduct(pr, po);
        double b2 = DotProduct(ps, po);

        // inverse matrix
        // (the matrix is symmetric, so we can use the cofactor matrix)
        double det11 = t00 * t22 - t02 * t02;
        double det12 = t00 * t12 - t01 * t02;
        double det22 = t00 * t11 - t01 * t01;
        double i00 = det00 / det;
        double i01 = -det01 / det;
        double i02 = det02 / det;
        double i11 = det11 / det;
        double i12 = -det12 / det;
        double i22 = det22 / det;

        double a1 = i00 * b0 + i01 * b1 + i02 * b2;
        double a2 = i01 * b0 + i11 * b1 + i12 * b2;
        double a3 = i02 * b0 + i12 * b1 + i22 * b2;
        double a0 = 1.0 - a1 - a2 - a3;

        // check if the origin is very close to one of the faces of the
        // simplex. In this case, a 2D projection will be more accurate.
        if (a0 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            smp.Vrtx[0, 0] = smp.Vrtx[3, 0];
            smp.Vrtx[0, 1] = smp.Vrtx[3, 1];
            smp.Vrtx[0, 2] = smp.Vrtx[3, 2];
            smp.VrtxIdx[0, 0] = smp.VrtxIdx[3, 0];
            smp.VrtxIdx[0, 1] = smp.VrtxIdx[3, 1];
            W2D(bd1, bd2, smp);
        }
        else if (a1 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            smp.Vrtx[1, 0] = smp.Vrtx[3, 0];
            smp.Vrtx[1, 1] = smp.Vrtx[3, 1];
            smp.Vrtx[1, 2] = smp.Vrtx[3, 2];
            smp.VrtxIdx[1, 0] = smp.VrtxIdx[3, 0];
            smp.VrtxIdx[1, 1] = smp.VrtxIdx[3, 1];
            W2D(bd1, bd2, smp);
        }
        else if (a2 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            smp.Vrtx[2, 0] = smp.Vrtx[3, 0];
            smp.Vrtx[2, 1] = smp.Vrtx[3, 1];
            smp.Vrtx[2, 2] = smp.Vrtx[3, 2];
            smp.VrtxIdx[2, 0] = smp.VrtxIdx[3, 0];
            smp.VrtxIdx[2, 1] = smp.VrtxIdx[3, 1];
            W2D(bd1, bd2, smp);
        }
        else if (a3 < _gkEpsilon)
        {
            smp.NVrtx = 3;
            W2D(bd1, bd2, smp);
        }

        // Compute witness points
        // This is done by blending the original points using
        // the barycentric coordinates
        double[] w00 = bd1.Coord[smp.VrtxIdx[0, 0]];
        double[] w01 = bd2.Coord[smp.VrtxIdx[0, 1]];
        double[] w10 = bd1.Coord[smp.VrtxIdx[1, 0]];
        double[] w11 = bd2.Coord[smp.VrtxIdx[1, 1]];
        double[] w20 = bd1.Coord[smp.VrtxIdx[2, 0]];
        double[] w21 = bd2.Coord[smp.VrtxIdx[2, 1]];
        double[] w30 = bd1.Coord[smp.VrtxIdx[3, 0]];
        double[] w31 = bd2.Coord[smp.VrtxIdx[3, 1]];
        for (int t = 0; t < 3; t++)
        {
            smp.Witnesses[0, t] = w00[t] * a0 + w10[t] * a1 + w20[t] * a2 + w30[t] * a3;
            smp.Witnesses[1, t] = w01[t] * a0 + w11[t] * a1 + w21[t] * a2 + w31[t] * a3;
        }
    }
}

#include "NurbsIntersection.h"

#include "NurbsSurface.h"
#include "NurbsCurve.h"
#include "NurbsUtil.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <limits>

namespace {
    bool find_surface_intersection_points(const NurbsSurface& s1, const NurbsSurface& s2,
                                        std::vector<Point3>& intersectionPoints,
                                        std::vector<std::pair<double, double>>& params1,
                                        std::vector<std::pair<double, double>>& params2);

    bool connect_intersection_points(const std::vector<Point3>& points,
                                   const std::vector<std::pair<double, double>>& params1,
                                   const std::vector<std::pair<double, double>>& params2,
                                   NurbsIntersectionCurve& curve);

    double point_surface_distance(const Point3& point, const NurbsSurface& surface,
                                double& u, double& v, double tol = 1e-6);

    bool march_along_curve(const NurbsSurface& s1, const NurbsSurface& s2,
                          double u1_start, double v1_start, double u2_start, double v2_start,
                          std::vector<Point3>& points,
                          std::vector<std::pair<double, double>>& params1,
                          std::vector<std::pair<double, double>>& params2,
                          double step_size = 0.01, int max_steps = 1000);

    bool find_surface_intersection_points(const NurbsSurface& s1, const NurbsSurface& s2,
                                        std::vector<Point3>& intersectionPoints,
                                        std::vector<std::pair<double, double>>& params1,
                                        std::vector<std::pair<double, double>>& params2)
    {
        intersectionPoints.clear();
        params1.clear();
        params2.clear();

        // Simple grid sampling approach for finding intersection points
        // In a full implementation, this would use more sophisticated methods

        const int grid_size = 20; // Sample points per surface
        const double tol = 1e-6;

        for (int i = 0; i <= grid_size; ++i) {
            for (int j = 0; j <= grid_size; ++j) {
                double u1 = static_cast<double>(i) / grid_size;
                double v1 = static_cast<double>(j) / grid_size;

                Point3 p1;
                s1.evaluate_clamped(u1, v1, p1);

                // Project p1 onto s2
                double u2_proj, v2_proj;
                double dist = point_surface_distance(p1, s2, u2_proj, v2_proj, tol);

                if (dist < tol) {
                    // Check if this point is already found (within tolerance)
                    bool duplicate = false;
                    for (size_t k = 0; k < intersectionPoints.size(); ++k) {
                        if (intersectionPoints[k].distance_square(p1) < tol * tol) {
                            duplicate = true;
                            break;
                        }
                    }

                    if (!duplicate) {
                        intersectionPoints.push_back(p1);
                        params1.emplace_back(u1, v1);
                        params2.emplace_back(u2_proj, v2_proj);
                    }
                }
            }
        }

        return !intersectionPoints.empty();
    }

    bool connect_intersection_points(const std::vector<Point3>& points,
                                   const std::vector<std::pair<double, double>>& params1,
                                   const std::vector<std::pair<double, double>>& params2,
                                   NurbsIntersectionCurve& curve)
    {
        if (points.empty()) return false;

        curve.samples.clear();
        curve.closed = false;

        // For now, just add all points as samples
        // In a full implementation, this would sort and connect them properly
        for (size_t i = 0; i < points.size(); ++i) {
            NurbsIntersectionSample sample;
            sample.point = points[i];
            sample.uA = params1[i].first;
            sample.vA = params1[i].second;
            sample.uB = params2[i].first;
            sample.vB = params2[i].second;
            curve.samples.push_back(sample);
        }

        return true;
    }

    double point_surface_distance(const Point3& point, const NurbsSurface& surface,
                                double& u, double& v, double tol)
    {
        // Simple projection using surface evaluation
        // In a full implementation, this would use Newton iteration or similar

        double min_dist = std::numeric_limits<double>::max();
        double best_u = 0.5, best_v = 0.5;

        const int samples = 10;
        for (int i = 0; i <= samples; ++i) {
            for (int j = 0; j <= samples; ++j) {
                double test_u = static_cast<double>(i) / samples;
                double test_v = static_cast<double>(j) / samples;

                Point3 p;
                surface.evaluate_clamped(test_u, test_v, p);

                double dist = point.distance_square(p);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_u = test_u;
                    best_v = test_v;
                }
            }
        }

        u = best_u;
        v = best_v;
        return std::sqrt(min_dist);
    }

    bool march_along_curve(const NurbsSurface& s1, const NurbsSurface& s2,
                          double u1_start, double v1_start, double u2_start, double v2_start,
                          std::vector<Point3>& points,
                          std::vector<std::pair<double, double>>& params1,
                          std::vector<std::pair<double, double>>& params2,
                          double step_size, int max_steps)
    {
        // Basic marching algorithm - in a full implementation this would be more sophisticated
        double u1 = u1_start, v1 = v1_start;
        double u2 = u2_start, v2 = v2_start;

        points.push_back(Point3()); // Placeholder
        params1.emplace_back(u1, v1);
        params2.emplace_back(u2, v2);

        // This is a simplified implementation
        return true;
    }
}

// Evaluate residual F and Jacobian J for parameters p = [u1,v1,u2,v2]
static void eval_F_and_J(const NurbsSurface& s1, const NurbsSurface& s2,
                         const double p[4], double F[3], double J[3][4])
{
    double u1 = p[0];
    double v1 = p[1];
    double u2 = p[2];
    double v2 = p[3];

    Point3 P1; Point3 d1u, d1v, d1uu, d1uv, d1vv;
    s1.evaluate_derivatives(u1, v1, d1u, d1v, d1uu, d1uv, d1vv);
    s1.evaluate(u1, v1, P1);

    Point3 P2; Point3 d2u, d2v, d2uu, d2uv, d2vv;
    s2.evaluate_derivatives(u2, v2, d2u, d2v, d2uu, d2uv, d2vv);
    s2.evaluate(u2, v2, P2);

    Point3 R = P1 - P2;
    F[0] = R.x(); F[1] = R.y(); F[2] = R.z();

    // J = [dS1/du dS1/dv -dS2/du -dS2/dv]
    J[0][0] = d1u.x(); J[0][1] = d1v.x(); J[0][2] = -d2u.x(); J[0][3] = -d2v.x();
    J[1][0] = d1u.y(); J[1][1] = d1v.y(); J[1][2] = -d2u.y(); J[1][3] = -d2v.y();
    J[2][0] = d1u.z(); J[2][1] = d1v.z(); J[2][2] = -d2u.z(); J[2][3] = -d2v.z();
}

// Solve small linear system A x = b for 4x4 A using Gaussian elimination
static bool solve4x4(double A[4][4], double b[4], double xOut[4])
{
    const int N = 4;
    double M[4][5];
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j) M[i][j] = A[i][j];
        M[i][4] = b[i];
    }

    for (int i = 0; i < N; ++i)
    {
        int piv = i;
        for (int r = i+1; r < N; ++r)
            if (fabs(M[r][i]) > fabs(M[piv][i])) piv = r;
        if (fabs(M[piv][i]) < 1e-15) return false;
        if (piv != i) for (int c = i; c < N+1; ++c) std::swap(M[i][c], M[piv][c]);

        double div = M[i][i];
        for (int c = i; c < N+1; ++c) M[i][c] /= div;
        for (int r = 0; r < N; ++r) if (r != i)
        {
            double mul = M[r][i];
            for (int c = i; c < N+1; ++c) M[r][c] -= mul * M[i][c];
        }
    }

    for (int i = 0; i < N; ++i) xOut[i] = M[i][N];
    return true;
}

// Levenberg-Marquardt local refine for seed parameters
static bool lm_refine_seed(const NurbsSurface& s1, const NurbsSurface& s2, double p[4], const IntersectionOptions& opt)
{
    double lambda = opt.initialDamping;
    double F[3]; double J[3][4];

    for (int iter = 0; iter < opt.maxNewtonIterations; ++iter)
    {
        eval_F_and_J(s1, s2, p, F, J);
        double res2 = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
        if (res2 <= opt.geomTol * opt.geomTol) return true;

        // Build normal equations A = J^T J, g = -J^T F
        double A[4][4] = {0}; double g[4] = {0};
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < 3; ++k) sum += J[k][i] * J[k][j];
                A[i][j] = sum;
            }
        for (int i = 0; i < 4; ++i)
        {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) sum += J[k][i] * F[k];
            g[i] = -sum;
        }

        // Damping
        for (int i = 0; i < 4; ++i) A[i][i] += lambda;

        double delta[4];
        if (!solve4x4(A, g, delta))
        {
            lambda *= 10.0;
            if (lambda > 1e12) return false;
            continue;
        }

        // try step
        double pTrial[4];
        for (int i = 0; i < 4; ++i) pTrial[i] = p[i] + delta[i];
        // clamp to [0,1]
        for (int i = 0; i < 4; ++i) if (pTrial[i] < 0.) pTrial[i] = 0.; else if (pTrial[i] > 1.) pTrial[i] = 1.;

        double Ftrial[3]; double Jtrial[3][4];
        eval_F_and_J(s1, s2, pTrial, Ftrial, Jtrial);
        double res2Trial = Ftrial[0]*Ftrial[0] + Ftrial[1]*Ftrial[1] + Ftrial[2]*Ftrial[2];

        if (res2Trial < res2)
        {
            // accept
            for (int i = 0; i < 4; ++i) p[i] = pTrial[i];
            lambda *= 0.1;
            if (lambda < 1e-16) lambda = 1e-16;
        }
        else
        {
            lambda *= 10.0;
            if (lambda > 1e12) return false;
        }
    }

    // final residual check
    eval_F_and_J(s1, s2, p, F, J);
    double res2 = F[0]*F[0] + F[1]*F[1] + F[2]*F[2];
    return res2 <= opt.geomTol * opt.geomTol;
}

// Main intersection function implementation (coarse seed detection + local refinement)
void compute_surface_intersection(const NurbsSurface& surface1, const NurbsSurface& surface2,
                                NurbsIntersectionResult& result,
                                const IntersectionOptions& options)
{
    result.curves.clear();
    result.hasPartialOverlap = false;

    const int resA = std::max(1, options.seedSamplingRes);
    const int resB = std::max(1, options.seedSamplingRes);

    struct Sample { Point3 p; double u; double v; };
    std::vector<Sample> samplesA; samplesA.reserve((resA+1)*(resA+1));
    std::vector<Sample> samplesB; samplesB.reserve((resB+1)*(resB+1));

    for (int i = 0; i <= resA; ++i)
    {
        double u = (double)i / resA;
        for (int j = 0; j <= resA; ++j)
        {
            double v = (double)j / resA;
            Point3 p; surface1.evaluate(u, v, p);
            samplesA.push_back({p,u,v});
        }
    }

    for (int i = 0; i <= resB; ++i)
    {
        double u = (double)i / resB;
        for (int j = 0; j <= resB; ++j)
        {
            double v = (double)j / resB;
            Point3 p; surface2.evaluate(u, v, p);
            samplesB.push_back({p,u,v});
        }
    }

    double seedTol = options.geomTol * options.seedTolScale;
    double seedTolSq = seedTol * seedTol;

    std::vector<std::array<double,4>> seeds;

    for (const auto& sa : samplesA)
    {
        double best = 1e300; int ibest = -1;
        for (int k = 0; k < (int)samplesB.size(); ++k)
        {
            double d2 = sa.p.distance_square(samplesB[k].p);
            if (d2 < best) { best = d2; ibest = k; }
        }
        if (best <= seedTolSq && ibest >= 0)
        {
            seeds.push_back({sa.u, sa.v, samplesB[ibest].u, samplesB[ibest].v});
        }
    }

    // refine unique seeds: simple dedup by proximity in param space
    // refine seeds and deduplicate based on refined parameters
    const double dupParamTol = 1e-6;
    std::vector<std::array<double,4>> refinedSeeds;
    for (const auto& s : seeds)
    {
        double p[4] = { s[0], s[1], s[2], s[3] };
        if (!lm_refine_seed(surface1, surface2, p, options))
            continue;

        bool dup = false;
        for (const auto& r : refinedSeeds)
        {
            double du = p[0]-r[0]; double dv = p[1]-r[1]; double du2 = p[2]-r[2]; double dv2 = p[3]-r[3];
            if (du*du + dv*dv + du2*du2 + dv2*dv2 < dupParamTol*dupParamTol) { dup = true; break; }
        }
        if (!dup) refinedSeeds.push_back({p[0],p[1],p[2],p[3]});
    }

    // For each refined seed, trace full curve
    for (const auto& seed : refinedSeeds)
    {
        double pseed[4] = { seed[0], seed[1], seed[2], seed[3] };
        Point3 pa; surface1.evaluate(pseed[0], pseed[1], pa);
        Point3 pb; surface2.evaluate(pseed[2], pseed[3], pb);
        Point3 seedPoint = (pa + pb) * 0.5;

        NurbsIntersectionCurve curve;
        // forward declaration for trace helper is below, so we call a local lambda instead
        auto trace_from_seed_fn = [&](const double seedP[4], const Point3& seedPoint,
                                      const IntersectionOptions& opt, NurbsIntersectionCurve& outCurve)
        {
            const double stepInit = opt.initialStep;
            const double minStep = opt.initialStep * 1e-3;
            const int maxSteps = 1000;

            auto trace_dir = [&](int sign)
            {
                double pcur[4] = { seedP[0], seedP[1], seedP[2], seedP[3] };
                Point3 Pcur = seedPoint;
                double step = stepInit;

                for (int iter = 0; iter < maxSteps; ++iter)
                {
                    Point3 d1u, d1v, d1uu, d1uv, d1vv;
                    surface1.evaluate_derivatives(pcur[0], pcur[1], d1u, d1v, d1uu, d1uv, d1vv);
                    Point3 n1 = d1u.cross_product(d1v);

                    Point3 d2u, d2v, d2uu, d2uv, d2vv;
                    surface2.evaluate_derivatives(pcur[2], pcur[3], d2u, d2v, d2uu, d2uv, d2vv);
                    Point3 n2 = d2u.cross_product(d2v);

                    Point3 t3 = n1.cross_product(n2);
                    double tn = t3.norm();
                    if (tn < 1e-12) break;
                    t3 = t3 / tn;

                    Point3 pPred3 = Pcur + t3 * (step * sign);

                    double u1p, v1p, u2p, v2p;
                    u1p = v1p = u2p = v2p = 0.5;
                    point_surface_distance(pPred3, surface1, u1p, v1p, opt.geomTol);
                    point_surface_distance(pPred3, surface2, u2p, v2p, opt.geomTol);

                    double pTrial[4] = { u1p, v1p, u2p, v2p };
                    bool ok = lm_refine_seed(surface1, surface2, pTrial, opt);
                    if (!ok)
                    {
                        step *= 0.5;
                        if (step < minStep) break;
                        continue;
                    }

                    Point3 pa2, pb2;
                    surface1.evaluate(pTrial[0], pTrial[1], pa2);
                    surface2.evaluate(pTrial[2], pTrial[3], pb2);
                    Point3 Pnew = (pa2 + pb2) * 0.5;

                    NurbsIntersectionSample sample;
                    sample.point = Pnew;
                    sample.uA = pTrial[0]; sample.vA = pTrial[1]; sample.uB = pTrial[2]; sample.vB = pTrial[3];
                    if (sign > 0)
                        outCurve.samples.push_back(sample);
                    else
                        outCurve.samples.insert(outCurve.samples.begin(), sample);

                    for (int i = 0; i < 4; ++i) pcur[i] = pTrial[i];
                    Pcur = Pnew;
                    step = std::min(step * 1.2, 0.1);
                }
            };

            outCurve.samples.clear();
            NurbsIntersectionSample seedS;
            seedS.point = seedPoint;
            seedS.uA = seedP[0]; seedS.vA = seedP[1]; seedS.uB = seedP[2]; seedS.vB = seedP[3];
            outCurve.samples.push_back(seedS);

            trace_dir(+1);
            trace_dir(-1);
        };

        trace_from_seed_fn(pseed, seedPoint, options, curve);
        curve.closed = false;
        if (curve.samples.size() >= 2)
            result.curves.push_back(curve);
    }

    // If no seeds found, leave curves empty; partial overlap detection kept false for now
    result.hasPartialOverlap = false;
}

// Trace curve from a refined seed using 3D tangent predictor and LM corrector
static void trace_from_seed(const NurbsSurface& s1, const NurbsSurface& s2,
                            const double seedP[4], const Point3& seedPoint,
                            const IntersectionOptions& opt, NurbsIntersectionCurve& outCurve)
{
    const double stepInit = opt.initialStep;
    const double minStep = opt.initialStep * 1e-3;
    const int maxSteps = 1000;

    auto trace_dir = [&](int sign)
    {
        double pcur[4] = { seedP[0], seedP[1], seedP[2], seedP[3] };
        Point3 Pcur = seedPoint;
        double step = stepInit;

        for (int iter = 0; iter < maxSteps; ++iter)
        {
            // compute normals
            Point3 d1u, d1v, d1uu, d1uv, d1vv;
            s1.evaluate_derivatives(pcur[0], pcur[1], d1u, d1v, d1uu, d1uv, d1vv);
            Point3 n1 = d1u.cross_product(d1v);

            Point3 d2u, d2v, d2uu, d2uv, d2vv;
            s2.evaluate_derivatives(pcur[2], pcur[3], d2u, d2v, d2uu, d2uv, d2vv);
            Point3 n2 = d2u.cross_product(d2v);

            Point3 t3 = n1.cross_product(n2);
            double tn = t3.norm();
            if (tn < 1e-12) break;
            t3 = t3 / tn;

            Point3 pPred3 = Pcur + t3 * (step * sign);

            // project predicted 3D point onto both surfaces (coarse)
            double u1p, v1p, u2p, v2p;
            u1p = v1p = u2p = v2p = 0.5;
            point_surface_distance(pPred3, s1, u1p, v1p, opt.geomTol);
            point_surface_distance(pPred3, s2, u2p, v2p, opt.geomTol);

            double pTrial[4] = { u1p, v1p, u2p, v2p };
            bool ok = lm_refine_seed(s1, s2, pTrial, opt);
            if (!ok)
            {
                // reduce step and retry
                step *= 0.5;
                if (step < minStep) break;
                continue;
            }

            // refined point
            Point3 pa, pb;
            s1.evaluate(pTrial[0], pTrial[1], pa);
            s2.evaluate(pTrial[2], pTrial[3], pb);
            Point3 Pnew = (pa + pb) * 0.5;

            // append sample
            NurbsIntersectionSample sample;
            sample.point = Pnew;
            sample.uA = pTrial[0]; sample.vA = pTrial[1]; sample.uB = pTrial[2]; sample.vB = pTrial[3];
            if (sign > 0)
                outCurve.samples.push_back(sample);
            else
                outCurve.samples.insert(outCurve.samples.begin(), sample);

            // advance
            for (int i = 0; i < 4; ++i) pcur[i] = pTrial[i];
            Pcur = Pnew;
            // increase step mildly
            step = std::min(step * 1.2, 0.1);
        }
    };

    // initialize curve with seed
    NurbsIntersectionSample seedS;
    seedS.point = seedPoint;
    seedS.uA = seedP[0]; seedS.vA = seedP[1]; seedS.uB = seedP[2]; seedS.vB = seedP[3];
    outCurve.samples.clear();
    outCurve.samples.push_back(seedS);

    // trace forward and backward
    trace_dir(+1);
    trace_dir(-1);
}
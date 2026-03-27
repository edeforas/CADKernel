#include "MeshBoolean.h"
#include "MeshUtil.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>

namespace
{
	Point3 triangle_center(const Triangle3& t)
	{
		return (t.p1() + t.p2() + t.p3()) / 3.;
	}

	void mesh_extent(const Mesh& m, Point3& pMin, Point3& pMax)
	{
		assert(m.nb_vertices() > 0);
		m.get_vertex(0, pMin);
		pMax = pMin;

		for (int i = 1; i < m.nb_vertices(); ++i)
		{
			Point3 p;
			m.get_vertex(i, p);

			if (p.x() < pMin.x()) pMin.x() = p.x();
			if (p.y() < pMin.y()) pMin.y() = p.y();
			if (p.z() < pMin.z()) pMin.z() = p.z();

			if (p.x() > pMax.x()) pMax.x() = p.x();
			if (p.y() > pMax.y()) pMax.y() = p.y();
			if (p.z() > pMax.z()) pMax.z() = p.z();
		}
	}

	double mesh_diagonal(const Mesh& m)
	{
		if (m.nb_vertices() == 0)
			return 1.;

		Point3 pMin, pMax;
		mesh_extent(m, pMin, pMax);
		return (pMax - pMin).norm();
	}

	bool point_on_mesh_surface(const Mesh& m, const Point3& p, double epsilon)
	{
		for (int i = 0; i < m.nb_triangles(); ++i)
		{
			if (m.is_triangle_unlinked(i))
				continue;

			Triangle3 t;
			m.get_triangle(i, t);
			Plane3 plane(t);
			if (std::fabs(plane.distance_to(p)) > epsilon)
				continue;

			if (t.contains(p))
				return true;
		}

		return false;
	}

	bool point_inside_closed_mesh_with_dir(const Mesh& m, const Point3& p, const Point3& dir, double rayLength, double epsilon, bool& bAmbiguous)
	{
		bAmbiguous = false;
		Point3 pEnd = p + dir * rayLength;
		Segment3 ray(p, pEnd);
		std::vector<double> hitDistances;

		for (int i = 0; i < m.nb_triangles(); ++i)
		{
			if (m.is_triangle_unlinked(i))
				continue;

			Triangle3 t;
			m.get_triangle(i, t);

			Point3 hit;
			if (!t.intersect_with(ray, hit))
				continue;

			double distance = (hit - p).norm();
			if (distance <= epsilon)
			{
				bAmbiguous = true;
				return false;
			}

			hitDistances.push_back(distance);
		}

		if (hitDistances.empty())
			return false;

		std::sort(hitDistances.begin(), hitDistances.end());
		int iUniqueHits = 0;
		for (size_t i = 0; i < hitDistances.size(); ++i)
		{
			if ((i == 0) || (std::fabs(hitDistances[i] - hitDistances[i - 1]) > epsilon))
				iUniqueHits++;
		}

		return (iUniqueHits % 2) == 1;
	}

	bool point_inside_closed_mesh(const Mesh& m, const Point3& p)
	{
		const double diagonal = mesh_diagonal(m);
		const double epsilon = std::max(1.e-7, diagonal * 1.e-7);

		if (point_on_mesh_surface(m, p, epsilon))
			return true;

		double rayLength = (diagonal > 1.e-8 ? diagonal : 1.) * 10.;

		std::vector<Point3> rayDirs;
		rayDirs.push_back(Point3(1., 0.3125, 0.125).normalized());
		rayDirs.push_back(Point3(0.125, 1., 0.4375).normalized());
		rayDirs.push_back(Point3(0.3125, 0.0625, 1.).normalized());

		for (const Point3& dir : rayDirs)
		{
			bool bAmbiguous = false;
			bool bInside = point_inside_closed_mesh_with_dir(m, p, dir, rayLength, epsilon, bAmbiguous);
			if (!bAmbiguous)
				return bInside;
		}

		return false;
	}

	void add_triangle_to_mesh(const Triangle3& t, Mesh& out)
	{
		out.add_triangle(t.p1(), t.p2(), t.p3());
	}

	void cleanup_boolean_mesh(Mesh& m)
	{
		MeshUtil::merge_vertices(m, 1.e-6);
		MeshUtil::remove_empty_triangles(m, 1.e-8);

		std::set<std::vector<int>> seenTriangles;
		for (int i = 0; i < m.nb_triangles(); ++i)
		{
			if (m.is_triangle_unlinked(i))
				continue;

			int i1, i2, i3;
			m.get_triangle_vertices(i, i1, i2, i3);

			std::vector<int> key;
			key.push_back(i1);
			key.push_back(i2);
			key.push_back(i3);
			std::sort(key.begin(), key.end());

			if (seenTriangles.find(key) != seenTriangles.end())
			{
				m.unlink_triangle(i);
				continue;
			}

			seenTriangles.insert(key);
		}
	}

	void remove_frontier_triangles(Mesh& out, const Mesh& A, const Mesh& B)
	{
		const double diagonal = std::max(mesh_diagonal(A), mesh_diagonal(B));
		const double epsilon = std::max(1.e-7, diagonal * 1.e-5);

		for (int i = 0; i < out.nb_triangles(); ++i)
		{
			if (out.is_triangle_unlinked(i))
				continue;

			Triangle3 t;
			out.get_triangle(i, t);
			Point3 c = triangle_center(t);
			Point3 n = t.normal();

			bool bFrontierA = false;
			bool bFrontierB = false;
			if (n.norm_square() > 1.e-20)
			{
				const Point3 cPlus = c + n * epsilon;
				const Point3 cMinus = c - n * epsilon;

				const bool bAPlus = point_inside_closed_mesh(A, cPlus);
				const bool bAMinus = point_inside_closed_mesh(A, cMinus);
				const bool bBPlus = point_inside_closed_mesh(B, cPlus);
				const bool bBMinus = point_inside_closed_mesh(B, cMinus);

				bFrontierA = (bAPlus != bAMinus);
				bFrontierB = (bBPlus != bBMinus);
			}

			if ((point_on_mesh_surface(A, c, epsilon) && point_on_mesh_surface(B, c, epsilon)) || (bFrontierA && bFrontierB))
				out.unlink_triangle(i);
		}
	}

	void remove_triangles_outside_reference(Mesh& out, const Mesh& reference)
	{
		for (int i = 0; i < out.nb_triangles(); ++i)
		{
			if (out.is_triangle_unlinked(i))
				continue;

			Triangle3 t;
			out.get_triangle(i, t);
			if (!point_inside_closed_mesh(reference, triangle_center(t)))
				out.unlink_triangle(i);
		}
	}

	void split_mesh_by_other(Mesh& meshToSplit, const Mesh& splitterMesh)
	{
		const int iNbPasses = 2;
		for (int iPass = 0; iPass < iNbPasses; ++iPass)
		{
			for (int j = 0; j < splitterMesh.nb_triangles(); ++j)
			{
				if (splitterMesh.is_triangle_unlinked(j))
					continue;

				Triangle3 tSplitter;
				splitterMesh.get_triangle(j, tSplitter);
				BoundingBox3 bboxSplitter(tSplitter);

				int iNbTrianglesSnapshot = meshToSplit.nb_triangles();
				for (int iTriangle = 0; iTriangle < iNbTrianglesSnapshot; ++iTriangle)
				{
					if (meshToSplit.is_triangle_unlinked(iTriangle))
						continue;

					Triangle3 tCurrent;
					meshToSplit.get_triangle(iTriangle, tCurrent);
					BoundingBox3 bboxCurrent(tCurrent);
					if (!bboxCurrent.intersect_with(bboxSplitter))
						continue;

					meshToSplit.split_triangle(iTriangle, tSplitter);
				}
			}
		}
	}

	void add_mesh_with_flipped_triangles(const Mesh& in, Mesh& out)
	{
		for (int i = 0; i < in.nb_triangles(); ++i)
		{
			if (in.is_triangle_unlinked(i))
				continue;

			Triangle3 t;
			in.get_triangle(i, t);
			out.add_triangle(t.p1(), t.p3(), t.p2());
		}
	}

	void classify_triangles(Mesh& meshToClassify, const Mesh& meshReference, Mesh& outInside)
	{
		const double diagonal = mesh_diagonal(meshReference);
		const double surfaceEpsilon = std::max(1.e-7, diagonal * 1.e-6);

		for (int i = 0; i < meshToClassify.nb_triangles(); ++i)
		{
			if (meshToClassify.is_triangle_unlinked(i))
				continue;

			Triangle3 t;
			meshToClassify.get_triangle(i, t);

			const Point3 center = triangle_center(t);
			bool bInside = point_inside_closed_mesh(meshReference, center);

			if (!bInside)
			{
				if (point_on_mesh_surface(meshReference, center, surfaceEpsilon))
				{
					bInside = true;
				}
				else
				{
					Point3 normal = t.normal();
					if (normal.norm_square() > 1.e-20)
					{
						const Point3 p1 = center + normal * surfaceEpsilon;
						const Point3 p2 = center - normal * surfaceEpsilon;
						const bool b1 = point_inside_closed_mesh(meshReference, p1);
						const bool b2 = point_inside_closed_mesh(meshReference, p2);

						if (b1 != b2)
							bInside = true;
					}
				}
			}

			if (!bInside)
				continue;

			meshToClassify.unlink_triangle(i);
			add_triangle_to_mesh(t, outInside);
		}
	}
}

MeshBoolean::MeshBoolean()
{
}

void MeshBoolean::compute_split(const Mesh& A, const Mesh& B, Mesh& Aoutside, Mesh& Boutside, Mesh& AInB, Mesh& BInA)
{
	AInB.clear();
	BInA.clear();

	Aoutside = A;
	Boutside = B;

	for (int i = 0; i < A.nb_triangles(); i++)
	{
		if (A.is_triangle_unlinked(i))
			continue;

		Triangle3 tA;
		A.get_triangle(i, tA);
		BoundingBox3 bboxA(tA);

		for (int j = 0; j < B.nb_triangles(); j++)
		{
			if (B.is_triangle_unlinked(j))
				continue;

			Triangle3 tB;
			B.get_triangle(j, tB);
			BoundingBox3 bboxB(tB);

			if (bboxA.intersect_with(bboxB) == false)
				continue;

			Aoutside.split_triangle(i, tB);
			Boutside.split_triangle(j, tA);
		}
	}

	classify_triangles(Aoutside, B, AInB);
	classify_triangles(Boutside, A, BInA);
}

void MeshBoolean::compute_intersection(const Mesh& A, const Mesh& B, Mesh& intersection)
{
	Mesh Aoutside, Boutside, AInB, BInA;
	compute_split(A, B, Aoutside, Boutside, AInB, BInA);

	intersection.clear();
	intersection.add_mesh(AInB);
	intersection.add_mesh(BInA);
	remove_frontier_triangles(intersection, A, B);
	cleanup_boolean_mesh(intersection);
}

void MeshBoolean::compute_union(const Mesh& A, const Mesh& B, Mesh& meshUnion)
{
	Mesh Aoutside, Boutside, AInB, BInA;
	compute_split(A, B, Aoutside, Boutside, AInB, BInA);

	meshUnion.clear();
	meshUnion.add_mesh(Aoutside);
	meshUnion.add_mesh(Boutside);
	remove_frontier_triangles(meshUnion, A, B);
	cleanup_boolean_mesh(meshUnion);
}

void MeshBoolean::compute_difference(const Mesh& A, const Mesh& B, Mesh& meshDifference)
{
	Mesh Aoutside, Boutside, AInB, BInA;
	compute_split(A, B, Aoutside, Boutside, AInB, BInA);

	meshDifference.clear();
	meshDifference.add_mesh(Aoutside);
	add_mesh_with_flipped_triangles(BInA, meshDifference);
	remove_frontier_triangles(meshDifference, A, B);
	remove_triangles_outside_reference(meshDifference, A);
	cleanup_boolean_mesh(meshDifference);
}
///////////////////////////////////////////////////////////////////////////
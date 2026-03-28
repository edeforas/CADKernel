#include "StepFile.h"

#include "NurbsSolid.h"
#include "NurbsSurface.h"
#include "NurbsCurve.h"
#include "NurbsTrimmedSurface.h"
#include "Geometry.h"
#include "NurbsUtil.h"

#include <algorithm>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <typeinfo>

namespace
{
	double sanitize_double(double value)
	{
		if (!std::isfinite(value))
			return 0.;
		return value;
	}

	std::vector<double> sanitize_weights_step(const std::vector<double>& input, int expectedSize)
	{
		std::vector<double> weights;
		weights.reserve(expectedSize);

		if ((int)input.size() != expectedSize)
		{
			weights.assign(expectedSize, 1.);
			return weights;
		}

		for (int i = 0; i < expectedSize; ++i)
		{
			double w = std::fabs(sanitize_double(input[i]));
			if (w < 1.e-12)
				w = 1.e-12;
			weights.push_back(w);
		}

		return weights;
	}
}


StepWriter::StepWriter() :
	_iItemIndent(10),
	_iGeomContextId(-1),
	_iProductDefShapeId(-1),
	_repKind(REP_NONE)
{
}

StepWriter::~StepWriter()
{
	close();
}

bool StepWriter::is_open()
{
	return _f.is_open();
}

void StepWriter::open(const string& filename)
{
	close();

	_sNameFile = filename;
	_f.open(filename.c_str(), std::ios::out | std::ios::trunc);
	if (!_f.is_open())
		return;

	_repKind = REP_NONE;
	_repItems.clear();

	_f << std::setprecision(15);
	write_header();
}

void StepWriter::close()
{
	if (!_f.is_open())
		return;

	flush_representation();
	write_footer();
	_f.close();
}

int StepWriter::next_id()
{
	return _iItemIndent++;
}



void StepWriter::write_header()
{
	_f << "ISO-10303-21;" << endl;
	_f << "HEADER;" << endl;
	_f << "FILE_DESCRIPTION(('NURBS export from CADKernel'),'2;1');" << endl;
	_f << "FILE_NAME('" << _sNameFile << "','2026-03-02T00:00:00',('CADKernel'),('CADKernel'),'CADKernel StepWriter','CADKernel','');" << endl;
	_f << "FILE_SCHEMA(('AUTOMOTIVE_DESIGN_CC2'));" << endl;
	_f << "ENDSEC;" << endl;
	_f << "DATA;" << endl;

	_iItemIndent = 10;

	const int originId = next_id();
	_f << "#" << originId << "=CARTESIAN_POINT('',(0.,0.,0.));" << endl;

	const int axisZId = next_id();
	_f << "#" << axisZId << "=DIRECTION('',(0.,0.,1.));" << endl;

	const int axisXId = next_id();
	_f << "#" << axisXId << "=DIRECTION('',(1.,0.,0.));" << endl;

	const int placementId = next_id();
	_f << "#" << placementId << "=AXIS2_PLACEMENT_3D('',#" << originId << ",#" << axisZId << ",#" << axisXId << ");" << endl;

	const int lengthUnitId = next_id();
	_f << "#" << lengthUnitId << "=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));" << endl;

	const int planeAngleUnitId = next_id();
	_f << "#" << planeAngleUnitId << "=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));" << endl;

	const int solidAngleUnitId = next_id();
	_f << "#" << solidAngleUnitId << "=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());" << endl;

	const int uncertaintyId = next_id();
	_f << "#" << uncertaintyId << "=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-6),#" << lengthUnitId << ",'distance_accuracy_value','confusion accuracy');" << endl;

	_iGeomContextId = next_id();
	_f << "#" << _iGeomContextId << "=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#" << uncertaintyId << "))GLOBAL_UNIT_ASSIGNED_CONTEXT((#" << lengthUnitId << ",#" << planeAngleUnitId << ",#" << solidAngleUnitId << "))REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'));" << endl;

	const int appContextId = next_id();
	_f << "#" << appContextId << "=APPLICATION_CONTEXT('automotive_design');" << endl;

	const int appProtocolId = next_id();
	_f << "#" << appProtocolId << "=APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#" << appContextId << ");" << endl;

	const int productContextId = next_id();
	_f << "#" << productContextId << "=PRODUCT_CONTEXT('',#" << appContextId << ",'mechanical');" << endl;

	const int productId = next_id();
	_f << "#" << productId << "=PRODUCT('CADKernelModel','CADKernelModel','',(#" << productContextId << "));" << endl;

	const int productFormationId = next_id();
	_f << "#" << productFormationId << "=PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#" << productId << ",.NOT_KNOWN.);" << endl;

	const int productDefContextId = next_id();
	_f << "#" << productDefContextId << "=PRODUCT_DEFINITION_CONTEXT('part definition',#" << appContextId << ",'design');" << endl;

	const int productDefId = next_id();
	_f << "#" << productDefId << "=PRODUCT_DEFINITION('','',#" << productFormationId << ",#" << productDefContextId << ");" << endl;

	_iProductDefShapeId = next_id();
	_f << "#" << _iProductDefShapeId << "=PRODUCT_DEFINITION_SHAPE('','',#" << productDefId << ");" << endl;
}

void StepWriter::write_footer()
{
	_f << "ENDSEC;" << endl;
	_f << "END-ISO-10303-21;" << endl;
}

void StepWriter::queue_representation_item(int itemId, RepresentationKind kind)
{
	if (itemId <= 0)
		return;

	if (_repItems.empty())
		_repKind = kind;
	else if (_repKind != kind)
		_repKind = REP_SHAPE;

	_repItems.push_back(itemId);
}

void StepWriter::flush_representation()
{
	if (!_f.is_open())
		return;
	if (_repItems.empty())
		return;
	if ((_iGeomContextId < 0) || (_iProductDefShapeId < 0))
		return;

	const int repId = next_id();
	if (_repKind == REP_ADVANCED_BREP)
		_f << "#" << repId << "=ADVANCED_BREP_SHAPE_REPRESENTATION('',(";
	else if (_repKind == REP_MANIFOLD_SURFACE)
		_f << "#" << repId << "=MANIFOLD_SURFACE_SHAPE_REPRESENTATION('',(";
	else
		_f << "#" << repId << "=SHAPE_REPRESENTATION('',(";

	for (int i = 0; i < (int)_repItems.size(); ++i)
	{
		_f << "#" << _repItems[i];
		if (i + 1 < (int)_repItems.size())
			_f << ",";
	}
	_f << "),#" << _iGeomContextId << ");" << endl;

	const int sdrId = next_id();
	_f << "#" << sdrId << "=SHAPE_DEFINITION_REPRESENTATION(#" << _iProductDefShapeId << ",#" << repId << ");" << endl;

	_repItems.clear();
	_repKind = REP_NONE;
}

void StepWriter::write_cartesian_point(const Point3& p)
{
	if (!_f.is_open())
		return;

	const int id = next_id();
	_f << std::uppercase << "#" << id << "=CARTESIAN_POINT('',("
		<< sanitize_double(p.x()) << ","
		<< sanitize_double(p.y()) << ","
		<< sanitize_double(p.z()) << "));" << endl;
}

int StepWriter::write_surface_entity(const NurbsSurface& n)
{
	if (!_f.is_open())
		return -1;

	int nU = n.nb_points_u();
	int nV = n.nb_points_v();
	std::vector<Point3> points = n.points();
	if ((nU <= 0) || (nV <= 0) || ((int)points.size() != nU * nV))
		return -1;

	int degreeU = n.degree_u();
	int degreeV = n.degree_v();
	if (degreeU < 0) degreeU = 0;
	if (degreeV < 0) degreeV = 0;
	if (degreeU >= nU) degreeU = nU - 1;
	if (degreeV >= nV) degreeV = nV - 1;

	std::vector<double> weights = sanitize_weights_step(n.weights(), nU * nV);
	const NurbsUtil::KnotAnalysis uDataOrig = NurbsUtil::analyze_knots(n.knots_u(), degreeU, nU);
	const NurbsUtil::KnotAnalysis vDataOrig = NurbsUtil::analyze_knots(n.knots_v(), degreeV, nV);
	if (uDataOrig.unique_knots.empty() || vDataOrig.unique_knots.empty())
		return -1;

	// Handle closed surfaces by duplicating control points
	NurbsUtil::KnotAnalysis uData = uDataOrig;
	NurbsUtil::KnotAnalysis vData = vDataOrig;
	bool closedU = n.is_closed_u();
	bool closedV = n.is_closed_v();

	if (closedU && !closedV) {
		// Adjust knots to non-periodic
		if (uData.unique_knots.back() < 1.0 - 1e-6) {
			uData.unique_knots.push_back(1.0);
			uData.multiplicities.push_back(degreeU + 1);
		}
		uData.multiplicities[0] = degreeU + 1;
		uData.multiplicities.back() = degreeU + 1;
		closedU = false;
	}

	// For closedV, keep as is for now, since torus worked with it

	const int firstPointId = _iItemIndent;
	for (const auto& p : points)
		write_cartesian_point(p);

	const int surfaceId = next_id();
	_f << "#" << surfaceId << "=(BOUNDED_SURFACE()B_SPLINE_SURFACE(" << degreeU << "," << degreeV << ",(" << endl;

	for (int v = 0; v < nV; ++v)
	{
		_f << "(";
		for (int u = 0; u < nU; ++u)
		{
			_f << "#" << (firstPointId + v * nU + u);
			if (u + 1 < nU)
				_f << ",";
		}
		_f << ")";
		if (v + 1 < nV)
			_f << ",";
		_f << endl;
	}

	_f << "),.UNSPECIFIED.," << (closedU ? ".T." : ".F.") << "," << (closedV ? ".T." : ".F.") << ",.F.)" << endl;

	_f << "B_SPLINE_SURFACE_WITH_KNOTS((";
	for (int i = 0; i < (int)uData.multiplicities.size(); ++i)
	{
		_f << uData.multiplicities[i];
		if (i + 1 < (int)uData.multiplicities.size())
			_f << ",";
	}

	_f << "),(";
	for (int i = 0; i < (int)vData.multiplicities.size(); ++i)
	{
		_f << vData.multiplicities[i];
		if (i + 1 < (int)vData.multiplicities.size())
			_f << ",";
	}

	_f << "),(";
	for (int i = 0; i < (int)uData.unique_knots.size(); ++i)
	{
		_f << sanitize_double(uData.unique_knots[i]);
		if (i + 1 < (int)uData.unique_knots.size())
			_f << ",";
	}

	_f << "),(";
	for (int i = 0; i < (int)vData.unique_knots.size(); ++i)
	{
		_f << sanitize_double(vData.unique_knots[i]);
		if (i + 1 < (int)vData.unique_knots.size())
			_f << ",";
	}
	_f << "),.UNSPECIFIED.)" << endl;
	_f << "GEOMETRIC_REPRESENTATION_ITEM()" << endl;

	_f << "RATIONAL_B_SPLINE_SURFACE((" << endl;
	for (int v = 0; v < nV; ++v)
	{
		_f << "(";
		for (int u = 0; u < nU; ++u)
		{
			_f << sanitize_double(weights[v * nU + u]);
			if (u + 1 < nU)
				_f << ",";
		}
		_f << ")";
		if (v + 1 < nV)
			_f << ",";
		_f << endl;
	}

	_f << "))REPRESENTATION_ITEM('')SURFACE());" << endl;
	return surfaceId;
}

int StepWriter::write_curve_entity(const NurbsCurve& n)
{
	if (!_f.is_open())
		return -1;

	int nP = n.nb_points();
	const std::vector<Point3>& points = n.points();
	if ((nP <= 0) || ((int)points.size() != nP))
		return -1;

	int degree = n.degree();
	if (degree < 0) degree = 0;
	if (degree >= nP) degree = nP - 1;

	const std::vector<double> weights = sanitize_weights_step(n.weights(), nP);
	const NurbsUtil::KnotAnalysis knotData = NurbsUtil::analyze_knots(n.knots(), degree, nP);
	if (knotData.unique_knots.empty())
		return -1;

	const int firstPointId = _iItemIndent;
	for (const auto& p : points)
		write_cartesian_point(p);

	const int curveId = next_id();
	_f << "#" << curveId << "=(BOUNDED_CURVE()B_SPLINE_CURVE(" << degree << ",(" << endl;

	for (int i = 0; i < nP; ++i)
	{
		_f << "#" << (firstPointId + i);
		if (i + 1 < nP)
			_f << ",";
	}
	_f << "),.UNSPECIFIED.," << (n.is_closed() ? ".T." : ".F.") << ",.F.)" << endl;

	_f << "B_SPLINE_CURVE_WITH_KNOTS((";
	for (int i = 0; i < (int)knotData.multiplicities.size(); ++i)
	{
		_f << knotData.multiplicities[i];
		if (i + 1 < (int)knotData.multiplicities.size())
			_f << ",";
	}

	_f << "),(";
	for (int i = 0; i < (int)knotData.unique_knots.size(); ++i)
	{
		_f << sanitize_double(knotData.unique_knots[i]);
		if (i + 1 < (int)knotData.unique_knots.size())
			_f << ",";
	}
	_f << "),.UNSPECIFIED.)" << endl;
	_f << "GEOMETRIC_REPRESENTATION_ITEM()" << endl;

	_f << "RATIONAL_B_SPLINE_CURVE((";
	for (int i = 0; i < nP; ++i)
	{
		_f << sanitize_double(weights[i]);
		if (i + 1 < nP)
			_f << ",";
	}
	_f << "))REPRESENTATION_ITEM('')CURVE());" << endl;
	return curveId;
}

int StepWriter::write_advanced_face(const NurbsSurface& n)
{
	const int surfaceId = write_surface_entity(n);
	if (surfaceId < 0)
		return -1;

	const int faceId = next_id();
	const NurbsTrimmedSurface* trimmed = dynamic_cast<const NurbsTrimmedSurface*>(&n);
	if (trimmed && !trimmed->trim_loops().empty()) {
		std::vector<int> boundIds;
		bool hasOuter = false;
		for (const auto& loop : trimmed->trim_loops())
		{
			if (loop.points.size() < 3)
				continue;

			std::vector<Point3> loop3d;
			loop3d.reserve(loop.points.size());
			for (const auto& uv : loop.points)
			{
				Point3 p;
				n.evaluate(uv.u, uv.v, p);
				loop3d.push_back(p);
			}

			const bool bHole = loop.hole || hasOuter;
			const int boundId = write_trim_loop_bound(loop3d, bHole);
			if (boundId > 0)
			{
				boundIds.push_back(boundId);
				if (!bHole)
					hasOuter = true;
			}
		}
		_f << "#" << faceId << "=ADVANCED_FACE('',(";
		for (int i = 0; i < (int)boundIds.size(); ++i)
		{
			_f << "#" << boundIds[i];
			if (i + 1 < (int)boundIds.size())
				_f << ",";
		}
		_f << "),#" << surfaceId << ",.T.);" << endl;
	} else {
		_f << "#" << faceId << "=ADVANCED_FACE('',(),#" << surfaceId << ",.T.);" << endl;
	}
	return faceId;
}

namespace
{
	bool almost_same_point(const Point3& a, const Point3& b)
	{
		return (a - b).norm_square() <= 1.e-24;
	}
}

int StepWriter::write_trim_loop_bound(const std::vector<Point3>& rawLoop, bool bHole)
{
	if (!_f.is_open())
		return -1;

	std::vector<Point3> loop;
	loop.reserve(rawLoop.size());
	for (int i = 0; i < (int)rawLoop.size(); ++i)
	{
		if (!loop.empty() && almost_same_point(loop.back(), rawLoop[i]))
			continue;
		loop.push_back(rawLoop[i]);
	}

	if (loop.size() < 3)
		return -1;

	if (!almost_same_point(loop.front(), loop.back()))
		loop.push_back(loop.front());

	if (loop.size() < 4)
		return -1;

	const int vertexCount = (int)loop.size() - 1;
	if (vertexCount < 3)
		return -1;

	std::vector<int> vertexPointIds(vertexCount, -1);
	for (int i = 0; i < vertexCount; ++i)
	{
		const int pId = next_id();
		vertexPointIds[i] = pId;
		_f << "#" << pId << "=CARTESIAN_POINT('',("
			<< sanitize_double(loop[i].x()) << ","
			<< sanitize_double(loop[i].y()) << ","
			<< sanitize_double(loop[i].z()) << "));" << endl;
	}

	std::vector<int> vertexIds(vertexCount, -1);
	for (int i = 0; i < vertexCount; ++i)
	{
		const int vId = next_id();
		vertexIds[i] = vId;
		_f << "#" << vId << "=VERTEX_POINT('',#" << vertexPointIds[i] << ");" << endl;
	}

	std::vector<int> orientedEdgeIds;
	orientedEdgeIds.reserve(vertexCount);

	for (int i = 0; i < vertexCount; ++i)
	{
		const int iNext = (i + 1) % vertexCount;
		Point3 d = loop[iNext] - loop[i];
		double dNorm = d.norm();
		if (dNorm <= 1.e-12)
			continue;
		d /= dNorm;

		const int dirId = next_id();
		_f << "#" << dirId << "=DIRECTION('',("
			<< sanitize_double(d.x()) << ","
			<< sanitize_double(d.y()) << ","
			<< sanitize_double(d.z()) << "));" << endl;

		const int vecId = next_id();
		_f << "#" << vecId << "=VECTOR('',#" << dirId << ",1.);" << endl;

		const int lineId = next_id();
		_f << "#" << lineId << "=LINE('',#" << vertexPointIds[i] << ",#" << vecId << ");" << endl;

		const int edgeCurveId = next_id();
		_f << "#" << edgeCurveId << "=EDGE_CURVE('',#" << vertexIds[i] << ",#" << vertexIds[iNext] << ",#" << lineId << ",.T.);" << endl;

		const int orientedEdgeId = next_id();
		_f << "#" << orientedEdgeId << "=ORIENTED_EDGE('',*,*,#" << edgeCurveId << ",.T.);" << endl;

		orientedEdgeIds.push_back(orientedEdgeId);
	}

	if (orientedEdgeIds.size() < 3)
		return -1;

	const int loopId = next_id();
	_f << "#" << loopId << "=EDGE_LOOP('',(";
	for (int i = 0; i < (int)orientedEdgeIds.size(); ++i)
	{
		_f << "#" << orientedEdgeIds[i];
		if (i + 1 < (int)orientedEdgeIds.size())
			_f << ",";
	}
	_f << "));" << endl;

	const int boundId = next_id();
	if (bHole)
		_f << "#" << boundId << "=FACE_BOUND('',#" << loopId << ",.T.);" << endl;
	else
		_f << "#" << boundId << "=FACE_OUTER_BOUND('',#" << loopId << ",.T.);" << endl;

	return boundId;
}
////////////////////////////////////////////////////////////////////////////////////////////////
void StepWriter::write(const NurbsCurve& n)
{
	if (!_f.is_open())
		return;
	if ((_iGeomContextId < 0) || (_iProductDefShapeId < 0))
		return;

	const int curveId = write_curve_entity(n);
	if (curveId < 0)
		return;

	queue_representation_item(curveId, REP_SHAPE);
}
////////////////////////////////////////////////////////////////////////////////////////////////
void StepWriter::write(const NurbsSolid& n)
{
	if (!_f.is_open())
		return;
	if ((_iGeomContextId < 0) || (_iProductDefShapeId < 0))
		return;

	std::vector<int> faceIds;
	faceIds.reserve(n.surfaces().size());
	for (const auto& s : n.surfaces())
	{
		const int id = write_advanced_face(s);
		if (id > 0)
			faceIds.push_back(id);
	}

	if (faceIds.empty())
		return;

	const int shellId = next_id();
	_f << "#" << shellId << "=CLOSED_SHELL('',(";
	for (int i = 0; i < (int)faceIds.size(); ++i)
	{
		_f << "#" << faceIds[i];
		if (i + 1 < (int)faceIds.size())
			_f << ",";
	}
	_f << "));" << endl;

	const int brepId = next_id();
	_f << "#" << brepId << "=MANIFOLD_SOLID_BREP('',#" << shellId << ");" << endl;
	queue_representation_item(brepId, REP_ADVANCED_BREP);
}
////////////////////////////////////////////////////////////////////////////////////////////////
void StepWriter::write(const NurbsSurface& n)
{
	if (!_f.is_open())
		return;
	if ((_iGeomContextId < 0) || (_iProductDefShapeId < 0))
		return;

	const int surfaceId = write_surface_entity(n);
	if (surfaceId < 0)
		return;
	queue_representation_item(surfaceId, REP_SHAPE);
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
namespace
{
	class StepTextParser
	{
	public:
		StepTextParser(const std::string& text, size_t start) : _text(text), _i(start) {}

		void skip_spaces()
		{
			while (_i < _text.size() && std::isspace(static_cast<unsigned char>(_text[_i])))
				++_i;
		}

		bool expect(char c)
		{
			skip_spaces();
			if (_i >= _text.size() || _text[_i] != c)
				return false;
			++_i;
			return true;
		}

		bool parse_int(int& value)
		{
			skip_spaces();
			char* endPtr = NULL;
			long v = std::strtol(_text.c_str() + _i, &endPtr, 10);
			if (endPtr == _text.c_str() + _i)
				return false;
			_i = static_cast<size_t>(endPtr - _text.c_str());
			value = static_cast<int>(v);
			return true;
		}

		bool parse_double(double& value)
		{
			skip_spaces();
			char* endPtr = NULL;
			double v = std::strtod(_text.c_str() + _i, &endPtr);
			if (endPtr == _text.c_str() + _i)
				return false;
			_i = static_cast<size_t>(endPtr - _text.c_str());
			value = v;
			return true;
		}

		bool parse_ref(int& id)
		{
			skip_spaces();
			if (_i >= _text.size() || _text[_i] != '#')
				return false;
			++_i;
			return parse_int(id);
		}

		bool parse_step_bool(bool& value)
		{
			skip_spaces();
			if ((_i + 2) >= _text.size() || _text[_i] != '.' || _text[_i + 2] != '.')
				return false;

			if (_text[_i + 1] == 'T')
				value = true;
			else if (_text[_i + 1] == 'F')
				value = false;
			else
				return false;

			_i += 3;
			return true;
		}

		bool consume_until(char c)
		{
			size_t pos = _text.find(c, _i);
			if (pos == std::string::npos)
				return false;
			_i = pos;
			return true;
		}

		bool parse_int_list(std::vector<int>& values)
		{
			values.clear();
			if (!expect('('))
				return false;

			while (true)
			{
				int v = 0;
				if (!parse_int(v))
					return false;
				values.push_back(v);

				skip_spaces();
				if (_i >= _text.size())
					return false;
				if (_text[_i] == ',')
				{
					++_i;
					continue;
				}
				if (_text[_i] == ')')
				{
					++_i;
					break;
				}
				return false;
			}

			return true;
		}

		bool parse_ref_list(std::vector<int>& values)
		{
			values.clear();
			if (!expect('('))
				return false;

			while (true)
			{
				int id = 0;
				if (!parse_ref(id))
					return false;
				values.push_back(id);

				skip_spaces();
				if (_i >= _text.size())
					return false;
				if (_text[_i] == ',')
				{
					++_i;
					continue;
				}
				if (_text[_i] == ')')
				{
					++_i;
					break;
				}
				return false;
			}

			return true;
		}

		bool parse_double_list(std::vector<double>& values)
		{
			values.clear();
			if (!expect('('))
				return false;

			while (true)
			{
				double v = 0.;
				if (!parse_double(v))
					return false;
				values.push_back(v);

				skip_spaces();
				if (_i >= _text.size())
					return false;
				if (_text[_i] == ',')
				{
					++_i;
					continue;
				}
				if (_text[_i] == ')')
				{
					++_i;
					break;
				}
				return false;
			}

			return true;
		}

		bool parse_ref_grid(std::vector<std::vector<int>>& rows)
		{
			rows.clear();
			if (!expect('('))
				return false;

			while (true)
			{
				if (!expect('('))
					return false;

				std::vector<int> row;
				while (true)
				{
					int id = 0;
					if (!parse_ref(id))
						return false;
					row.push_back(id);

					skip_spaces();
					if (_i >= _text.size())
						return false;
					if (_text[_i] == ',')
					{
						++_i;
						continue;
					}
					if (_text[_i] == ')')
					{
						++_i;
						break;
					}
					return false;
				}

				rows.push_back(row);

				skip_spaces();
				if (_i >= _text.size())
					return false;
				if (_text[_i] == ',')
				{
					++_i;
					continue;
				}
				if (_text[_i] == ')')
				{
					++_i;
					break;
				}
				return false;
			}

			return true;
		}

		bool parse_double_grid(std::vector<std::vector<double>>& rows)
		{
			rows.clear();
			if (!expect('('))
				return false;

			while (true)
			{
				if (!expect('('))
					return false;

				std::vector<double> row;
				while (true)
				{
					double value = 0.;
					if (!parse_double(value))
						return false;
					row.push_back(value);

					skip_spaces();
					if (_i >= _text.size())
						return false;
					if (_text[_i] == ',')
					{
						++_i;
						continue;
					}
					if (_text[_i] == ')')
					{
						++_i;
						break;
					}
					return false;
				}

				rows.push_back(row);

				skip_spaces();
				if (_i >= _text.size())
					return false;
				if (_text[_i] == ',')
				{
					++_i;
					continue;
				}
				if (_text[_i] == ')')
				{
					++_i;
					break;
				}
				return false;
			}

			return true;
		}

	private:
		const std::string& _text;
		size_t _i;
	};

	std::vector<double> expand_knots(const std::vector<int>& multiplicities, const std::vector<double>& uniqueKnots)
	{
		std::vector<double> knots;
		if (multiplicities.size() != uniqueKnots.size())
			return knots;

		for (size_t i = 0; i < multiplicities.size(); ++i)
			for (int k = 0; k < multiplicities[i]; ++k)
				knots.push_back(uniqueKnots[i]);

		return knots;
	}

	bool parse_cartesian_points(const std::string& content, std::map<int, Point3>& pointsById)
	{
		pointsById.clear();

		size_t pos = 0;
		while (true)
		{
			size_t pointPos = content.find("=CARTESIAN_POINT", pos);
			if (pointPos == std::string::npos)
				break;

			size_t hashPos = content.rfind('#', pointPos);
			if (hashPos == std::string::npos)
				return false;

			char* endId = NULL;
			long id = std::strtol(content.c_str() + hashPos + 1, &endId, 10);
			if (endId == content.c_str() + hashPos + 1)
				return false;

			size_t tupleStart = content.find('(', pointPos);
			tupleStart = content.find('(', tupleStart + 1);
			if (tupleStart == std::string::npos)
				return false;

			size_t tupleEnd = content.find(')', tupleStart + 1);
			if (tupleEnd == std::string::npos)
				return false;

			std::string tuple = content.substr(tupleStart + 1, tupleEnd - tupleStart - 1);
			std::stringstream ss(tuple);
			std::string sx, sy, sz;
			if (!std::getline(ss, sx, ',') || !std::getline(ss, sy, ',') || !std::getline(ss, sz, ','))
				return false;

			double x = std::strtod(sx.c_str(), NULL);
			double y = std::strtod(sy.c_str(), NULL);
			double z = std::strtod(sz.c_str(), NULL);
			pointsById[(int)id] = Point3(x, y, z);

			pos = tupleEnd + 1;
		}

		return !pointsById.empty();
	}

	bool parse_entities(const std::string& content, std::map<int, std::string>& entities)
	{
		entities.clear();

		size_t pos = 0;
		while (true)
		{
			size_t hashPos = content.find('#', pos);
			if (hashPos == std::string::npos)
				break;

			char* endId = NULL;
			long id = std::strtol(content.c_str() + hashPos + 1, &endId, 10);
			if (endId == content.c_str() + hashPos + 1)
			{
				pos = hashPos + 1;
				continue;
			}

			size_t eqPos = content.find('=', static_cast<size_t>(endId - content.c_str()));
			if (eqPos == std::string::npos)
				return false;

			size_t semiPos = content.find(';', eqPos + 1);
			if (semiPos == std::string::npos)
				return false;

			entities[(int)id] = content.substr(eqPos + 1, semiPos - eqPos - 1);
			pos = semiPos + 1;
		}

		return !entities.empty();
	}

	bool parse_surface_entity(const std::string& entity, const std::map<int, Point3>& pointsById, NurbsSurface& n)
	{
		const size_t bsplinePos = entity.find("B_SPLINE_SURFACE(");
		if (bsplinePos == std::string::npos)
			return false;

		StepTextParser pSurface(entity, bsplinePos + std::string("B_SPLINE_SURFACE").size());
		if (!pSurface.expect('('))
			return false;

		int degreeU = 0;
		int degreeV = 0;
		if (!pSurface.parse_int(degreeU))
			return false;
		if (!pSurface.expect(','))
			return false;
		if (!pSurface.parse_int(degreeV))
			return false;
		if (!pSurface.expect(','))
			return false;

		std::vector<std::vector<int>> pointIdRows;
		if (!pSurface.parse_ref_grid(pointIdRows))
			return false;
		if (pointIdRows.empty() || pointIdRows[0].empty())
			return false;

		if (!pSurface.expect(','))
			return false;
		if (!pSurface.consume_until(','))
			return false;
		if (!pSurface.expect(','))
			return false;

		bool closedU = false;
		bool closedV = false;
		if (!pSurface.parse_step_bool(closedU))
			return false;
		if (!pSurface.expect(','))
			return false;
		if (!pSurface.parse_step_bool(closedV))
			return false;

		const size_t knotsPos = entity.find("B_SPLINE_SURFACE_WITH_KNOTS(", bsplinePos);
		if (knotsPos == std::string::npos)
			return false;

		StepTextParser pKnots(entity, knotsPos + std::string("B_SPLINE_SURFACE_WITH_KNOTS").size());
		if (!pKnots.expect('('))
			return false;

		std::vector<int> multU;
		std::vector<int> multV;
		std::vector<double> uniqueU;
		std::vector<double> uniqueV;
		if (!pKnots.parse_int_list(multU))
			return false;
		if (!pKnots.expect(','))
			return false;
		if (!pKnots.parse_int_list(multV))
			return false;
		if (!pKnots.expect(','))
			return false;
		if (!pKnots.parse_double_list(uniqueU))
			return false;
		if (!pKnots.expect(','))
			return false;
		if (!pKnots.parse_double_list(uniqueV))
			return false;

		const size_t weightsPos = entity.find("RATIONAL_B_SPLINE_SURFACE(", knotsPos);
		if (weightsPos == std::string::npos)
			return false;

		StepTextParser pWeights(entity, weightsPos + std::string("RATIONAL_B_SPLINE_SURFACE").size());
		if (!pWeights.expect('('))
			return false;

		std::vector<std::vector<double>> weightRows;
		if (!pWeights.parse_double_grid(weightRows))
			return false;

		const int nV = (int)pointIdRows.size();
		const int nU = (int)pointIdRows[0].size();
		for (int v = 1; v < nV; ++v)
			if ((int)pointIdRows[v].size() != nU)
				return false;

		if ((int)weightRows.size() != nV)
			return false;
		for (int v = 0; v < nV; ++v)
			if ((int)weightRows[v].size() != nU)
				return false;

		std::vector<Point3> points;
		points.reserve(nU * nV);
		std::vector<double> weights;
		weights.reserve(nU * nV);

		for (int v = 0; v < nV; ++v)
		{
			for (int u = 0; u < nU; ++u)
			{
				const int id = pointIdRows[v][u];
				std::map<int, Point3>::const_iterator it = pointsById.find(id);
				if (it == pointsById.end())
					return false;

				points.push_back(it->second);
				weights.push_back(weightRows[v][u]);
			}
		}

		std::vector<double> knotsU = expand_knots(multU, uniqueU);
		std::vector<double> knotsV = expand_knots(multV, uniqueV);
		if ((int)knotsU.size() != nU + degreeU + 1)
			return false;
		if ((int)knotsV.size() != nV + degreeV + 1)
			return false;

		n.clear();
		n.set_degree(degreeU, degreeV);
		n.set_points(points, nU, nV);
		n.set_knots_u(knotsU);
		n.set_knots_v(knotsV);
		n.set_weights(weights);
		n.set_closed_u(closedU);
		n.set_closed_v(closedV);

		return true;
	}

	bool parse_curve_entity(const std::string& entity, const std::map<int, Point3>& pointsById, NurbsCurve& n)
	{
		const size_t bsplinePos = entity.find("B_SPLINE_CURVE(");
		if (bsplinePos == std::string::npos)
			return false;

		StepTextParser pCurve(entity, bsplinePos + std::string("B_SPLINE_CURVE").size());
		if (!pCurve.expect('('))
			return false;

		int degree = 0;
		if (!pCurve.parse_int(degree))
			return false;
		if (!pCurve.expect(','))
			return false;

		std::vector<int> pointIds;
		if (!pCurve.parse_ref_list(pointIds))
			return false;
		if (pointIds.empty())
			return false;

		if (!pCurve.expect(','))
			return false;
		if (!pCurve.consume_until(','))
			return false;
		if (!pCurve.expect(','))
			return false;

		bool closed = false;
		if (!pCurve.parse_step_bool(closed))
			return false;

		const size_t knotsPos = entity.find("B_SPLINE_CURVE_WITH_KNOTS(", bsplinePos);
		if (knotsPos == std::string::npos)
			return false;

		StepTextParser pKnots(entity, knotsPos + std::string("B_SPLINE_CURVE_WITH_KNOTS").size());
		if (!pKnots.expect('('))
			return false;

		std::vector<int> mult;
		std::vector<double> unique;
		if (!pKnots.parse_int_list(mult))
			return false;
		if (!pKnots.expect(','))
			return false;
		if (!pKnots.parse_double_list(unique))
			return false;

		const size_t weightsPos = entity.find("RATIONAL_B_SPLINE_CURVE(", knotsPos);
		if (weightsPos == std::string::npos)
			return false;

		StepTextParser pWeights(entity, weightsPos + std::string("RATIONAL_B_SPLINE_CURVE").size());
		if (!pWeights.expect('('))
			return false;

		std::vector<double> weights;
		if (!pWeights.parse_double_list(weights))
			return false;

		const int nP = (int)pointIds.size();
		if ((int)weights.size() != nP)
			return false;

		std::vector<Point3> points;
		points.reserve(nP);

		for (int i = 0; i < nP; ++i)
		{
			const int id = pointIds[i];
			std::map<int, Point3>::const_iterator it = pointsById.find(id);
			if (it == pointsById.end())
				return false;

			points.push_back(it->second);
		}

		std::vector<double> knots = expand_knots(mult, unique);
		if ((int)knots.size() != nP + degree + 1)
			return false;

		n.clear();
		n.set_degree(degree);
		n.set_points(points);
		n.set_knots(knots);
		n.set_weights(weights);

		return true;
	}
}

StepReader::StepReader()
{
}

StepReader::~StepReader()
{
}

void StepReader::open(const string& filename)
{
	std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
	if (!in.is_open())
		throw std::runtime_error("Cannot open file: " + filename);

	std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	if (content.find("ISO-10303-21;") == std::string::npos)
		throw std::runtime_error("Not a valid STEP file: " + filename);

	_pointsById.clear();
	if (!parse_cartesian_points(content, _pointsById))
		throw std::runtime_error("Failed to parse cartesian points");

	_entities.clear();
	if (!parse_entities(content, _entities))
		throw std::runtime_error("Failed to parse entities");

	_currentEntity = _entities.begin();
}

bool StepReader::read(NurbsSurface& n)
{
	for (; _currentEntity != _entities.end(); ++_currentEntity)
	{
		if (_currentEntity->second.find("B_SPLINE_SURFACE(") != std::string::npos)
		{
			bool result = parse_surface_entity(_currentEntity->second, _pointsById, n);
			++_currentEntity;
			return result;
		}
	}
	return false;
}

bool StepReader::read(NurbsSolid& n)
{
	for (; _currentEntity != _entities.end(); ++_currentEntity)
	{
		const std::string& entity = _currentEntity->second;
		const size_t brepPos = entity.find("MANIFOLD_SOLID_BREP(");
		if (brepPos != std::string::npos)
		{
			int shellId = -1;
			StepTextParser pBrep(entity, brepPos + std::string("MANIFOLD_SOLID_BREP").size());
			if (!pBrep.expect('('))
				continue;
			if (!pBrep.consume_until(','))
				continue;
			if (!pBrep.expect(','))
				continue;
			if (!pBrep.parse_ref(shellId))
				continue;

			auto itShell = _entities.find(shellId);
			if (itShell == _entities.end())
				continue;

			const std::string& shellEntity = itShell->second;
			const size_t shellPos = shellEntity.find("CLOSED_SHELL(");
			if (shellPos == std::string::npos)
				continue;

			StepTextParser pShell(shellEntity, shellPos + std::string("CLOSED_SHELL").size());
			if (!pShell.expect('('))
				continue;
			if (!pShell.consume_until(','))
				continue;
			if (!pShell.expect(','))
				continue;

			std::vector<int> faceIds;
			if (!pShell.parse_ref_list(faceIds))
				continue;
			if (faceIds.empty())
				continue;

			n.clear();
			for (size_t i = 0; i < faceIds.size(); ++i)
			{
				auto itFace = _entities.find(faceIds[i]);
				if (itFace == _entities.end())
					continue;

				const std::string& faceEntity = itFace->second;
				const size_t facePos = faceEntity.find("ADVANCED_FACE(");
				if (facePos == std::string::npos)
					continue;

				StepTextParser pFace(faceEntity, facePos + std::string("ADVANCED_FACE").size());
				if (!pFace.expect('('))
					continue;
				if (!pFace.consume_until(','))
					continue;
				if (!pFace.expect(','))
					continue;
				if (!pFace.consume_until(','))
					continue;
				if (!pFace.expect(','))
					continue;

				int surfaceId = -1;
				if (!pFace.parse_ref(surfaceId))
					continue;

				auto itSurface = _entities.find(surfaceId);
				if (itSurface == _entities.end())
					continue;

				NurbsSurface surface;
				if (!parse_surface_entity(itSurface->second, _pointsById, surface))
					continue;

				n.add_surface(surface);
			}

			if (!n.surfaces().empty())
			{
				++_currentEntity;
				return true;
			}
		}
	}
	return false;
}

bool StepReader::read(NurbsCurve& n)
{
	for (; _currentEntity != _entities.end(); ++_currentEntity)
	{
		if (_currentEntity->second.find("B_SPLINE_CURVE(") != std::string::npos)
		{
			bool result = parse_curve_entity(_currentEntity->second, _pointsById, n);
			++_currentEntity;
			return result;
		}
	}
	return false;
}

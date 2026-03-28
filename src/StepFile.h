#ifndef _StepFile_
#define _StepFile_

#include <string>
#include <vector>
#include <fstream>
using namespace std;

class NurbsSurface;
class NurbsSolid;
class NurbsCurve;
class NurbsTrimmedSurface;
class Point3;

class StepWriter
{
public:
	StepWriter();
	virtual ~StepWriter();

	void open(const string& filename);
	bool is_open();
	void close();

	void write(const NurbsSurface& n);
	void write(const NurbsSolid& n);
	void write(const NurbsCurve& n);

private:
	enum RepresentationKind
	{
		REP_NONE,
		REP_SHAPE,
		REP_ADVANCED_BREP,
		REP_MANIFOLD_SURFACE
	};

	int next_id();
	void write_header();
	void write_footer();
	int write_surface_entity(const NurbsSurface& n);
	int write_curve_entity(const NurbsCurve& n);
	int write_advanced_face(const NurbsSurface& n);
	int write_trim_loop_bound(const vector<Point3>& loop, bool bHole);
	void write_cartesian_point(const Point3& p);
	void queue_representation_item(int itemId, RepresentationKind kind);
	void flush_representation();

	ofstream _f;
	string _sNameFile;
	int _iItemIndent;
	int _iGeomContextId;
	int _iProductDefShapeId;
	RepresentationKind _repKind;
	vector<int> _repItems;
};

class StepReader
{
public:
	StepReader();
	virtual ~StepReader();

	void open(const string& filename);
	bool read(NurbsSurface& n);
	bool read(NurbsSolid& n);
	bool read(NurbsCurve& n);

private:
	std::map<int, Point3> _pointsById;
	std::map<int, std::string> _entities;
	std::map<int, std::string>::const_iterator _currentEntity;
};

#endif
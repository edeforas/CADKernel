#include "StepFile.h"

#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <cctype>
#include <iterator>
#include <set>
#include <string>
using namespace std;

void test_bool(bool b, const string& sMessage = "")
{
	if (!b)
	{
		cerr << "Test Error: " << sMessage.c_str() << endl;
		exit(-1);
	}
}

string read_file_text(const string& filename)
{
	ifstream in(filename.c_str(), ios::in | ios::binary);
	test_bool(in.is_open(), "could not open file for integrity check: " + filename);

	return string((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
}

int count_occurrences(const string& text, const string& token)
{
	if (token.empty())
		return 0;

	int count = 0;
	size_t pos = 0;
	while ((pos = text.find(token, pos)) != string::npos)
	{
		++count;
		pos += token.size();
	}

	return count;
}

void test_file_not_empty(const string& filename)
{
	ifstream in(filename.c_str(), ios::binary | ios::ate);
	test_bool(in.is_open(), "could not open exported file: " + filename);
	test_bool(in.tellg() > 0, "exported STEP file is empty: " + filename);
}

set<int> extract_defined_ids(const string& s)
{
	set<int> ids;
	for (size_t i = 0; i < s.size(); ++i)
	{
		if (s[i] != '#')
			continue;

		size_t j = i + 1;
		if (j >= s.size() || !isdigit(static_cast<unsigned char>(s[j])))
			continue;

		int id = 0;
		while (j < s.size() && isdigit(static_cast<unsigned char>(s[j])))
		{
			id = id * 10 + (s[j] - '0');
			++j;
		}

		size_t k = j;
		while (k < s.size() && isspace(static_cast<unsigned char>(s[k])))
			++k;

		if (k < s.size() && s[k] == '=')
			ids.insert(id);
	}

	return ids;
}

set<int> extract_referenced_ids(const string& s)
{
	set<int> ids;
	for (size_t i = 0; i < s.size(); ++i)
	{
		if (s[i] != '#')
			continue;

		size_t j = i + 1;
		if (j >= s.size() || !isdigit(static_cast<unsigned char>(s[j])))
			continue;

		int id = 0;
		while (j < s.size() && isdigit(static_cast<unsigned char>(s[j])))
		{
			id = id * 10 + (s[j] - '0');
			++j;
		}

		if (j < s.size())
		{
			char c = s[j];
			if (c == ')' || c == ',' || c == ';' || isspace(static_cast<unsigned char>(c)))
				ids.insert(id);
		}
	}

	return ids;
}

void test_step_integrity(const string& filename)
{
	const string content = read_file_text(filename);

	test_bool(content.find("ISO-10303-21;") != string::npos, "missing STEP header in " + filename);
	test_bool(content.find("HEADER;") != string::npos, "missing HEADER section in " + filename);
	test_bool(content.find("DATA;") != string::npos, "missing DATA section in " + filename);
	test_bool(content.find("END-ISO-10303-21;") != string::npos, "missing STEP footer in " + filename);
	test_bool(content.find("B_SPLINE_SURFACE_WITH_KNOTS") != string::npos, "missing B-spline knots entity in " + filename);
	test_bool(content.find("RATIONAL_B_SPLINE_SURFACE") != string::npos, "missing rational B-spline entity in " + filename);
	const bool hasSurfaceRep = (content.find("SHAPE_REPRESENTATION") != string::npos) ||
		(content.find("ADVANCED_BREP_SHAPE_REPRESENTATION") != string::npos) ||
		(content.find("MANIFOLD_SURFACE_SHAPE_REPRESENTATION") != string::npos);
	test_bool(hasSurfaceRep, "missing STEP representation entity in " + filename);
	test_bool(content.find("SHAPE_REPRESENTATION_RELATIONSHIP('','',#14,#") == string::npos, "invalid SHAPE_REPRESENTATION_RELATIONSHIP type in " + filename);

	const set<int> defined_ids = extract_defined_ids(content);
	const set<int> referenced_ids = extract_referenced_ids(content);

	test_bool(!defined_ids.empty(), "no STEP entity definitions found in " + filename);
	for (set<int>::const_iterator it = referenced_ids.begin(); it != referenced_ids.end(); ++it)
		test_bool(defined_ids.find(*it) != defined_ids.end(), "unresolved STEP entity reference in " + filename);
}

void test_step_export_torus()
{
	NurbsSurface ns;
	NurbsFactory::create_torus(50, 20, ns);

	StepWriter sw;
	sw.open("test_torus.step");
	test_bool(sw.is_open(), "step writer should open file");
	sw.write(ns);
	sw.close();

	std::ifstream in("test_torus.step");
	test_bool(in.is_open(), "step file should exist");
	test_file_not_empty("test_torus.step");
	test_step_integrity("test_torus.step");
}

void test_step_export_with_invalid_knots_weights()
{
	NurbsSurface ns;
	std::vector<Point3> points = {
		Point3(0., 0., 0.), Point3(1., 0., 0.),
		Point3(0., 1., 0.), Point3(1., 1., 0.)
	};

	ns.set_degree(3, 3);
	ns.set_points(points, 2, 2);
	ns.set_knots_u({ 1., 0., -10. });
	ns.set_knots_v({ 42. });
	ns.set_weights({ 0., -2., std::numeric_limits<double>::quiet_NaN() });

	StepWriter sw;
	sw.open("test_invalid_surface.step");
	test_bool(sw.is_open(), "step writer should open invalid-surface file");
	sw.write(ns);
	sw.close();

	std::ifstream in("test_invalid_surface.step");
	test_bool(in.is_open(), "invalid surface step file should exist");
	test_file_not_empty("test_invalid_surface.step");
	test_step_integrity("test_invalid_surface.step");
}

void test_step_read_torus()
{
	NurbsSurface loaded;
	StepReader reader;
	reader.open("test_torus.step");
	test_bool(reader.read(loaded), "step reader should load torus file");
	test_bool(!reader.read(loaded), "should not read a second surface");

	test_bool(loaded.nb_points_u() > 1, "loaded torus should have at least 2 control points in U");
	test_bool(loaded.nb_points_v() > 1, "loaded torus should have at least 2 control points in V");
	test_bool((int)loaded.points().size() == loaded.nb_points_u() * loaded.nb_points_v(), "loaded torus control points size mismatch");
	test_bool((int)loaded.weights().size() == loaded.nb_points_u() * loaded.nb_points_v(), "loaded torus weights size mismatch");
	test_bool((int)loaded.knots_u().size() == loaded.nb_points_u() + loaded.degree_u() + 1, "loaded torus U knots size mismatch");
	test_bool((int)loaded.knots_v().size() == loaded.nb_points_v() + loaded.degree_v() + 1, "loaded torus V knots size mismatch");
}

void test_step_read_solid_box()
{
	NurbsSolid source;
	NurbsFactory::create_box(10, 20, 30, source);
	test_bool(!source.surfaces().empty(), "source solid box should have surfaces");

	StepWriter sw;
	sw.open("test_box.step");
	test_bool(sw.is_open(), "step writer should open solid box file");
	sw.write(source);
	sw.close();

	NurbsSolid loaded;
	StepReader reader;
	reader.open("test_box.step");
	test_bool(reader.read(loaded), "step reader should load solid box file");
	test_bool(!reader.read(loaded), "should not read a second solid");
	test_bool(loaded.surfaces().size() == source.surfaces().size(), "loaded solid box surface count mismatch");

	for (size_t i = 0; i < loaded.surfaces().size(); ++i)
	{
		const NurbsSurface& s = loaded.surfaces()[i];
		test_bool(s.nb_points_u() > 1, "loaded solid face should have at least 2 control points in U");
		test_bool(s.nb_points_v() > 1, "loaded solid face should have at least 2 control points in V");
		test_bool((int)s.points().size() == s.nb_points_u() * s.nb_points_v(), "loaded solid face points size mismatch");
		test_bool((int)s.weights().size() == s.nb_points_u() * s.nb_points_v(), "loaded solid face weights size mismatch");
		test_bool((int)s.knots_u().size() == s.nb_points_u() + s.degree_u() + 1, "loaded solid face U knots size mismatch");
		test_bool((int)s.knots_v().size() == s.nb_points_v() + s.degree_v() + 1, "loaded solid face V knots size mismatch");
	}
}

void test_step_read_solid_torus()
{
	NurbsSolid source;
	NurbsFactory::create_torus(30, 10, source);
	test_bool(!source.surfaces().empty(), "source solid torus should have surfaces");

	StepWriter sw;
	sw.open("test_solid_torus.step");
	test_bool(sw.is_open(), "step writer should open solid torus file");
	sw.write(source);
	sw.close();

	NurbsSolid loaded;
	StepReader reader;
	reader.open("test_solid_torus.step");
	test_bool(reader.read(loaded), "step reader should load solid torus file");
	test_bool(!reader.read(loaded), "should not read a second solid");
	test_bool(loaded.surfaces().size() == source.surfaces().size(), "loaded solid torus surface count mismatch");

	for (size_t i = 0; i < loaded.surfaces().size(); ++i)
	{
		const NurbsSurface& s = loaded.surfaces()[i];
		test_bool(s.nb_points_u() > 1, "loaded torus face should have at least 2 control points in U");
		test_bool(s.nb_points_v() > 1, "loaded torus face should have at least 2 control points in V");
		test_bool((int)s.points().size() == s.nb_points_u() * s.nb_points_v(), "loaded torus face points size mismatch");
		test_bool((int)s.weights().size() == s.nb_points_u() * s.nb_points_v(), "loaded torus face weights size mismatch");
		test_bool((int)s.knots_u().size() == s.nb_points_u() + s.degree_u() + 1, "loaded torus face U knots size mismatch");
		test_bool((int)s.knots_v().size() == s.nb_points_v() + s.degree_v() + 1, "loaded torus face V knots size mismatch");
	}

	// re write to compare
	StepWriter swout;
	swout.open("test_solid_torus_2ndread.step");
	test_bool(swout.is_open(), "step writer should open file");
	swout.write(loaded);
	swout.close();
}

void test_step_export_multiple_solids_single_representation()
{
	NurbsSolid torus;
	NurbsFactory::create_torus(10, 3, torus);

	NurbsSolid box;
	NurbsFactory::create_box(10, 8, 6, box);

	NurbsSolid sphere;
	NurbsFactory::create_sphere(10, sphere);

	StepWriter sw;
	sw.open("test_multi_solid.step");
	test_bool(sw.is_open(), "step writer should open multi solid file");
	sw.write(torus);
	sw.write(box);
	sw.write(sphere);
	sw.close();

	const string content = read_file_text("test_multi_solid.step");
	test_bool(count_occurrences(content, "MANIFOLD_SOLID_BREP") == 3, "multi solid export should contain three solid breps");
	test_bool(count_occurrences(content, "ADVANCED_BREP_SHAPE_REPRESENTATION") == 1, "multi solid export should use one shared brep representation");
	test_bool(count_occurrences(content, "SHAPE_DEFINITION_REPRESENTATION") == 1, "multi solid export should use one shape definition representation");
}

void test_step_read_multiple_solids()
{
	StepReader reader;
	reader.open("test_multi_solid.step");

	NurbsSolid solid1, solid2, solid3;
	test_bool(reader.read(solid1), "should read first solid (torus)");
	test_bool(reader.read(solid2), "should read second solid (box)");
	test_bool(reader.read(solid3), "should read third solid (sphere)");
	test_bool(!reader.read(solid1), "should not read fourth solid");

	// Basic checks
	test_bool(solid1.surfaces().size() > 0, "first solid should have surfaces");
	test_bool(solid2.surfaces().size() > 0, "second solid should have surfaces");
	test_bool(solid3.surfaces().size() > 0, "third solid should have surfaces");
}

///////////////////////////////////////////////////////////////////////////
int main()
{
	test_step_export_torus();
	test_step_export_with_invalid_knots_weights();
	test_step_read_torus();
	test_step_read_solid_box();
	test_step_read_solid_torus();
	test_step_export_multiple_solids_single_representation();
	test_step_read_multiple_solids();

	cout << "Test Finished.";
	return 0;
}
///////////////////////////////////////////////////////////////////////////
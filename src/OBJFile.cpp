#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <cassert>
using namespace std;

#include "OBJFile.h"

namespace OBJFile
{
	bool save(const string& filename, const Mesh& to_mesh)
	{
		OBJWriter ow;
		ow.open(filename);
		if (!ow.is_open())
			return false;

		ow.write(to_mesh);
		ow.close();
		return true;
	}
	///////////////////////////////////////////////////////////////////////////
	// todo add vertex aliasing on load
	bool load(const string& filename, MeshSolid& solid)
	{
		ifstream f(filename);
		if (!f)
			return false;

		solid.clear();
		Mesh m;
		while (!f.eof())
		{
			string sType, sLine;
			getline(f, sLine);

			if (sLine.size() == 0)
				continue;

			stringstream ss(sLine);
			ss >> sType;

			if (sType=="g")
			{			
				if(!m.empty())
					solid.add_surface(m);

				m.clear();
				string sGroupName;
				ss >> sGroupName;
			}

			// read vertice
			if (sType == "v")
			{
				double x, y, z;
				ss >> x >> y >> z;
				m.add_vertex(Point3(x, y, z));
			}

			// read facet
			if (sType == "f")
			{
				int i1, i2, i3;
				ss >> i1 >> i2 >> i3;

				m.add_triangle(i1 - 1, i2 - 1, i3 - 1);
			}
		}

		if(!m.empty())
			solid.add_surface(m);
		return true;
	}
	///////////////////////////////////////////////////////////////////////////

	bool load(const string& filename, Mesh& mesh)
	{
		MeshSolid solid;
		if (!load(filename, solid))
			return false;

		solid.to_mesh(mesh);
		return true;
	}

	///////////////////////////////////////////////////////////////////////////

}

OBJWriter::OBJWriter()
{
	_verticesCount = 0;
}

OBJWriter::~OBJWriter()
{
	if (_f.is_open())
		_f.close();
}

void OBJWriter::open(const string& filename)
{
	_f.open(filename);
	_iGroup = 0;
}

bool OBJWriter::is_open()
{
	return _f.is_open();
}

void OBJWriter::close()
{
	_f.close();
}

void OBJWriter::write(const Mesh& to_mesh)
{
	assert(_f.is_open());

	double colorR = 0., colorG = 0., colorB = 0.;
	bool bHasColor = to_mesh.get_color() != -1;
	if (bHasColor)
	{
		colorR = (to_mesh.get_color() >> 16) / 255.;
		colorG = ((to_mesh.get_color() & 0xFF00) >> 8) / 255.;
		colorB = (to_mesh.get_color() & 0xFF) / 255.;
	}

	// write vertices
	Point3 vertex;
	for (int i = 0; i < to_mesh.nb_vertices(); i++)
	{
		to_mesh.get_vertex(i, vertex);
		if (bHasColor)
			// save color as in: http://paulbourke.net/dataformats/obj/colour.html
			_f << "v " << vertex._x << " " << vertex._y << " " << vertex._z << " " << colorR << " " << colorG << " " << colorB << endl;
		else
			_f << "v " << vertex._x << " " << vertex._y << " " << vertex._z << endl;
	}
	_f << endl;

	// write triangle vertices
	int iVertex1 = 0, iVertex2 = 0, iVertex3 = 0;

	for (int i = 0; i < to_mesh.nb_triangles(); i++)
	{
		if (to_mesh.is_triangle_unlinked(i))
			continue;

		to_mesh.get_triangle_vertices(i, iVertex1, iVertex2, iVertex3);
		_f << "f " << iVertex1 + 1 + _verticesCount << " " << iVertex2 + 1 + _verticesCount << " " << iVertex3 + 1 + _verticesCount << endl;
	}

	_verticesCount += to_mesh.nb_vertices();
}

void OBJWriter::write(const vector<Point3>& polyline)
{
	assert(_f.is_open());
	
	// write vertices
	Point3 vertex;
	for (int i = 0; i < polyline.size(); i++)
	{
		_f << "v " << polyline[i]._x << " " << polyline[i]._y << " " << polyline[i]._z << endl;
	}
	_f << endl;

	//write polyline
	_f << "l ";
	for (int i = 0; i < polyline.size(); i++)
	{
		_f << i + _verticesCount << " ";
	}

	_verticesCount += polyline.size();
}

void OBJWriter::write(const MeshSolid& solid)
{
	assert(_f.is_open());

	for (int iFace = 0; iFace < solid.nb_surfaces(); iFace++)
	{
		// write each face of the solid as a group
		_f << "g " << "face_" <<iFace+_iGroup << endl;
		write(solid.surface(iFace));
	}

	_iGroup += solid.nb_surfaces();
}

bool OBJFile::save(const string& filename, const MeshSolid& solid)
{
	OBJWriter ow;
	ow.open(filename);
	if (!ow.is_open())
		return false;

	ow.write(solid);
	ow.close();
	return true;
}


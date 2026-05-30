#ifndef _OBJFile_
#define _OBJFile_

#include <string>
#include <vector>
#include <fstream>
using namespace std;

#include "Mesh.h"
#include "MeshSolid.h"

class OBJWriter
{
public:
	OBJWriter();
	virtual ~OBJWriter();

	void open(const string& filename);
	bool is_open();
	void close();

	void write(const Mesh& to_mesh);
	void write(const MeshSolid& solid);
	void write(const vector<Point3>& polyline);
	
private:
	ofstream _f;
	int _iGroup;
	int _verticesCount;
};

//helper functions to load and save OBJ files
namespace OBJFile
{
	bool save(const string& filename, const Mesh& to_mesh);
	bool save(const string& filename, const MeshSolid& solid);

	bool load(const string& filename, MeshSolid& solid);
	bool load(const string& filename, Mesh& mesh);
}


#endif
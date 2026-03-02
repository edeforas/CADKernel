#include "Face.h"
#include <cassert>

///////////////////////////////////////////////////////////////////////////
Face::Face()
{
	_iColor = 0x808080; //grey
}

Face::Face(const Face& f) :
	_mesh(f._mesh),
	_iColor(f._iColor)
{
}

Face::Face(const Mesh& m) :
	_mesh(m),
	_iColor(0x808080) //grey
{
}

Face::~Face()
{
}

Face& Face::operator=(const Face& other)
{
	if (this == &other)
		return *this;

	_mesh = other._mesh;
	_iColor = other._iColor;
	return *this;
}

void Face::set_color(int iColor)
{
	_iColor = iColor;
}

Mesh& Face::mesh()
{
	return _mesh;
}
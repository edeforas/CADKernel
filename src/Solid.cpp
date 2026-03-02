#include "Solid.h"
#include <cassert>

///////////////////////////////////////////////////////////////////////////
Solid::Solid()
{
}

Solid::Solid(const Solid& f) :
	_mesh(f._mesh),
	_faces(f._faces)
{
}

Solid::Solid(const Mesh& m) :
	_mesh(m)
{
}

Solid::~Solid()
{
}

Solid& Solid::operator=(const Solid& other)
{
	if (this == &other)
		return *this;

	_mesh = other._mesh;
	_faces = other._faces;
	return *this;
}


Mesh& Solid::mesh()
{
	return _mesh;
}


void Solid::add_face(const Face& f)
{
	_faces.push_back(f);
}
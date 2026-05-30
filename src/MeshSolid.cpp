#include "MeshSolid.h"
#include <cassert>
#include "Mesh.h"

MeshSolid::MeshSolid()
{ }

MeshSolid::~MeshSolid()
{ }

void MeshSolid::add_surface(Mesh& pSurface)
{
    _vSurfaces.push_back(pSurface);
}

 Mesh& MeshSolid::surface(int iSurface) 
{
    assert(iSurface >= 0 && iSurface < (int)_vSurfaces.size());

    return _vSurfaces[iSurface];
}

const Mesh& MeshSolid::surface(int iSurface) const
{
    assert(iSurface >= 0 && iSurface < (int)_vSurfaces.size());

    return _vSurfaces[iSurface];
}

int MeshSolid::nb_surfaces() const
{
    return (int)_vSurfaces.size();
}

bool MeshSolid::empty() const
{
    return _vSurfaces.empty();
}

void MeshSolid::apply_transform(const Transform& t)
{
    for (Mesh& surface : _vSurfaces)
    {
        surface.apply_transform(t);
    }
}

void MeshSolid::append(const MeshSolid& src)
{
    for (const Mesh& surface : src._vSurfaces)
    {
        _vSurfaces.push_back(surface);
    }
}

void MeshSolid::to_mesh(Mesh& m) const
{
    for (const Mesh& surface : _vSurfaces)
    {
        m.add_mesh(surface);
    }
}

void MeshSolid::clear()
{
    _vSurfaces.clear();
}
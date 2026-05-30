#ifndef _MeshSolid_
#define _MeshSolid_

#include <vector>
using namespace std;
class Mesh;

#include "Transform.h"

class MeshSolid
{
public:
    MeshSolid();
    virtual ~MeshSolid();
	void append(const MeshSolid& src);

    void add_surface(Mesh& pSurface);  
    int nb_surfaces() const;
	bool empty() const;
    
    Mesh& surface(int iSurface);
    const Mesh& surface(int iSurface) const;

    void apply_transform(const Transform& t);

    void to_mesh(Mesh& m) const;
    
    void  clear();
private:
    vector<Mesh> _vSurfaces;
};

#endif
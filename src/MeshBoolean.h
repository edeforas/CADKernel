#ifndef _MeshBoolean_
#define _MeshBoolean_

#include "Mesh.h"

///////////////////////////////////////////////////////////////////////////
class MeshBoolean
{
public:
    MeshBoolean();
    
    //split mesh into boolean Shells, WIP
    void compute_split(const Mesh& A, const Mesh& B, Mesh& Aonly, Mesh& Bonly, Mesh& AInB, Mesh& BInA);

	//compute mesh intersection shell A∩B
	void compute_intersection(const Mesh& A, const Mesh& B, Mesh& intersection);

    //compute mesh union shell A∪B
    void compute_union(const Mesh& A, const Mesh& B, Mesh& meshUnion);

    //compute mesh difference shell A-B
    void compute_difference(const Mesh& A, const Mesh& B, Mesh& meshDifference);

};
///////////////////////////////////////////////////////////////////////////

#endif
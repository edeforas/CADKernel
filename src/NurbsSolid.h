#ifndef NurbsSolid_
#define NurbsSolid_

#include "NurbsTrimmedSurface.h"

#include <vector>
class Transform;

///////////////////////////////////////////////////////////////////////////
class NurbsSolid
{
public:
	NurbsSolid();
	virtual ~NurbsSolid();
	void clear();

	void add_surface(const NurbsSurface& ns);
	std::vector<NurbsTrimmedSurface>& surfaces();
	const std::vector<NurbsTrimmedSurface>& surfaces() const;
	bool empty() const;

	void apply_transform(const Transform& t);

	void append(const NurbsSolid& src);

private:
	std::vector<NurbsTrimmedSurface> _surfaces;
};
///////////////////////////////////////////////////////////////////////////

#endif
#ifndef NurbsSolid_
#define NurbsSolid_

#include "NurbsSurface.h"

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
	std::vector<NurbsSurface>& surfaces();
	const std::vector<NurbsSurface>& surfaces() const;

	void apply_transform(const Transform& t);

	void append(const NurbsSolid& src);

private:
	std::vector<NurbsSurface> _surfaces;
};
///////////////////////////////////////////////////////////////////////////

#endif
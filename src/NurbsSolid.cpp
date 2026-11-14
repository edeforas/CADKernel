#include "NurbsSolid.h"

#include "Transform.h"

///////////////////////////////////////////////////////////////////////////
NurbsSolid::NurbsSolid()
{
	clear();
}

NurbsSolid::~NurbsSolid()
{
}

void NurbsSolid::clear()
{
	_surfaces.clear();
}

std::vector<NurbsTrimmedSurface>& NurbsSolid::surfaces()
{
	return _surfaces;
}

const std::vector<NurbsTrimmedSurface>& NurbsSolid::surfaces() const
{
	return _surfaces;
}

void NurbsSolid::add_surface(const NurbsSurface& ns)
{
	const NurbsTrimmedSurface* trimmedSurface = ns.trimming();
	if (trimmedSurface)
		_surfaces.push_back(*trimmedSurface);
	else
		_surfaces.emplace_back(ns);
}


bool NurbsSolid::empty() const
{
	return _surfaces.empty();
}

///////////////////////////////////////////////////////////////////////////
void NurbsSolid::apply_transform(const Transform& t)
{
	for (auto& f : surfaces())
		for (auto& p : f.points())
			t.apply(p);
}
///////////////////////////////////////////////////////////////////////////
void NurbsSolid::append(const NurbsSolid& src)
{
	for (const auto& s : src.surfaces())
		add_surface(s);
}

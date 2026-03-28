#ifndef NurbsTrimmedSurface_
#define NurbsTrimmedSurface_

#include "NurbsSurface.h"

#include <vector>

struct NurbsUvPoint
{
	double u;
	double v;
};

struct NurbsTrimLoop
{
	std::vector<NurbsUvPoint> points;
	bool hole;
};

class NurbsTrimmedSurface : public NurbsSurface
{
public:
	NurbsTrimmedSurface();
	NurbsTrimmedSurface(const NurbsSurface& n);
	NurbsTrimmedSurface(const NurbsTrimmedSurface& n);
	virtual ~NurbsTrimmedSurface();
	NurbsTrimmedSurface& operator=(const NurbsTrimmedSurface& other);

	void add_outer_loop(const std::vector<NurbsUvPoint>& loopPoints);
	void add_inner_loop(const std::vector<NurbsUvPoint>& loopPoints);

	const std::vector<NurbsTrimLoop>& trim_loops() const;
	virtual bool is_trimmed() const override;
	virtual NurbsTrimmedSurface* trimming() override;
	virtual const NurbsTrimmedSurface* trimming() const override;


private:
	void add_loop(const std::vector<NurbsUvPoint>& loopPoints, bool bHole);
	static bool normalize_loop(std::vector<NurbsUvPoint>& loopPoints, bool bHole);

	std::vector<NurbsTrimLoop> _trimLoops;
};

#endif

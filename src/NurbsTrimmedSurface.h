#ifndef NurbsTrimmedSurface_
#define NurbsTrimmedSurface_

#include "NurbsSurface.h"

#include <vector>

class NurbsCurve;

struct NurbsUvPoint
{
	double u;
	double v;
};

class NurbsTrimLoop
{
public:
	NurbsTrimLoop();
	NurbsTrimLoop(const NurbsTrimLoop& other) ;
	NurbsTrimLoop& operator=(const NurbsTrimLoop& other);

	~NurbsTrimLoop();

	bool has_curve() const;
	
	// Get points, sampling the curve if needed
	const std::vector<NurbsUvPoint>& get_points(int numSamples = 100) const;

	std::vector<NurbsUvPoint> points;  // For backward compatibility and when curve is not used
	NurbsCurve* curve;  // Smooth curve representation (takes ownership if set)
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
	void add_full_outer_loop();

	void add_outer_loop(const NurbsCurve& curve);
	void add_inner_loop(const NurbsCurve& curve);

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

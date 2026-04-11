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
	std::vector<NurbsUvPoint> points;  // For backward compatibility and when curve is not used
	NurbsCurve* curve;  // Smooth curve representation (takes ownership if set)
	bool hole;
	
	NurbsTrimLoop() : curve(nullptr), hole(false) {}
	NurbsTrimLoop(const NurbsTrimLoop& other) 
		: points(other.points), hole(other.hole), curve(nullptr)
	{
		if (other.curve)
			curve = new NurbsCurve(*other.curve);
	}
	NurbsTrimLoop& operator=(const NurbsTrimLoop& other)
	{
		if (this != &other)
		{
			points = other.points;
			hole = other.hole;
			delete curve;
			curve = nullptr;
			if (other.curve)
				curve = new NurbsCurve(*other.curve);
		}
		return *this;
	}
	~NurbsTrimLoop()
	{
		delete curve;
	}
	
	bool has_curve() const { return curve != nullptr; }
	
	// Get points, sampling the curve if needed
	const std::vector<NurbsUvPoint>& get_points(int numSamples = 100) const
	{
		if (points.empty() && curve != nullptr)
		{
			// Lazy sampling: generate points from curve on first access
			const_cast<std::vector<NurbsUvPoint>&>(points).clear();
			
			const std::vector<double>& knots = curve->knots();
			double uMin = knots.front();
			double uMax = knots.back();
			
			for (int i = 0; i < numSamples; ++i)
			{
				double t = uMin + (uMax - uMin) * i / numSamples;
				Point3 p;
				curve->evaluate(t, p);
				const_cast<std::vector<NurbsUvPoint>&>(points).push_back({ p.x(), p.y() });
			}
		}
		return points;
	}
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

#include <iostream>
using namespace std;

#include "NurbsFactory.h"
#include "NurbsSolid.h"
#include "NurbsSurface.h"
#include "NurbsUtil.h"
#include "NurbsFillet.h"
#include "OBJFile.h"
#include "StepFile.h"

////////////////////////////////////////////////////////////////////////////////
int main()
{
	cout << "Creating a box" << endl;
	NurbsSolid box;
	NurbsFactory::create_box(10., 10., 10., box);

	cout << "Box has " << box.surfaces().size() << " surfaces" << endl;

	cout << "Creating fillet on box edge" << endl;
	NurbsFillet fillet;
	NurbsSurface filletSurface;
	// Try different surface pairs
	for (int i = 0; i < (int)box.surfaces().size(); ++i)
	{
		for (int j = i + 1; j < (int)box.surfaces().size(); ++j)
		{
			cout << "Trying surfaces " << i << " and " << j << endl;
			bool success = fillet.create_fillet_or_chamfer_on_solid(box, i, j, true, 1.0, filletSurface); // true for fillet
			if (success)
			{
				cout << "Fillet created successfully between surfaces " << i << " and " << j << endl;

				// Save the original box
				Mesh m;
				NurbsUtil::to_mesh(box, m, 8);
				OBJFile::save("sample_nurbs_fillet_box.obj", m);

				// Save the fillet surface
				Mesh mFillet;
				NurbsUtil::to_mesh(filletSurface, mFillet, 8);
				OBJFile::save("sample_nurbs_fillet_surface.obj", mFillet);

				// Save to STEP
				StepWriter sw;
				sw.open("sample_nurbs_fillet.step");
				sw.write(box);
				sw.write(filletSurface);
				sw.close();

				return 0; // Exit after first success
			}
		}
	}

	cout << "No suitable edge found for fillet" << endl;
	return 1;
}
////////////////////////////////////////////////////////////////////////////////
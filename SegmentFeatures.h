#include <vector>

#include "opencv2/core/core_c.h"
#include "opencv2/core/core.hpp"
#include <opencv2/features2d/features2d.hpp>

using namespace std;
using namespace cv;

struct SegmentFeatures{

	double r; //ratio: (segment area)/(whole image) ratio
	vector<double> elongation; // min & max elongation lambda1, lambda2
	//Point mass; // center of mass
	int segmentLocation;
	vector<double> AverageFilterOutputs;

};
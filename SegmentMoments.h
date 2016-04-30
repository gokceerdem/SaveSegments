#include <vector>

#include "opencv2/core/core_c.h"
#include "opencv2/core/core.hpp"
#include <opencv2/features2d/features2d.hpp>

using namespace std;
using namespace cv;

struct SegmentMoments{
	double m00, m20, m02, m11;
	
	//vector<double> segmentmoments; 
	
};
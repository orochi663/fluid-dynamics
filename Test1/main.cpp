#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include "FluidCanvas.h"

/**************************************************************************/
int main(int argc, char** argv)
{
    FluidCanvas lCanvas(cv::Vec3i(64, 64, 1), 0.1, 1.0, 1.0);
    return 0;
}

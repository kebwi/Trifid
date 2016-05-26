#ifndef __DEBAYER_PROCS__
#define __DEBAYER_PROCS__

#include "ImageReal.h"
#include "Image3ChnlReal.h"

void bayer(Image3ChnlReal rgbImage, ImageReal &grayImage);
void deBayer(ImageReal &grayImage, Image3ChnlReal &rgbImage);

#endif

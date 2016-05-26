#ifndef __FOURIER_TRANSFORM__
#define __FOURIER_TRANSFORM__

#include "FloatType.h"
#include "ImageReal.h"
#include "Image3ChnlReal.h"

enum GainType { GT_NONE, GT_LINEAR, GT_SQRT3, GT_LOG10, GT_LOGLOG10 };

bool assertValidDims(int width, int height);
int findNextPow2(int val);

void generatePowerSpectrum(const ImageReal &ft, ImageReal &ps, GainType gainType, bool shiftCenter);

bool doFourierAnalysis(const Float_t *grayImage, Float_t *fourierTransform, int width, int height);
bool doFourierAnalysis(const ImageReal &grayImage, ImageReal &fourierTransform);

//This function assumes that the image and the filterFT already have the same dimensions
bool multiplyFourierTransforms(Float_t *ft1, Float_t *ft2,
	int width, int height, Float_t *prod);

//This function assumes that the image and the filterFT already have the same dimensions
bool multiplyFourierTransformsAndInvert(Float_t *ft1, Float_t *ft2,
	int width, int height, Float_t *convolvedImg);

//This function assumes that the image and the filterFT already have the same dimensions
bool convolveImageWithFilter(Float_t *img, Float_t *filterFT,
	int width, int height, Float_t *imgConv);

bool embedImageInPow2Bounds(const ImageReal &img, int newWidth, int newHeight, ImageReal &imgOut, int &offsetX, int &offsetY);
bool embedImageInPow2Bounds(const Image3ChnlReal &img, int newWidth, int newHeight, Image3ChnlReal &imgOut, int &offsetX, int &offsetY);
void cropImage(const ImageReal &img, int newWidth, int newHeight, int offsetX, int offsetY, ImageReal &imgOut);
void cropImage(const Image3ChnlReal &img, int newWidth, int newHeight, int offsetX, int offsetY, Image3ChnlReal &imgOut);

void shiftImage(const ImageReal &img, int shiftX, int shiftY, ImageReal &imgOut);
void shiftImage(const Image3ChnlReal &img, int shiftX, int shiftY, Image3ChnlReal &imgOut);

bool convolveImages(const Float_t *img1, const Float_t *img2, Float_t *imgConv);
bool convolveImages(const ImageReal &img1, const ImageReal &img2, ImageReal &imgConv);

void fourierTest();

#endif

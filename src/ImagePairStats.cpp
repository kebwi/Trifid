#include "ImagePairStats.h"
#include "ImageReal.h"
#include "Image3ChnlReal.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <assert.h>

using namespace std;

#define SQR(x) ((x) * (x))

//******************************************************************************
//Extern Globals
extern int gBitDepth;

//******************************************************************************
//Global Declarations
static const int skEdgeSkip = 0;

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

Float_t calcLabDeltaE(const Float_t l1, const Float_t a1, const Float_t b1,
						const Float_t l2, const Float_t a2, const Float_t b2)
{
	return sqrt(SQR(l1 - l2) + SQR(a1 - a2) + SQR(b1 - b2));
}

Float_t calcImgLabDeltaE(const Image3ChnlReal &img1, const Image3ChnlReal &img2)
{
	assert(img1.width() > 0 && img1.height() > 0);
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	
	Float_t deltaEmean = 0;
	int width = img1.width(), height = img1.height();
	int numNonborderPixels = 0;
	for (int y = skEdgeSkip; y < height - skEdgeSkip; y++)
		for (int x = skEdgeSkip; x < width - skEdgeSkip; x++)
		{
			deltaEmean += calcLabDeltaE((*img1[0])(x, y), (*img1[1])(x, y), (*img1[2])(x, y),
										(*img2[0])(x, y), (*img2[1])(x, y), (*img2[2])(x, y));
			numNonborderPixels++;
		}
	
	deltaEmean /= (Float_t)numNonborderPixels;
		
	return deltaEmean;
}

Float_t calc1ChnlME(const ImageReal &img1, const ImageReal &img2)
{
	assert(img1.width() > 0 && img1.height() > 0);
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	
	Float_t mse = 0;
	
	int width = img1.width(), height = img1.height();
	int numNonborderPixels = 0;
	for (int y = skEdgeSkip; y < height - skEdgeSkip; y++)
		for (int x = skEdgeSkip; x < width - skEdgeSkip; x++)
		{
			mse += fabs(img2(x, y) - img1(x, y));
			numNonborderPixels++;
		}
	
	assert(mse >= 0);
	
	mse /= (Float_t)numNonborderPixels;
	
	return mse;
}

Float_t calc1ChnlMSE(const ImageReal &img1, const ImageReal &img2)
{
	assert(img1.width() > 0 && img1.height() > 0);
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	
	Float_t mse = 0;
	
	int width = img1.width(), height = img1.height();
	int numNonborderPixels = 0;
	for (int y = skEdgeSkip; y < height - skEdgeSkip; y++)
		for (int x = skEdgeSkip; x < width - skEdgeSkip; x++)
		{
			mse += SQR(img2(x, y) - img1(x, y));
			numNonborderPixels++;
		}
	
	assert(mse >= 0);
	
	mse /= (Float_t)numNonborderPixels;
	
	return mse;
}

Float_t calcRGB_ME(const Image3ChnlReal &img1, const Image3ChnlReal &img2)
{
	assert(img1.width() > 0 && img1.height() > 0);
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	
	Float_t meanErr = 0;
	
	int width = img1.width(), height = img1.height();
	int numRelevantPixelComponents = 0;
	for (int c = 0; c < 3; c++)
		for (int y = skEdgeSkip; y < height - skEdgeSkip; y++)
			for (int x = skEdgeSkip; x < width - skEdgeSkip; x++)
			{
				meanErr += fabs((*img2[c])(x, y) - (*img1[c])(x, y));
				numRelevantPixelComponents++;
			}
	
	assert(meanErr >= 0);
	
	meanErr /= numRelevantPixelComponents;
	
	return meanErr;
}

Float_t calcRGB_MSE(const Image3ChnlReal &img1, const Image3ChnlReal &img2)
{
	assert(img1.width() > 0 && img1.height() > 0);
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	
	Float_t mse = 0;
	
	int width = img1.width(), height = img1.height();
	int numRelevantPixelComponents = 0;
	for (int c = 0; c < 3; c++)
		for (int y = skEdgeSkip; y < height - skEdgeSkip; y++)
			for (int x = skEdgeSkip; x < width - skEdgeSkip; x++)
			{
				mse += SQR((*img2[c])(x, y) - (*img1[c])(x, y));
				numRelevantPixelComponents++;
			}
	
	assert(mse >= 0);
	
	mse /= numRelevantPixelComponents;
	
	return mse;
}

void calc1ChnlMSE_PSNR(const ImageReal &img1, const ImageReal &img2,
						Float_t &mse, Float_t &psnr)
{
	assert(img1.width() > 0 && img1.height() > 0);
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	
	mse = calc1ChnlMSE(img1, img2);
	
	assert(mse > 0);
	
	//It is assumed that ImageReals are in the range [0-1] and therefore maxPixelValue is 1.0
	Float_t maxPixelValue = 1.0;
	psnr = 20.0 * log10(maxPixelValue / sqrt(mse));
}

void calcRGB_MSE_PSNR(const Image3ChnlReal &img1, const Image3ChnlReal &img2,
						Float_t &mse, Float_t &psnr)
{
	assert(img1.width() > 0 && img1.height() > 0);
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	
	mse = calcRGB_MSE(img1, img2);
	
	assert(mse > 0);
	
	//It is assumed that ImageReals are in the range [0-1] and therefore maxPixelValue is 1.0
	Float_t maxPixelValue = 1.0;
	psnr = 20.0 * log10(maxPixelValue / sqrt(mse));
}

Float_t calculate1ChnlPsnrScore(const ImageReal &imgOri, const ImageReal &img)
{
	assert(imgOri.width() == img.width() && imgOri.height() == img.height());
	
	Float_t mse, psnr;
	calc1ChnlMSE_PSNR(imgOri, img, mse, psnr);
	return psnr;
}

Float_t calculateRgbPsnrScore(const Image3ChnlReal &imgOri, const Image3ChnlReal &img)
{
	assert(imgOri.width() == img.width() && imgOri.height() == img.height());
	
	Float_t mse, psnr;
	calcRGB_MSE_PSNR(imgOri, img, mse, psnr);
	return psnr;
}

Float_t calculateLabScore(const Image3ChnlReal &imgOri, const Image3ChnlReal &img)
{
	assert(imgOri.width() == img.width() && imgOri.height() == img.height());
	
	//We must invert the sign of the lab delta error so that higher lab delta errors are better.
	//This way, psnr and lab delta error can be compared in a similar fashion, i.e.,
	//that higher scores will be considered better scores.
	
	return -calcImgLabDeltaE(imgOri, img);
}

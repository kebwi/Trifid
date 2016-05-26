#include "FourierTransform.h"
#include "ImageProcs.h"
#include "main.h"
#include "Random.h"	//Debug
#include "fftext.h"
#include "fftlib.h"
#include "fft2d.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>

using namespace std;

//******************************************************************************
//Extern Globals

//******************************************************************************
//Global Declarations
static Float_t sInvertedGain = 1;
static const double e = 2.71828182845904523536;

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

//Only accept powers of two excluding greater than or equal to 16
bool assertValidDims(int width, int height)
{
	int numOnes = 0;
	for (int bitPos = 0; bitPos < sizeof(int) * 8; bitPos++)
		if (width & (1 << bitPos))
			++numOnes;
	if (width < 16 || numOnes != 1)
		return false;
	
	numOnes = 0;
	for (int bitPos = 0; bitPos < sizeof(int) * 8; bitPos++)
		if (height & (1 << bitPos))
			++numOnes;
	if (height < 16 || numOnes != 1)
		return false;
	
	return true;
}

int findNextPow2(int val)
{
	for (int i = 2; i <= 65536; i<<=1)
		if (val <= i)
			return i;
	return -1;
}

void generatePowerSpectrum(const ImageReal &ft, ImageReal &ps, GainType gainType, bool shiftCenter)
{
	long dim = ft.width();
	assert(ft.height() == dim);
	long halfDim = dim / 2;
	int totalPixels = ft.totalPixels();
	
	ps.resize(ft.width(), ft.height());
	
	ImageReal psDCinCorner(ft.width(), ft.height());
	
	//Main body
	for (int y = 0; y < dim; y++)
		for (int x = 2; x < dim; x += 2)
			psDCinCorner(x / 2, y) = pow(ft(x, y), 2) + pow(ft(x + 1, y), 2);
	
	//Upper left corner (DC term)
	psDCinCorner(0, 0) = pow(ft(0, 0), 2) + pow(ft(1, 0), 2);
	
	//Left edge (DC stuff)
	for (int y = 1; y < halfDim; y++)
		psDCinCorner(0, y) = psDCinCorner(0, dim - y) = pow(ft(0, y), 2) + pow(ft(1, y), 2);
	psDCinCorner(0, halfDim) = pow(ft(1, 0), 2);
	
	//Center column (Nyquist)
	for (int y = 1; y < halfDim; y++)
		psDCinCorner(halfDim, y) = psDCinCorner(halfDim, dim - y) =
			pow(ft(0, halfDim + y), 2) + pow(ft(1, halfDim + y), 2);
	psDCinCorner(halfDim, 0) = pow(ft(0, halfDim), 2);
	psDCinCorner(halfDim, halfDim) = pow(ft(1, halfDim), 2);
	
	//Reflection of main body into right half of image
	for (int y = 1; y < dim; y++)
		for (int x = 1; x < halfDim; x++)
			psDCinCorner(halfDim + x, y) = psDCinCorner(halfDim - x, dim - y);
	
	//Reflection along top edge (DC stuff)
	for (int x = 1; x <= halfDim; x++)
		psDCinCorner(dim - x, 0) = psDCinCorner(x, 0);
	
	for (int i = 0; i < totalPixels; i++)
		assert(psDCinCorner[i] >= 0);
	
	if (gainType != GT_NONE)
	{
		//Gain
		for (int i = 0; i < totalPixels; i++)
			switch (gainType)
			{
				case GT_LINEAR:
					//psArr[y * dim + x] /= sInvertedGain;
					psDCinCorner[i] /= sInvertedGain;
					break;
				case GT_SQRT3:
					//psArr[y * dim + x] = sqrt(sqrt(sqrt(psArr[y * dim + x])));
					psDCinCorner[i] = sqrt(sqrt(sqrt(psDCinCorner[i])));
					break;
				case GT_LOG10:
					//psArr[y * dim + x] = log10(psArr[y * dim + x]);
					if (psDCinCorner[i] > 0)
						psDCinCorner[i] = log10(psDCinCorner[i]);
					else psDCinCorner[i] = 0;
					break;
				case GT_LOGLOG10:
					//psArr[y * dim + x] = log10(log10(psArr[y * dim + x]));
					if (psDCinCorner[i] > 0)
						psDCinCorner[i] = log10(log10(psDCinCorner[i]));
					else psDCinCorner[i] = 0;
					break;
				case GT_NONE:
					//I just added this case statement to make a warning go away
					assert(false);
			}
		
		for (int i = 0; i < totalPixels; i++)
			assert(!isnan(psDCinCorner[i]));
	}
	
	//Floor nonnegative
	for (int i = 0; i < totalPixels; i++)
		if (psDCinCorner[i] < 0)
			psDCinCorner[i] = 0;
	
	//Normalize
	if (gainType != GT_LINEAR)
	{
		Float_t maxVal = 0;
		for (int i = 0; i < totalPixels; i++)
			if (psDCinCorner[i] > maxVal)
				maxVal = psDCinCorner[i];
		assert(maxVal > 0);
		for (int i = 0; i < totalPixels; i++)
			psDCinCorner[i] /= maxVal;
	}
	
	//Floor/Ceiling to [0-1]
	for (int i = 0; i < totalPixels; i++)
		if (psDCinCorner[i] < 0)
			psDCinCorner[i] = 0;
		else if (psDCinCorner[i] > 1.0)
			psDCinCorner[i] = 1.0;
	
	for (int i = 0; i < totalPixels; i++)
		assert(!isnan(psDCinCorner[i]));
	
	//180 degree phase shift (put DC term in center instead of upper left corner)
	//This is wasteful, the power spectrum should be produced into this arrangement in the first place
	if (shiftCenter)
		for (int y = 0; y < dim; y++)
			for (int x = 0; x < dim; x++)
				ps(x, y) = psDCinCorner((x + halfDim) % dim, (y + halfDim) % dim);
	else
		for (int y = 0; y < dim; y++)
			for (int x = 0; x < dim; x++)
				ps(x, y) = psDCinCorner(x, y);
}

//Input and output can be the same
bool doFourierAnalysis(const Float_t *grayImage, Float_t *fourierTransform, int width, int height)
{
	if (!assertValidDims(width, height))
		return false;
	
	int err = fft2dInit(log2(height), log2(width));
	if (err != 0)
		return false;
	
	if (!fourierTransform)
		fourierTransform = new Float_t(width * height * sizeof(Float_t));
	if (fourierTransform != grayImage)
		memcpy(fourierTransform, grayImage, width * height * sizeof(Float_t));
	rfft2d(fourierTransform, log2(height), log2(width));
	
	return true;
}

bool doFourierAnalysis(const ImageReal &grayImage, ImageReal &fourierTransform)
{
	int width = grayImage.width(), height = grayImage.height();
	
	if (!assertValidDims(width, height))
		return false;
	
	int err = fft2dInit(log2(height), log2(width));
	if (err != 0)
		return false;
	
	fourierTransform = grayImage;
	rfft2d(fourierTransform.img(), log2(height), log2(width));
	
	return true;
}

//This function assumes that the image and the filterFT already have the same dimensions
bool multiplyFourierTransforms(Float_t *ft1, Float_t *ft2,
	int width, int height, Float_t *prod)
{
	if (!assertValidDims(width, height))
		return false;
	
	if (!prod)
		prod = new Float_t[width * height];
	rspect2dprod(ft1, ft2, prod, height, width);
	
	return true;
}

//This function assumes that the image and the filterFT already have the same dimensions
bool multiplyFourierTransformsAndInvert(Float_t *ft1, Float_t *ft2,
	int width, int height, Float_t *convolvedImg)
{
	if (!multiplyFourierTransforms(ft1, ft2, width, height, convolvedImg))
		return false;
	
	int err = fft2dInit(log2(height), log2(width));
	if (err != 0)
		return false;
	
	//Invert the Fourier Transform
	rifft2d(convolvedImg, log2(height), log2(width));
	
	return true;
}

//This function assumes that the image and the filterFT already have the same dimensions
bool convolveImageWithFilter(Float_t *img, Float_t *filterFT,
	int width, int height, Float_t *imgConv)
{
	if (!assertValidDims(width, height))
		return false;
	
	//Get the Fourier Transform of the image
	Float_t *imgFT = new Float_t[width * height];
	assert(doFourierAnalysis(img, imgFT, width, height));
	
	//Multiplty the Fourier Transforms and invert the FT
	return multiplyFourierTransformsAndInvert(imgFT, filterFT, width, height, imgConv);
}

bool embedImageInPow2Bounds(const ImageReal &img, int newWidth, int newHeight, ImageReal &imgOut, int &offsetX, int &offsetY)
{
	if (!assertValidDims(newWidth, newHeight))
		return false;
	
	imgOut.resize(newWidth, newHeight);
	int width = img.width(), height = img.height();
	offsetX = (newWidth - width) / 2;
	offsetY = (newHeight - height) / 2;
	
	//By default, an embedded odd-length signal will be positioned such that the center pixel
	//resides at lower right, but I prefer that it be at upper left, so shift it by one.
	if (width % 2 == 1)
		offsetX++;
	if (height % 2 == 1)
		offsetY++;
	
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			imgOut(x + offsetX, y + offsetY) = img(x, y);
	
	return true;
}

bool embedImageInPow2Bounds(const Image3ChnlReal &img, int newWidth, int newHeight, Image3ChnlReal &imgOut, int &offsetX, int &offsetY)
{
	if (!assertValidDims(newWidth, newHeight))
		return false;
	
	imgOut.resize(newWidth, newHeight);
	int width = img.width(), height = img.height();
	offsetX = (newWidth - width) / 2;
	offsetY = (newHeight - height) / 2;
	
	//By default, an embedded odd-length signal will be positioned such that the center pixel
	//resides at lower right, but I prefer that it be at upper left, so shift it by one.
	if (width % 2 == 1)
		offsetX++;
	if (height % 2 == 1)
		offsetY++;
	
	for (int c = 0; c < 3; c++)
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)
				(*imgOut[c])(x + offsetX, y + offsetY) = (*img[c])(x, y);
	
	return true;
}

void cropImage(const ImageReal &img, int newWidth, int newHeight, int offsetX, int offsetY, ImageReal &imgOut)
{
	int width = img.width(), height = img.height();
	assert(offsetX + newWidth <= width && offsetY + newHeight <= height);
	
	imgOut.resize(newWidth, newHeight);
	
	for (int y = 0; y < newHeight; y++)
		for (int x = 0; x < newWidth; x++)
			imgOut(x, y) = img(x + offsetX, y + offsetY);
}

void cropImage(const Image3ChnlReal &img, int newWidth, int newHeight, int offsetX, int offsetY, Image3ChnlReal &imgOut)
{
	int width = img.width(), height = img.height();
	assert(offsetX + newWidth <= width && offsetY + newHeight <= height);
	
	imgOut.resize(newWidth, newHeight);
	
	for (int c = 0; c < 3; c++)
		for (int y = 0; y < newHeight; y++)
			for (int x = 0; x < newWidth; x++)
				(*imgOut[c])(x, y) = (*img[c])(x + offsetX, y + offsetY);
}

void shiftImage(const ImageReal &img, int shiftX, int shiftY, ImageReal &imgOut)
{
	assert(&img != &imgOut);
	
	if (shiftX < 0)
		shiftX += img.width();
	shiftX %= img.width();
	if (shiftY < 0)
		shiftY += img.height();
	shiftY %= img.height();
	
	imgOut.resize(img.width(), img.height());
	int width = img.width(), height = img.height();
	
	int xx, yy;
	for (int y = 0; y < height; y++)
	{
		yy = (y + shiftY) % height;
		for (int x = 0; x < width; x++)
		{
			xx = (x + shiftX) % width;
			imgOut(x, y) = img(xx, yy);
		}
	}
}

void shiftImage(const Image3ChnlReal &img, int shiftX, int shiftY, Image3ChnlReal &imgOut)
{
	assert(&img != &imgOut);
	
	if (shiftX < 0)
		shiftX += img.width();
	shiftX %= img.width();
	if (shiftY < 0)
		shiftY += img.height();
	shiftY %= img.height();
	
	imgOut.resize(img.width(), img.height());
	int width = img.width(), height = img.height();
	
	int xx, yy;
	for (int c = 0; c < 3; c++)
		for (int y = 0; y < height; y++)
		{
			yy = (y + shiftY) % height;
			for (int x = 0; x < width; x++)
			{
				xx = (x + shiftX) % width;
				(*imgOut[c])(x, y) = (*img[c])(xx, yy);
			}
		}
}

bool convolveImages(const ImageReal &img1, const ImageReal &img2, ImageReal &imgConv)
{
	int maxWidth = max(img1.width(), img2.width());
	int maxHeight = max(img1.height(), img2.height());
	int maxHalfWidth = maxWidth / 2, maxHalfHeight = maxHeight / 2;
	
	if (!assertValidDims(maxWidth, maxHeight))
		return false;
	
	ImageReal img1Temp;
	int offsetX, offsetY;
	assert(embedImageInPow2Bounds(img1, maxWidth, maxHeight, img1Temp, offsetX, offsetY));
	ImageReal img2Temp;
	assert(embedImageInPow2Bounds(img2, maxWidth, maxHeight, img2Temp, offsetX, offsetY));
	
	//Half shift the second image
	ImageReal img2Shifted(maxWidth, maxHeight);
	for (int y = 0; y < maxHeight; y++)
		for (int x = 0; x < maxWidth; x++)
			img2Shifted(x, y) =
				img2Temp((x + maxHalfWidth) % maxWidth, (y + maxHalfHeight) % maxHeight);
	
	//Get the Fourier Transform of the images
	ImageReal img1Ft(maxWidth, maxHeight);
	assert(doFourierAnalysis(img1Temp, img1Ft));
	ImageReal img2Ft(maxWidth, maxHeight);
	assert(doFourierAnalysis(img2Shifted, img2Ft));
	
	//Multiplty the Fourier Transforms
	imgConv.resize(maxWidth, maxHeight);
	rspect2dprod(img1Ft.img(), img2Ft.img(), imgConv.img(), maxHeight, maxWidth);
	
	int err = fft2dInit(log2(maxHeight), log2(maxWidth));
	if (err != 0)
		return false;
	
	//Invert the Fourier Transform
	rifft2d(imgConv.img(), log2(maxHeight), log2(maxWidth));
	
	return true;
}

void fourierTest()
{
	SeedRandom(0);
	
	cout << "Float type sizes on this machine:" << endl;
	cout << "    float:       " << sizeof(float) << endl;
	cout << "    double:      " << sizeof(double) << endl;
	cout << "    long double: " << sizeof(long double) << endl;
	cout << endl;
	
	cout << "Fourier Transform correctness test:" << endl;
	cout << "This test being performed with a floating point type of size: " << sizeof(Float_t) << endl;
	cout << "See FloatType.h for information on how to set the floating point type." << endl;
	cout << endl;
	
	Float_t a[16];
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
			a[j * 4 + i] = RandZeroFloat(1.0);
	
	cout.precision(8);
	cout << "Input array and copy/paste for matlab script:" << endl;
	cout << "a = [ ";
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++)
		{
			cout << a[j * 4 + i];
			if (i != 4 - 1)
				cout << ", ";
			else if (j != 4 - 1)
				cout << "; ";
			else cout << " ]" << endl;
		}
	cout << "b = fft2(a)" << endl;
	cout << "c = ifft2(b)" << endl;
	cout << endl;
	
	cout << "Original signal:" << endl;
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
			cout << setw(15) << a[j * 4 + i];
		cout << endl;
	}
	cout << endl;
	
	int err = fft2dInit(log2(4), log2(4));
	if (err != 0)
		cout << err << endl;
	
	rfft2d(a, log2(4), log2(4));
	
	cout << "Fourier Transform:" << endl;
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
			cout << setw(15) << a[j * 4 + i];
		cout << endl;
	}
	cout << endl;
	
	rifft2d(a, log2(4), log2(4));
	
	cout << "Inverse Fourier Transform (should match original signal above):" << endl;
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
			cout << setw(15) << a[j * 4 + i];
		cout << endl;
	}
	cout << endl;
	
	//=======================================================================
	
	cout << "Fourier Transform timing test" << endl;
	cout << "This test being performed with a floating point type of size: " << sizeof(Float_t) << endl;
	cout << "See FloatType.h for information on how to set the floating point type." << endl;
	cout << endl;
	
	int dim = 512;
	Float_t *a2 = new Float_t[dim * dim];
	for (int j = 0; j < dim; j++)
		for (int i = 0; i < dim; i++)
			a2[j * dim + i] = RandZeroFloat(1.0);
	
	err = fft2dInit(log2(dim), log2(dim));
	if (err != 0)
		cout << err << endl;
	
	Float_t meanTime_s = 0;
	int numPasses = 100;
	for (int i = 0; i < numPasses; i++)
	{
		unsigned long startTime = getCurrentTime_ms();
		rfft2d(a2, log2(dim), log2(dim));
		unsigned long elapsedTime = getCurrentTime_ms() - startTime;
		meanTime_s += (Float_t)elapsedTime / CURRENT_TIME_PER_SEC_F;
	}
	meanTime_s /= (Float_t)numPasses;
	cout << "Mean time to perform 2D Fourier Transform on a " << dim << " x " << dim << " image: "
		<< meanTime_s << " seconds" << endl << endl;
	
	delete [] a2;
}

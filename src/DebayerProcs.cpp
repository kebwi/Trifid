#include "DebayerProcs.h"
#include "Bayer.h"
#include "main.h"
#include "ImageReal.h"
#include "Image3ChnlReal.h"
#include "ImageProcs.h"
#include "SteerablePyramid.h"
#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;

//******************************************************************************
//Extern Globals

//******************************************************************************
//Global Declarations
static int sBitDepth = 16;	//Only used by the hue Debayer procedure

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

void nn_procOnePixel(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	int x, int y)
{
	Color color = determineBayerColor(x, y);
	int width = grayImage.width(), height = grayImage.height();
	
	Float_t bayerValue = 0, bayerValueRight = 0, bayerValueDown = 0, bayerValueRightDown = 0;
	bool rightEdge = false, bottomEdge = false;
	
	bayerValue = grayImage(x, y);
	if (x + 1 < width)
		bayerValueRight = grayImage(x + 1, y);
	else rightEdge = true;
	if (y + 1 < height)
		bayerValueDown = grayImage(x, y + 1);
	else bottomEdge = true;
	if (x + 1 < width && y + 1 < height)
		bayerValueRightDown = grayImage(x + 1, y + 1);
	
	//Fill in the output
	switch (color)
	{
		case R:
			rgbImage(R, x, y) = bayerValue;
			if (!rightEdge)
			{
				rgbImage(G, x, y) = bayerValueRight;
				if (!bottomEdge)
					rgbImage(B, x, y) = bayerValueRightDown;
			}
			break;
		case G:
			{
				Color colorRight;
				if (x < width - 1)
					colorRight = determineBayerColor(x + 1, y);
				else colorRight = determineBayerColor(x - 1, y);
				assert(colorRight != G);
				if (colorRight == R)
				{
					if (!rightEdge)
						rgbImage(R, x, y) = bayerValueRight;
					rgbImage(G, x, y) = bayerValue;
					if (!bottomEdge)
						rgbImage(B, x, y) = bayerValueDown;
				}
				else //if (colorRight == B)
				{
					if (!bottomEdge)
						rgbImage(R, x, y) = bayerValueDown;
					rgbImage(G, x, y) = bayerValue;
					if (!rightEdge)
						rgbImage(B, x, y) = bayerValueRight;
				}
			}
			break;
		case B:
			if (!rightEdge)
			{
				if (!bottomEdge)
					rgbImage(R, x, y) = bayerValueRightDown;
				rgbImage(G, x, y) = bayerValueRight;
			}
			rgbImage(B, x, y) = bayerValue;
			break;
	}
}

void nn_deBayer(ImageReal &grayImage, Image3ChnlReal &rgbImage)
{
	int width = grayImage.width(), height = grayImage.height();
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			nn_procOnePixel(grayImage, rgbImage, x, y);
}

//======================================================================
#pragma mark -

//The two missing colors are always returned in the order red/green, green/blue, or red/blue
void bilinear_calcMissingColors(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	int x, int y, Float_t &color1, Float_t &color2)
{
	Color color = determineBayerColor(x, y);
	int width = grayImage.width(), height = grayImage.height();
	
	Float_t color1L = 0, color2L = 0;
	color1 = color2 = 0;
	int divisor;
	
	if (color == R)
	{
		//Determine green component
		divisor = 0;
		if (y > 0)
		{
			color1L += grayImage(x, y - 1);
			divisor++;
		}
		if (x < width - 1)
		{
			color1L += grayImage(x + 1, y);
			divisor++;
		}
		if (y < height - 1)
		{
			color1L += grayImage(x, y + 1);
			divisor++;
		}
		if (x > 0)
		{
			color1L += grayImage(x - 1, y);
			divisor++;
		}
		assert(divisor > 0);
		color1 = color1L / (Float_t)divisor;
		
		//Determine blue component
		divisor = 0;
		if (x > 0 && y > 0)
		{
			color2L += grayImage(x - 1, y - 1);
			divisor++;
		}
		if (x < width - 1 && y > 0)
		{
			color2L += grayImage(x + 1, y - 1);
			divisor++;
		}
		if (x < width - 1 && y < height - 1)
		{
			color2L += grayImage(x + 1, y + 1);
			divisor++;
		}
		if (x > 0 && y < height - 1)
		{
			color2L += grayImage(x - 1, y + 1);
			divisor++;
		}
		assert(divisor > 0);
		color2 = color2L / (Float_t)divisor;
	}
	else if (color == G)
	{
		Color colorRight;
		if (x < width - 1)
			colorRight = determineBayerColor(x + 1, y);
		else colorRight = determineBayerColor(x - 1, y);
		assert(colorRight != G);
		if (colorRight == R)
		{
			//Determine red component
			divisor = 0;
			if (x > 0)
			{
				color1L += grayImage(x - 1, y);
				divisor++;
			}
			if (x < width - 1)
			{
				color1L += grayImage(x + 1, y);
				divisor++;
			}
			assert(divisor > 0);
			color1 = color1L / (Float_t)divisor;
			
			//Determine blue component
			divisor = 0;
			if (y > 0)
			{
				color2L += grayImage(x, y - 1);
				divisor++;
			}
			if (y < height - 1)
			{
				color2L += grayImage(x, y + 1);
				divisor++;
			}
			assert(divisor > 0);
			color2 = color2L / (Float_t)divisor;
		}
		else //if (colorRight == B)
		{
			//Determine red component
			divisor = 0;
			if (y > 0)
			{
				color1L += grayImage(x, y - 1);
				divisor++;
			}
			if (y < height - 1)
			{
				color1L += grayImage(x, y + 1);
				divisor++;
			}
			assert(divisor > 0);
			color1 = color1L / (Float_t)divisor;
			
			//Determine blue component
			divisor = 0;
			if (x > 0)
			{
				color2L += grayImage(x - 1, y);
				divisor++;
			}
			if (x < width - 1)
			{
				color2L += grayImage(x + 1, y);
				divisor++;
			}
			assert(divisor > 0);
			color2 = color2L / (Float_t)divisor;
		}
	}
	else //if (color == B)
	{
		//Determine red component
		divisor = 0;
		if (x > 0 && y > 0)
		{
			color1L += grayImage(x - 1, y - 1);
			divisor++;
		}
		if (x < width - 1 && y > 0)
		{
			color1L += grayImage(x + 1, y - 1);
			divisor++;
		}
		if (x < width - 1 && y < height - 1)
		{
			color1L += grayImage(x + 1, y + 1);
			divisor++;
		}
		if (x > 0 && y < height - 1)
		{
			color1L += grayImage(x - 1, y + 1);
			divisor++;
		}
		assert(divisor > 0);
		color1 = color1L / (Float_t)divisor;
		
		//Determine green component
		divisor = 0;
		if (y > 0)
		{
			color2L += grayImage(x, y - 1);
			divisor++;
		}
		if (x < width - 1)
		{
			color2L += grayImage(x + 1, y);
			divisor++;
		}
		if (y < height - 1)
		{
			color2L += grayImage(x, y + 1);
			divisor++;
		}
		if (x > 0)
		{
			color2L += grayImage(x - 1, y);
			divisor++;
		}
		assert(divisor > 0);
		color2 = color2L / (Float_t)divisor;
	}
}

void bilinear_procOnePixel(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	int x, int y)
{
	Color color = determineBayerColor(x, y);
	
	Float_t color1 = 0, color2 = 0;
	bilinear_calcMissingColors(grayImage, rgbImage, x, y, color1, color2);
	
	//Fill in the output
	Float_t bayerValue = grayImage(x, y);
	switch (color)
	{
		case R:
			rgbImage(R, x, y) = bayerValue;
			rgbImage(G, x, y) = color1;
			rgbImage(B, x, y) = color2;
			break;
		case G:
			rgbImage(R, x, y) = color1;
			rgbImage(G, x, y) = bayerValue;
			rgbImage(B, x, y) = color2;
			break;
		case B:
			rgbImage(R, x, y) = color1;
			rgbImage(G, x, y) = color2;
			rgbImage(B, x, y) = bayerValue;
			break;
	}
}

void bilinear_deBayer(ImageReal &grayImage, Image3ChnlReal &rgbImage)
{
	int width = grayImage.width(), height = grayImage.height();
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			bilinear_procOnePixel(grayImage, rgbImage, x, y);
}

//======================================================================
#pragma mark -

void hue_constructGreenChannel(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	int x, int y)
{
	assert(determineBayerColor(x, y) != G);
	
	int width = grayImage.width(), height = grayImage.height();
	
	Float_t colorL = 0;
	int divisor = 0;
	if (y > 0)
	{
		colorL += grayImage(x, y - 1);
		divisor++;
	}
	if (x < width - 1)
	{
		colorL += grayImage(x + 1, y);
		divisor++;
	}
	if (y < height - 1)
	{
		colorL += grayImage(x, y + 1);
		divisor++;
	}
	if (x > 0)
	{
		colorL += grayImage(x - 1, y);
		divisor++;
	}
	assert(divisor > 0);
	
	rgbImage(G, x, y) = (Float_t)colorL / (Float_t)divisor;
}

void hue_constructGreenChannelEdgeSensitive(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	int x, int y)
{
	assert(determineBayerColor(x, y) != G);
	
	//Use bilinear interpolation along the edges of the image
	if (x == 0 || y == 0 || x == grayImage.width() - 1 || y == grayImage.height() - 1)
		return hue_constructGreenChannel(grayImage, rgbImage, x, y);
	
	//Use only two of the surrounding four green pixels to interpolate,
	//either the two horizontal neighbors or the two vertical neighbors,
	//depending on which pair preserves the edge best, i.e., interpolate
	//along the edge instead of across it.  If there is no edge here,
	//use all four pixels, i.e., bilinear interpolation.
	
	Float_t horDiff = grayImage(x - 1, y) - grayImage(x + 1, y);
	if (horDiff < 0)
		horDiff = -horDiff;
	Float_t verDiff = grayImage(x, y - 1) - grayImage(x, y + 1);
	if (verDiff < 0)
		verDiff = -verDiff;
	
	//If the horizontal and vertical differences are similar, then assume no edge and bilinear interpolate.
	//For unnoisy lena128, edgeRatioThreshold=.25 yields the best PSNR and edgeRatioThreshold=.75 yields
	//the best lab delta error.  Nonedge hue yields better PSNR than .75, but worst lab delta error than
	//either .25 or .75.  Overall, I like the appearance .75 best. Not sure if this will hold for other
	//images or for noisy bayered images.
	Float_t edgeRatioThreshold = .75;
	if ((horDiff < verDiff && horDiff > verDiff * edgeRatioThreshold) ||
		(verDiff < horDiff && verDiff > horDiff * edgeRatioThreshold))
		return hue_constructGreenChannel(grayImage, rgbImage, x, y);
	
	Float_t colorL = 0;
	if (horDiff < verDiff)	//The edge is vertical, so interpolate horizontally
	{
		colorL += grayImage(x - 1, y);
		colorL += grayImage(x + 1, y);
	}
	else if (verDiff < horDiff)	//The edge is horizontal, so interpolate vertically
	{
		colorL += grayImage(x, y - 1);
		colorL += grayImage(x, y + 1);
	}
	else return hue_constructGreenChannel(grayImage, rgbImage, x, y);
	
	rgbImage(G, x, y) = colorL / 2.0;
}

void hue_constructRedOrBlueAtBlueOrRed(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	ImageReal &rgbImageRorB, int x, int y)
{
	if (&rgbImageRorB == rgbImage[0])
		assert(determineBayerColor(x, y) == B);
	else assert(determineBayerColor(x, y) == R);
	
	int width = grayImage.width(), height = grayImage.height();
	
	Float_t color1L = rgbImage(G, x, y);
	Float_t color2L = 0;
	int divisor = 0;
	if (x > 0 && y > 0)
	{
		color2L += (Float_t)grayImage(x - 1, y - 1) /
					(Float_t)rgbImage(G, x - 1, y - 1);
		divisor++;
	}
	if (x < width - 1 && y > 0)
	{
		color2L += (Float_t)grayImage(x + 1, y - 1) /
					(Float_t)rgbImage(G, x + 1, y - 1);
		divisor++;
	}
	if (x < width - 1 && y < height - 1)
	{
		color2L += (Float_t)grayImage(x + 1, y + 1) /
					(Float_t)rgbImage(G, x + 1, y + 1);
		divisor++;
	}
	if (x > 0 && y < height - 1)
	{
		color2L += (Float_t)grayImage(x - 1, y + 1) /
					(Float_t)rgbImage(G, x - 1, y + 1);
		divisor++;
	}
	
	assert(divisor > 0);
	color1L = ((Float_t)color1L / (Float_t)divisor) * color2L;
	if (color1L > pow(2, sBitDepth) - 1)
		color1L = pow(2, sBitDepth) - 1;
	
	rgbImageRorB(x, y) = color1L;
}

void hue_constructHorizontalRedOrBlueAtGreen(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	ImageReal &rgbImageRorB, int x, int y)
{
	assert(determineBayerColor(x, y) == G);
	if (&rgbImageRorB == rgbImage[0] && x > 0)
		assert(determineBayerColor(x - 1, y) == R);
	else if (&rgbImageRorB == rgbImage[2] && x > 0)
		assert(determineBayerColor(x - 1, y) == B);
	
	int width = grayImage.width();
	
	Float_t color1L = rgbImage(G, x, y);
	Float_t color2L = 0;
	int divisor = 0;
	
	if (x > 0)
	{
		color2L += (Float_t)grayImage(x - 1, y) /
					(Float_t)rgbImage(G, x - 1, y);
		divisor++;
	}
	if (x < width - 1)
	{
		color2L += (Float_t)grayImage(x + 1, y) /
					(Float_t)rgbImage(G, x + 1, y);
		divisor++;
	}
	
	assert(divisor > 0);
	color1L = ((Float_t)color1L / (Float_t)divisor) * color2L;
	if (color1L > pow(2, sBitDepth) - 1)
		color1L = pow(2, sBitDepth) - 1;
	
	rgbImageRorB(x, y) = color1L;
}

void hue_constructVerticalRedOrBlueAtGreen(
	ImageReal &grayImage, Image3ChnlReal &rgbImage,
	ImageReal &rgbImageRorB, int x, int y)
{
	assert(determineBayerColor(x, y) == G);
	if (&rgbImageRorB == rgbImage[0] && y > 0)
		assert(determineBayerColor(x, y - 1) == R);
	else if (&rgbImageRorB == rgbImage[2] && y > 0)
		assert(determineBayerColor(x, y - 1) == B);
	
	int width = grayImage.width();
	
	Float_t color1L = rgbImage(G, x, y);
	Float_t color2L = 0;
	int divisor = 0;
	
	if (x > 0)
	{
		color2L += (Float_t)grayImage(x, y - 1) /
					(Float_t)rgbImage(G, x, y - 1);
		divisor++;
	}
	if (x < width - 1)
	{
		color2L += (Float_t)grayImage(x, y + 1) /
					(Float_t)rgbImage(G, x, y + 1);
		divisor++;
	}
	
	assert(divisor > 0);
	color1L = ((Float_t)color1L / (Float_t)divisor) * color2L;
	if (color1L > pow(2, sBitDepth) - 1)
		color1L = pow(2, sBitDepth) - 1;
	
	rgbImageRorB(x, y) = color1L;
}

void hue_deBayer(ImageReal &grayImage, Image3ChnlReal &rgbImage)
{
	int startX, startY;
	int width = grayImage.width(), height = grayImage.height();
	
	//Transfer the green channel from Bayer to output
	startX = (determineBayerColor(0, 0) == G) ? 0 : 1;
	for (int y = 0; y < height; y++)
		for (int x = (y % 2 == 0 ? startX : (startX + 1) % 2); x < width; x += 2)
			rgbImage(G, x, y) = grayImage(x, y);
	
	//Transfer the red channel from Bayer to output
	determineBlockCoord(R, startX, startY);
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			rgbImage(R, x, y) = grayImage(x, y);
	
	//Transfer the blue channel from Bayer to output
	determineBlockCoord(B, startX, startY);
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			rgbImage(B, x, y) = grayImage(x, y);
	
	//Use bilinear or edge-sensitive interpolation to construct the incomplete green channel
	startX = (determineBayerColor(0, 0) == G) ? 1 : 0;
	for (int y = 0; y < height; y++)
		for (int x = (y % 2 == 0 ? startX : (startX + 1) % 2); x < width; x += 2)
			//hue_constructGreenChannel(grayImage, rgbImage, x, y);
			hue_constructGreenChannelEdgeSensitive(grayImage, rgbImage, x, y);
	
	//Treating the completed green channel as a luminance channel,
	//assume the red and blue channels are chrominance channels and construct them.
	
	//Calculate red pixels at Bayer blue locations
	determineBlockCoord(B, startX, startY);
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			hue_constructRedOrBlueAtBlueOrRed(grayImage, rgbImage, *rgbImage[0], x, y);
	
	//Calculate blue pixels at Bayer red locations
	determineBlockCoord(R, startX, startY);
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			hue_constructRedOrBlueAtBlueOrRed(grayImage, rgbImage, *rgbImage[2], x, y);
	
	//Calculate red pixels at Bayer green locations
	determineBlockCoord(R, startX, startY);
	startX = ++startX % 2;
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			hue_constructHorizontalRedOrBlueAtGreen(grayImage, rgbImage, *rgbImage[0], x, y);
	
	determineBlockCoord(R, startX, startY);
	startY = ++startY % 2;
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			hue_constructVerticalRedOrBlueAtGreen(grayImage, rgbImage, *rgbImage[0], x, y);
	
	//Calculate blue pixels at Bayer green locations
	determineBlockCoord(B, startX, startY);
	startX = ++startX % 2;
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			hue_constructHorizontalRedOrBlueAtGreen(grayImage, rgbImage, *rgbImage[2], x, y);
	
	determineBlockCoord(B, startX, startY);
	startY = ++startY % 2;
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			hue_constructVerticalRedOrBlueAtGreen(grayImage, rgbImage, *rgbImage[2], x, y);
}

//======================================================================
#pragma mark -

void altProj_repairGreen(ImageReal &grayImage, Image3ChnlReal &rgbImage)
{
	int width = grayImage.width(), height = grayImage.height();
	int rx, ry, g1x, g2y, bx, by;
	int redArrayPos = 0, greenArrayPos = 0, blueArrayPos = 0, secondGreenUnusedPos = 0;
	determineBlockColorOrder(rx, ry, g1x, g2y, bx, by, redArrayPos, greenArrayPos, blueArrayPos, secondGreenUnusedPos);
	
	//Generate downsampled images of red, blue, and two greens, one that overlaps red and one that overlaps blue
	ImageReal downsampledR(width, height), downsampledB(width, height);
	ImageReal downsampledGR(width, height), downsampledGB(width, height);
	
	grayImage.downsample(downsampledR, rx, ry);
	grayImage.downsample(downsampledB, bx, by);
	rgbImage[1]->downsample(downsampledGR, rx, ry);
	rgbImage[1]->downsample(downsampledGB, bx, by);
	
	int halfWidth = width / 2, halfHeight = height / 2;
	assert(downsampledR.width() == halfWidth && downsampledR.height() == halfHeight);
	
	//Perform a wavelet analysis of red, blue, and each of the interpolated greens that overlapped red and blue
	SteerablePyramid waveletR, waveletB, waveletGR, waveletGB;		int arrayOffset = 1;
	//Daubechies waveletR, waveletB, waveletGR, waveletGB;		int arrayOffset = 0;
	
	waveletR.DoFullWaveletAnalysis(downsampledR.img(), halfWidth, halfHeight, true);
	cout << ".";	cout.flush();
	waveletB.DoFullWaveletAnalysis(downsampledB.img(), halfWidth, halfHeight, true);
	cout << ".";	cout.flush();
	waveletGR.DoFullWaveletAnalysis(downsampledGR.img(), halfWidth, halfHeight, true);
	cout << ".";	cout.flush();
	waveletGB.DoFullWaveletAnalysis(downsampledGB.img(), halfWidth, halfHeight, true);
	cout << ".";	cout.flush();
	
	int numLevels = waveletGR.NumLevels();
	
	//Find the correlation between green and red and between green and blue on a per level basis
	vector<Float_t> greenRedCorrs, greenBlueCorrs;	//One correlation per level
	waveletGR.findCorrelations(waveletR, greenRedCorrs);
	waveletGB.findCorrelations(waveletB, greenBlueCorrs);
	/*
	cout << "Green/Red correlations: ";
	for (int i = 0; i < greenRedCorrs.size(); i++)
		cout << setw(15) << greenRedCorrs[i];
	cout << endl;
	cout << "Green/Blue correlations:";
	for (int i = 0; i < greenBlueCorrs.size(); i++)
		cout << setw(15) << greenBlueCorrs[i];
	cout << endl;
	*/
	vector<Float_t> means, stdDevs;
	vector<Float_t> thresholdsR(numLevels), thresholdsB(numLevels);
	
	waveletGR.CalcStats(means, stdDevs);
	for (int i = 0; i < numLevels; i++)
		if (i == 0)
			thresholdsR[i] = (((greenRedCorrs[i] + 1) / 2) * means[i + arrayOffset]);
		else thresholdsR[i] = -1;
	
	waveletGB.CalcStats(means, stdDevs);
	for (int i = 0; i < numLevels; i++)
		if (i == 0)
			thresholdsB[i] = (((greenBlueCorrs[i] + 1) / 2) * means[i + arrayOffset]);
		else thresholdsB[i] = -1;
	
	//for (int i = 0; i < numLevels; i++)
	//	cout << "Red/Blue thresholds: " << thresholdsR[i] << "    " << thresholdsB[i] << endl;
	
	//Replace the high frequency coefficients of green with those of red and blue
	waveletGR.CopyCoefficients(waveletR, thresholdsR);
	waveletGB.CopyCoefficients(waveletB, thresholdsB);
	
	//Reconstruct green.  Note that this will only affect the interpolated green pixels.
	waveletGR.DoFullWaveletSynthesis(downsampledGR.img());
	cout << ".";	cout.flush();
	waveletGB.DoFullWaveletSynthesis(downsampledGB.img());
	cout << ".";	cout.flush();
	
	downsampledGR.upsample(*rgbImage[1], rx, ry);
	downsampledGB.upsample(*rgbImage[1], bx, by);
}

bool altProj_repairRedBlue(SteerablePyramid &steerPyrR, const SteerablePyramid &steerPyrG, SteerablePyramid &steerPyrB,
							ImageReal &r, ImageReal &b,
							const ImageReal &rRealOri, const ImageReal &bRealOri)
{
	int width = r.width(), height = r.height();
	assert(b.width() == width && b.height() == height);
	
	int numLevels = steerPyrG.NumLevels();
	
	//Perform a wavelet analysis of red and blue
	steerPyrR.DoFullWaveletAnalysis(r.img(), width, height, true);
	cout << ".";	cout.flush();
	steerPyrB.DoFullWaveletAnalysis(b.img(), width, height, true);
	cout << ".";	cout.flush();
	
	//Compare the high frequency components of red and blue with green.
	//Use the similarity/correlation between high frequency red/green and blue/green
	//to derive thresholds to be used below.
	vector<Float_t> greenRedCorrs, greenBlueCorrs;	//One correlation per level
	steerPyrG.findCorrelations(steerPyrR, greenRedCorrs);
	steerPyrG.findCorrelations(steerPyrB, greenBlueCorrs);
	/*
	cout << "Green/Red correlations: ";
	for (int i = 0; i < greenRedCorrs.size(); i++)
		cout << setw(15) << greenRedCorrs[i];
	cout << endl;
	cout << "Green/Blue correlations:";
	for (int i = 0; i < greenBlueCorrs.size(); i++)
		cout << setw(15) << greenBlueCorrs[i];
	cout << endl;
	*/
	int arrayOffset = 1;	//1 for Steerable Pyramid, 0 for Daubechies
	
	vector<Float_t> means, stdDevs;
	steerPyrG.CalcStats(means, stdDevs);
	assert(means.size() == steerPyrG.NumLevels() + arrayOffset);
	
	vector<Float_t> thresholdsR(numLevels), thresholdsB(numLevels);
	for (int i = 0; i < numLevels; i++)
	{
		if (i == 0)
		{
			thresholdsR[i] = (((greenRedCorrs[i] + 1) / 2) * means[i + 1]);
			thresholdsB[i] = (((greenBlueCorrs[i] + 1) / 2) * means[i + 1]);
		}
		else thresholdsR[i] = thresholdsB[i] = -1;
	}
	
	//for (int i = 0; i < numLevels; i++)
	//	cout << "Red/Blue thresholds: " << thresholdsR[i] << "    " << thresholdsB[i] << endl;
	
	//Push the high frequency coefficients of red and blue toward those of green
	steerPyrR.CopyCoefficients(steerPyrG, thresholdsR);
	steerPyrB.CopyCoefficients(steerPyrG, thresholdsB);
	
	//Reconstruct red and blue
	steerPyrR.DoFullWaveletSynthesis(r.img());
	cout << ".";	cout.flush();
	steerPyrB.DoFullWaveletSynthesis(b.img());
	cout << ".";	cout.flush();
	
	//Copy original red and blue pixels into the reconstruction
	int startX, startY;
	determineBlockCoord(R, startX, startY);
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			r(x, y) = rRealOri(x, y);
	
	determineBlockCoord(B, startX, startY);
	for (int y = startY; y < height; y += 2)
		for (int x = startX; x < width; x += 2)
			b(x, y) = bRealOri(x, y);
	
	//Determine if convergence has been attained
	return true;
}

/*
NOTES:
The following setting are configurable, with options labeled in square brackets
and used as references in notes at the bottom:
-- Initial green interpolation:
	-- Bilinear [A0]
	-- Edge-sensitive [A1]
-- Wavelet:
	-- Daubechies 9/7 [B0]
	-- Steerable pyramid 4 orientations [B1]
-- Number of levels down (starting from the highest frequency coefficients) to be transferred from one channel to another:
	-- One level [C0]
	-- Two levels [C1]
	-- All levels [C2]
-- When transferring from red/blue to green, threshold should be:
	-- Zero, as in original alt proj publication [D0]
	-- Correlation-sensitive, much like green to red/blue is done [D1]
-- Number of red/blue update iterations:
	-- 1 [E0]
	-- 2 [E1]

Results for Lena 16 bit with no noise:
  -- [A0 B1 C0 D1 E1]:
	-- Yields the lowest Lab Delta E
	-- Yields the lowest scielab error
	-- Has the second best subjective appearance
  -- [A1 B1 C0 D1 E1]:
	-- Yields the highest PSNR in all channels and in RGB overall
	-- Yields the second lowest scielab error
	-- Has the best subjective appearance
*/
void altProj_deBayer(ImageReal &grayImage, Image3ChnlReal &rgbImage)
{
	cout << "DeBayer alternating projections (A-E)  ";	cout.flush();
	
	int width = grayImage.width(), height = grayImage.height();
	
	//Start by using bilinear interpolation to reconstruct all three incomplete channels
	cout << "a";	cout.flush();
	bilinear_deBayer(grayImage, rgbImage);
	cout << "A";	cout.flush();
	
	//Alternatively, we can use edge-sensitive interpolation to construct the incomplete green channel.
	//Doing this effectively replaces the interpolated green pixels from the bilinear interpolation above.
	//A more efficient approach would not perform both methods on the green channel.
	//In one experiment I got a more accurate green reconstruction using bilinear however.
	
	cout << "b";	cout.flush();
	int startX = (determineBayerColor(0, 0) == G) ? 1 : 0;
	for (int y = 0; y < height; y++)
		for (int x = (y % 2 == 0 ? startX : (startX + 1) % 2); x < width; x += 2)
			hue_constructGreenChannelEdgeSensitive(grayImage, rgbImage, x, y);
	cout << "B";	cout.flush();
	
	//Repair the interpolated green pixels
	cout << "c";	cout.flush();
	altProj_repairGreen(grayImage, rgbImage);
	cout << "C";	cout.flush();
	
	//Generate wavelet decompositions for all three channels and analyze green
	ImageReal rReal = *rgbImage[0], gReal = *rgbImage[1], bReal = *rgbImage[2];
	
	SteerablePyramid steerPyrR, steerPyrG, steerPyrB;
	cout << "d";	cout.flush();
	assert(rgbImage[1]->width() == width && rgbImage[1]->height() == height);
	steerPyrG.DoFullWaveletAnalysis(gReal.img(), width, height, true);
	cout << "D";	cout.flush();
	
	//Loop until a stopping criterion is satisifed
	int numIterations = 0;
	while (true)
	{
		cout << "e";	cout.flush();
		bool converged = altProj_repairRedBlue(steerPyrR, steerPyrG, steerPyrB,
												rReal, bReal,
												*rgbImage[0], *rgbImage[2]);
		
		//if (converged)
		//	break;
		if (++numIterations == 2)	//2 yields highest performance on unnoisy lena. 
			break;
	}
	cout << "E";	cout.flush();
	
	*rgbImage[0] = rReal;
	*rgbImage[2] = bReal;
	
	cout << endl;
}

//======================================================================
#pragma mark -

void bayer(Image3ChnlReal rgbImage, ImageReal &grayImage)
{
	cout << endl << "Bayerizing..." << endl;
	
	int width = rgbImage.width(), height = rgbImage.height();
	
	if (grayImage.width() != rgbImage.width() || grayImage.height() != rgbImage.height())
		die("Bayer output isn't same size as bayer input");
	
	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
		{
			Color color = determineBayerColor(x, y);
			switch (color)
			{
				case R:
					grayImage(x, y) = rgbImage(R, x, y);
					break;
				case G:
					grayImage(x, y) = rgbImage(G, x, y);
					break;
				case B:
					grayImage(x, y) = rgbImage(B, x, y);
					break;
			}
		}
	
	grayImage.setFilename(rgbImage.filename());
	grayImage.appendFilename("_b");
}

void deBayer(ImageReal &grayImage, Image3ChnlReal &rgbImage)
{
	cout << endl << "DeBayerizing..." << endl;
	
	if (rgbImage.width() != grayImage.width() || rgbImage.height() != grayImage.height())
		die("DeBayer output isn't same size as deBayer input");
	
	string suffix;
	switch (gDeBayerMethod)
	{
		case DB_NEAREST_NEIGHBOR:
			nn_deBayer(grayImage, rgbImage);
			suffix = "-nn";
			break;
		case DB_BILINEAR:
			bilinear_deBayer(grayImage, rgbImage);
			suffix = "-bl";
			break;
		case DB_SMOOTH_HUE:
			hue_deBayer(grayImage, rgbImage);
			suffix = "-hue";
			break;
		case DB_ALTERNATING_PROJECTIONS:
			altProj_deBayer(grayImage, rgbImage);
			suffix = "-ap";
			break;
	}
	
	rgbImage.setFilename(grayImage.filename());
	rgbImage.appendFilename("_db");
	rgbImage.appendFilename(suffix);
}

#include "ImageProcs.h"
#include "main.h"
#include "Random.h"
#include "ImageReal.h"
#include "Image3ChnlReal.h"
#include "SteerablePyramid.h"
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>	//debug

using namespace std;

//******************************************************************************
//Extern Globals
extern bool gVerbose;
extern bool gGenerateIntermediateImages;

//******************************************************************************
//Global Declarations
static int sInputHistogram[256] = { 0 };
static double sMeanNoise[256] = { 0 };	//Mean noise applied per input value
static double sMeanNoiseProporation[256] = { 0 };	//Mean noise relative to nonnoised output
static double sMeanOutputValue[256] = { 0 };	//Output value relative to input value
static double sMeanError[256] = { 0 };	//Diff between final value and input value

static const Float_t kAutoGainConvergence = .001;

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

void rgb2xyz(Float_t r, Float_t g, Float_t b,
			Float_t &x, Float_t &y, Float_t &z)
{
	assert(r >= 0 && r <= 1);
	assert(g >= 0 && g <= 1);
	assert(b >= 0 && b <= 1);
	
	//Good resources:
	//http://www.aim-dtp.net/aim/technology/color-spaces_faq/color-space.faq.txt
	//http://www.brucelindbloom.com/
	
	//----------------------------------------------------------------------
	
	//Taken from http://www.cs.rit.edu/~ncs/color/t_convert.html
	/*
	static const Float_t kM[3][3] =	{
											{ 0.412453, 0.357580, 0.180423 },
											{ 0.212671, 0.715160, 0.072169 },
											{ 0.019334, 0.119193, 0.950227 }
										};
	
	x = r * kM[0][0] + g * kM[0][1] + b * kM[0][2];
	y = r * kM[1][0] + g * kM[1][1] + b * kM[1][2];
	z = r * kM[2][0] + g * kM[2][1] + b * kM[2][2];
	*/
	//----------------------------------------------------------------------
	
	//Taken from http://www.easyrgb.com/math.php
	//and from http://en.wikipedia.org/wiki/SRGB_color_space
	
	/*
	//Convert srgb to linear rgb (if necessary)
	//http://www.sjbrown.co.uk/?article=gamma
	Float_t alpha = .055, gamma = 2.4;
	if (r > 0.04045)
		r = pow((( r + alpha) / (1 + alpha)), gamma);
	else r /= 12.92;
	if (g > 0.04045)
		g = pow((( g + alpha) / (1 + alpha)), gamma);
	else g /= 12.92;
	if (b > 0.04045)
		b = pow((( b + alpha) / (1 + alpha)), gamma);
	else b /= 12.92;
	*/
	//What the hell are these x100 multiplications for?
	/*
	r *= 100;
	g *= 100;
	b *= 100;
	*/
	//Linear rgb to xyz
	//http://www.easyrgb.com/math.php
	//x = r * 0.4124 + g * 0.3576 + b * 0.1805;
	//y = r * 0.2126 + g * 0.7152 + b * 0.0722;
	//z = r * 0.0193 + g * 0.1192 + b * 0.9505;
	
	//More precise values from http://www.mediamacros.com/item/item-1006687658/
	x = r * 0.412453 + g * 0.357580 + b * 0.180423;
	y = r * 0.212671 + g * 0.715160 + b * 0.072169;
	z = r * 0.019334 + g * 0.119193 + b * 0.950227;
}

void xyz2rgb(Float_t x, Float_t y, Float_t z,
			Float_t &r, Float_t &g, Float_t &b)
{
	r = x *  3.240479 + y * -1.537150 + z * -0.498535;
	g = x * -0.969256 + y *  1.875992 + z *  0.041556;
	b = x *  0.055648 + y * -0.204043 + z *  1.057311;
	
	if (r < 0)
		r = 0;
	else if (r > 1)
		r = 1;
	
	if (g < 0)
		g = 0;
	else if (g > 1)
		g = 1;
	
	if (b < 0)
		b = 0;
	else if (b > 1)
		b = 1;
}

void xyz2cielab(Float_t x, Float_t y, Float_t z,
				Float_t &l, Float_t &a, Float_t &b)
{
	//assert(x >= 0 && x <= 1);
	//assert(y >= 0 && y <= 1);
	//assert(z >= 0 && z <= 1);
	
	//----------------------------------------------------------------------
	
	//Taken from http://www.easyrgb.com/math.php
	//and http://en.wikipedia.org/wiki/SRGB_color_space
	//and http://www.mediamacros.com/item/item-1006687658/
	//and http://www.cs.rit.edu/~ncs/color/t_convert.html
	
	//Two degree D65 white point tristimulus values
	x /= 95.047;
	y /= 100.0;
	z /= 108.883;
	
	if (x > 0.008856)
		x = pow(x, 1.0 / 3.0);
	else x = (7.787 * x) + (16.0 / 116.0);
	if (y > 0.008856)
		y = pow(y, 1.0 / 3.0);
	else y = (7.787 * y) + (16.0 / 116.0);
	if (z > 0.008856)
		z = pow(z, 1.0 / 3.0);
	else z = (7.787 * z) + (16.0 / 116.0);
	
	//This resource ( http://www.cs.rit.edu/~ncs/color/t_convert.html ) computes l like this:
	if (y > 0.008856)
		l = (116.0 * y) - 16.0;
	else l = 903.3 * y;
	//l = (116 * y) - 16;
	a = 500.0 * (x - y);
	b = 200.0 * (y - z);
}

void cielab2xyz(Float_t l, Float_t a, Float_t b,
				Float_t &x, Float_t &y, Float_t &z)
{
	y = (l + 16.0) / 116.0;
	x = a / 500.0 + y;
	z = y - b / 200.0;
	
	if (pow(y, 3) > .008856)
		y = pow(y, 3);
	else y = (y - 16.0 / 116.0) / 7.787;
	if (pow(x, 3) > .008856)
		x = pow(x, 3);
	else x = (x - 16.0 / 116.0) / 7.787;
	if (pow(z, 3) > .008856)
		z = pow(z, 3);
	else z = (z - 16.0 / 116.0) / 7.787;
	
	x *= 95.047;
	y *= 100.0;
	z *= 108.883;
}

void testRGBXYZLABRoutines()
{
	Float_t minorCutoff = .000001;
	Float_t majorCutoff = .001;
	
	cout.precision(8);
	
	cout << "RGB coordinates specified as input, then run through RGB->XYZ->LAB->XYZ->RGB." << endl
		<< "Final XYZ and RGB difs are presented and compared to originals with minor and major cutoffs: " << endl
		<< minorCutoff << " and " << majorCutoff << "." << endl << endl;
	
	cout << "Display legend:" << endl;
	cout	<< setw(12) << "R in"
			<< setw(12) << "G in"
			<< setw(12) << "B in" << "    "
			<< setw(12) << "X out"
			<< setw(12) << "Y out"
			<< setw(12) << "Z out" << "    " << endl
			<< setw(12) << "R2 out"
			<< setw(12) << "G2 out"
			<< setw(12) << "B2 out" << "    "
			<< setw(12) << "X2 out"
			<< setw(12) << "Y2 out"
			<< setw(12) << "Z2 out" << endl
			<< setw(12) << "|X2 - X|"
			<< setw(12) << "|Y2 - Y|"
			<< setw(12) << "|Z2 - Z|" << "    "
			<< setw(12) << "|R2 - R|"
			<< setw(12) << "|G2 - G|"
			<< setw(12) << "|B2 - B|" << endl << endl;
			
	for (Float_t r = 0; r <= 1.0; r += .25)
		for (Float_t g = 0; g <= 1.0; g += .25)
			for (Float_t b = 0; b <= 1.0; b += .25)
			{
				Float_t x, y, z;
				Float_t l, a, bb;
				Float_t x2, y2, z2;
				Float_t r2, g2, b2;
				
				rgb2xyz(r, g, b, x, y, z);
				xyz2cielab(x, y, z, l, a, bb);
				cielab2xyz(l, a, bb, x2, y2, z2);
				xyz2rgb(x2, y2, z2, r2, g2, b2);
				
				cout	<< setw(12) << r
						<< setw(12) << g
						<< setw(12) << b << "    "
						<< setw(12) << x
						<< setw(12) << y
						<< setw(12) << z << endl
						<< setw(12) << r2
						<< setw(12) << g2
						<< setw(12) << b2 << "    "
						<< setw(12) << x2
						<< setw(12) << y2
						<< setw(12) << z2 << endl
						<< setw(12) << fabs(x2 - x)
						<< setw(12) << fabs(y2 - y)
						<< setw(12) << fabs(z2 - z) << "    "
						<< setw(12) << fabs(r2 - r)
						<< setw(12) << fabs(g2 - g)
						<< setw(12) << fabs(b2 - b) << endl;
				
				if (fabs(r2 - r) > minorCutoff || fabs(g2 - g) > minorCutoff || fabs(b2 - b) > minorCutoff ||
					fabs(x2 - x) > minorCutoff || fabs(y2 - y) > minorCutoff || fabs(z2 - z) > minorCutoff)
					cout << "    -------- minor error in previous test --------" << endl;
				
				if (fabs(r2 - r) > majorCutoff || fabs(g2 - g) > majorCutoff || fabs(b2 - b) > majorCutoff ||
					fabs(x2 - x) > majorCutoff || fabs(y2 - y) > majorCutoff || fabs(z2 - z) > majorCutoff)
					cout << "    ******** MAJOR ERROR IN PREVIOUS TEST ********" << endl;
				
				cout << endl;
			}
	
	exit(0);
}

void printColor(Float_t r, Float_t g, Float_t b)
{
	assert(r >= 0 && r <= 1 &&
			g >= 0 && g <= 1 &&
			b >= 0 && b <= 1);
	
	Float_t x, y, z, l, a, b2;
	
	rgb2xyz(r, g, b, x, y, z);
	xyz2cielab(x, y, z, l, a, b2);
	
	cout << "RGB:";
	cout << setw(15) << r;
	cout << setw(15) << g;
	cout << setw(15) << b << endl;
	cout << "XYZ:";
	cout << setw(15) << x;
	cout << setw(15) << y;
	cout << setw(15) << z << endl;
	cout << "LAB:";
	cout << setw(15) << l;
	cout << setw(15) << a;
	cout << setw(15) << b2 << endl << endl;
}

void rgbImgToLabImg(Float_t maxVal, const Image3ChnlReal &rgbImg, Image3ChnlReal &labImg)
{
	labImg.resize(rgbImg.width(), rgbImg.height());
	
	Float_t r, g, b;
	Float_t x, y, z;
	Float_t l, a, b2;
	int totalPixels = rgbImg.totalPixels();
	for (int i = 0; i < totalPixels; i++)
	{
		r = (*rgbImg[0])[i] / maxVal;
		g = (*rgbImg[1])[i] / maxVal;
		b = (*rgbImg[2])[i] / maxVal;
		
		rgb2xyz(r, g, b, x, y, z);
		xyz2cielab(x, y, z, l, a, b2);
		/*
		if (l < 0 || a < 0 || b2 < 0)
			cout << setw(10) << r
				<< setw(10) << g
				<< setw(10) << b
				<< setw(15) << x
				<< setw(10) << y
				<< setw(10) << z
				<< setw(15) << l
				<< setw(10) << a
				<< setw(10) << b2 << endl;
		*/
		(*labImg[0])[i] = l;
		(*labImg[1])[i] = a;
		(*labImg[2])[i] = b2;
		/*
		if (i % 128 == 0)
		{
			cout << "RGB:";
			cout << setw(15) << r;
			cout << setw(15) << g;
			cout << setw(15) << b << endl;
			cout << "XYZ:";
			cout << setw(15) << x;
			cout << setw(15) << y;
			cout << setw(15) << z << endl;
			cout << "LAB:";
			cout << setw(15) << l;
			cout << setw(15) << a;
			cout << setw(15) << b2 << endl << endl;
		}
		*/
	}
	
	labImg.setFilename(rgbImg.filename());
	labImg.appendFilename("_lab");
}

void labImgToRgbImg(Float_t maxVal, const Image3ChnlReal &labImg, Image3ChnlReal &rgbImg)
{
	rgbImg.resize(labImg.width(), labImg.height());
	
	Float_t r, g, b;
	Float_t x, y, z;
	Float_t l, a, b2;
	int totalPixels = labImg.totalPixels();
	for (int i = 0; i < totalPixels; i++)
	{
		l = (*labImg[0])[i];
		a = (*labImg[1])[i];
		b2 = (*labImg[2])[i];
		
		cielab2xyz(l, a, b2, x, y, z);
		xyz2rgb(x, y, z, r, g, b);
		
		(*rgbImg[0])[i] = r * maxVal;
		(*rgbImg[1])[i] = g * maxVal;
		(*rgbImg[2])[i] = b * maxVal;
	}
	
	rgbImg.setFilename(labImg.filename());
	rgbImg.appendFilename("_rgb");
}

Float_t calcNonoutlierMeanPixelValue(ImageReal &grayImage, Float_t outlierSigmas)
{
	int totalPixels = grayImage.totalPixels();
	
	//Find the mean pixel value
	Float_t meanVal = 0;
	for (int i = 0; i < totalPixels; i++)
		meanVal += grayImage[i];
	meanVal /= totalPixels;
	
	//Find the variance of pixel value
	Float_t variance = 0;
	for (int i = 0; i < totalPixels; i++)
		variance += pow(grayImage[i] - meanVal, 2);
	variance /= (totalPixels * 3);
	
	//Find the standard deviation of pixel value
	Float_t stddev = sqrt(variance);
	
	//Define outliers as some number of standard deviations
	Float_t minOutlier = meanVal - stddev * outlierSigmas;
	Float_t maxOutlier = meanVal + stddev * outlierSigmas;
	
	//Find the mean pixel value of nonoutliers
	Float_t nonOutlierMeanVal = 0;
	int numNonOutliers = 0;
		for (int i = 0; i < totalPixels; i++)
			if (grayImage[i] >= minOutlier && grayImage[i] <= maxOutlier)
			{
				nonOutlierMeanVal += grayImage[i];
				numNonOutliers++;
			}
	nonOutlierMeanVal /= numNonOutliers;
	
	return nonOutlierMeanVal;
}

Float_t calcNonoutlierMeanPixelValue(Image3ChnlReal &rgbImage, Float_t outlierSigmas)
{
	int totalPixels = rgbImage.totalPixels();
	
	//Find the mean pixel value
	Float_t meanVal = 0;
	for (int c = 0; c < 3; c++)
		for (int i = 0; i < totalPixels; i++)
			meanVal += (*rgbImage[c])[i];
	meanVal /= (totalPixels * 3);
	
	//Find the variance of pixel value
	Float_t variance = 0;
	for (int c = 0; c < 3; c++)
		for (int i = 0; i < totalPixels; i++)
			variance += pow((*rgbImage[c])[i] - meanVal, 2);
	variance /= (totalPixels * 3);
	
	//Find the standard deviation of pixel value
	Float_t stddev = sqrt(variance);
	
	//Define outliers as some number of standard deviations
	Float_t minOutlier = meanVal - stddev * outlierSigmas;
	Float_t maxOutlier = meanVal + stddev * outlierSigmas;
	
	//Find the mean pixel value of nonoutliers
	Float_t nonOutlierMeanVal = 0;
	int numNonOutliers = 0;
	for (int c = 0; c < 3; c++)
		for (int i = 0; i < totalPixels; i++)
			if ((*rgbImage[c])[i] >= minOutlier && (*rgbImage[c])[i] <= maxOutlier)
			{
				nonOutlierMeanVal += (*rgbImage[c])[i];
				numNonOutliers++;
			}
	nonOutlierMeanVal /= numNonOutliers;
	
	return nonOutlierMeanVal;
}

//======================================================================
#pragma mark -

Float_t luxToPPFD(
	Float_t lux,
	Float_t lensFocalLength,	//Unused in the current PPFD calculation, remnant of an earlier version
	Float_t lensFnumber,
	Float_t pixelPitchMicrons,
	Float_t fullWellCapacity
)
{
	cout.precision(15);
	
	cout << setw(40) << left << "  lux:" << lux << endl;
	Float_t luxAfterLensTransmission = lux * .9;
	cout << setw(40) << left << "  luxAfterLensTransmission:" << luxAfterLensTransmission << endl;
	cout << setw(40) << left << "  lensFnumber:" << lensFnumber << endl;
	Float_t luxAtImagePlane = luxAfterLensTransmission / pow(lensFnumber, 2);
	cout << setw(40) << left << "  luxAtImagePlane:" << luxAtImagePlane << endl;
	Float_t photoSensitivePixelPitchMicrons = pixelPitchMicrons - .25;	// half of a .5u width dead zone between pixels
	cout << setw(40) << left << "  photoSensitivePixelPitchMicrons:" << photoSensitivePixelPitchMicrons << endl;
	Float_t photoSensitivePixelAreaM = pow((photoSensitivePixelPitchMicrons / 1000000), 2);
	cout << setw(40) << left << "  photoSensitivePixelAreaM:" << photoSensitivePixelAreaM << endl;
	Float_t lumensPerPixel = luxAtImagePlane * photoSensitivePixelAreaM;
	cout << setw(40) << left << "  lumensPerPixel:" << lumensPerPixel << endl;
	Float_t ppfdPerPixel = lumensPerPixel * .0185;
	cout << setw(40) << left << "  ppfdPerPixel:" << ppfdPerPixel << endl;
	Float_t photonsPerPixelPerSec = ppfdPerPixel * 6.0221415e17;
	cout << setw(40) << left << "  photonsPerPixelPerSec:" << photonsPerPixelPerSec << endl;
	
	Float_t saturationTime = fullWellCapacity / photonsPerPixelPerSec;
	cout << setw(40) << left << "  saturationTime at 100% reflectance:" << saturationTime << endl;
	
	cout.precision(6);
	return photonsPerPixelPerSec;
}

//A detailed description of the noise model is available at the top of the function body
void addSensorNoise(ImageReal &grayImage, bool applyBayer, NoiseParam noiseParams)
{
	/*
	Many references were used in the construction of the following model.  However, one stands out:
		00276126.pdf (http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=276126)
		Radiometric CCD Camera Calibration and Noise Estimation,
		Glenn E. Healey and Raghava Kondepudy,
		IEEE Transactions on Pattern Analysis and Machine Intelligence,
		Vol 16., No. 3, pp. 267-276, Mar 1994.
		
		Another very helpful reference was:
		Analog & Digital Signal Processing,
		H. Baher,
		John Wiley & Sons, 2nd ed., 2001.
		
		See the README for additional references.
	
	------------------------------------------------------------------------------------------------
	
	A full CCD/CMOS sensor noise model is quite complicated.
	The following is a list and description of most of the contributing factors, as described in Healey and Kondepudy.
	
	a,b					Pixel location in col, row.
	m					Temperature.  Affects d and h (not sure about g).
	T					Exposure time.  Affects x, d, h, and g.
	B(a,b)				Irradiance.  The amount of energy/light incident at a pixel (generally a function of wavelength, but that is ignored here).
	I(a,b,T)			Ideal sensor value (this is what we truly want to know/read/sense).
	K(a,b)				Bias (quantum efficiency, vignetting, and dust).
							Scalar (multiplicative) per-pixel sensitivity, best represented as strictly < 1 overall (it can only decrease response, not increase it).
							Quantum efficiency is a per-pixel ratio of electron-flux to photon-flux strictly < 1.  Note this is generally highly wavelength dependent as well.
							Quantum efficiency distribution is probably normal around the CCD's overall efficiency (e.g. 50%-90% for some CCDs).
							Vignetting represents a "brighter" response at the center of the CCD due to a lens effect with a presumably Gaussian drop-off.
							Flecks of dust on the optics appear as dim fuzzy circular (or annulus in the case of telescopes with a secondary mirror) spots on the CCD.
	d(m,T)				Dark current.
							Additive accumulation of thermally-generated free electrons.
							Positive (can only increase the perceived value).
							Can vary slightly between pixels, thus d(a,b,T,m), but this is ignored by this model.
							Note that hot pixels are basically pixel-dependent dark current and are represented separately in this model.
							I am assuming this is Poisson, and therefore, normally distributed, but I don't know the mean or variance dependencies.
	g(a,b,T)			Amplifier glow.
							Additive.
							Positive (can only increase the perceived value).
							Negligable for all but very long exposures, thus, ignored by this model.
							For most sensors, amp glow is a (presumably circular) gradient originating at one corner or just off-frame from one corner of the sensor.
							The gradient dropoff function is currently unknown, probably 1/r^2 I'm guessing, perhaps 1/r (depending on the degree to which spherical radiance is confined to the plane of the CCD circuitry by reflection from adjacent surfaces).
							Not sure if this depends on m, temperature.
							Not sure if this increases with A, amplifier intensity/gain.
	A					Amplifier gain (amp glow)
							Scalar (multiplicative).
	h(a,b,m,T,A)		Hot pixel.
							Additive charge leakage.
							Positive (can only increase the perceived value).
							Only notable in a smattering of pixels, tens per CCD, gets worse with temp and exposure.
							Scales with the amplifier gain A, not sure how.
							Grows much faster than dark current relative to time (saturates in seconds instead of minutes).
							I am assuming this is Poisson, and therefore, normally distributed, but I don't know the mean or variance dependencies.
	Av					Variance of amplifier noise (constant)
							Additive.
							Normal with 0 mean.
							Variance is independent of sensor signal.
							Amplifier noise dominates shot noise at low signal levels, thereby determining the noise floor.
	v1 = A^2*(I*K+d+g)	Variance of signal dependent noise (basically shot noise)
							Healey & Kondepudy, Eq. 26, with the addition of amp glow g.
							Additive.
							Variance increases as I increases, possibly v should equal I because it is Poisson, I'm not sure.
							One ref states that the uncertainty in signal level equals the square root of the number of photons collected.
							Strictly speaking, shot noise can only be expressed in integer terms, i.e., in discrete photons, which doesn't jive with this program's real-valued pixels.
							Normally distributed with mean 0 and variance v.
	v2 = A^2 * Av		Variance of signal independent noise (basically amplifier noise)
							Healey & Kondepudy, Eq. 27, dismissing the quantization error variance (q^2/12) shown in that equation.
							Amplifier noise dominates shot noise at low signal levels, thereby determining the noise floor.
							I think the variance of h should be included here somehow, not sure how though.
	v3 = v2 + v1		Variance of noise
							Healey & Kondepudy, Eq. 28.
	n					Noise, zero mean normal with variance v3
	y = A*(I*K+d+g+h)+n	Sensor output
							Healey & Kondepudy, Eq. 20, with the addition of amp glow g and hot pixel h.
							I am not absolutely certain that h is being applied correctly, but I think this makes sense.
	
	------------------------------------------------------------------------------------------------
	
	Slightly different attempt, matching variable names to Healey & Kondepudy as much as possible:
	
	un(r)				Uniformally distributed random variable over range [-.5r,.5r], with implicit variance (r^2 / 12)
	nm(v)				Normally distributed random variable with variance v
	a,b					Pixel location
	m					Temperature scalar
	T					Exposure time
	B(a,b)				Irradiance (implicitly a function of time, but not exposure time)
	I(a,b) =			Ideal value (the number of photons/electrons should have ideally been detected/generated)
		B * T
	K(a,b)				Bias (quantum efficiency, vignetting, dust)
	d(m,T)				Dark current
							A more sophisticated model would use d(a,b,T,m), but this model ignores location-dependency of dark current
	g(a,b,T)			Amplifier glow
	h(a,b,m,T,A)		Hot pixel
	A					Amplifier gain
	u(a,b) =			Expected value without noise
		A*(I*K+d+g+h)		Healey & Kondepudy, Eq. 22, with the addition of amp glow g and hot pixel h
	UNUSED:
	Ns(a,b) =			Zero mean Poisson shot noise with variance dependent on K*I, d and g
		nm(I*K+g+h)			Healey & Kondepudy, Eq. 2 and following paragraph
	Nrv					Nr variance
							Unsure how this is determined, whether it varies with A, anything...
	UNUSED:
	Nr =				Zero mean amplifier noise independent of the number of collected electrons
		nm(Nrv)				Healey & Kondepudy, Eq. 3 and preceeding paragraph
	Niv(a,b) =			Ni variance
		A^2*(I*K+d+g)		Healey & Kondepudy, Eq. 26
	UNUSED:
	Ni(a,b) =			Part of noise dependent on number of collected electrons
		A * Ns				Healey & Kondepudy, Eq. 24
	fw =				Full well capacity of a pixel (number of electrons it can hold)
							See note in definition of rn
	rn					Root mean square CCD read noise (expected number of electrons with zero illuminance, zero exposure time, zero dark current)
							See note in definition of ccdDR
	sensorSNR =			CCD's signal-to-noise ratio
		fw / rn				Common definition of signal-to-noise ratio (max recordable signal without clipping over noise floor)
	ccdDR =				CCD's dynamic range (in dB)
		(20 * log(SNR))		Common definition of dynamic range w.r.t. SNR
							The ADC should use enough bits to encode this dynamic range
	V =					Analog voltage (microvolts) passed to the analog/digital converter (through the amplifier first)
							See note in definition of q
	bt =				Number of quantization bits
							See note in definition of q
	q =					Quantization step
		1 / 2^bt			Baher, p. 417, Eq. 11.2.  See also relationship between bt and q at the end of Eq. 11.3 in determining quantization variance.
							According to Healey & Kondepudy, Eq. 6 and following paragraph, bt and q should be chosen so that (V <= (2^bt - .5)q)
		1 / (2^bt - 1)		My variation of Baher's equation which seems obviously flawed to me.  See note in the code below where q is calculated.
	qc =				Quantization variance
		q^2 / 12			Healey & Kondepudy, Eq. 6 and following paragraph, and Baher, p. 417, Eq. 11.13
							See also definiition of un above
	Pc =				Peak power of a digital coder
		q^2*2^(2*bt-3)		Baher, p. 419, Eq. 11.12
	Rc =				Coding dynamic range of a digital coder
		3*(2^bt-1)			Baher, p. 419, Eq. 11.13
		10*log(Pc/qc)		Baher, p. 419, Eq. 11.14
		6.02*bt+1.76dB		Baher, p. 419, Eq. 11.14 (8-bit yields 50dB)
							See note in definition of Rcs
	Rcs	=				Coding dynamic range of a digital coder of a signal that is scaled to the ADC's max power
		6.02*bt-3.1876dB	Baher, p. 420, Eq. 11.20 (8-bit yields 45dB, 16-bit yields 100dB)
							I *think* the important point of all of this is that bt (and consequently q) must be chosen such that the ADC's dynamic range (Rcs) meets or exceeds the CCD's dynamic range (ccdDR)
	Nq =				Zero mean quantization noise, uniform distribution over the the range [-.5q,.5q]]
		un(q)				Healey & Kondepudy, Eq. 6 and following paragraph
	Ncv =				Nc variance
		A^2*Nrv+(q^2/12)	Healey & Kondepudy, Eq. 27
	UNUSED:
	Nc =				Part of noise independent of number of collected electrons
		A * Nr + Nq			Healey & Kondepudy, Eq. 25
	Nv(a,b) =			N variance (total noise variance)
		Niv + Ncv			Healey & Kondepudy, Eq. 28
	UNUSED:
	N(a,b) =			Zero mean noise
		Ni + Nc				Healey & Kondepudy, Eq. 23
	N(a,b) =			Zero mean noise
		nm(Nv)				Healey & Kondepudy, Eq. 28 and preceeding paragraph
	y(a,b) =			Sensor output
		u + N				Healey & Kondepudy, Eq. 21
	
	------------------------------------------------------------------------------------------------
	
	Forms of noise not included this model:
		Blooming occurs when a saturated pixel overflows into nearby pixels, often more in one cardinal direction than the other, depending on the CCD's physical architecture.
		There is a quantization error described as N_Q in Healey & Kondepudy, which this model disregards.
		Transfer error is an inefficiency in transfer from the CCD to the amplifier.  It is generally negligable, i.e., most sensors and amplifers have a very efficient transter.
		Analog-to-digital converter (ADC) noise is 0 mean additive and uniform (not normal), independent of measured signal.
	*/
	
	cout << "Adding grayscale sensor noise..." << endl;
	
	for (int i = 0; i < 256; i++)
	{
		sInputHistogram[i] = 0;
		sMeanNoise[i] = 0;
		sMeanNoiseProporation[i] = 0;
		sMeanOutputValue[i] = 0;
		sMeanError[i] = 0;
	}
	
	noiseParams.selfCheck();
	
	//cout << "Noise params received by addSensorNoise():" << endl;
	//noiseParams.dump();
	
	//It is important to *not* normalize the image first since this gray image may be a RGB
	//layer of an RGB image and we don't want to normalize the three layers independently, i.e., color shift.
	
	//Calculate the image-plane luminance as PPFD in photons/s/px
	cout << "Converting scene min/max luminance to photon flux..." << endl;
	Float_t maxPPFD = luxToPPFD(noiseParams.lux_, 0, noiseParams.focalRatio_, noiseParams.pixelPitchMicrons_, noiseParams.fullWellCapacity_);
	
	//Calculate the reflectance range
	Float_t reflectanceRange = noiseParams.maxReflectance_ - noiseParams.minReflectance_;
	
	//Debug
	Float_t fullReflectanceSaturationTimeS = noiseParams.fullWellCapacity_ / maxPPFD;
	Float_t maxQE = max(max(noiseParams.quantumEfficiency_[0], noiseParams.quantumEfficiency_[1]), noiseParams.quantumEfficiency_[2]);
	Float_t maxReflectanceQESaturationTimeS = noiseParams.fullWellCapacity_ / (maxPPFD * noiseParams.maxReflectance_ * maxQE);
	
	//Calculate the temperature scalar for dark current and hot pixels
	const Float_t kRoomTemperatureC = 20;
	const Float_t kDarkCurrentTempScalarRate = 6.0;	//Every diff in temp by this much is worth a factor of 2 in dark current rate.
	Float_t temperatureScalar = pow(2.0, fabs(noiseParams.temperatureC_ - kRoomTemperatureC) / kDarkCurrentTempScalarRate);
	assert(temperatureScalar > 0);
	if (noiseParams.temperatureC_ < kRoomTemperatureC)
		temperatureScalar = 1.0 / temperatureScalar;
	
	//Calculate the dark current electrons per pixel at the indicated temperature
	Float_t meanDarkCurrentElectronsPerPix =
		noiseParams.darkCurrentRoomTempElectronsPerPixPerSec_
		* temperatureScalar
		* noiseParams.exposureSecs_;
	
	//Convert some values to electron counts
	Float_t meanAmpGlowElectronsPerPix = noiseParams.ampGlowElectronsPerPixPerSec_ * noiseParams.exposureSecs_;
	Float_t ampVarFactorScaled = noiseParams.ampVarFactor_;// * noiseParams.fullWellCapacity_;	//Is this correct?
	
	//Calculate the hot pixel severity factor
	Float_t hotPixelFactorScaled =
		temperatureScalar
		* noiseParams.exposureSecs_
		//* noiseParams.ampGain_	//Not sure whether this should go here
		* noiseParams.hotPixelFactor_
		* noiseParams.fullWellCapacity_;	//Is this right?
	
	//I have absolutely no idea how the variance of the read noise is determined or how it tends to vary in real CCDs.
	//I will assume the standard deviation equals the mean.
	Float_t readNoiseVar = pow(noiseParams.readNoise_, 2);
	
	//Determine the ADC's minimum required number of bits.
	//Is this all correct?
	Float_t sensorSNR = noiseParams.fullWellCapacity_ / noiseParams.readNoise_;
	Float_t sensorDRdB = 20.0 * log10(sensorSNR);
	
	Float_t adcRequiredBits = 0;
	while (true)
	{
		Float_t adcDRdB = 6.02 * ++adcRequiredBits - 3.1876;
		if (adcDRdB >= sensorDRdB)
			break;
	}
	
	//Determine the quantized range in discrete units
	int maxQuantizedValue = pow(2, adcRequiredBits) - 1;
	
	//Set the quantization step, q.
	//Notice the subtraction by one (q = 1 / (2^bt - 1)), not shown in Baher p.417 eq. 11.2.
	//This subtraction seems correct to me since the resulting histogram will have 2^bt bins.
	//Baher's equation (q = 1 / 2^bt) yields a histogram with (2^bt + 1) bins!
	Float_t quantizationStep = 1.0 / (pow(2.0, adcRequiredBits) - 1);	//This step size is based on a range [0,1]
	Float_t quantizationStepScaled = quantizationStep * (maxQuantizedValue + 1);	//I'm *relatively* sure this is correct
	Float_t quantizationStepScaledHalved = quantizationStepScaled / 2.0;	//Used for actual quantizing (rounding)
	
	//A^2
	Float_t ampGainSquared = pow(noiseParams.ampGain_, 2);
	
	//Final Ncv.
	//Notice that quantization error isn't included here.  While it is often included in similar
	//equations in the literature in an *analysis* of error, e.g., Healey & Kondepudy, Eq. 27,
	//it shouldn't be explicitly added as a form of noise variance here.  Such noise will implicitly
	//occur during the quantization process later.  It hardly matters though because as far as
	//I can tell, read noise always vastly dominates quantization error, and since one is derived
	//from the other, I cannot think of a combination of settings in which this would not occur.
	//In other words, what's the point of being concerned about quantization error at all; what difference does it make?
	Float_t signalIndepVar =
		ampVarFactorScaled * ampGainSquared * (
			+ readNoiseVar	//I *think* this is where this goes
		);	//Excluded -> + pow(quantizationStepScaled, 2.0) / 12.0;
	
	cout << endl;
	cout << "Quantization step and ADC bits w.r.t. CCD dynamic range..." << endl;
	cout << "  maxQuantizedValue: " << maxQuantizedValue << endl;
	cout << "  quantizationStep: " << quantizationStep << endl;
	cout << "  quantizationStepScaled: " << quantizationStepScaled << endl;
	cout << "  quantizationStepScaled^2 / 12: " << (pow(quantizationStepScaled, 2.0) / 12.0) << endl;
	cout << "  CCD SNR: " << sensorSNR << "  and DR: " << sensorDRdB << " dB" << endl;
	cout << "  ADC required bits: " << adcRequiredBits << "  yielding ADC DR of " << 6.02 * adcRequiredBits - 3.1876 << " dB" << endl;
	cout << "  ADC value range: " << "[0," << maxQuantizedValue << "]" << endl;
	
	if (adcRequiredBits < 10)
	{
		stringstream msgSS;
		msgSS << "ADC required bits (" << adcRequiredBits << "b) is unrealistically low (below 10b)";
		warning(msgSS.str().c_str());
	}
	if (adcRequiredBits > 16)
	{
		stringstream msgSS;
		msgSS << "ADC required bits (" << adcRequiredBits << "b) is unrealistically high (above 16b)";
		warning(msgSS.str().c_str());
	}
	
	//Check that the well density is within a realistic range 
	Float_t wellDensity = noiseParams.fullWellCapacity_ / pow(noiseParams.pixelPitchMicrons_, 2);
	if (wellDensity < 400)
	{
		stringstream msgSS;
		msgSS << "Well density (" << wellDensity << "e/u^2) is unrealistically low (below 400e/u^2)";
		warning(msgSS.str().c_str());
	}
	if (wellDensity > 1800)
	{
		stringstream msgSS;
		msgSS << "Well density (" << wellDensity << "e/u^2) is unrealistically high (above 1800e/u^2)";
		warning(msgSS.str().c_str());
	}
	
	cout << endl;
	cout << "Settings:" << endl;
	//cout << setw(35) << left << "  minPPFD: " << minPPFD << endl;
	cout << setw(40) << left << "  maxPPFD: " << maxPPFD << endl;
	cout << setw(40) << left << "  wellDensity: " << wellDensity << endl;//(wellDensity > 0 ? wellDensity : "N/A, pixel pitch required") << endl;
	cout << setw(40) << left << "  fullReflectanceSaturationTimeS: " << fullReflectanceSaturationTimeS << endl;
	cout << setw(40) << left << "  maxReflectanceQESaturationTimeS: " << maxReflectanceQESaturationTimeS << endl;
	cout << setw(40) << left << "  temperatureScalar: " << temperatureScalar << endl;
	cout << setw(40) << left << "  meanDarkCurrentElectronsPerPix: " << meanDarkCurrentElectronsPerPix << endl;
	cout << setw(40) << left << "  hotPixelFactorScaled: " << hotPixelFactorScaled << endl;
	cout << setw(40) << left << "  meanAmpGlowElectronsPerPix: " << meanAmpGlowElectronsPerPix << endl;
	cout << setw(40) << left << "  ampGain/squared: " << noiseParams.ampGain_ << "  " << ampGainSquared << endl;
	cout << setw(40) << left << "  readNoise/squared: " << noiseParams.readNoise_ << "  " << readNoiseVar << endl;
	cout << setw(40) << left << "  quantizationStepScaled^2 / 12: " << (pow(quantizationStepScaled, 2.0) / 12.0) << endl;
	cout << setw(40) << left << "  signalIndepVar: " << signalIndepVar << endl;
	
	if (noiseParams.exposureSecs_ > maxReflectanceQESaturationTimeS)
	{
		stringstream msgSS;
		msgSS << "Exposure time (" << noiseParams.exposureSecs_ << "s) exceeds saturation time (" << maxReflectanceQESaturationTimeS << "s)";
		warning(msgSS.str().c_str());
	}
	
	//Debug
	Float_t a[20] = { 0 };
	
	Float_t meanSignalDepVar = 0;
	
	//Apply the noise model on a per-pixel basis
	int totalPixels = grayImage.totalPixels();
	int numSaturatedPixelsA = 0, numSaturatedPixelsB = 0, numSaturatedPixelsC = 0, numSaturatedPixelsD = 0;
	for (int i = 0; i < totalPixels; i++)
	{
		//Get the Bayer color and the quantum efficiency for that color
		Color bayerColor = R;
		Float_t quantumEfficiency = 1.0;
		if (applyBayer)
		{
			int x = i % grayImage.width();
			int y = i / grayImage.width();
			bayerColor = determineBayerColor(x, y);
			quantumEfficiency = noiseParams.quantumEfficiency_[bayerColor];
		}
		
		//Add the input value to the input histogram
		int inputVal = grayImage[i] * 256;
		if (inputVal > 255)
			inputVal = 255;
		sInputHistogram[inputVal]++;
		
		//Get the scene reflectance for one pixel (and one channel ala Bayer) and convert to PPFD,
		//i.e., photons/electrons detected/generated per sec, B.
		Float_t onePixelPPFD = maxPPFD * (reflectanceRange * grayImage[i] + noiseParams.minReflectance_);
		//Float_t onePixelPPFD = maxPPFD * grayImage[i];
		
		//Convert the flux density (photons/s/px) to actual flux (actual number of photons received by this pixel)
		//by exposing the pixel for some period of time.  Note that photon flux is discrete (integer).
		int onePixelPPF = roundToNearestInt(onePixelPPFD * noiseParams.exposureSecs_, 0, 0);//noiseParams.fullWellCapacity_);
		
		//Convert the flux (number of photons) to the electron count, I.
		//Note that quantum efficiency is ignored here.  It is applied through K, the bias, if applied at all.
		//Note that electron count is discrete (integer).
		int idealElectronCount = onePixelPPF;	//Note, quantum efficiency is not modeled here.
		assert(idealElectronCount >= 0);
		
		//Determine the bias, K.
		//Note that there is no site-to-site variation, e.g., vignetting or dust, in this model,
		//but there is color-dependent quantum efficiency.
		Float_t bias = quantumEfficiency;
		
		//Determine the hot pixel value
		Float_t hotPixelVal =
			noiseParams.hotPixelScalars_
			? (*noiseParams.hotPixelScalars_)[i] * hotPixelFactorScaled
			: 0;
		
		//Calculate the actual number of electrons collected by the pixel.
		//Note that this is where blooming would be applied.
		//Note that electron count is discrete (integer).
		int actualElectronCount = roundToNearestInt(
			idealElectronCount * bias + meanDarkCurrentElectronsPerPix + hotPixelVal + meanAmpGlowElectronsPerPix,
			0, 0);	//No clipping inside the rounding function
		assert(actualElectronCount >= 0);
		if (actualElectronCount > noiseParams.fullWellCapacity_)
		{
			//This is where blooming would be applied
			numSaturatedPixelsA++;
			actualElectronCount = noiseParams.fullWellCapacity_;
		}
		
		//Calculate the expected value without noise, u.
		//I suspect that perhaps this shouldn't be clipped to full well capacity.
		//By the time we are applying the amplifier gain, we are off the pixel.
		Float_t expectedVal = noiseParams.ampGain_ * actualElectronCount;
		
		//Calculate the signal dependent variance, Niv (exclude hotPixelVal).
		Float_t actualSignalElectronCount = min(
			idealElectronCount * bias + meanDarkCurrentElectronsPerPix + meanAmpGlowElectronsPerPix,
			(Float_t)noiseParams.fullWellCapacity_);
		Float_t signalDepVar = ampGainSquared * ampVarFactorScaled * actualSignalElectronCount;
		
		//Calculate total noise variance, Nv.
		Float_t noiseVar = signalDepVar + signalIndepVar;
		
		//Generate the noise, N.
		Float_t noise = RandGenNorm(0, sqrt(noiseVar));
		
		//Generate the sensor output, y.
		//By the time we are applying the amplifier gain, we are off the pixel.
		Float_t noisedVal = roundToNearestInt(expectedVal + noise, 0, 0/*noiseParams.fullWellCapacity_*/);
		if (noisedVal < 0)
			noisedVal = 0;
		assert(noisedVal >= 0);
		if (noisedVal > noiseParams.fullWellCapacity_)
			numSaturatedPixelsB++;
		
		//Quantize the value.
		//Is this definitely the correct place to do this, effectively quantizing *before* deBayering?
		//I suspect there are two quantization stages.
		//One here, from the CCD, and one later, after deBayering, for digital image storage.  Is this true?
		//Frequently, the CCD uses more bits than the final image, for example.
		//Note that the later stage is already implicit in this program in ImageInterface::writeFile().
		noisedVal /= noiseParams.fullWellCapacity_;	//Normalize noisedVal to [0,1]
		noisedVal *= maxQuantizedValue;	//Scale noisedVal to [0,maxQuantizedValue]
		noisedVal =
			((unsigned long)((noisedVal + quantizationStepScaledHalved)
				/ quantizationStepScaled))
			* quantizationStepScaled;
		//noisedVal = roundToNearestInt(noisedVal, 0, 0/*noiseParams.fullWellCapacity_*/);
		assert(maxQuantizedValue >= 0);
		if (noisedVal > maxQuantizedValue)
		{
			numSaturatedPixelsC++;
			noisedVal = maxQuantizedValue;
		}
		
		//Debug
		int idx = -1;
		a[++idx] += onePixelPPFD;
		a[++idx] += onePixelPPF;
		a[++idx] += idealElectronCount;
		a[++idx] += actualElectronCount;
		a[++idx] += expectedVal;
		a[++idx] += hotPixelVal;
		a[++idx] += meanDarkCurrentElectronsPerPix + hotPixelVal + meanAmpGlowElectronsPerPix;
		a[++idx] += signalDepVar;
		a[++idx] += sqrt(noiseVar);
		a[++idx] += fabs(noise);
		a[++idx] += noisedVal;
		
		//Keep track of the proportion of the output signal that is noise
		sMeanNoiseProporation[inputVal] +=
			fabs(noise) /
				(noiseParams.ampGain_
				* (idealElectronCount * bias + meanDarkCurrentElectronsPerPix + hotPixelVal + meanAmpGlowElectronsPerPix));
		sMeanNoise[inputVal] += fabs(noise);
		
		//Scale back to the range [0-1] for storage in the ImageInterface object
		noisedVal /= maxQuantizedValue;//noiseParams.fullWellCapacity_;
		assert(noisedVal >= 0 && noisedVal <= 1.0);
		
		//Correct the color shift that resulted from the quantum efficiency.
		//Is this really the best way to fix this?  Is this how a real camera would deal with it?
		//I think this is basically the white-balancing problem.
		/*
		if (applyBayer)
		{
			Float_t quantumEfficiencyColorUnshift = 1.0 / quantumEfficiency;
			noisedVal *= quantumEfficiencyColorUnshift;
		}
		*/
		//Clip the final value
		if (noisedVal > 1.0)
			noisedVal = 1.0;
		
		int finalVal = round(noisedVal * 256);
		if (finalVal > 255)
		{
			finalVal = 255;
			numSaturatedPixelsD++;
		}
		sMeanOutputValue[inputVal] += finalVal;
		sMeanError[inputVal] += abs(finalVal - inputVal);
		
		//Store the final value in the output image
		grayImage[i] = noisedVal;
	}
	
	for (int i = 0; i < 20; i++)
		a[i] /= totalPixels;
	int idx = -1;
	cout << endl;
	++idx; cout << setw(40) << left << "Mean ppfd: "					<< a[idx] << " of range [" << maxPPFD * noiseParams.minReflectance_ << ", " << maxPPFD * noiseParams.maxReflectance_<< "]" << endl;
	++idx; cout << setw(40) << left << "Mean ppf: "					<< a[idx] << " (" << (a[idx] / noiseParams.fullWellCapacity_ * 100) << "% of full well capacity)" << endl;
	++idx; cout << setw(40) << left << "Mean ideal electron count: "	<< a[idx] << endl;
	++idx; cout << setw(40) << left << "Mean actual electron count: "	<< a[idx] << endl;
	++idx; cout << setw(40) << left << "Mean expected value: "			<< a[idx] << endl;
	++idx; cout << setw(40) << left << "Mean hot pixel: "				<< a[idx] << endl;
	++idx; cout << setw(40) << left << "Mean glow: "					<< a[idx] << " (" << (a[idx] / noiseParams.fullWellCapacity_ * 100) << "% of full well capacity)" << endl;
	++idx; cout << setw(40) << left << "Mean noise sig dep var: "		<< a[idx] << endl;
	++idx; cout << setw(40) << left << "Mean noise std dev: "			<< a[idx] << endl;
	++idx; cout << setw(40) << left << "Mean |noise|: "				<< a[idx] << " (" << (a[idx] / noiseParams.fullWellCapacity_ * 100) << "% of full well capacity)" << endl;
	++idx; cout << setw(40) << left << "Mean noised val: "				<< a[idx] << " of range " << noiseParams.fullWellCapacity_ << endl;
	cout << setw(40) << left << "Num saturated pixels A: "		<< numSaturatedPixelsA << endl;
	cout << setw(40) << left << "Num saturated pixels B: "		<< numSaturatedPixelsB << endl;
	cout << setw(40) << left << "Num saturated pixels C: "		<< numSaturatedPixelsC << endl;
	cout << setw(40) << left << "Num saturated pixels D: "		<< numSaturatedPixelsD << endl;
	
	grayImage.appendFilename("_n-");
	grayImage.appendFilename(noiseParams.label_);
	
	for (int i = 0; i < 256; i++)
		if (sInputHistogram[i] > 0)
		{
			sMeanNoise[i] /= sInputHistogram[i];
			sMeanNoiseProporation[i] /= sInputHistogram[i];
			sMeanOutputValue[i] /= sInputHistogram[i];
			sMeanError[i] /= sInputHistogram[i];
		}
		else assert(sMeanNoise[i] == 0 && sMeanNoiseProporation[i] == 0 && sMeanOutputValue[i] == 0 && sMeanError[i] == 0);
}

void addSensorNoiseWithAutoISO(ImageReal &grayImage, bool applyBayer, NoiseParam noiseParams)
{
	cout << "Adding grayscale sensor noise with autoISO..." << endl;
	
	//Check to see if the image is normalized
	int totalPixels = grayImage.totalPixels();
	Float_t maxVal = 0;
	for (int i = 0; i < totalPixels; i++)
	{
		assert(grayImage[i] >= 0 && grayImage[i] <= 1.0);
		if (grayImage[i] > maxVal)
			maxVal = grayImage[i];
	}
	if (maxVal != 1.0)
		cout << "  ***  NOTE  ***  This input image is not normalized yet" << endl;
	
	cout << "Performing auto-ISO search for proper amplifier gain to achieve good saturation..." << endl;
	
	//Perform "auto-ISO" by altering the amplifier gain until image saturation is achieved which matches the input image.
	//This is not how a digital camera would perform auto-ISO (there is no "input image").
	//A real camera would aim for some rational indication of good saturation,
	//e.g., a mean pixel value of X some fraction of 1.0, or perhaps some ratio Y of clipped pixels to total pixels, something like that.
	//Using the input image as the target saturation ensures that all output images have the same saturation and that they correspond to the input image,
	//which then enables subsequent comparisons across a set of images (PSNR, Lab Delta E, etc.).
	Float_t nonOutlierMeanPixelValOri = calcNonoutlierMeanPixelValue(grayImage, 3);
	cout << "  Input image mean intensity (auto-gain target intensity): " << nonOutlierMeanPixelValOri << endl << endl;
	ImageReal grayImageOri = grayImage;
	Float_t maxPPFD = luxToPPFD(noiseParams.lux_, 0, noiseParams.focalRatio_, noiseParams.pixelPitchMicrons_, noiseParams.fullWellCapacity_);
	noiseParams.ampGain_ = noiseParams.fullWellCapacity_ / (maxPPFD * noiseParams.exposureSecs_);
	Float_t stepSize = noiseParams.ampGain_ * .5;	//Beginning of binary search (well, approximately a binary search at any rate)
	vector<pair<Float_t, Float_t> > ampGainHistory;
	Float_t nonOutlierMeanPixelVal = 0;
	int lastStepDir = 0;
	while (true)
	{
		cout << endl << "--------------------------------------------------------------------------------" << endl;
		
		int historyIdx = 0;
		for (historyIdx = 0; historyIdx < ampGainHistory.size(); historyIdx++)
			if (ampGainHistory[historyIdx].first == noiseParams.ampGain_)
				break;
		
		int numSaturatedPixels = totalPixels;
		if (historyIdx == ampGainHistory.size())
		{
			//Restore the original image
			grayImage = grayImageOri;
			
			//Noisify the image
			cout << endl << setw(40) << left << "  Noisifying with ags: " << noiseParams.ampGain_ << endl;
			addSensorNoise(grayImage, applyBayer, noiseParams);
			
			//Determine the saturation
			nonOutlierMeanPixelVal = calcNonoutlierMeanPixelValue(grayImage, 3);
			cout << endl << "    nonOutlierMeanPixelVal: " << nonOutlierMeanPixelVal << endl;
		}
		else	//This amp gain has already been tried
		{
			maxVal = ampGainHistory[historyIdx].second;
			cout << "  Already tried amp gain " << noiseParams.ampGain_ << ".  maxVal = " << maxVal << endl;
		}
		
		if (nonOutlierMeanPixelVal >= nonOutlierMeanPixelValOri - kAutoGainConvergence
			&& nonOutlierMeanPixelVal <= nonOutlierMeanPixelValOri + kAutoGainConvergence)	//Acceptable saturation
		{
			assert(historyIdx == ampGainHistory.size());
			cout << "    Good saturation achieved" << endl;
			break;
		}
		
		if (stepSize < .001)
		{
			cout << "Bailing on binary search, step size is too small." << endl;
			break;
		}
		
		ampGainHistory.push_back(pair<Float_t, Float_t>(noiseParams.ampGain_, maxVal));
		
		//cout << "    lastStepDir:" << setw(21) << lastStepDir << endl;
		if (nonOutlierMeanPixelVal < nonOutlierMeanPixelValOri - kAutoGainConvergence)	//Increase gain
		{
			if (lastStepDir == -1)
				stepSize *= .5;
			//cout << "    Increasing ags from " << setw(13) << ampGainScalar;
			noiseParams.ampGain_ += stepSize;
			lastStepDir = 1;
			//cout << " to " << setw(10) << ampGainScalar << " with ss " << setw(10) << stepSize << endl;
		}
		else	//Decrease gain
		{
			if (lastStepDir == 1)
				stepSize *= .5;
			//cout << "    Decreasing ags from " << setw(13) << ampGainScalar;
			noiseParams.ampGain_ -= stepSize;
			if (noiseParams.ampGain_ < 1)
			{
				noiseParams.ampGain_ = 1;
				stepSize *= .5;
				//cout << "    Increasing ags from " << setw(13) << ampGainScalar;
				noiseParams.ampGain_ += stepSize;
				lastStepDir = 1;
			}
			lastStepDir = -1;
			//cout << " to " << setw(10) << ampGainScalar << " with ss " << setw(10) << stepSize << endl;
		}
		//cout << "    lastStepDir:" << setw(21) << lastStepDir << endl;
		
		assert(noiseParams.ampGain_ >= 1);//0);
	}
	
	cout << "Final ags: " << noiseParams.ampGain_ << endl;
	/*
	ofstream ofs;
	ofs.open("noiseData.txt");
	for (int i = 0; i < 256; i++)
		ofs << i << '\t' << sMeanNoise[i] << '\t' << sMeanNoiseProporation[i] << '\t'
			<< sMeanOutputValue[i] << '\t' << sMeanError[i] << endl;
	ofs.close();
	*/
}

void addSensorNoise(Image3ChnlReal &rgbImage, bool applyBayer, NoiseParam noiseParams)
{
	cout << "Adding color sensor noise..." << endl;
	
	//Check to see if the image is normalized
	int totalPixels = rgbImage.totalPixels();
	Float_t maxVal = 0;
	for (int c = 0; c < 3; c++)
		for (int i = 0; i < totalPixels; i++)
		{
			assert((*rgbImage[c])[i] >= 0 && (*rgbImage[c])[i] <= 1.0);
			if ((*rgbImage[c])[i] > maxVal)
				maxVal = (*rgbImage[c])[i];
		}
	if (maxVal != 1.0)
		cout << "***  NOTE  ***  This input image is not normalized yet" << endl;
	
	//Noisify the image
	cout << endl << setw(40) << left << "  Noisifying with ags: " << noiseParams.ampGain_ << endl;
	for (int c = 0; c < 3; c++)
		addSensorNoise(*rgbImage[c], applyBayer, noiseParams);
	
	rgbImage.appendFilename("_n-");
	rgbImage.appendFilename(noiseParams.label_);
}

void addSensorNoiseWithAutoISO(Image3ChnlReal &rgbImage, bool applyBayer, NoiseParam noiseParams)
{
	cout << "Adding color sensor noise with autoISO..." << endl;
	
	//Check to see if the image is normalized
	int totalPixels = rgbImage.totalPixels();
	Float_t maxVal = 0;
	for (int c = 0; c < 3; c++)
		for (int i = 0; i < totalPixels; i++)
		{
			assert((*rgbImage[c])[i] >= 0 && (*rgbImage[c])[i] <= 1.0);
			if ((*rgbImage[c])[i] > maxVal)
				maxVal = (*rgbImage[c])[i];
		}
	if (maxVal != 1.0)
		cout << "  ***  NOTE  ***  This input image is not normalized yet" << endl;
	
	cout << "Performing auto-ISO search for proper amplifier gain to achieve good saturation..." << endl;
	
	//Perform "auto-ISO" by altering the amplifier gain until image saturation is achieved which matches the input image.
	//This is not how a digital camera would perform auto-ISO (there is no "input image").
	//A real camera would aim for some rational indication of good saturation,
	//e.g., a mean pixel value of X some fraction of 1.0, or perhaps some ratio Y of clipped pixels to total pixels, something like that.
	//Using the input image as the target saturation ensures that all output images have the same saturation and that they correspond to the input image,
	//which then enables subsequent comparisons across a set of images (PSNR, Lab Delta E, etc.).
	Float_t nonOutlierMeanPixelValOri = calcNonoutlierMeanPixelValue(rgbImage, 3);
	cout << "  Input image mean intensity (auto-gain target intensity): " << nonOutlierMeanPixelValOri << endl << endl;
	Image3ChnlReal rgbImageOri = rgbImage;
	Float_t maxPPFD = luxToPPFD(noiseParams.lux_, 0, noiseParams.focalRatio_, noiseParams.pixelPitchMicrons_, noiseParams.fullWellCapacity_);
	noiseParams.ampGain_ = noiseParams.fullWellCapacity_ / (maxPPFD * noiseParams.exposureSecs_);
	Float_t stepSize = noiseParams.ampGain_ * .5;	//Beginning of binary search (well, approximately a binary search at any rate)
	vector<pair<Float_t, Float_t> > ampGainHistory;
	Float_t nonOutlierMeanPixelVal = 0;
	int lastStepDir = 0;
	while (true)
	{
		cout << endl << "--------------------------------------------------------------------------------" << endl;
		
		int historyIdx = 0;
		for (historyIdx = 0; historyIdx < ampGainHistory.size(); historyIdx++)
			if (ampGainHistory[historyIdx].first == noiseParams.ampGain_)
				break;
		
		int numSaturatedPixels = totalPixels;
		if (historyIdx == ampGainHistory.size())
		{
			//Restore the original image
			rgbImage = rgbImageOri;
			
			//Noisify the image
			cout << endl << setw(40) << left << "  Noisifying with ags: " << noiseParams.ampGain_ << endl;
			for (int c = 0; c < 3; c++)
				addSensorNoise(*rgbImage[c], applyBayer, noiseParams);
			
			//Determine the saturation
			nonOutlierMeanPixelVal = calcNonoutlierMeanPixelValue(rgbImage, 3);
			cout << endl << "    nonOutlierMeanPixelVal: " << nonOutlierMeanPixelVal << endl;
		}
		else	//This amp gain has already been tried
		{
			maxVal = ampGainHistory[historyIdx].second;
			cout << "  Already tried amp gain " << noiseParams.ampGain_ << ".  maxVal = " << maxVal << endl;
		}
		
		if (nonOutlierMeanPixelVal >= nonOutlierMeanPixelValOri - kAutoGainConvergence
			&& nonOutlierMeanPixelVal <= nonOutlierMeanPixelValOri + kAutoGainConvergence)	//Acceptable saturation
		{
			assert(historyIdx == ampGainHistory.size());
			cout << "    Good saturation achieved" << endl;
			break;
		}
		
		if (stepSize < .001)
		{
			cout << "Bailing on binary search, step size is too small." << endl;
			break;
		}
		
		ampGainHistory.push_back(pair<Float_t, Float_t>(noiseParams.ampGain_, maxVal));
		
		//cout << "    lastStepDir:" << setw(21) << lastStepDir << endl;
		if (nonOutlierMeanPixelVal < nonOutlierMeanPixelValOri - kAutoGainConvergence)	//Increase gain
		{
			if (lastStepDir == -1)
				stepSize *= .5;
			//cout << "    Increasing ags from " << setw(13) << ampGainScalar;
			noiseParams.ampGain_ += stepSize;
			lastStepDir = 1;
			//cout << " to " << setw(10) << ampGainScalar << " with ss " << setw(10) << stepSize << endl;
		}
		else	//Decrease gain
		{
			if (lastStepDir == 1)
				stepSize *= .5;
			//cout << "    Decreasing ags from " << setw(13) << ampGainScalar;
			noiseParams.ampGain_ -= stepSize;
			if (noiseParams.ampGain_ < 1)
			{
				noiseParams.ampGain_ = 1;
				cout << "Bailing on binary search, auto-gain descended below 1.0." << endl;
				break;
			}
			lastStepDir = -1;
			//cout << " to " << setw(10) << ampGainScalar << " with ss " << setw(10) << stepSize << endl;
		}
		//cout << "    lastStepDir:" << setw(21) << lastStepDir << endl;
		
		assert(noiseParams.ampGain_ >= 1);//0);
	}
	
	cout << "Final ags: " << noiseParams.ampGain_ << endl;
	
	rgbImage.appendFilename("_n-");
	rgbImage.appendFilename(noiseParams.label_);
}

//======================================================================
#pragma mark -

void denoiseGrayImage(string filename, string manualShrinkageThresholdStr)
{
	cout << "Denoising gray image..." << endl;
	
	//Open the image
	ImageReal grayImage;
	grayImage.readFile(filename);
	
	//Construct a steerable pyramid of the image
	SteerablePyramid steerPyr;
	steerPyr.DoFullWaveletAnalysis(grayImage.img(), grayImage.width(), grayImage.height(), true);
	
	if (gGenerateIntermediateImages)
		steerPyr.GeneratePGMimages(grayImage.filename().c_str());
	
	int numLevels = steerPyr.NumLevels();
	cout << "Steerable pyramid will be built with " << numLevels << " levels" << endl;
	if (gManualShrinkageThresholds[0].size() != numLevels + 1)
		gManualShrinkageThresholds[0].resize(numLevels + 1);
	cout << "Thresholds:" << endl;
	cout << "  H:  " << setw(10) << gManualShrinkageThresholds[0][0] << endl;
	cout << "  B:  ";
	for (int i = 1; i < gManualShrinkageThresholds[0].size(); i++)
		cout << setw(10) << gManualShrinkageThresholds[0][i];
	cout << endl << endl;
	
	//Shrink the coefficients
	steerPyr.SoftThreshold(gManualShrinkageThresholds[0]);
	
	if (gGenerateIntermediateImages)
	{
		string s = grayImage.filename() + "_sk";
		steerPyr.GeneratePGMimages(s.c_str());
	}
	
	//Synthesize the steerable pyramid
	ImageReal grayImageShrunk(grayImage.width(), grayImage.height());
	steerPyr.DoFullWaveletSynthesis(grayImageShrunk.img());
	/*
	stringstream denoiseParamSS;
	denoiseParamSS << "_";
	for (int c = 0; c < 3; c++)
	{
		if (c == 0)
			denoiseParamSS << "R";
		else if (c == 1)
			denoiseParamSS << "G";
		else denoiseParamSS << "B";
		
		denoiseParamSS << gManualShrinkageThresholds[c][0];
		for (int i = 1; i < gManualShrinkageThresholds[c].size(); i++)
			denoiseParamSS << "," << gManualShrinkageThresholds[c][i];
	}
	*/
	grayImageShrunk.setFilename(grayImage.filename());
	grayImageShrunk.appendFilename("_dn");
	grayImageShrunk.appendFilename(manualShrinkageThresholdStr);
	grayImageShrunk.writeFile(1);
}

void denoiseRGBImage(string filename, string manualShrinkageThresholdStr)
{
	cout << "Denoising color image..." << endl;
	
	//Open the image
	Image3ChnlReal rgbImage;
	rgbImage.readFile(filename);
	
	//Construct a steerable pyramid of the image
	SteerablePyramid steerPyr[3];
	for (int c = 0; c < 3; c++)
		steerPyr[c].DoFullWaveletAnalysis(rgbImage[c]->img(), rgbImage[c]->width(), rgbImage[c]->height(), true);
	
	if (gGenerateIntermediateImages)
	{
		for (int c = 0; c < 3; c++)
		{
			string s = rgbImage.filename();
			if (c == 0)
				s += "_R";
			else if (c == 1)
				s += "_G";
			else s += "_B";
			steerPyr[c].GeneratePGMimages(s.c_str());
		}
	}
	
	int numLevels = steerPyr[0].NumLevels();
	cout << "Steerable pyramid will be built with " << numLevels << " levels" << endl;
	if (gManualShrinkageThresholds[0].size() != numLevels + 1)
		for (int c = 0; c < 3; c++)
			gManualShrinkageThresholds[c].resize(numLevels + 1);
	cout << "Thresholds:" << endl;
	for (int c = 0; c < 3; c++)
	{
		if (c == 0)
			cout << "  Red:" << endl;
		else if (c == 1)
			cout << "  Green:" << endl;
		else cout << "  Blue:" << endl;
		
		cout << "    H:  " << setw(10) << gManualShrinkageThresholds[c][0] << endl;
		cout << "    B:  ";
		for (int i = 1; i < gManualShrinkageThresholds[c].size(); i++)
			cout << setw(10) << gManualShrinkageThresholds[c][i];
		cout << endl;
	}
	cout << endl;
	
	//Shrink the coefficients
	for (int c = 0; c < 3; c++)
		steerPyr[c].SoftThreshold(gManualShrinkageThresholds[c]);
	
	if (gGenerateIntermediateImages)
	{
		for (int c = 0; c < 3; c++)
		{
			string s = rgbImage.filename();
			if (c == 0)
				s += "_R_sk";
			else if (c == 1)
				s += "_G_sk";
			else s += "_B_sk";
			steerPyr[c].GeneratePGMimages(s.c_str());
		}
	}
	
	//Synthesize the steerable pyramid
	Image3ChnlReal rgbImageShrunk(rgbImage.width(), rgbImage.height());
	for (int c = 0; c < 3; c++)
		steerPyr[c].DoFullWaveletSynthesis(rgbImageShrunk[c]->img());
	/*
	stringstream denoiseParamSS;
	denoiseParamSS << "_";
	for (int c = 0; c < 3; c++)
	{
		if (c == 0)
			denoiseParamSS << "R";
		else if (c == 1)
			denoiseParamSS << "G";
		else denoiseParamSS << "B";
		
		denoiseParamSS << gManualShrinkageThresholds[c][0];
		for (int i = 1; i < gManualShrinkageThresholds[c].size(); i++)
			denoiseParamSS << "," << gManualShrinkageThresholds[c][i];
	}
	*/
	rgbImageShrunk.setFilename(rgbImage.filename());
	rgbImageShrunk.appendFilename("_dn");
	rgbImageShrunk.appendFilename(manualShrinkageThresholdStr);
	rgbImageShrunk.writeFile(1);
}

void diffImages(const ImageReal &img1, const ImageReal &img2, ImageReal &diff)
{
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	int width = img1.width(), height = img1.height();
	diff.resize(width, height);
	int totalPixels = img1.totalPixels();
	assert(totalPixels = width * height);
	
	Float_t maxDif = 0;
	for (int i = 0; i < totalPixels; i++)
	{
		diff[i] = fabs(img1[i] - img2[i]);
		if (diff[i] > maxDif)
			maxDif = diff[i];
	}
	
	//cout << "Max dif = " << maxDif << "        x 255 = " << maxDif * 255.0 << "        x 65535 = " << maxDif * 65535.0 << endl;
	
	diff.setFilename(img1.filename());
	diff.appendFilename("_dif");
}

void diffImages(const Image3ChnlReal &img1, const Image3ChnlReal &img2, Image3ChnlReal &diff)
{
	assert(img1.width() == img2.width() && img1.height() == img2.height());
	int width = img1.width(), height = img1.height();
	diff.resize(width, height);
	int totalPixels = img1.totalPixels();
	assert(totalPixels = width * height);
	
	for (int c = 0; c < 3; c++)
		diffImages(*img1[c], *img2[c], *diff[c]);
	
	diff.setFilename(img1.filename());
	diff.appendFilename("_dif");
}

#ifndef __IMAGE_PROCS__
#define __IMAGE_PROCS__

#include "FloatType.h"
#include "ImageReal.h"
#include "Image3ChnlReal.h"
#include "Bayer.h"
#include <string>
#include <iostream>
#include <iomanip>

struct NoiseParam
{
	//Filename suffix label
	std::string label_;
	
	//--------------------------------------------------------------------------------
	//Environmental Properties
	
	//Scene illumination (lux)
	Float_t lux_;
	
	//Scene reflectance corresponding input pixels with minimum (0) and maximum (1) value
	Float_t minReflectance_;
	Float_t maxReflectance_;
	
	//Dark current rate and hot pixel severity will double with every 6 degrees C
	Float_t temperatureC_;
	
	//--------------------------------------------------------------------------------
	//Camera Properties
	
	//Lens/telescope focal length
	//Float_t focalLengthMM_;
	
	//Pixel pitch (size along one dimension).  Assuming a square pixel, surface area will then be the square of this parameter of course.
	Float_t pixelPitchMicrons_;
	
	//Most CCDs hold between 700 to 1100 electrons per square micron, although this ought to additionally vary with well depth.
	//Typical small webcams have 5.6um pixels, so about 30,000 full well capacity.  Digital cameras, especially high-end cameras
	//have much larger pixels both in surface area and depth.
	unsigned long fullWellCapacity_;
	
	//Read noise.  Mean electroncs generated per pixel from circuit noise.
	Float_t readNoise_;
	
	//Bias (quantum efficiency, vignetting, and dust scalars).
	//Quantum efficiency is a per-pixel ratio of electron-flux to photon-flux strictly < 1.  Note this is generally highly wavelength dependent as well.
	//Vignetting represents a "brighter" response at the center of the CCD due to a lens effect with a presumably Gaussian drop-off.
	//Flecks of dust on the optics appear as dim fuzzy circular (or annulus in the case of telescopes with a secondary mirror) spots on the CCD.
	//ImageReal *biases_;
	
	//How noisy is the amplifier
	Float_t ampVarFactor_;
	
	//Dark current
	Float_t darkCurrentRoomTempElectronsPerPixPerSec_;
	
	//How severe is hot pixel noise
	Float_t hotPixelFactor_;
	
	//A smattering of bright pixels across the CCD
	ImageReal *hotPixelScalars_;
	
	//Quantum efficiency at three "wavelengths", actually, three bands, approximately 400-500nm, 500-600nm, and 600-700nm
	Float_t quantumEfficiency_[3];
	
	//How quickly does radial amp glow (from the corner of the CCD) grow w.r.t. exposure time.
	//This is a single value, but ideally it should be tapered at 1/radius from one corner of the CCD.
	Float_t ampGlowElectronsPerPixPerSec_;
	
	//--------------------------------------------------------------------------------
	//Camera settings
	
	Float_t exposureSecs_;
	Float_t ampGain_;
	Float_t focalRatio_;
	
	//--------------------------------------------------------------------------------
	
	NoiseParam() :
		label_(""),
		lux_(0),
		minReflectance_(0),
		maxReflectance_(0),
		temperatureC_(0),
		//focalLengthMM_(0),
		pixelPitchMicrons_(0),
		fullWellCapacity_(0),
		readNoise_(0),
		//biases_(NULL),
		ampVarFactor_(0),
		darkCurrentRoomTempElectronsPerPixPerSec_(0),
		hotPixelFactor_(0),
		hotPixelScalars_(NULL),
		ampGlowElectronsPerPixPerSec_(0),
		exposureSecs_(0),
		ampGain_(0),
		focalRatio_(0)
	{
		quantumEfficiency_[R] = quantumEfficiency_[G] = quantumEfficiency_[B] = 0;
	}
	
	NoiseParam(const NoiseParam &np) :
		label_(""),
		lux_(0),
		minReflectance_(0),
		maxReflectance_(0),
		temperatureC_(0),
		//focalLengthMM_(0),
		pixelPitchMicrons_(0),
		fullWellCapacity_(0),
		readNoise_(0),
		//biases_(NULL),
		ampVarFactor_(0),
		darkCurrentRoomTempElectronsPerPixPerSec_(0),
		hotPixelFactor_(0),
		hotPixelScalars_(NULL),
		ampGlowElectronsPerPixPerSec_(0),
		exposureSecs_(0),
		ampGain_(0),
		focalRatio_(0)
	{
		quantumEfficiency_[R] = quantumEfficiency_[G] = quantumEfficiency_[B] = 0;
		
		*this = np;
		
		selfCheck();
	}
	
	~NoiseParam()
	{
		//delete biases_;
		delete hotPixelScalars_;
	}
	
	const NoiseParam& operator=(const NoiseParam &np)
	{
		if (&np == this)
			return *this;
		
		set(
			np.label_.c_str(),
			np.lux_,
			np.minReflectance_,
			np.maxReflectance_,
			np.temperatureC_,
			//np.focalLengthMM_,
			np.pixelPitchMicrons_,
			np.fullWellCapacity_,
			np.readNoise_,
			np.ampVarFactor_,
			np.darkCurrentRoomTempElectronsPerPixPerSec_,
			np.hotPixelFactor_,
			np.hotPixelScalars_,
			np.quantumEfficiency_,
			np.ampGlowElectronsPerPixPerSec_,
			np.exposureSecs_,
			np.ampGain_,
			np.focalRatio_
		);
		
		selfCheck();
		
		return *this;
	}
	
	void set(
		const char *label,
		Float_t lux,
		Float_t minReflectance,
		Float_t maxReflectance,
		Float_t temperatureC,
		//Float_t focalLengthMM,
		Float_t pixelPitchMicrons,
		unsigned long fullWellCapacity,
		Float_t readNoise,
		//ImageReal *bias,
		Float_t ampVarFactor,
		Float_t darkCurrentRoomTempElectronsPerPixPerSec,
		Float_t hotPixelFactor,
		ImageReal *hotPixels,
		const Float_t quantumEfficiency[3],
		Float_t ampGlowElectronsPerPixPerSec,
		Float_t exposureSecs,
		Float_t ampGain,
		Float_t lensFnumber
	)
	{
		label_ = label;
		lux_ = lux;
		minReflectance_ = minReflectance;
		maxReflectance_ = maxReflectance;
		temperatureC_ = temperatureC;
		//focalLengthMM_ = focalLengthMM;
		pixelPitchMicrons_ = pixelPitchMicrons;
		fullWellCapacity_ = fullWellCapacity;
		readNoise_ = readNoise;
		//if (biases_)
		//	*biases_ = *bias;
		//else biases_ = new ImageReal(*bias);
		ampVarFactor_ = ampVarFactor;
		darkCurrentRoomTempElectronsPerPixPerSec_ = darkCurrentRoomTempElectronsPerPixPerSec;
		hotPixelFactor_ = hotPixelFactor;
		delete hotPixelScalars_;
		hotPixelScalars_ = NULL;
		if (hotPixels)
			hotPixelScalars_ = new ImageReal(*hotPixels);
		quantumEfficiency_[R] = quantumEfficiency[R];
		quantumEfficiency_[G] = quantumEfficiency[G];
		quantumEfficiency_[B] = quantumEfficiency[B];
		ampGlowElectronsPerPixPerSec_ = ampGlowElectronsPerPixPerSec;
		exposureSecs_ = exposureSecs;
		ampGain_ = ampGain;
		focalRatio_ = lensFnumber;
		
		selfCheck();
	}
	
	void selfCheck()
	{
		assert(label_.length() > 0);
		assert(lux_ > 0);
		assert(minReflectance_ >= 0);
		assert(maxReflectance_ > minReflectance_ && maxReflectance_ <= 1);
		assert(temperatureC_ >= -100 && temperatureC_ <= 50);
		//assert(focalLengthMM_ > 0);
		assert(pixelPitchMicrons_ > 0);
		assert(fullWellCapacity_ >= 1);
		assert(readNoise_ > 0 && readNoise_ < fullWellCapacity_);
		//Assert that all bias pixels are >= 0 and <= 1
		assert(ampVarFactor_ >= 0);
		assert(darkCurrentRoomTempElectronsPerPixPerSec_ >= 0);
		assert(hotPixelFactor_ >= 0);
		assert(quantumEfficiency_[R] >= 0 && quantumEfficiency_[R] <= 1);
		assert(quantumEfficiency_[G] >= 0 && quantumEfficiency_[G] <= 1);
		assert(quantumEfficiency_[B] >= 0 && quantumEfficiency_[B] <= 1);
		assert(ampGlowElectronsPerPixPerSec_ >= 0);
		assert(exposureSecs_ > 0);
		assert(ampGain_ >= 0 || ampGain_ == -1);
		assert(focalRatio_ > 0);
	}
	
	void dump()
	{
		std::cout << "NoiseParam:" << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "label_:" << label_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "lux_:" << lux_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "minReflectance_:" << minReflectance_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "maxReflectance_:" << maxReflectance_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "temperatureC_:" << temperatureC_ << std::endl;
		//std::cout << "  " << std::setw(45) << std::left << "focalLengthMM_:" << focalLengthMM_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "pixelPitchMicrons_:" << pixelPitchMicrons_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "fullWellCapacity_:" << fullWellCapacity_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "readNoise_:" << readNoise_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "ampVarFactor_:" << ampVarFactor_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "darkCurrentRoomTempElectronsPerPixPerSec_:" << darkCurrentRoomTempElectronsPerPixPerSec_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "hotPixelFactor_:" << hotPixelFactor_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "quantumEfficiency_:" << quantumEfficiency_[R] << "  " << quantumEfficiency_[G] << "  " << quantumEfficiency_[B] << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "ampGlowElectronsPerPixPerSec_:" << ampGlowElectronsPerPixPerSec_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "exposureSecs_:"		<< exposureSecs_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "ampGain_:" << ampGain_ << std::endl;
		std::cout << "  " << std::setw(45) << std::left << "focalRatio_:" << focalRatio_ << std::endl;
	}
};

extern int gBitDepth;

//======================================================================

void rgb2xyz(Float_t r, Float_t g, Float_t b,
			Float_t &x, Float_t &y, Float_t &z);

void xyz2rgb(Float_t x, Float_t y, Float_t z,
			Float_t &r, Float_t &g, Float_t &b);

void xyz2cielab(Float_t x, Float_t y, Float_t z,
				Float_t &l, Float_t &a, Float_t &b);

void cielab2xyz(Float_t l, Float_t a, Float_t b,
				Float_t &x, Float_t &y, Float_t &z);

void testRGBXYZLABRoutines();

void printColor(Float_t r, Float_t g, Float_t b);

void rgbImgToLabImg(Float_t maxVal, const Image3ChnlReal &rgbImg, Image3ChnlReal &labImg);
void labImgToRgbImg(Float_t maxVal, const Image3ChnlReal &labImg, Image3ChnlReal &rgbImg);

//======================================================================
#pragma mark -

//A detailed description of the noise model is available at the top of the function body
void addSensorNoise(ImageReal &grayImage, bool applyBayer, NoiseParam noiseParams);
void addSensorNoiseWithAutoISO(ImageReal &grayImage, bool applyBayer, NoiseParam noiseParams);
void addSensorNoise(Image3ChnlReal &rgbImage, bool applyBayer, NoiseParam noiseParams);
void addSensorNoiseWithAutoISO(Image3ChnlReal &rgbImage, bool applyBayer, NoiseParam noiseParams);

void denoiseGrayImage(std::string filename, std::string manualShrinkageThresholdStr);
void denoiseRGBImage(std::string filename, std::string manualShrinkageThresholdStr);

void diffImages(const ImageReal &img1, const ImageReal &img2, ImageReal &diff);
void diffImages(const Image3ChnlReal &img1, const Image3ChnlReal &img2, Image3ChnlReal &diff);

#endif

#include "SteerablePyramid.h"
#include "FourierTransform.h"
#include "ImageReal.h"	//DEBUG
#include "main.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

//******************************************************************************
//Extern Globals

//******************************************************************************
//Global Declarations

static const string skFilterSetFilename = "steerablePyramidFilterSet.txt";	//THIS FILE BETTER BE THERE!
static int skFilterDims[7] = { 0 };	//Effectively const, after initialization

//An array of filters
//array[filterId][filterCoefficient]
static Float_t *skFilters[7] = { NULL };	//Effectively const, after initialization
static Float_t *skFiltersInv[7] = { NULL };	//Effectively const, after initialization

//This struct must be defined below the declaration of skFilterDims
struct FilterFT
{
	int width_, height_;
	
	//Seven arrays of vectors.
	//Each vector is of levels and stores arrays of coefficients.
	vector<Float_t*> filtersFTs[7];
	vector<Float_t*> filtersInvFTs[7];
	
	FilterFT() : width_(0), height_(0)
	{}
	
	FilterFT(const FilterFT &fft)
	{
		*this = fft;
	}
	
	~FilterFT()
	
	{
		for (int i = 0; i < 7; i++)
			for (int j = 0; j < filtersFTs[i].size(); j++)
			{
				assert(filtersFTs[i][j]);
				delete [] filtersFTs[i][j];
				delete [] filtersInvFTs[i][j];
			}
	}
	
	const FilterFT& operator=(const FilterFT &fft)
	{
		if (&fft == this)
			return *this;
		
		width_ = fft.width_;
		height_ = fft.height_;
		
		for (int i = 0; i < 7; i++)
		{
			filtersFTs[i].resize(fft.filtersFTs[i].size());
			filtersInvFTs[i].resize(fft.filtersInvFTs[i].size());
			for (int j = 0; j < filtersFTs[i].size(); j++)
			{
				int totalPixels = (width_ * height_) / pow(4, j);
				filtersFTs[i][j] = new Float_t[totalPixels];
				filtersInvFTs[i][j] = new Float_t[totalPixels];
				memcpy(filtersFTs[i][j], fft.filtersFTs[i][j], totalPixels * sizeof(Float_t));
				memcpy(filtersInvFTs[i][j], fft.filtersInvFTs[i][j], totalPixels * sizeof(Float_t));
			}
		}
		
		return *this;
	}
	
	void init()
	{
		assert(filtersFTs[0].size() == 0);
		
		int numLevels = 0;
		int shorterDim = min(width_, height_);
		while (shorterDim > skFilterDims[2])	//Limit is determined by the dimensions of the L1 filter (probably 17x17)
		{
			numLevels++;
			shorterDim /= 2;
		}
		
		for (int i = 0; i < 7; i++)
		{
			filtersFTs[i].resize(numLevels);
			filtersInvFTs[i].resize(numLevels);
			for (int j = 0; j < numLevels; j++)
			{
				int totalPixels = (width_ * height_) / pow(4, j);
				filtersFTs[i][j] = new Float_t[totalPixels];
				filtersInvFTs[i][j] = new Float_t[totalPixels];
				
				//Make absolutely certain the space is cleared to zeroes since much of it won't be overwritten during initialization.
				memset(filtersFTs[i][j], 0, totalPixels * sizeof(Float_t));
				memset(filtersInvFTs[i][j], 0, totalPixels * sizeof(Float_t));
			}
		}
	}
	
	void resize(int width, int height)
	{
		width_ = width;
		height_ = height;
		init();
	}
};
static vector<FilterFT> skFiltersFTs;	//Effectively const, after initialization

//Profiling timers, milliseconds, see main.h and getCurrentTime_ms()
static unsigned long gT1, gT2, gT3, gT4, gT5;

Float_t *SteerablePyramid::scratch_2D_A_ = NULL;
Float_t *SteerablePyramid::scratch_2D_B_ = NULL;
int SteerablePyramid::scratchSize_ = 0;

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

SteerablePyramid::SteerablePyramid() :
	h0Coeffs_(NULL),
	b1Coeffs_(NULL), b2Coeffs_(NULL), b3Coeffs_(NULL), b4Coeffs_(NULL), l1Coeffs_(NULL),
	width_(0), height_(0), numLevels_(0),
	yLookup_(NULL),
	useFTconvolution_(false)
{
	InitFilters();
	CreateScratch();
}

SteerablePyramid::SteerablePyramid(const SteerablePyramid &sp) :
	h0Coeffs_(NULL),
	b1Coeffs_(NULL), b2Coeffs_(NULL), b3Coeffs_(NULL), b4Coeffs_(NULL), l1Coeffs_(NULL),
	width_(0), height_(0), numLevels_(0),
	yLookup_(NULL),
	useFTconvolution_(false)
{
	InitFilters();
	CreateScratch();
	
	*this = sp;
}

SteerablePyramid::~SteerablePyramid()
{
	Destroy();
}

#pragma mark -

const SteerablePyramid& SteerablePyramid::operator=(const SteerablePyramid& sp)
{
	if (&sp == this)
		return *this;
	
	Init(sp.width_, sp.height_);
	
	useFTconvolution_ = sp.useFTconvolution_;
	
	int totalPixels = width_ * height_;
	memcpy(h0Coeffs_, sp.h0Coeffs_, totalPixels * sizeof(Float_t));
	
	for (int i = 0; i < numLevels_; i++)
	{
		int totalPixels = (width_ * height_) / (int)(pow(4, i));
		memcpy(b1Coeffs_[i], sp.b1Coeffs_[i], totalPixels * sizeof(Float_t));
		memcpy(b2Coeffs_[i], sp.b2Coeffs_[i], totalPixels * sizeof(Float_t));
		memcpy(b3Coeffs_[i], sp.b3Coeffs_[i], totalPixels * sizeof(Float_t));
		memcpy(b4Coeffs_[i], sp.b4Coeffs_[i], totalPixels * sizeof(Float_t));
		memcpy(l1Coeffs_[i], sp.l1Coeffs_[i], (totalPixels / 4) * sizeof(Float_t));
	}
	
	return *this;
}

bool SteerablePyramid::operator==(const SteerablePyramid& sp) const
{
	if (&sp == this)
		return true;
	
	if (width_ != sp.width_ || height_ != sp.height_ || numLevels_ != sp.numLevels_)
		return false;
	
	for (int lvl = 0; lvl < numLevels_; lvl++)
	{
		int totalPixels = (width_ * height_) / (int)(pow(4, lvl));
		if (lvl == 0)
		{
			for (int i = 0; i < totalPixels; i++)
				if (h0Coeffs_[i] != sp.h0Coeffs_[i])
					return false;
		}
		
		for (int i = 0; i < totalPixels; i++)
		{
			if (b1Coeffs_[lvl][i] != sp.b1Coeffs_[lvl][i])
				return false;
			if (b2Coeffs_[lvl][i] != sp.b2Coeffs_[lvl][i])
				return false;
			if (b3Coeffs_[lvl][i] != sp.b3Coeffs_[lvl][i])
				return false;
			if (b4Coeffs_[lvl][i] != sp.b4Coeffs_[lvl][i])
				return false;
		}
		
		totalPixels /= 4;
		for (int i = 0; i < totalPixels; i++)
			if (l1Coeffs_[lvl][i] != sp.l1Coeffs_[lvl][i])
				return false;
	}
	
	return true;
}

//static
void SteerablePyramid::readFilterSet()
{
	ifstream ifs;
	ifs.open(skFilterSetFilename.c_str());
	
	if (!ifs)
	{
		string errMsg = "Could not find input file '" + skFilterSetFilename + "'.  Please make sure it is in the same directory as the executable.";
		die(errMsg.c_str());
	}
	
	for (int filter = 0; filter < 7; filter++)
	{
		string label;
		int filterDim;
		ifs >> label >> filterDim;
		
		assert(label == "H0" || label == "L0" || label == "L1" || label == "B1" || label == "B2" || label == "B3" || label == "B4");
		assert(filterDim > 0 && filterDim <= 25 && filterDim % 2 == 1);
		
		skFilterDims[filter] = filterDim;
		int numTaps = filterDim * filterDim;
		skFilters[filter] = new Float_t[numTaps];
		skFiltersInv[filter] = new Float_t[numTaps];
		
		for (int i = 0; i < numTaps; i++)
		{
			ifs >> skFilters[filter][i];
			assert(ifs);
		}
	}
	assert(ifs);
	ifs.get();
	assert(!ifs);
	
	ifs.close();
}

//Init Filters by dividing by 10000 and by creating inverse filters
//static
void SteerablePyramid::InitFilters()
{
	//Only init the filters once
	static bool alreadyInited = false;
	if (alreadyInited)
		return;
	alreadyInited = true;
	
	readFilterSet();
	
	//Normalize the filters
	//This is definitely required for the first filter set, which is defined using integers
	//I am unsure why the third filter set used an L0 and L1 that sum to 2 instead of 1.  Perhaps normalizing them is not correct.
	/*
	for (int filter = 0; filter < 7; filter++)
	{
		Float_t divisor = 0;
		int numTaps = skFilterDims[filter] * skFilterDims[filter];
		for (int i = 0; i < numTaps; i++)
			divisor += fabs(skFilters[filter][i]);
		cout << "Filter " << filter << " has " << numTaps << " taps with sum: " << divisor << endl;
		for (int i = 0; i < numTaps; i++)
			skFilters[filter][i] /= divisor;
	}
	cout << endl;
	*/
	if (skFilterSetFilename == "filterSet1.txt")
	{
		//The filters are integers with power 10000.
		//Therefore they must be scaled by 1/10000 to achieve unit power.
		//Note, however, that L0 and L1 do not quite sum to 10000.
		//Therefore, they are scaled slightly differently.
		for (int filter = 0; filter < 7; filter++)
		{
			for (int i = 0; i < 289; i++)
			{
				if (filter == 1)
					skFilters[filter][i] /= 9996.0;
				else if (filter == 2)
					skFilters[filter][i] /= 10005.0;
				else skFilters[filter][i] /= 10000.0;
			}
		}
	}
	if (skFilterSetFilename == "steerablePyramidFilterSet.txt")
	{
		//H0 sums to 0
		//L0 and L1 sum to 2
		//B1-4 sum to 0
		//To get correct reconstruction, it is necessary to scale just L1 to unit power but to leave L0 at power 2.  Why?
		for (int i = 0; i < 289; i++)
			skFilters[2][i] /= 2.0;
	}
	
	//Define the inverted filters
	for (int filter = 0; filter < 7; filter++)
	{
		int numTaps = skFilterDims[filter] * skFilterDims[filter];
		for (int i = 0; i < numTaps; i++)
			skFiltersInv[filter][i] = skFilters[filter][(numTaps - 1) - i];
	}
}

int SteerablePyramid::getFilterFTIndex()
{
	assert(skFilterDims[0] > 0);
	
	//Search the currently defined filter FTs for one with the appropriate dimensions.
	//If one is found, return its index.
	for (int i = 0; i < skFiltersFTs.size(); i++)
		if (skFiltersFTs[i].width_ == width_ && skFiltersFTs[i].height_ == height_)
			return i;
	
	//No filter FT with the appropriate dimensions was found, so add a new one and return its index.
	FilterFT filterFT;
	skFiltersFTs.push_back(filterFT);
	skFiltersFTs.back().resize(width_, height_);
	
	Float_t *shiftBuffer = new Float_t[width_ * height_];
	for (int i = 0; i < 7; i++)
		for (int lvl = 0 ; lvl < numLevels_; lvl++)
		{	
			int filterDim = skFilterDims[i];
			
			int horDim = width_ / pow(2, lvl);
			int verDim = height_ / pow(2, lvl);
			assert(horDim >= filterDim && verDim >= filterDim);
			
			memset(shiftBuffer, 0, horDim * verDim * sizeof(Float_t));
			
			int xOffset = (horDim - filterDim) / 2;
			int yOffset = (verDim - filterDim) / 2;
			
			//By default, an embedded odd-length signal will be positioned such that the center pixel
			//resides at lower right, but I prefer that it be at upper left, so shift it by one.
			assert(filterDim % 2 == 1);	//All the filters should have odd dimensions
			if (filterDim % 2 == 1)
			{
				xOffset++;
				yOffset++;
			}
			
			//Copy the filter to the temp buffer 
			for (int y = 0; y < filterDim; y++)
				for (int x = 0; x < filterDim; x++)
					shiftBuffer[(y + yOffset) * horDim + (x + xOffset)] =
						skFilters[i][y * filterDim + x];
			
			//Shift (with wrap) the center to the corners
			for (int y = 0; y < verDim; y++)
				for (int x = 0; x < horDim; x++)
					skFiltersFTs.back().filtersFTs[i][lvl][y * horDim + x] =
						shiftBuffer[((y + verDim / 2) % verDim) * horDim + (x + horDim / 2) % horDim];
			
			//Do the Fourier analysis
			assert(doFourierAnalysis(skFiltersFTs.back().filtersFTs[i][lvl],
									skFiltersFTs.back().filtersFTs[i][lvl], horDim, verDim));
			
			//Copy the inverse filter to the temp buffer 
			for (int y = 0; y < filterDim; y++)
				for (int x = 0; x < filterDim; x++)
					shiftBuffer[(y + yOffset) * horDim + (x + xOffset)] =
						skFiltersInv[i][y * filterDim + x];
			
			//Shift (with wrap) the center to the corners
			for (int y = 0; y < verDim; y++)
				for (int x = 0; x < horDim; x++)
					skFiltersFTs.back().filtersInvFTs[i][lvl][y * horDim + x] =
						shiftBuffer[((y + verDim / 2) % verDim) * horDim + (x + horDim / 2) % horDim];
			
			//Do the Fourier analysis
			assert(doFourierAnalysis(skFiltersFTs.back().filtersInvFTs[i][lvl],
									skFiltersFTs.back().filtersInvFTs[i][lvl], horDim, verDim));
		}
	delete [] shiftBuffer;
	
	return skFiltersFTs.size() - 1;
}

void SteerablePyramid::CreateScratch()
{
	//Note that scratch_2D_A_ will be memory leaked when the last instance of SteerablePyramid is deleted!!!
	//Reasonably simple to fix, but this is just a dinky little experimental program.  Tsk tsk.
	if (scratchSize_ < width_ * height_)
	{
		if (scratch_2D_A_)
			delete [] scratch_2D_A_;
		if (scratch_2D_B_)
			delete [] scratch_2D_B_;
		
		scratchSize_ = width_ * height_;
		
		scratch_2D_A_ = new Float_t[scratchSize_];
		assert(scratch_2D_A_);
		scratch_2D_B_ = new Float_t[scratchSize_];
		assert(scratch_2D_B_);
	}
}

//static
void SteerablePyramid::GenerateFilterPGMimages(const char* filenameHeader, bool plotLog)
{
	InitFilters();
	
	Float_t maxAbsTap = 0;
	for (int i = 0; i < 7; i++)
	{
		int numTaps = skFilterDims[i] * skFilterDims[i];
		for (int j = 0; j < numTaps; j++)
			if (fabs(skFilters[i][j]) > maxAbsTap)
				maxAbsTap = fabs(skFilters[i][j]);
	}
	Float_t scalar = 127.0 / maxAbsTap;
	//Float_t scalar = 127.0 / (plotLog ? log(5896.0) : 5896.0);
	
	for (int i = 0; i < 7; i++)
	{
		string filename(filenameHeader);
		
		Float_t *filter = skFilters[i];
		switch (i)
		{
			case 0: filename += "_H0"; break;
			case 1: filename += "_L0"; break;
			case 2: filename += "_L1"; break;
			case 3: filename += "_B1"; break;
			case 4: filename += "_B2"; break;
			case 5: filename += "_B3"; break;
			case 6: filename += "_B4"; break;
		}
		
		filename += ".pgm";
		
		ofstream ofs;;
		ofs.open(filename.c_str());
		assert(ofs);
		
		//Write pgm header
		ofs << "P5" << endl
			<< "#." << endl
			<< skFilterDims[i] << " " << skFilterDims[i] << endl
			<< 255 << endl;
		
		Float_t pixelF;
		int pixel;
		unsigned char pixelUC;
		int numTaps = skFilterDims[i] * skFilterDims[i];
		for (int j = 0; j < numTaps; j++)
		{
			Float_t tap = filter[j];
			
			pixelF = tap * scalar;
			assert(pixelF >= -127 && pixelF <= 127);
			pixelF += 128;
			pixel = pixelF;
			
			cout << setw(4) << pixel;
			if (j % skFilterDims[i] == skFilterDims[i] - 1)
				cout << endl;
			
			assert(pixel >= 0 && pixel <= 255);
			
			pixelUC = pixel;
			ofs.put(pixelUC);
		}
		cout << endl << endl;
		
		ofs.close();
	}
}

//static
void SteerablePyramid::GetFilters(vector<vector<Float_t> > &filters)
{
	InitFilters();
	
	filters.clear();
	filters.resize(7);
	
	for (int i = 0; i < 7; i++)
	{
		filters[i].resize(skFilterDims[i] * skFilterDims[i]);
		for (int j = 0; j < skFilterDims[i] * skFilterDims[i]; j++)
			filters[i][j] = skFilters[i][j];
	}
}

void SteerablePyramid::Destroy()
{
	if (yLookup_)
		delete [] yLookup_;
	yLookup_ = NULL;
	
	if (h0Coeffs_)
		delete [] h0Coeffs_;
	h0Coeffs_ = NULL;
	
	if (b1Coeffs_)
	{
		for (int i = 0; i < numLevels_; i++)
			delete [] b1Coeffs_[i];
		delete [] b1Coeffs_;
	}
	b1Coeffs_ = NULL;
	
	if (b2Coeffs_)
	{
		for (int i = 0; i < numLevels_; i++)
			delete [] b2Coeffs_[i];
		delete [] b2Coeffs_;
	}
	b2Coeffs_ = NULL;
	
	if (b3Coeffs_)
	{
		for (int i = 0; i < numLevels_; i++)
			delete [] b3Coeffs_[i];
		delete [] b3Coeffs_;
	}
	b3Coeffs_ = NULL;
	
	if (b4Coeffs_)
	{
		for (int i = 0; i < numLevels_; i++)
			delete [] b4Coeffs_[i];
		delete [] b4Coeffs_;
	}
	b4Coeffs_ = NULL;
	
	if (l1Coeffs_)
	{
		for (int i = 0; i < numLevels_; i++)
			delete [] l1Coeffs_[i];
		delete [] l1Coeffs_;
	}
	l1Coeffs_ = NULL;
}

void SteerablePyramid::Init(int width, int height)
{
	if (width_ == width && height_ == height)
		return;
	width_ = width;
	height_ = height;
	
	CreateScratch();
	
	Destroy();
	
	DetermineNumLevels();
	
	yLookup_ = new int[height_];
	
	int totalPixels = width_ * height_;
	h0Coeffs_ = new Float_t[totalPixels];
	
	b1Coeffs_ = new Float_t*[numLevels_];
	for (int i = 0; i < numLevels_; i++)
		b1Coeffs_[i] = new Float_t[totalPixels / (int)(pow(4, i))];
	
	b2Coeffs_ = new Float_t*[numLevels_];
	for (int i = 0; i < numLevels_; i++)
		b2Coeffs_[i] = new Float_t[totalPixels / (int)(pow(4, i))];
	
	b3Coeffs_ = new Float_t*[numLevels_];
	for (int i = 0; i < numLevels_; i++)
		b3Coeffs_[i] = new Float_t[totalPixels / (int)(pow(4, i))];
	
	b4Coeffs_ = new Float_t*[numLevels_];
	for (int i = 0; i < numLevels_; i++)
		b4Coeffs_[i] = new Float_t[totalPixels / (int)(pow(4, i))];
	
	l1Coeffs_ = new Float_t*[numLevels_];
	for (int i = 0; i < numLevels_; i++)
		l1Coeffs_[i] = new Float_t[totalPixels / (int)(pow(4, i + 1))];
}

void SteerablePyramid::DetermineNumLevels()
{
	numLevels_ = 0;
	int shorterDim = min(width_, height_);
	while (shorterDim > skFilterDims[2])	//Limit is determined by the dimensions of the L1 filter (probably 17x17)
	{
		numLevels_++;
		shorterDim /= 2;
	}
}

const int SteerablePyramid::NumLevels() const
{
	return numLevels_;
}

void SteerablePyramid::SetAllZero()
{
	for (int i = 0; i < numLevels_; i++)
	{
		int horDim = width_ / pow(2, i);
		int verDim = height_ / pow(2, i);
		int totalPixels = horDim * verDim;
		int totalBytes = totalPixels * sizeof(Float_t);
		
		if (i == 0)
			memset(h0Coeffs_, 0, totalBytes);
		memset(b1Coeffs_[i], 0, totalBytes);
		memset(b2Coeffs_[i], 0, totalBytes);
		memset(b3Coeffs_[i], 0, totalBytes);
		memset(b4Coeffs_[i], 0, totalBytes);
		memset(l1Coeffs_[i], 0, totalBytes / 4);
	}
}

#pragma mark -

void SteerablePyramid::ConvolveImage(const Float_t* image, int width, int height, Float_t kernel[289], Float_t* result, bool initValues)
{
	//Should be using the other version of this function
	assert(false);
	
	Float_t* resultPtr = result;
	if (initValues)
		memset(result, 0, width * height * sizeof(Float_t));
	
	int widthMinus8 = width - 8, widthMinus16 = width - 16;
	int heightMinus8 = height - 8;
	int y, x;
	int xx, yy;
	int ky, kx;
	
	Float_t* kernelPtr;
	const Float_t* imagePtr;
	Float_t sum;
	
	unsigned long t = getCurrentTime_ms();
	//Do the interior first, then the border
	for (y = 8; y < heightMinus8; y++)
	{
		resultPtr = &result[yLookup_[y] + 8];
		for (x = 0; x < widthMinus16; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			sum = 0;
			for (ky = -8; ky <= 8; ky++)
			{
				imagePtr = image + yLookup_[y + ky] + x;
				//for (kx = 0; kx < 17; kx++)
				//	sum += *imagePtr++ * *kernelPtr++;
				
				//Unroll the loop.  This will be done by compiler optimizations anyway,
				//but I like to leave the optimizations off to ease debugging.  Also,
				//remove the postincrements, which can be costly.
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
				sum += *imagePtr * *kernelPtr;    ++imagePtr;    ++kernelPtr;
			}
			*resultPtr++ += sum;
		}
	}
	gT2 += getCurrentTime_ms() - t;
	
	t = getCurrentTime_ms();
	//Now the border
	//Horizontal stripes top and bottom
	for (y = 0; y < 8; y++)
	{
		resultPtr = &result[y * width + 8];
		for (x = 8; x < widthMinus8; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -8; ky <= 8; ky++)
			{
				yy = y + ky;
				
				if (yy < 0)
					yy = -yy;//yy += height;
				else if (yy >= height)
					assert(false);//yy = (height - 1) - (yy - height);//yy -= height;
					
				for (kx = -8; kx <= 8; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
				}
			}
			resultPtr++;
		}
	}
				
	for (y = heightMinus8; y < height; y++)
	{
		resultPtr = &result[y * width + 8];
		for (x = 8; x < widthMinus8; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -8; ky <= 8; ky++)
			{
				yy = y + ky;
				
				if (yy < 0)
					assert(false);//yy = -yy;//yy += height;
				else if (yy >= height)
					yy = (height - 1) - (yy - height);//yy -= height;
				
				for (kx = -8; kx <= 8; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
				}
			}
			resultPtr++;
		}
	}
	
	//Vertical stripes left and right
	for (y = 8; y < heightMinus8; y++)
	{
		resultPtr = &result[y * width];
		for (x = 0; x < 8; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -8; ky <= 8; ky++)
			{
				yy = y + ky;
				for (kx = -8; kx <= 8; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					
					if (xx < 0)
						xx = -xx;//xx += width;
					else if (xx >= width)
						assert(false);//xx = (width - 1) - (xx - width);//xx -= width;
					
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
				}
			}
			resultPtr++;
		}
				
		resultPtr = &result[y * width] + widthMinus8;
		for (x = widthMinus8; x < width; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -8; ky <= 8; ky++)
			{
				yy = y + ky;
				for (kx = -8; kx <= 8; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					
					if (xx < 0)
						assert(false);//xx = -xx;//xx += width;
					else if (xx >= width)
						xx = (width - 1) - (xx - width);//xx -= width;
					
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
				}
			}
			resultPtr++;
		}
	}
	
	//Corners
	for (y = 0; y < 8; y++)
		for (x = 0; x < 8; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -8; ky <= 8; ky++)
				for (kx = -8; kx <= 8; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					yy = y + ky;
					
					if (xx < 0)
						xx = -xx;//xx += width;
					else if (xx >= width)
						assert(false);//xx = (width - 1) - (xx - width);//xx -= width;
					
					if (yy < 0)
						yy = -yy;//yy += height;
					else if (yy >= height)
						assert(false);//yy = (height - 1) - (yy - height);//yy -= height;
					
					result[y * width + x] += image[yLookup_[yy] + xx] * *kernelPtr;
					
					//Find the offset pixel location
					xx = widthMinus8 + x + kx;
					yy = y + ky;
					
					if (xx < 0)
						assert(false);//xx = -xx;//xx += width;
					else if (xx >= width)
						xx = (width - 1) - (xx - width);//xx -= width;
					
					if (yy < 0)
						yy = -yy;//yy += height;
					else if (yy >= height)
						assert(false);//yy = (height - 1) - (yy - height);//yy -= height;
					
					result[y * width + (widthMinus8 + x)] += image[yLookup_[yy] + xx] * *kernelPtr;
					
					//Find the offset pixel location
					xx = x + kx;
					yy = heightMinus8 + y + ky;
					
					if (xx < 0)
						xx = -xx;//xx += width;
					else if (xx >= width)
						assert(false);//xx = (width - 1) - (xx - width);//xx -= width;
					
					if (yy < 0)
						assert(false);//yy = -yy;//yy += height;
					else if (yy >= height)
						yy = (height - 1) - (yy - height);//yy -= height;
					
					result[(heightMinus8 + y) * width + x] += image[yLookup_[yy] + xx] * *kernelPtr;
					
					//Find the offset pixel location
					xx = widthMinus8 + x + kx;
					yy = heightMinus8 + y + ky;
					
					if (xx < 0)
						assert(false);//xx = -xx;//xx += width;
					else if (xx >= width)
						xx = (width - 1) - (xx - width);//xx -= width;
					
					if (yy < 0)
						assert(false);//yy = -yy;//yy += height;
					else if (yy >= height)
						yy = (height - 1) - (yy - height);//yy -= height;
					
					result[(heightMinus8 + y) * width + (widthMinus8 + x)] += image[yLookup_[yy] + xx] * *kernelPtr;
					
					kernelPtr++;
				}
		}
	gT3 += getCurrentTime_ms() - t;
}

void SteerablePyramid::ConvolveImage(const Float_t* image, int width, int height, Float_t kernel[], int kernelDim, Float_t* result, bool initValues)
{
	Float_t* resultPtr = result;
	if (initValues)
		memset(result, 0, width * height * sizeof(Float_t));
	
	assert(kernelDim % 2 == 1);	//Kernel dim should be odd
	int kernelHalfDim = kernelDim / 2;
	int widthMinusKernelHalfDim = width - kernelHalfDim, widthMinusKernelDimMinusOne = width - (kernelDim - 1);
	int heightMinusKernelHalfDim = height - kernelHalfDim;
	int y, x;
	int xx, yy;
	int ky, kx;
	
	Float_t* kernelPtr;
	const Float_t* imagePtr;
	Float_t sum;
	
	//Convolve the interior
	unsigned long t = getCurrentTime_ms();
	for (y = kernelHalfDim; y < heightMinusKernelHalfDim; y++)
	{
		resultPtr = &result[yLookup_[y] + kernelHalfDim];
		for (x = 0; x < widthMinusKernelDimMinusOne; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			sum = 0;
			for (ky = -kernelHalfDim; ky <= kernelHalfDim; ky++)
			{
				imagePtr = image + yLookup_[y + ky] + x;
				for (kx = 0; kx < kernelDim; kx++)
				{
					sum += *imagePtr++ * *kernelPtr++;
					
					//if (kernelDim == 17 && width == 32 && y == kernelHalfDim && x == 0)
					//	cout << setw(6) << x + kx << setw(4) << y + ky;
				}
				//if (kernelDim == 17 && width == 32 && y == kernelHalfDim && x == 0)
				//	cout << endl;
			}
			*resultPtr++ += sum;
		}
	}
	gT2 += getCurrentTime_ms() - t;
	
	//Convolve the border
	t = getCurrentTime_ms();
	
	//Horizontal stripes top and bottom
	for (y = 0; y < kernelHalfDim; y++)
	{
		resultPtr = &result[yLookup_[y] + kernelHalfDim];
		for (x = kernelHalfDim; x < widthMinusKernelHalfDim; x++, resultPtr++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -kernelHalfDim; ky <= kernelHalfDim; ky++)
			{
				yy = y + ky;
				
				if (yy < 0)
					yy = handleNegEdge(yy, height);
				else if (yy >= height)
					assert(false);
				
				if (yy == -1)
				{
					kernelPtr += kernelDim;
					continue;
				}
				
				for (kx = -kernelHalfDim; kx <= kernelHalfDim; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					assert(xx >= 0 && xx < width);
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
					
					//if (kernelDim == 17 && width == 32 && y == 0 && x == kernelHalfDim)
					//	cout << setw(6) << xx << setw(4) << yy;
				}
				//if (kernelDim == 17 && width == 32 && y == 0 && x == kernelHalfDim)
				//	cout << endl;
			}
		}
	}
				
	for (y = heightMinusKernelHalfDim; y < height; y++)
	{
		resultPtr = &result[yLookup_[y] + kernelHalfDim];
		for (x = kernelHalfDim; x < widthMinusKernelHalfDim; x++, resultPtr++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -kernelHalfDim; ky <= kernelHalfDim; ky++)
			{
				yy = y + ky;
				
				if (yy < 0)
					assert(false);
				else if (yy >= height)
					yy = handlePosEdge(yy, height);
				
				if (yy == -1)
				{
					kernelPtr += kernelDim;
					continue;
				}
				
				for (kx = -kernelHalfDim; kx <= kernelHalfDim; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					assert(xx >= 0 && xx < width);
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
				}
			}
		}
	}
	
	//Vertical stripes left and right
	for (y = kernelHalfDim; y < heightMinusKernelHalfDim; y++)
	{
		resultPtr = &result[y * width];
		for (x = 0; x < kernelHalfDim; x++, resultPtr++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -kernelHalfDim; ky <= kernelHalfDim; ky++)
			{
				yy = y + ky;
				for (kx = -kernelHalfDim; kx <= kernelHalfDim; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					
					if (xx < 0)
						xx = handleNegEdge(xx, width);
					else if (xx >= width)
						assert(false);
					
					if (xx == -1)
					{
						kernelPtr++;
						continue;
					}
					
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
				}
			}
		}
				
		resultPtr = &result[y * width] + widthMinusKernelHalfDim;
		for (x = widthMinusKernelHalfDim; x < width; x++, resultPtr++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -kernelHalfDim; ky <= kernelHalfDim; ky++)
			{
				yy = y + ky;
				for (kx = -kernelHalfDim; kx <= kernelHalfDim; kx++)
				{
					//Find the offset pixel location
					xx = x + kx;
					
					if (xx < 0)
						assert(false);
					else if (xx >= width)
						xx = handlePosEdge(xx, width);
					
					if (xx == -1)
					{
						kernelPtr++;
						continue;
					}
					
					*resultPtr += image[yLookup_[yy] + xx] * *kernelPtr++;
				}
			}
		}
	}
	
	//Corners
	for (y = 0; y < kernelHalfDim; y++)
		for (x = 0; x < kernelHalfDim; x++)
		{
			//For each component in the kernel
			kernelPtr = kernel;
			for (ky = -kernelHalfDim; ky <= kernelHalfDim; ky++)
				for (kx = -kernelHalfDim; kx <= kernelHalfDim; kx++, kernelPtr++)
				{
					//Find the offset pixel location
					xx = x + kx;
					yy = y + ky;
					
					if (xx < 0)
						xx = handleNegEdge(xx, width);
					else if (xx >= width)
						assert(false);
					
					if (yy < 0)
						yy = handleNegEdge(yy, height);
					else if (yy >= height)
						assert(false);
					
					if (xx != -1 && yy != -1)
						result[y * width + x] += image[yLookup_[yy] + xx] * *kernelPtr;
					
					//Find the offset pixel location
					xx = widthMinusKernelHalfDim + x + kx;
					yy = y + ky;
					
					if (xx < 0)
						assert(false);
					else if (xx >= width)
						xx = handlePosEdge(xx, width);
					
					if (yy < 0)
						yy = handleNegEdge(yy, height);
					else if (yy >= height)
						assert(false);
					
					if (xx != -1 && yy != -1)
						result[y * width + (widthMinusKernelHalfDim + x)] += image[yLookup_[yy] + xx] * *kernelPtr;
					
					//Find the offset pixel location
					xx = x + kx;
					yy = heightMinusKernelHalfDim + y + ky;
					
					if (xx < 0)
						xx = handleNegEdge(xx, width);
					else if (xx >= width)
						assert(false);
					
					if (yy < 0)
						assert(false);
					else if (yy >= height)
						yy = handlePosEdge(yy, height);
					
					if (xx != -1 && yy != -1)
						result[(heightMinusKernelHalfDim + y) * width + x] += image[yLookup_[yy] + xx] * *kernelPtr;
					
					//Find the offset pixel location
					xx = widthMinusKernelHalfDim + x + kx;
					yy = heightMinusKernelHalfDim + y + ky;
					
					if (xx < 0)
						assert(false);
					else if (xx >= width)
						xx = handlePosEdge(xx, width);
					
					if (yy < 0)
						assert(false);
					else if (yy >= height)
						yy = handlePosEdge(yy, height);
					
					if (xx != -1 && yy != -1)
						result[(heightMinusKernelHalfDim + y) * width + (widthMinusKernelHalfDim + x)] += image[yLookup_[yy] + xx] * *kernelPtr;
				}
		}
	gT3 += getCurrentTime_ms() - t;
}

#pragma mark -

void SteerablePyramid::DoOneStepWaveletAnalysis(Float_t* pixelValues, int level)
{
	int horDim = width_ / pow(2, level);
	int verDim = height_ / pow(2, level);
	
	for (int i = 0; i < verDim; i++)
		yLookup_[i] = i * horDim;
	
	Float_t *l0OriginalCoeffs = NULL;
	
	const Float_t *pixelPtr;
	Float_t *coeffPtr;
	
	if (level == 0)
	{
		ConvolveImage(pixelValues, horDim, verDim, skFilters[0], skFilterDims[0], h0Coeffs_, true);
		
		l0OriginalCoeffs = new Float_t[horDim * verDim];
		assert(l0OriginalCoeffs);
		ConvolveImage(pixelValues, horDim, verDim, skFilters[1], skFilterDims[1], l0OriginalCoeffs, true);
	}
	else l0OriginalCoeffs = l1Coeffs_[level - 1];
	
	ConvolveImage((const Float_t*)l0OriginalCoeffs, horDim, verDim, skFilters[3], skFilterDims[3], b1Coeffs_[level], true);
	ConvolveImage((const Float_t*)l0OriginalCoeffs, horDim, verDim, skFilters[4], skFilterDims[4], b2Coeffs_[level], true);
	ConvolveImage((const Float_t*)l0OriginalCoeffs, horDim, verDim, skFilters[5], skFilterDims[5], b3Coeffs_[level], true);
	ConvolveImage((const Float_t*)l0OriginalCoeffs, horDim, verDim, skFilters[6], skFilterDims[6], b4Coeffs_[level], true);
	ConvolveImage((const Float_t*)l0OriginalCoeffs, horDim, verDim, skFilters[2], skFilterDims[2], pixelValues, true);
	
	coeffPtr = l1Coeffs_[level];
	for (int y = 0; y < verDim / 2; y++)
	{
		pixelPtr = &pixelValues[(y * 2) * horDim];
		for (int x = 0; x < horDim / 2; x++)
		{
			*coeffPtr++ = *pixelPtr;
			pixelPtr += 2;
		}
	}
	
	if (level == 0)
		delete [] l0OriginalCoeffs;
}

void SteerablePyramid::DoOneStepWaveletAnalysisUsingFT(Float_t* pixelValues, int level)
{
	int horDim = width_ / pow(2, level);
	int verDim = height_ / pow(2, level);
	
	int filterIndex = getFilterFTIndex();
	
	for (int i = 0; i < verDim; i++)
		yLookup_[i] = i * horDim;
	
	Float_t *l0OriginalCoeffs = NULL;
	
	if (level == 0)
	{
		assert(doFourierAnalysis(pixelValues, scratch_2D_B_, horDim, verDim));
		assert(multiplyFourierTransformsAndInvert(scratch_2D_B_, skFiltersFTs[filterIndex].filtersFTs[0][level], horDim, verDim, h0Coeffs_));
		
		l0OriginalCoeffs = new Float_t[horDim * verDim];
		assert(l0OriginalCoeffs);
		assert(multiplyFourierTransformsAndInvert(scratch_2D_B_, skFiltersFTs[filterIndex].filtersFTs[1][level], horDim, verDim, l0OriginalCoeffs));
	}
	else l0OriginalCoeffs = l1Coeffs_[level - 1];
	
	assert(doFourierAnalysis(l0OriginalCoeffs, scratch_2D_B_, horDim, verDim));
	
	assert(multiplyFourierTransformsAndInvert(scratch_2D_B_, skFiltersFTs[filterIndex].filtersFTs[3][level], horDim, verDim, b1Coeffs_[level]));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_B_, skFiltersFTs[filterIndex].filtersFTs[4][level], horDim, verDim, b2Coeffs_[level]));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_B_, skFiltersFTs[filterIndex].filtersFTs[5][level], horDim, verDim, b3Coeffs_[level]));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_B_, skFiltersFTs[filterIndex].filtersFTs[6][level], horDim, verDim, b4Coeffs_[level]));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_B_, skFiltersFTs[filterIndex].filtersFTs[2][level], horDim, verDim, pixelValues));
	
	const Float_t *pixelPtr;
	Float_t *coeffPtr = l1Coeffs_[level];
	for (int y = 0; y < verDim / 2; y++)
	{
		pixelPtr = &pixelValues[(y * 2) * horDim];
		for (int x = 0; x < horDim / 2; x++)
		{
			*coeffPtr++ = *pixelPtr;
			pixelPtr += 2;
		}
	}
	
	if (level == 0)
		delete [] l0OriginalCoeffs;
}

void SteerablePyramid::DoFullWaveletAnalysis(const Float_t* pixelValues, int width, int height, bool useFTconvolution)
{
	gT1 = getCurrentTime_ms();
	gT2 = gT3 = gT4 = gT5 = 0;
	
	Init(width, height);
	
	useFTconvolution_ = useFTconvolution;
	
	if (!useFTconvolution_)
		cout << endl << "SteerablePyramid Fourier Transform convolution disabled" << endl;
	
	int totalPixels = width_ * height_;
	/*
	Float_t minPixel = 999999, maxPixel = -999999;
	for (int i = 0; i < totalPixels; i++)
	{
		if (pixelValues[i] < minPixel)
			minPixel = pixelValues[i];
		if (pixelValues[i] > maxPixel)
			maxPixel = pixelValues[i];
	}
	cout << "Analysis Min/Max Pixels: " << minPixel << "    " << maxPixel << endl;
	*/
	memcpy(scratch_2D_A_, pixelValues, width_ * height_ * sizeof(Float_t));
	
	//cout << "  Performing analysis with " << numLevels_ << " levels" <<  endl;
	for (int i = 0; i < numLevels_; i++)
	{
		//cout << "    Analyzing level " << i << "..." << endl;
		if (!useFTconvolution_)
			DoOneStepWaveletAnalysis(scratch_2D_A_, i);
		else DoOneStepWaveletAnalysisUsingFT(scratch_2D_A_, i);
	}
	
	//cout << "  Done" << endl;
	
	gT1 = getCurrentTime_ms() - gT1;
	Float_t t1Secs = (Float_t)gT1 / CURRENT_TIME_PER_SEC_F;
	Float_t t2Secs = (Float_t)gT2 / CURRENT_TIME_PER_SEC_F;
	Float_t t3Secs = (Float_t)gT3 / CURRENT_TIME_PER_SEC_F;
	Float_t t4Secs = (Float_t)gT4 / CURRENT_TIME_PER_SEC_F;
	Float_t t5Secs = (Float_t)gT5 / CURRENT_TIME_PER_SEC_F;
	Float_t t2Fraction = (Float_t)gT2 / (Float_t)gT1;
	Float_t t3Fraction = (Float_t)gT3 / (Float_t)gT1;
	Float_t t4Fraction = (Float_t)gT4 / (Float_t)gT1;
	Float_t t5Fraction = (Float_t)gT5 / (Float_t)gT1;
	/*
	cout << "  Times:" << endl;
	cout << "    Total: " << t1Secs << endl;
	cout << "    T2: ";
	cout << setw(10) << t2Secs << "   Ratio: ";
	cout << setw(10) << t2Fraction << endl;
	cout << "    T3: ";
	cout << setw(10) << t3Secs << "   Ratio: ";
	cout << setw(10) << t3Fraction << endl;
	*/
	//cout << "  At2 " << t2Fraction << "  At3 " << t3Fraction << "  ";
}

void SteerablePyramid::DoOneStepWaveletSynthesis(int level, Float_t* pixelValues, Float_t *upsampledThumbnail)
{
	int horDim = width_ / pow(2, level);
	int verDim = height_ / pow(2, level);
	int totalPixels = horDim * verDim;
	
	for (int i = 0; i < verDim; i++)
		yLookup_[i] = i * horDim;
	
	Float_t *finalImage = NULL;
	Float_t *pixelPtr = NULL;
	Float_t *coeffPtr = NULL;
	
	if (level == 0)
		finalImage = new Float_t[totalPixels];
	
	//Is this necessary or did new clear out the allocated memory?
	memset(upsampledThumbnail, 0, totalPixels * sizeof(Float_t));
	
	//Upsample the growing thumbnail from the base thumbnail (lowest level) or the previous iteration's thumbnail
	
	if (level == numLevels_ - 1)	//Lowest level, upsample actual stored thumbnail
		pixelPtr = l1Coeffs_[level];
	else	//Intermediate level, upsample growing thumbnail from previous iteration
		pixelPtr = pixelValues;
	
	for (int y = 0; y < verDim / 2; y++)
	{
		coeffPtr = &upsampledThumbnail[(y * 2) * horDim];
		for (int x = 0; x < horDim / 2; x++)
		{
			*coeffPtr = *pixelPtr++ * 4;
			coeffPtr += 2;
		}
	}
	
	ConvolveImage(upsampledThumbnail, horDim, verDim, skFiltersInv[2], skFilterDims[2], pixelValues, true);
	ConvolveImage(b1Coeffs_[level], horDim, verDim, skFiltersInv[3], skFilterDims[3], pixelValues, false);
	ConvolveImage(b2Coeffs_[level], horDim, verDim, skFiltersInv[4], skFilterDims[4], pixelValues, false);
	ConvolveImage(b3Coeffs_[level], horDim, verDim, skFiltersInv[5], skFilterDims[5], pixelValues, false);
	ConvolveImage(b4Coeffs_[level], horDim, verDim, skFiltersInv[6], skFilterDims[6], pixelValues, false);
	
	if (level == 0)
	{
		ConvolveImage(pixelValues, horDim, verDim, skFiltersInv[1], skFilterDims[1], finalImage, true);
		ConvolveImage(h0Coeffs_, horDim, verDim, skFiltersInv[0], skFilterDims[0], finalImage, false);
	}
	
	if (level == 0)
		memcpy(pixelValues, finalImage, totalPixels * sizeof(Float_t));
	
	if (level == 0)
		delete [] finalImage;
}

void SteerablePyramid::DoOneStepWaveletSynthesisUsingFT(int level, Float_t* pixelValues, Float_t *upsampledThumbnail)
{
	int horDim = width_ / pow(2, level);
	int verDim = height_ / pow(2, level);
	int totalPixels = horDim * verDim;
	
	int filterIndex = getFilterFTIndex();
	
	for (int i = 0; i < verDim; i++)
		yLookup_[i] = i * horDim;
	
	Float_t *finalImage = NULL;
	Float_t *pixelPtr = NULL;
	Float_t *coeffPtr = NULL;
	
	if (level == 0)
		finalImage = new Float_t[totalPixels];
	
	//Is this necessary or did new clear out the allocated memory?
	memset(upsampledThumbnail, 0, totalPixels * sizeof(Float_t));
	
	//Upsample the growing thumbnail from the base thumbnail (lowest level) or the previous iteration's thumbnail
	
	if (level == numLevels_ - 1)	//Lowest level, upsample actual stored thumbnail
		pixelPtr = l1Coeffs_[level];
	else	//Intermediate level, upsample growing thumbnail from previous iteration
		pixelPtr = pixelValues;
	
	for (int y = 0; y < verDim / 2; y++)
	{
		coeffPtr = &upsampledThumbnail[(y * 2) * horDim];
		for (int x = 0; x < horDim / 2; x++)
		{
			*coeffPtr = *pixelPtr++ * 4;
			coeffPtr += 2;
		}
	}
	
	assert(doFourierAnalysis(upsampledThumbnail, scratch_2D_A_, horDim, verDim));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_A_, skFiltersFTs[filterIndex].filtersInvFTs[2][level], horDim, verDim, pixelValues));
	
	assert(doFourierAnalysis(b1Coeffs_[level], scratch_2D_A_, horDim, verDim));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_A_, skFiltersFTs[filterIndex].filtersInvFTs[3][level], horDim, verDim, scratch_2D_B_));
	for (int i = 0; i < totalPixels; i++)
		pixelValues[i] += scratch_2D_B_[i];
	
	assert(doFourierAnalysis(b2Coeffs_[level], scratch_2D_A_, horDim, verDim));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_A_, skFiltersFTs[filterIndex].filtersInvFTs[4][level], horDim, verDim, scratch_2D_B_));
	for (int i = 0; i < totalPixels; i++)
		pixelValues[i] += scratch_2D_B_[i];
	
	assert(doFourierAnalysis(b3Coeffs_[level], scratch_2D_A_, horDim, verDim));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_A_, skFiltersFTs[filterIndex].filtersInvFTs[5][level], horDim, verDim, scratch_2D_B_));
	for (int i = 0; i < totalPixels; i++)
		pixelValues[i] += scratch_2D_B_[i];
	
	assert(doFourierAnalysis(b4Coeffs_[level], scratch_2D_A_, horDim, verDim));
	assert(multiplyFourierTransformsAndInvert(scratch_2D_A_, skFiltersFTs[filterIndex].filtersInvFTs[6][level], horDim, verDim, scratch_2D_B_));
	for (int i = 0; i < totalPixels; i++)
		pixelValues[i] += scratch_2D_B_[i];
	
	if (level == 0)
	{
		assert(doFourierAnalysis(pixelValues, scratch_2D_A_, horDim, verDim));
		assert(multiplyFourierTransformsAndInvert(scratch_2D_A_, skFiltersFTs[filterIndex].filtersInvFTs[1][level], horDim, verDim, finalImage));
		
		assert(doFourierAnalysis(h0Coeffs_, scratch_2D_A_, horDim, verDim));
		assert(multiplyFourierTransformsAndInvert(scratch_2D_A_, skFiltersFTs[filterIndex].filtersInvFTs[0][level], horDim, verDim, scratch_2D_B_));
		for (int i = 0; i < totalPixels; i++)
			finalImage[i] += scratch_2D_B_[i];
	}
	
	if (level == 0)
		memcpy(pixelValues, finalImage, totalPixels * sizeof(Float_t));
	
	if (level == 0)
		delete [] finalImage;
}

void SteerablePyramid::DoFullWaveletSynthesis(Float_t* pixelValues)
{
	assert(width_ > 0 && height_ > 0);
	
	gT1 = getCurrentTime_ms();
	gT2 = gT3 = gT4 = gT5 = 0;
	
	int totalPixels = width_ * height_;
	
	//Allocate some scratch storage outside the levels loop so it doesn't have to be allocated each time through the loop
	gT4 = getCurrentTime_ms();
	Float_t *upsampledThumbnail = new Float_t[totalPixels];
	gT4 = getCurrentTime_ms() - gT4;
	
	//It shouldn't be necessary to erase the output image here.
	//memset(pixelValues, 0, totalPixels * sizeof(Float_t));
	
	//cout << "  Performing synthesis with " << numLevels_ << " levels" <<  endl;
	for (int i = numLevels_ - 1; i >= 0; i--)
	{
		//cout << "    Synthesizing level " << i << "..." << endl;
		if (!useFTconvolution_)
			DoOneStepWaveletSynthesis(i, pixelValues, upsampledThumbnail);
		else DoOneStepWaveletSynthesisUsingFT(i, pixelValues, upsampledThumbnail);
	}
	/*
	Float_t minPixel = 999999, maxPixel = -999999;
	for (int i = 0; i < totalPixels; i++)
	{
		if (pixelValues[i] < minPixel)
			minPixel = pixelValues[i];
		if (pixelValues[i] > maxPixel)
			maxPixel = pixelValues[i];
	}
	cout << "Synthesis Min/Max Pixels: " << minPixel << "    " << maxPixel << endl;
	*/
	//cout << "  Clipping image" << endl;
	for (int i = 0; i < totalPixels; i++)
		if (pixelValues[i] < 0)
			pixelValues[i] = 0;
		else if (pixelValues[i] > 1.0)
			pixelValues[i] = 1.0;
	/*
	maxPixel = -999999;
	for (int i = 0; i < totalPixels; i++)
		if (pixelValues[i] > maxPixel)
			maxPixel = pixelValues[i];
	string filename = "aaaUpsampleFilteredFinalClipped.pgm";
	GeneratePGMimage(pixelValues, width_, height_, 255.0 / maxPixel, 0, filename);
	*/
	//cout << "  Cleaning up" << endl;
	gT5 = getCurrentTime_ms();
	delete [] upsampledThumbnail;
	gT5 = getCurrentTime_ms() - gT5;
	
	//cout << "  Done" << endl;
	
	gT1 = getCurrentTime_ms() - gT1;
	Float_t t1Secs = (Float_t)gT1 / CURRENT_TIME_PER_SEC_F;
	Float_t t2Secs = (Float_t)gT2 / CURRENT_TIME_PER_SEC_F;
	Float_t t3Secs = (Float_t)gT3 / CURRENT_TIME_PER_SEC_F;
	Float_t t4Secs = (Float_t)gT4 / CURRENT_TIME_PER_SEC_F;
	Float_t t5Secs = (Float_t)gT5 / CURRENT_TIME_PER_SEC_F;
	Float_t t2Fraction = (Float_t)gT2 / (Float_t)gT1;
	Float_t t3Fraction = (Float_t)gT3 / (Float_t)gT1;
	Float_t t4Fraction = (Float_t)gT4 / (Float_t)gT1;
	Float_t t5Fraction = (Float_t)gT5 / (Float_t)gT1;
	/*
	cout << "  Times:" << endl;
	cout << "    Total: " << t1Secs << endl;
	cout << "    T2: ";
	cout << setw(10) << t2Secs << "   Ratio: ";
	cout << setw(10) << t2Fraction << endl;
	cout << "    T3: ";
	cout << setw(10) << t3Secs << "   Ratio: ";
	cout << setw(10) << t3Fraction << endl;
	cout << "    T4: ";
	cout << setw(10) << t4Secs << "   Ratio: ";
	cout << setw(10) << t4Fraction << endl;
	cout << "    T5: ";
	cout << setw(10) << t5Secs << "   Ratio: ";
	cout << setw(10) << t5Fraction << endl;
	*/
	//cout << "  St2 " << t2Fraction << "  St3 " << t3Fraction << "  ";
}

#pragma mark -

//returned vector will be (numLevels_+1) long.
//The first element will contain the h0 mean.
//The remaining elements will contain b1-4 means (so B1-4 level 0 mean is in the second position).
void SteerablePyramid::CalcStats(vector<Float_t> &means, vector<Float_t> &stdDevs) const
{
	assert(width_ > 0 && height_ > 0);
	
	//Init the vectors
	means.resize(numLevels_ + 1);
	stdDevs.resize(numLevels_ + 1);
	for (int i = 0; i <= numLevels_; i++)
		means[i] = stdDevs[i] = 0;
	
	//Get the H0 means
	int totalPixels = width_ * height_;
	for (int i = 0; i < totalPixels; i++)
		means[0] += fabs(h0Coeffs_[i]);
	means[0] /= totalPixels;
	
	//Get the B means
	for (int i = 0; i < numLevels_; i++)
	{
		totalPixels = (width_ * height_) / (int)(pow(4, i));
		for (int j = 0; j < totalPixels; j++)
		{
			means[i + 1] += fabs(b1Coeffs_[i][j]);
			means[i + 1] += fabs(b2Coeffs_[i][j]);
			means[i + 1] += fabs(b3Coeffs_[i][j]);
			means[i + 1] += fabs(b4Coeffs_[i][j]);
		}
		means[i + 1] /= (totalPixels * 4);
	}
	
	//Get the H0 std devs
	for (int i = 0; i < totalPixels; i++)
		stdDevs[0] += fabs(h0Coeffs_[i] - means[0]);
	stdDevs[0] /= totalPixels;
	
	//Get the B std devs
	for (int i = 0; i < numLevels_; i++)
	{
		totalPixels = (width_ * height_) / (int)(pow(4, i));
		for (int j = 0; j < totalPixels; j++)
		{
			stdDevs[i + 1] += fabs(b1Coeffs_[i][j] - means[i]);
			stdDevs[i + 1] += fabs(b2Coeffs_[i][j] - means[i]);
			stdDevs[i + 1] += fabs(b3Coeffs_[i][j] - means[i]);
			stdDevs[i + 1] += fabs(b4Coeffs_[i][j] - means[i]);
		}
		stdDevs[i + 1] /= (totalPixels * 4);
	}
}

void SteerablePyramid::SoftThreshold(vector<Float_t> thresholds)
{
	assert(width_ > 0 && height_ > 0);
	
	assert(thresholds.size() == numLevels_ + 1);
	for (int i = 0; i < thresholds.size(); i++)
		assert(thresholds[0] >= 0);
	
	int totalPixels = width_ * height_;
	if (thresholds[0] > 0)
		for (int i = 0; i < totalPixels; i++)
		{
			if (fabs(h0Coeffs_[i]) < thresholds[0])
				h0Coeffs_[i] = 0;
			else if (h0Coeffs_[i] > 0)
				h0Coeffs_[i] -= thresholds[0];
			else if (h0Coeffs_[i] < 0)
				h0Coeffs_[i] += thresholds[0];
		}
	
	for (int i = 0; i < numLevels_; i++)
	{
		Float_t threshold = thresholds[i + 1];
		
		if (threshold == 0)
			continue;
		
		int totalPixels = (width_ * height_) / (int)(pow(4, i));
		for (int j = 0; j < totalPixels; j++)
		{
			if (fabs(b1Coeffs_[i][j]) < threshold)
				b1Coeffs_[i][j] = 0;
			else if (b1Coeffs_[i][j] > 0)
				b1Coeffs_[i][j] -= threshold;
			else if (b1Coeffs_[i][j] < 0)
				b1Coeffs_[i][j] += threshold;
				
			if (fabs(b2Coeffs_[i][j]) < threshold)
				b2Coeffs_[i][j] = 0;
			else if (b2Coeffs_[i][j] > 0)
				b2Coeffs_[i][j] -= threshold;
			else if (b2Coeffs_[i][j] < 0)
				b2Coeffs_[i][j] += threshold;
				
			if (fabs(b3Coeffs_[i][j]) < threshold)
				b3Coeffs_[i][j] = 0;
			else if (b3Coeffs_[i][j] > 0)
				b3Coeffs_[i][j] -= threshold;
			else if (b3Coeffs_[i][j] < 0)
				b3Coeffs_[i][j] += threshold;
				
			if (fabs(b4Coeffs_[i][j]) < threshold)
				b4Coeffs_[i][j] = 0;
			else if (b4Coeffs_[i][j] > 0)
				b4Coeffs_[i][j] -= threshold;
			else if (b4Coeffs_[i][j] < 0)
				b4Coeffs_[i][j] += threshold;
		}
	}
}

void SteerablePyramid::findCorrelations(const SteerablePyramid &sp, vector<Float_t> &lvlCorrs) const
{
	assert(width_ == sp.width_ && height_ == sp.height_ && numLevels_ == sp.numLevels_);
	
	lvlCorrs.resize(numLevels_);
	for (int level = 0; level < numLevels_; level++)
	{
		int totalPixels = (width_ * height_) / (int)(pow(4, level));
		int totalCoeffs = totalPixels * ((level == 0) ? 5 : 4);
		
		lvlCorrs[level] = 0;
		
		//This algorithm is from http://en.wikipedia.org/wiki/Correlation
		Float_t sum_sq_x = 0;
		Float_t sum_sq_y = 0;
		Float_t sum_coproduct = 0;
		Float_t mean_x = b1Coeffs_[level][0];
		Float_t mean_y = sp.b1Coeffs_[level][0];
		Float_t sweep, delta_x, delta_y;
		for (int i = 1; i < totalCoeffs; i++)
		{
			sweep = (i - 1.0) / i;
			
			if (i / totalPixels == 0)
			{
				delta_x = b1Coeffs_[level][i % totalPixels] - mean_x;
				delta_y = sp.b1Coeffs_[level][i % totalPixels] - mean_y;
			}
			else if (i / totalPixels == 1)
			{
				delta_x = b2Coeffs_[level][i % totalPixels] - mean_x;
				delta_y = sp.b2Coeffs_[level][i % totalPixels] - mean_y;
			}
			else if (i / totalPixels == 2)
			{
				delta_x = b3Coeffs_[level][i % totalPixels] - mean_x;
				delta_y = sp.b3Coeffs_[level][i % totalPixels] - mean_y;
			}
			else if (i / totalPixels == 3)
			{
				delta_x = b4Coeffs_[level][i % totalPixels] - mean_x;
				delta_y = sp.b4Coeffs_[level][i % totalPixels] - mean_y;
			}
			else if (i / totalPixels == 4)
			{
				delta_x = h0Coeffs_[i % totalPixels] - mean_x;
				delta_y = sp.h0Coeffs_[i % totalPixels] - mean_y;
			}
			else assert(false);
			
			sum_sq_x += delta_x * delta_x * sweep;
			sum_sq_y += delta_y * delta_y * sweep;
			sum_coproduct += delta_x * delta_y * sweep;
			mean_x += delta_x / i;
			mean_y += delta_y / i;
		}
		Float_t pop_sd_x = sqrt(sum_sq_x / totalPixels);
		Float_t pop_sd_y = sqrt(sum_sq_y / totalPixels);
		Float_t cov_x_y = sum_coproduct / totalPixels;
		Float_t denom = pop_sd_x * pop_sd_y;
		Float_t correlation = 1;
		if (denom != 0)
			correlation = cov_x_y / denom;
		else assert(cov_x_y == 0);	//So long as both numerator and denominator are 0, call the correlation 0/0, or 1...sorta kinda.
		
		if (isnan(correlation) || isinf(correlation) || correlation < -1 || correlation > 1)
			cout << "Correlation: " << correlation << endl;
		assert(!isnan(correlation));
		assert(!isinf(correlation));
		assert(correlation >= -1 && correlation <= 1);
		
		lvlCorrs[level] = correlation;
	}
}

void SteerablePyramid::CopyCoefficients(const SteerablePyramid& sp, vector<Float_t> thresholds)
{
	assert(width_ == sp.width_ && height_ == sp.height_ && numLevels_ == sp.numLevels_);
	assert(thresholds.size() == numLevels_);
	
	for (int level = 0; level < numLevels_; level++)
	{
		int totalPixels = (width_ * height_) / (int)(pow(4, level));
		Float_t threshold = thresholds[level];
		
		//cout << "Level/Threshold: " << level << "    " << threshold << endl;
		
		//A threshold of -1 is a flag to not copy this level
		if (threshold == -1)
			continue;
		
		//A lower threshold is more likely to cause a change and causes a more drastic change.
		//A threshold of zero changes all coefficients and sets them precisely to sp's coefficients.
		
		//Use threshold as a diff between the two coefficients
		for (int i = 0; i < totalPixels; i++)
		{
			if (level == 0)
			{
				if (h0Coeffs_[i] < sp.h0Coeffs_[i] - threshold)
					h0Coeffs_[i] = sp.h0Coeffs_[i] - threshold;
				else if (h0Coeffs_[i] > sp.h0Coeffs_[i] + threshold)
					h0Coeffs_[i] = sp.h0Coeffs_[i] + threshold;
			}
			
			if (b1Coeffs_[level][i] < sp.b1Coeffs_[level][i] - threshold)
				b1Coeffs_[level][i] = sp.b1Coeffs_[level][i] - threshold;
			else if (b1Coeffs_[level][i] > sp.b1Coeffs_[level][i] + threshold)
				b1Coeffs_[level][i] = sp.b1Coeffs_[level][i] + threshold;
			
			if (b2Coeffs_[level][i] < sp.b2Coeffs_[level][i] - threshold)
				b2Coeffs_[level][i] = sp.b2Coeffs_[level][i] - threshold;
			else if (b2Coeffs_[level][i] > sp.b2Coeffs_[level][i] + threshold)
				b2Coeffs_[level][i] = sp.b2Coeffs_[level][i] + threshold;
			
			if (b3Coeffs_[level][i] < sp.b3Coeffs_[level][i] - threshold)
				b3Coeffs_[level][i] = sp.b3Coeffs_[level][i] - threshold;
			else if (b3Coeffs_[level][i] > sp.b3Coeffs_[level][i] + threshold)
				b3Coeffs_[level][i] = sp.b3Coeffs_[level][i] + threshold;
			
			if (b4Coeffs_[level][i] < sp.b4Coeffs_[level][i] - threshold)
				b4Coeffs_[level][i] = sp.b4Coeffs_[level][i] - threshold;
			else if (b4Coeffs_[level][i] > sp.b4Coeffs_[level][i] + threshold)
				b4Coeffs_[level][i] = sp.b4Coeffs_[level][i] + threshold;
		}
		
		/*
		//Use threshold as a ratio between the two coefficients
		for (int i = 0; i < totalPixels; i++)
		{
			if (level == 0)
				if (fabs(h0Coeffs_[i]) < fabs(sp.h0Coeffs_[i]) * threshold)
					h0Coeffs_[i] = sp.h0Coeffs_[i];
			if (fabs(b1Coeffs_[level][i]) < fabs(sp.b1Coeffs_[level][i]) * threshold)
				b1Coeffs_[level][i] = sp.b1Coeffs_[level][i];
			if (fabs(b2Coeffs_[level][i]) < fabs(sp.b2Coeffs_[level][i]) * threshold)
				b2Coeffs_[level][i] = sp.b2Coeffs_[level][i];
			if (fabs(b3Coeffs_[level][i]) < fabs(sp.b3Coeffs_[level][i]) * threshold)
				b3Coeffs_[level][i] = sp.b3Coeffs_[level][i];
			if (fabs(b4Coeffs_[level][i]) < fabs(sp.b4Coeffs_[level][i]) * threshold)
				b4Coeffs_[level][i] = sp.b4Coeffs_[level][i];
		}
		*/
	}
}

void SteerablePyramid::FindStrongestCoefficients(vector<Float_t> &maxCoeffs, vector<Float_t> &maxThumbnailPixels)
{
	maxCoeffs.clear();
	maxThumbnailPixels.clear();
	
	//Find the highest magnitude coefficients
	for (int i = 0; i < 2/*numLevels_*/; i++)
	{
		Float_t maxCoeff = 0;
		int horDim = width_ / pow(2, i);
		int verDim = height_ / pow(2, i);
		int totalPixels = horDim * verDim;
		for (int j = 0; j < totalPixels; j++)
		{
			//if (i == 0 && fabs(h0Coeffs_[j]) > maxCoeff)
			//	maxCoeff = fabs(h0Coeffs_[j]);
			if (fabs(b1Coeffs_[i][j]) > maxCoeff)
				maxCoeff = fabs(b1Coeffs_[i][j]);
			if (fabs(b2Coeffs_[i][j]) > maxCoeff)
				maxCoeff = fabs(b2Coeffs_[i][j]);
			if (fabs(b3Coeffs_[i][j]) > maxCoeff)
				maxCoeff = fabs(b3Coeffs_[i][j]);
			if (fabs(b4Coeffs_[i][j]) > maxCoeff)
				maxCoeff = fabs(b4Coeffs_[i][j]);
		}
		maxCoeffs.push_back(maxCoeff);
	}
	
	//Find the highest thumbnail values
	for (int i = 0; i < numLevels_; i++)
	{
		Float_t maxThumbnailPixel = 0;
		int horDim = width_ / pow(2, i + 1);
		int verDim = height_ / pow(2, i + 1);
		int totalPixels = horDim * verDim;
		for (int j = 0; j < totalPixels; j++)
			if (fabs(l1Coeffs_[i][j]) > maxThumbnailPixel)
				maxThumbnailPixel = fabs(l1Coeffs_[i][j]);
		maxThumbnailPixels.push_back(maxThumbnailPixel);
	}
}

void SteerablePyramid::GeneratePGMimages(const char* filenameHeader)
{
	vector<Float_t> maxCoeffs, maxThumbnailPixels;
	FindStrongestCoefficients(maxCoeffs, maxThumbnailPixels);
	Float_t maxCoeff = 0;
	for (int i = 0; i < maxCoeffs.size(); i++)
		if (maxCoeffs[i] > maxCoeff)
			maxCoeff = maxCoeffs[i];
	Float_t maxThumbnailPixel = 0;
	for (int i = 0; i < maxThumbnailPixels.size(); i++)
		if (maxThumbnailPixels[i] > maxThumbnailPixel)
			maxThumbnailPixel = maxThumbnailPixels[i];
	
	//maxCoeff = .5;//еее
	Float_t coeffScalar = 127.0 / maxCoeff, thumbnailScalar = 255.0 / maxThumbnailPixel;
	
	cout << "maxCoeff/maxThumbnailPixel: " << maxCoeff << "    " << maxThumbnailPixel << endl;
	cout << "coeffScalar/thumnailScalar: " << coeffScalar << "    " << thumbnailScalar << endl;
	
	//Write each set of coefficients to a pgm file
	for (int i = 0; i < numLevels_; i++)
	{
		int horDim = width_ / pow(2, i);
		int verDim = height_ / pow(2, i);
		
		stringstream levelSS;
		levelSS << i;
		string filename;
		
		if (i == 0)
		{
			//Add an 'a' to push h0 to the beginning of the alphebetically sorted list of files
			filename = string(filenameHeader) + levelSS.str() + "_ah0" + ".pgm";
			GeneratePGMimage(h0Coeffs_, horDim, verDim, coeffScalar, 128, filename);
		}
		filename = string(filenameHeader) + levelSS.str() + "_b1" + ".pgm";
		GeneratePGMimage(b1Coeffs_[i], horDim, verDim, coeffScalar, 128, filename);
		filename = string(filenameHeader) + levelSS.str() + "_b2" + ".pgm";
		GeneratePGMimage(b2Coeffs_[i], horDim, verDim, coeffScalar, 128, filename);
		filename = string(filenameHeader) + levelSS.str() + "_b3" + ".pgm";
		GeneratePGMimage(b3Coeffs_[i], horDim, verDim, coeffScalar, 128, filename);
		filename = string(filenameHeader) + levelSS.str() + "_b4" + ".pgm";
		GeneratePGMimage(b4Coeffs_[i], horDim, verDim, coeffScalar, 128, filename);
		filename = string(filenameHeader) + levelSS.str() + "_l1" + ".pgm";
		GeneratePGMimage(l1Coeffs_[i], horDim / 2, verDim / 2, thumbnailScalar, 0, filename);
	}
}

void SteerablePyramid::GenerateMaxAbsImages(const char* filenameHeader)
{
	for (int i = 0; i < numLevels_; i++)
	{
		int horDim = width_ / pow(2, i);
		int verDim = height_ / pow(2, i);
		int totalPixelsThisLevel = horDim * verDim;
		Float_t *img = new Float_t[totalPixelsThisLevel];
		
		for (int j = 0; j < totalPixelsThisLevel; j++)
		{
			Float_t meanVal = 0;
			meanVal += b1Coeffs_[i][j];
			meanVal += b2Coeffs_[i][j];
			meanVal += b3Coeffs_[i][j];
			meanVal += b4Coeffs_[i][j];
			meanVal /= 4;
			
			if (meanVal < 0)	//Min
			{
				img[j] = b1Coeffs_[i][j];
				if (b2Coeffs_[i][j] < img[j])
					img[j] = b2Coeffs_[i][j];
				if (b3Coeffs_[i][j] < img[j])
					img[j] = b3Coeffs_[i][j];
				if (b4Coeffs_[i][j] < img[j])
					img[j] = b4Coeffs_[i][j];
			}
			else if (meanVal > 0)	//Max
			{
				img[j] = b1Coeffs_[i][j];
				if (b2Coeffs_[i][j] > img[j])
					img[j] = b2Coeffs_[i][j];
				if (b3Coeffs_[i][j] > img[j])
					img[j] = b3Coeffs_[i][j];
				if (b4Coeffs_[i][j] > img[j])
					img[j] = b4Coeffs_[i][j];
				
				img[j] = -img[j];
			}
			else img[j] = 1;
		}
		
		Float_t maxCoeff = 0;
		for (int j = 0; j < totalPixelsThisLevel; j++)
			if (fabs(img[j]) > maxCoeff)
				maxCoeff = fabs(img[j]);
		
		Float_t coeffScalar = 127.0 / maxCoeff;
		
		string filename(filenameHeader);
		filename += "_";
		filename += '0' + i;
		filename += ".pgm";
		
		GeneratePGMimage(img, horDim, verDim, coeffScalar, 128, filename);
		
		delete [] img;
	}
}

void SteerablePyramid::GeneratePGMimage(Float_t *coeffs, int horDim, int verDim, Float_t scalar, Float_t shift, string filename)
{
	ofstream ofs;
	ofs.open(filename.c_str());
	assert(ofs);
	
	//Write pgm header
	ofs << "P5" << endl
		<< "#." << endl
		<< horDim << " " << verDim << endl
		<< 255 << endl;
	
	int totalPixels = horDim * verDim;
	Float_t pixelF;
	int pixelInt;
	unsigned char pixel;
	for (int i = 0; i < totalPixels; i++)
	{
		pixelF = *coeffs++ * scalar + shift;
		if ((int)(pixelF * 10.0) % 10 >= 5)
			pixelF++;
		pixelInt = pixelF;
		if (pixelInt < 0 || pixelInt > 255)
		{
			/*
			cout << "Out of bounds pixel: "
				<< horDim << "    " << verDim << "    "
				<< i % horDim << "    " << i / verDim << "    "
				<< *coeffs << "    " << pixelInt << endl;
			*/
			if (pixelInt < 0)
				pixelInt = 0;
			else pixelInt = 255;
		}
		pixel = pixelInt;
		ofs.put(pixel);
	}
	
	ofs.close();
}

void SteerablePyramid::GenerateHistograms(int numBins, char* filename)
{
	vector<Float_t> maxCoeffs, maxThumbnailPixels;
	FindStrongestCoefficients(maxCoeffs, maxThumbnailPixels);
	Float_t maxCoeff = 0;
	for (int i = 0; i < maxCoeffs.size(); i++)
		if (maxCoeffs[i] > maxCoeff)
			maxCoeff = maxCoeffs[i];
	Float_t maxThumbnailPixel = 0;
	for (int i = 0; i < maxThumbnailPixels.size(); i++)
		if (maxThumbnailPixels[i] > maxThumbnailPixel)
			maxThumbnailPixel = maxThumbnailPixels[i];
	
	Float_t binSpan = maxCoeff / numBins;
	
	//Generate coefficient histograms
	vector<vector<int> > coeffHistograms(numLevels_);
	for (int i = 0; i < numLevels_; i++)
	{
		coeffHistograms[i].resize(numBins);
		
		int horDim = width_ / pow(2, i);
		int verDim = height_ / pow(2, i);
		
		if (i == 0)
			AccumulateHistograms(h0Coeffs_, horDim, verDim, coeffHistograms[i], binSpan);
		AccumulateHistograms(b1Coeffs_[i], horDim, verDim, coeffHistograms[i], binSpan);
		AccumulateHistograms(b2Coeffs_[i], horDim, verDim, coeffHistograms[i], binSpan);
		AccumulateHistograms(b3Coeffs_[i], horDim, verDim, coeffHistograms[i], binSpan);
		AccumulateHistograms(b4Coeffs_[i], horDim, verDim, coeffHistograms[i], binSpan);
	}
	
	binSpan = maxThumbnailPixel / numBins;
	
	//Generate thumbnail histograms
	vector<vector<int> > thumbnailHistograms(numLevels_);
	for (int i = 0; i < numLevels_; i++)
	{
		thumbnailHistograms[i].resize(numBins);
		
		int horDim = width_ / pow(2, i + 1);
		int verDim = height_ / pow(2, i + 1);
		
		AccumulateHistograms(l1Coeffs_[i], horDim, verDim, thumbnailHistograms[i], binSpan);
	}
	
	//Dump the histograms to a tab delimited file
	ofstream ofs;
	ofs.open(filename);
	assert(ofs);
	
	cout << "Histograms:" << endl;
	for (int i = 0; i < numBins; i++)
	{
		for (int j = 0; j < numLevels_; j++)
		{
			ofs << coeffHistograms[j][i] << '\t';
			cout << coeffHistograms[j][i] << '\t';
		}
		for (int j = 0; j < numLevels_; j++)
		{
			ofs << thumbnailHistograms[j][i] << '\t';
			cout << thumbnailHistograms[j][i] << '\t';
		}
		ofs << endl;
		
		cout << endl;
	}
	cout << "Done" << endl;
	
	ofs.close();
}

void SteerablePyramid::AccumulateHistograms(Float_t *coeffs, int horDim, int verDim, vector<int> &histogram, Float_t binSpan)
{
	int numBins = histogram.size();
	
	int totalPixels = horDim * verDim;
	int bin;
	for (int i = 0; i < totalPixels; i++)
	{
		bin = fabs(coeffs[i]) / binSpan;
		if (bin >= numBins)
		{
			cout << "Out of bounds histogram: "
					<< fabs(coeffs[i]) << "   "
					<< binSpan << "   "
					<< fabs(coeffs[i]) / binSpan
					<< "   "
					<< bin << "   "
					<< numBins << endl;
			bin = numBins - 1;
		}
		histogram[bin]++;
	}
}

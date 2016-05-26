#ifndef __STEERABLE_PYRAMID__
#define __STEERABLE_PYRAMID__

#include "FloatType.h"
#include <vector>
#include <string>

class SteerablePyramid
{
	public:
		SteerablePyramid();
		SteerablePyramid(const SteerablePyramid &sp);
		~SteerablePyramid();
		
		const SteerablePyramid& operator=(const SteerablePyramid& sp);
		bool operator==(const SteerablePyramid& sp) const;
		
		static void GenerateFilterPGMimages(const char* filenameHeader, bool plotLog);
		static void GetFilters(std::vector<std::vector<Float_t> > &filters);
		
		const int NumLevels() const;
		void SetAllZero();
		void CalcStats(std::vector<Float_t> &means, std::vector<Float_t> &stdDevs) const;
		void SoftThreshold(std::vector<Float_t> thresholds);
		void findCorrelations(const SteerablePyramid &sp, std::vector<Float_t> &lvlCorrs) const;
		void CopyCoefficients(const SteerablePyramid& sp, std::vector<Float_t> thresholds);
		void GeneratePGMimages(const char* filenameHeader);
		void GenerateMaxAbsImages(const char* filenameHeader);
		void GenerateHistograms(int numBins, char* filename);
		
		void DoFullWaveletAnalysis(const Float_t* pixelValues, int width, int height, bool useFTconvolution);
		void DoFullWaveletSynthesis(Float_t* pixelValues);
		
	private:
		Float_t *h0Coeffs_, **b1Coeffs_, **b2Coeffs_, **b3Coeffs_, **b4Coeffs_, **l1Coeffs_;
		int width_, height_, numLevels_;
		int *yLookup_;
		bool useFTconvolution_;
		
		static Float_t *scratch_2D_A_, *scratch_2D_B_;	//Used for temp work.  Could be allocated and destroyed when needed, but that wastes a lot of time.
		static int scratchSize_;
		
		static void readFilterSet();
		static void InitFilters();
		int getFilterFTIndex();
		
		void CreateScratch();
		void Destroy();
		void Init(int width, int height);
		void DetermineNumLevels();
		
		//Wrap gives the least reconstruction error, not reflect.  Not sure why.
		static inline int handleNegEdge(int coord, int dim)
		{
			//return -1;                           //Black
			return coord + dim;                  //Wrap
			//return -coord;                       //Reflect
			
			//Note that reflect pivots around the center of the 0th and (dim-1)th pixels,
			//not at the edges of the image as might be expected.  This is necessary because,
			//when convolving an upsampled image, pixels only reside at even xy coordinates.
		}
		
		static inline int handlePosEdge(int coord, int dim)
		{
			//return -1;                           //Black
			return coord - dim;                  //Wrap
			//return (dim - 2) - (coord - dim);    //Reflect
		}
		
		void ConvolveImage(const Float_t* image, int width, int height, Float_t kernel[289], Float_t* result, bool initValues);
		void ConvolveImage(const Float_t* image, int width, int height, Float_t kernel[], int kernelDim, Float_t* result, bool initValues);
		
		void FindStrongestCoefficients(std::vector<Float_t> &maxCoeffs, std::vector<Float_t> &maxThumbnailPixels);
		void GeneratePGMimage(Float_t *coeffs, int horDim, int verDim, Float_t scalar, Float_t shift, std::string filename);
		void AccumulateHistograms(Float_t *coeffs, int horDim, int verDim, std::vector<int> &histogram, Float_t binSpan);
		
		void DoOneStepWaveletAnalysis(Float_t* pixelValues, int level);
		void DoOneStepWaveletAnalysisUsingFT(Float_t* pixelValues, int level);
		void DoOneStepWaveletSynthesis(int level, Float_t* pixelValues, Float_t *upsampledThumbnail);
		void DoOneStepWaveletSynthesisUsingFT(int level, Float_t* pixelValues, Float_t *upsampledThumbnail);
};

#endif

#ifndef __IMAGE_RGB_REAL__
#define __IMAGE_RGB_REAL__

#include "ImageInterface.h"
#include "ImageReal.h"
#include "Bayer.h"	//For Color enum
#include "FloatType.h"
#include <vector>
#include <string>
#include <assert.h>

class Image3ChnlReal : public ImageInterface
{
	public:
		Image3ChnlReal();
		Image3ChnlReal(const Image3ChnlReal& ir);
		Image3ChnlReal(int width, int height);
		
		~Image3ChnlReal();
		
		const Image3ChnlReal& operator=(const Image3ChnlReal &ir);
		bool operator==(const Image3ChnlReal &im) const;
		bool operator!=(const Image3ChnlReal &im) const;
		
		bool resize(int width, int height);
		
		Float_t& operator()(Color c, unsigned int x, unsigned int y);
		
		const Float_t& operator()(Color c, unsigned int x, unsigned int y) const;
			
		ImageReal* operator[](unsigned int i)
		{
			//assert(i < totalPixels_);
			return img_[i];
		}
				
		const ImageReal* operator[](unsigned int i) const
		{
			//assert(i < totalPixels_);
			return img_[i];
		}
		
		void scale(Float_t scalar);
		void normalize(Float_t maxVal, bool normalizeUpToo, bool independentChannels);
		void flipHorizontal() {}
		void flipVertical() {}
		void rot90CW() {}
		void rot180() {}
		void rot90CCW() {}
		
		void erase();
		
		bool readFile(std::string filename);
		std::string writeFile(int bytesPerSample) const;
		std::string writeFile(const char *filename, int bytesPerSample) const;
		
	private:
		std::vector<ImageReal*> img_;
};

#endif

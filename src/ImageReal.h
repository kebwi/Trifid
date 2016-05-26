#ifndef __IMAGE_REAL__
#define __IMAGE_REAL__

#include "ImageInterface.h"
#include "FloatType.h"
#include <string>
#include <assert.h>

class ImageReal : public ImageInterface
{
	public:
		ImageReal();
		ImageReal(const ImageReal &ir);
		ImageReal(int width, int height);
		
		virtual ~ImageReal();
		
		const ImageReal& operator=(const ImageReal &ir);
		bool operator==(const ImageReal &im) const;
		bool operator!=(const ImageReal &im) const;
		
		bool resize(int width, int height);
		
		inline Float_t& operator()(unsigned int x, unsigned int y)
		{
			//assert(x < width_ && y < height_);
			return img_[y * width_ + x];
		}
		
		inline const Float_t& operator()(unsigned int x, unsigned int y) const
		{
			//assert(x < width_ && y < height_);
			return img_[y * width_ + x];
		}
		
		Float_t& operator[](unsigned int i)
		{
			//assert(i < totalPixels_);
			return img_[i];
		}
		
		const Float_t& operator[](unsigned int i) const
		{
			//assert(i < totalPixels_);
			return img_[i];
		}
		
		inline Float_t *img()
		{
			return img_;
		}
		
		void scale(Float_t scalar);
		void normalize(Float_t maxVal, bool normalizeUpToo);
		void flipHorizontal();
		void flipVertical();
		void rot90CW();
		void rot180();
		void rot90CCW();
		
		void erase();
		void downsample(ImageReal &imageOut, int startX, int startY) const;
		void upsample(ImageReal &imageOut, int startX, int startY) const;
		
		bool readFile(std::string filename);
		std::string writeFile(int bytesPerSample) const;
		std::string writeFile(const char *filename, int bytesPerSample) const;
		
	private:
		Float_t *img_;
};

#endif

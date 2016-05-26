#ifndef __IMAGE_INTERFACE__
#define __IMAGE_INTERFACE__

#include "FloatType.h"
#include <string>

class ImageInterface
{
	public:
		ImageInterface();
		ImageInterface(const ImageInterface &ir);
		ImageInterface(int width, int height);
		ImageInterface(std::string filename);
		
		virtual ~ImageInterface();
		
		virtual const ImageInterface& operator=(const ImageInterface &ir);
		virtual bool operator==(const ImageInterface &im) const;
		virtual bool operator!=(const ImageInterface &im) const;
		
		static void stripString(std::string &str);
		
		//return false is no change is made (the image is already the new size)
		virtual bool resize(int width, int height);
		
		inline std::string filename() const
		{
			return filename_;
		}
		
		inline int width() const
		{
			return width_;
		}
		
		inline int height() const
		{
			return height_;
		}
		
		inline int totalPixels() const
		{
			return totalPixels_;
		}
		
		void setFilename(std::string filename);
		inline void setFilename(const char *filename)
		{
			setFilename(std::string(filename));
		}
		
		void appendFilename(std::string filename);
		inline void appendFilename(const char *filename)
		{
			appendFilename(std::string(filename));
		}
		
		virtual void scale(Float_t scalar) = 0;
		virtual void flipHorizontal() = 0;
		virtual void flipVertical() = 0;
		virtual void rot90CW() = 0;
		virtual void rot180() = 0;
		virtual void rot90CCW() = 0;
		
		virtual bool readFile(std::string filename) = 0;
		virtual std::string writeFile(int bytesPerSample) const = 0;
		virtual std::string writeFile(const char *filename, int bytesPerSample) const = 0;
		
	protected:
		std::string filename_;
		int width_, height_, totalPixels_;
};

#endif

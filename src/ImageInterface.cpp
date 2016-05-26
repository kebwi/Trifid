#include "ImageInterface.h"
#include <assert.h>

using namespace std;

ImageInterface::ImageInterface() :
	filename_(""),
	width_(0),
	height_(0),
	totalPixels_(0)
{
}

ImageInterface::ImageInterface(const ImageInterface &ir) :
	filename_(""),
	width_(0),
	height_(0),
	totalPixels_(0)
{
	*this = ir;
}

ImageInterface::ImageInterface(int width, int height) :
	filename_(""),
	width_(width),
	height_(height),
	totalPixels_(width * height)
{
}

ImageInterface::ImageInterface(string filename) :
	filename_(filename),
	width_(0),
	height_(0),
	totalPixels_(0)
{
}

ImageInterface::~ImageInterface()
{
}

const ImageInterface& ImageInterface::operator=(const ImageInterface &ir)
{
	if (&ir == this)
		return *this;
	
	filename_ = ir.filename_;
	resize(ir.width_, ir.height_);
	
	return *this;
}

bool ImageInterface::operator==(const ImageInterface &im) const
{
	if (&im == this)
		return true;
	
	//Notice, no check against the filename
	
	if (width_ != im.width_ || height_ != im.height_)
		return false;
	assert(totalPixels_ == im.totalPixels_);
	
	return true;
}

bool ImageInterface::operator!=(const ImageInterface &im) const
{
	return !(im == *this);
}

//Strip extensions prefixes with '.' and/or '_' off end of string
//static
void ImageInterface::stripString(string &str)
{
	if (str.length() >= 4)	//Strip extension (presumably ".pgm" or ".ppm")
	{
		if (str[str.length() - 4] == '.'
			&& (isalpha(str[str.length() - 3]) || isalpha(str[str.length() - 2]) || isalpha(str[str.length() - 1])))
			str = str.substr(0, str.length() - 4);
	}
	
	for (int i = 0; i < 2; i++)	//Do this twice, once for "_out" and once for "_ori"
	{
		if (str.length() >= 4)	//Strip "_out" and "_ori" suffixes
		{
			if (str.substr(str.length() - 4, 4) == "_out" ||
				str.substr(str.length() - 4, 4) == "_ori")
				str = str.substr(0, str.length() - 4);
		}
	}
}

void ImageInterface::setFilename(string filename)
{
	filename_ = filename;
	stripString(filename_);
}

void ImageInterface::appendFilename(string filenameTail)
{
	stripString(filename_);
	filename_ += string(filenameTail);
	stripString(filename_);
}

bool ImageInterface::resize(int width, int height)
{
	if (width_ == width && height_ == height)
		return false;
	
	width_ = width;
	height_ = height;
	totalPixels_ = width_ * height_;
	
	return true;
}

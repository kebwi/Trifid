#include "ImageReal.h"
#include "main.h"
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>

using namespace std;

//******************************************************************************
//Extern Globals

//******************************************************************************
//Global Declarations

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

ImageReal::ImageReal() :
	img_(NULL)
{
}

ImageReal::ImageReal(const ImageReal &ir) :
	img_(NULL)
{
	*this = ir;
}

ImageReal::ImageReal(int width, int height) :
	ImageInterface(width, height),
	img_(NULL)
{
	img_ = new Float_t[totalPixels_];
	erase();
}

ImageReal::~ImageReal()
{
	if (img_)
		delete [] img_;
}

const ImageReal& ImageReal::operator=(const ImageReal &ir)
{
	if (&ir == this)
		return *this;
	
	ImageInterface::operator=(ir);
	
	resize(ir.width_, ir.height_);
	memcpy(img_, ir.img_, totalPixels_ * sizeof(Float_t));
	
	return *this;
}

bool ImageReal::operator==(const ImageReal &im) const
{
	if (&im == this)
		return true;
	
	if (!(ImageInterface::operator==(im)))
		return false;
	
	return (memcmp(img_, im.img_, totalPixels_ * sizeof(Float_t)) == 0);
}

bool ImageReal::operator!=(const ImageReal &im) const
{
	return !(im == *this);
}

bool ImageReal::resize(int width, int height)
{
	if (!ImageInterface::resize(width, height))
	{
		erase();
		return false;
	}
	
	if (img_)
		delete [] img_;
	img_ = new Float_t[totalPixels_];
	erase();
	
	return true;
}

void ImageReal::scale(Float_t scalar)
{
	for (int i = 0; i < totalPixels_; i++)
		img_[i] *= scalar;
}

void ImageReal::normalize(Float_t maxVal, bool normalizeUpToo)
{
	//Normalize the absolute value, some "images" might contain negative pixels
	Float_t maxCurrentVal = 0;
	for (int i = 0; i < totalPixels_; i++)
		if (fabs(img_[i]) > maxCurrentVal)
			maxCurrentVal = fabs(img_[i]);
	
	if (maxCurrentVal < maxVal && !normalizeUpToo)
		return;
	if (maxCurrentVal == 0)
		return;
	
	Float_t scalar = maxVal / maxCurrentVal;
	scale(scalar);
	
	//This doesn't really have anything to do with normalization
	for (int i = 0; i < totalPixels_; i++)
		img_[i] = fabs(img_[i]);
}

void ImageReal::flipHorizontal()
{
	Float_t temp;
	for (int y = 0; y < height_; y++)
	{
		int halfWidth = width_ / 2;
		for (int x = 0; x < halfWidth; x++)
		{
			temp = (*this)(x, y);
			(*this)(x, y) = (*this)((width_ - 1) - x, y);
			(*this)((width_ - 1) - x, y) = temp;
		}
	}
}

void ImageReal::flipVertical()
{
	Float_t temp;
	int halfHeight = height_ / 2;
	for (int y = 0; y < halfHeight; y++)
		for (int x = 0; x < width_; x++)
		{
			temp = (*this)(x, y);
			(*this)(x, y) = (*this)(x, (height_ - 1) - y);
			(*this)(x, (height_ - 1) - y) = temp;
		}
}

void ImageReal::rot90CW()
{
	ImageReal temp(width_, height_);
	Float_t *tempImgPtr = temp.img_;
	for (int y = 0; y < height_; y++)
		for (int x = 0; x < width_; x++)
			*tempImgPtr++ = img_[((height_ - 1) - x) * width_ + y];
	*this = temp;
}

void ImageReal::rot180()
{
	Float_t temp;
	int halfTotalPixels = totalPixels_ / 2;
	for (int i = 0; i < halfTotalPixels; i++)
	{
		temp = img_[i];
		img_[i] = img_[(totalPixels_ - 1) - i];
		img_[(totalPixels_ - 1) - i] = temp;
	}
}

void ImageReal::rot90CCW()
{
	ImageReal temp(width_, height_);
	Float_t *tempImgPtr = temp.img_;
	for (int y = 0; y < height_; y++)
		for (int x = 0; x < width_; x++)
			*tempImgPtr++ = img_
				[x * width_ + ((width_ - 1) - y)];
	*this = temp;
}

void ImageReal::erase()
{
	memset(img_, 0, totalPixels_ * sizeof(Float_t));
	//for (int i = 0; i < totalPixels_; i++)
	//	img_[i] = 0;
}

void ImageReal::downsample(ImageReal &imageOut, int startX, int startY) const
{
	assert(width_ % 2 == 0 && height_ % 2 == 0);
	
	int	widthOut = width_ / 2, heightOut = height_ / 2;
	if (imageOut.width() != widthOut || imageOut.height() != heightOut)
		imageOut.resize(widthOut, heightOut);
	
	int imageOutIdx = 0;
	for (int y = startY; y < height_; y += 2)
		for (int x = startX; x < width_; x += 2, ++imageOutIdx)
			imageOut[imageOutIdx] = img_[y * width_ + x];
}

//Note that the skipped rows and columns in the upsampled output are specifically
//*NOT* cleared to 0.  This permits upsampling to be cumulatively performed for
//each of the four bayer positions into a single final image.
void ImageReal::upsample(ImageReal &imageOut, int startX, int startY) const
{
	assert(width_ % 2 == 0 && height_ % 2 == 0);
	
	int widthOut = width_ * 2, heightOut = height_ * 2;
	if (imageOut.width() != widthOut || imageOut.height() != heightOut)
		imageOut.resize(widthOut, heightOut);
	
	int imageInIdx = 0;
	for (int y = startY; y < heightOut; y += 2)
		for (int x = startX; x < widthOut; x += 2, ++imageInIdx)
			imageOut[y * widthOut + x] = img_[imageInIdx];
}

bool ImageReal::readFile(string filename)
{
	filename_ = filename;
	
	ifstream ifs;
	ifs.open(filename_.c_str());
	
	//In theory, this function could return false on failure instead of dieing.
	if (!ifs)
	{
		string error = "Input filename error (file not found perhaps): " + filename_;
		die(error.c_str());
	}
	
	//Read pgm header
	string magicNumber, note;
	int width, height, maxValue;
	ifs >> magicNumber;
	assert(ifs.get() == '\n');
	while (ifs.peek() != '\n')
		note += ifs.get();
	ifs >> width >> height >> maxValue;
	assert(ifs.get() == '\n');
	
	if (magicNumber == "P3")
		die("Input file is a ASCII PGM file, binary PGM files are required");
	if (magicNumber != "P5")
		die("Input file does not appear to be a PGM file");
	/*
	cout << "PGM file opened for reading: " << filename_ << endl;
	cout << "  PGM header data:" << endl;
	cout << "    PGM magic number (should be P5):  " << magicNumber << endl;
	cout << "    Note:                             " << note << endl;
	cout << "    Width, Height:                    " << width << "   " << height << endl;
	cout << "    Max value:                        " << maxValue << endl;
	*/
	//Resize the images
	resize(width, height);
	
	//Read input image
	unsigned short maxPixelValueRead = 0;
	Float_t *imgPtr = img_;
	
	if (maxValue <= 255)	//8 bit image
	{
		unsigned char pixelReadC;
		for (int i = 0; i < totalPixels_; i++)
		{
			assert(ifs);
			pixelReadC = ifs.get();
			if (pixelReadC > maxPixelValueRead)
				maxPixelValueRead = pixelReadC;
			*imgPtr++ = pixelReadC / 255.0;
			assert(ifs);
		}
	}
	else	//16 bit image
	{
		unsigned short pixelReadS;
		char *charPtr;
		for (int i = 0; i < totalPixels_; i++)
		{
			assert(ifs);
			charPtr = (char*)&pixelReadS;
			*charPtr++ = ifs.get();
			*charPtr = ifs.get();
			if (pixelReadS > maxPixelValueRead)
				maxPixelValueRead = pixelReadS;
			*imgPtr++ = pixelReadS / 65535.0;
			assert(ifs);
		}
	}
	ifs.get();
	assert(!ifs);
	ifs.close();
	
	//cout << "  Max pixel value read: " << maxPixelValueRead << endl << endl;
	
	return true;
}

string ImageReal::writeFile(int bytesPerSample) const
{
	return writeFile(filename_.c_str(), bytesPerSample);
}

string ImageReal::writeFile(const char *filename, int bytesPerSample) const
{
	assert(bytesPerSample == 1 || bytesPerSample == 2);
	
	string filenameStr(filename);
	stripString(filenameStr);
	
	int numPasses = (bytesPerSample == 1 ? 1 : 2);
	for (int pass = 0; pass < numPasses; pass++)
	{
		bool flipEndian = false;
		if (pass != 0)
			flipEndian = true;
		
		string filenameOut = filenameStr;
		if (bytesPerSample == 2)
			filenameOut += "_16b";
		if (flipEndian && bytesPerSample == 2)
			filenameOut += "_endianFlip";
		filenameOut += "_out.pgm";
		
		ofstream ofs;
		ofs.open(filenameOut.c_str());
		assert(ofs);
		
		//Write pgm header
		cout << "PGM file opened for writing: " << filenameOut << endl << endl;
		ofs << "P5" << endl
			<< "#." << endl
			<< width_ << " " << height_ << endl
			<< ((bytesPerSample == 1) ? 255 : 65535) << endl;
		
		if (bytesPerSample == 1)
		{
			cout << "Writing 8-bit PGM file..." << endl;
			
			char pixelC;
			for (int i = 0; i < totalPixels_; i++)
			{
				pixelC = roundToNearestInt(img_[i] * 255.0, 0, 255);
				ofs.put(pixelC);
			}
		}
		else    //bytesPerSample == 2
		{
			if (flipEndian)
				cout << "Writing Endian-flipped 16-bit PGM file..." << endl;
			else cout << "Writing 16-bit PGM file..." << endl;
			
			char *charPtr;
			unsigned short pixelS;
			for (int i = 0; i < totalPixels_; i++)
			{
				pixelS = roundToNearestInt(img_[i] * 65535.0, 0, 65535);
				
				charPtr = (char*)&pixelS;
				
				if (!flipEndian)
				{
					ofs.put(*charPtr++);
					ofs.put(*charPtr);
				}
				else
				{
					ofs.put(*++charPtr);
					ofs.put(*--charPtr);
				}
			}
		}
		
		ofs.close();
	}
	
	return string("");//filenameOut;
}

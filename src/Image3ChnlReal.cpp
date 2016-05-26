#include "Image3ChnlReal.h"
#include "ImageReal.h"
#include "main.h"
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
void die(char *msg);

Image3ChnlReal::Image3ChnlReal()
{
	img_.push_back(new ImageReal());
	img_.push_back(new ImageReal());
	img_.push_back(new ImageReal());
}

Image3ChnlReal::Image3ChnlReal(const Image3ChnlReal &ir)
{
	img_.push_back(new ImageReal());
	img_.push_back(new ImageReal());
	img_.push_back(new ImageReal());
	
	*this = ir;
}

Image3ChnlReal::Image3ChnlReal(int width, int height) :
	ImageInterface(width, height)
{
	img_.push_back(new ImageReal(width_, height_));
	img_.push_back(new ImageReal(width_, height_));
	img_.push_back(new ImageReal(width_, height_));
}

Image3ChnlReal::~Image3ChnlReal()
{
	assert(img_.size() == 3);
	delete img_[0];
	delete img_[1];
	delete img_[2];
}

const Image3ChnlReal& Image3ChnlReal::operator=(const Image3ChnlReal &ir)
{
	if (&ir == this)
		return *this;
	
	ImageInterface::operator=(ir);
	
	resize(ir.width_, ir.height_);
	*img_[0] = *ir.img_[0];
	*img_[1] = *ir.img_[1];
	*img_[2] = *ir.img_[2];
	
	return *this;
}

bool Image3ChnlReal::operator==(const Image3ChnlReal &im) const
{
	if (&im == this)
		return true;
	
	if (!(ImageInterface::operator==(im)))
		return false;
	
	if (*img_[0] != *im.img_[0])
		return false;
	if (*img_[1] != *im.img_[1])
		return false;
	if (*img_[2] != *im.img_[2])
		return false;
	
	return true;
}

bool Image3ChnlReal::operator!=(const Image3ChnlReal &im) const
{
	return !(im == *this);
}

bool Image3ChnlReal::resize(int width, int height)
{
	if (!ImageInterface::resize(width, height))
	{
		erase();
		return false;
	}
	
	assert(img_.size() == 3);
	
	img_[0]->resize(width_, height_);
	img_[1]->resize(width_, height_);
	img_[2]->resize(width_, height_);
	
	return true;
}

Float_t& Image3ChnlReal::operator()(Color c, unsigned int x, unsigned int y)
{
	//assert(x < width_ && y < height_);
	return (*(img_[(int)c]))[y * width_ + x];
}

const Float_t& Image3ChnlReal::operator()(Color c, unsigned int x, unsigned int y) const
{
	//assert(x < width_ && y < height_);
	return (*(img_[(int)c]))[y * width_ + x];
}

void Image3ChnlReal::scale(Float_t scalar)
{
	img_[0]->scale(scalar);
	img_[1]->scale(scalar);
	img_[2]->scale(scalar);
}

void Image3ChnlReal::normalize(Float_t maxVal, bool normalizeUpToo, bool independentChannels)
{
	if (independentChannels)
	{
		img_[0]->normalize(maxVal, normalizeUpToo);
		img_[1]->normalize(maxVal, normalizeUpToo);
		img_[2]->normalize(maxVal, normalizeUpToo);
	}
	else
	{
		Float_t maxCurrentVal = 0;
		for (int i = 0; i < totalPixels_; i++)
		{
			if ((*img_[0])[i] > maxCurrentVal)
				maxCurrentVal = (*img_[0])[i];
			if ((*img_[1])[i] > maxCurrentVal)
				maxCurrentVal = (*img_[1])[i];
			if ((*img_[2])[i] > maxCurrentVal)
				maxCurrentVal = (*img_[2])[i];
		}
		//cout << "maxCurrentVal before RGB normalization: " << maxCurrentVal << endl;
		/*
		int minCurrentVal = 999999;
		for (int i = 0; i < totalPixels_; i++)
		{
			if ((*img_[0])[i] < minCurrentVal)
				minCurrentVal = (*img_[0])[i];
			if ((*img_[1])[i] < minCurrentVal)
				minCurrentVal = (*img_[1])[i];
			if ((*img_[2])[i] < minCurrentVal)
				minCurrentVal = (*img_[2])[i];
		}
		cout << "minCurrentVal before RGB normalization: " << minCurrentVal << endl;
		*/
		if (maxCurrentVal < maxVal && !normalizeUpToo)
			return;
		if (maxCurrentVal == 0)
			return;
		
		Float_t scalar = maxVal / maxCurrentVal;
		
		scale(scalar);
		/*
		maxCurrentVal = 0;
		for (int i = 0; i < totalPixels_; i++)
		{
			if ((*img_[0])[i] > maxCurrentVal)
				maxCurrentVal = (*img_[0])[i];
			if ((*img_[1])[i] > maxCurrentVal)
				maxCurrentVal = (*img_[1])[i];
			if ((*img_[2])[i] > maxCurrentVal)
				maxCurrentVal = (*img_[2])[i];
		}
		cout << "maxCurrentVal after RGB normalization: " << maxCurrentVal << endl;
		minCurrentVal = 999999;
		for (int i = 0; i < totalPixels_; i++)
		{
			if ((*img_[0])[i] < minCurrentVal)
				minCurrentVal = (*img_[0])[i];
			if ((*img_[1])[i] < minCurrentVal)
				minCurrentVal = (*img_[1])[i];
			if ((*img_[2])[i] < minCurrentVal)
				minCurrentVal = (*img_[2])[i];
		}
		cout << "minCurrentVal after RGB normalization: " << minCurrentVal << endl;
		*/
	}
}

void Image3ChnlReal::erase()
{
	img_[0]->erase();
	img_[1]->erase();
	img_[2]->erase();
}

bool Image3ChnlReal::readFile(string filename)
{
	filename_ = filename;
	
	ifstream ifs;
	ifs.open(filename_.c_str());
	if (!ifs)
	{
		string error = "Input filename error (file not found perhaps): " + filename_;
		die(error.c_str());
	}
	
	//Read ppm header
	string magicNumber, note;
	int width, height, maxValue;
	ifs >> magicNumber;
	assert(ifs.get() == '\n');
	while (ifs.peek() != '\n')
		note += ifs.get();
	ifs >> width >> height >> maxValue;
	assert(ifs.get() == '\n');
	
	if (magicNumber == "P3")
		die("Input file is a ASCII PPM file, binary PPM files are required");
	if (magicNumber != "P6")
		die("Input file does not appear to be a PPM file");
	/*
	cout << "PPM file read: " << filename_ << endl;
	cout << "  PPM header data:" << endl;
	cout << "    PPM magic number (should be P6):  " << magicNumber << endl;
	cout << "    Note:                             " << note << endl;
	cout << "    Width, Height:                    " << width << "   " << height << endl;
	cout << "    Max value:                        " << maxValue << endl;
	*/
	resize(width, height);
	
	//Read input image
	unsigned short maxPixelValueRead = 0;
	
	if (maxValue <= 255)	//8 bit image
	{
		unsigned char pixelReadC;
		for (int i = 0; i < totalPixels_; i++)
		{
			assert(ifs);
			
			for (int j = 0; j < 3; j++)
			{
				pixelReadC = ifs.get();
				if (pixelReadC > maxPixelValueRead)
					maxPixelValueRead = pixelReadC;
				(*img_[j])[i] = pixelReadC / 255.0;
			}
			
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
			
			for (int j = 0; j < 3; j++)
			{
				charPtr = (char*)&pixelReadS;
				*charPtr++ = ifs.get();
				*charPtr = ifs.get();
				if (pixelReadS > maxPixelValueRead)
					maxPixelValueRead = pixelReadS;
				(*img_[j])[i] = pixelReadS / 65535.0;
			}
			
			assert(ifs);
		}
	}
	ifs.get();
	assert(!ifs);
	ifs.close();
	
	//cout << "  Max pixel value read: " << maxPixelValueRead << endl << endl;
	
	//In theory, this function could return false on failure instead of dieing.
	return true;
}

string Image3ChnlReal::writeFile(int bytesPerSample) const
{
	return writeFile(filename_.c_str(), bytesPerSample);
}

string Image3ChnlReal::writeFile(const char *filename, int bytesPerSample) const
{
	assert(bytesPerSample == 1 || bytesPerSample == 2);
	assert(img_.size() == 3);
	
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
		filenameOut += "_out.ppm";
		
		ofstream ofs;
		ofs.open(filenameOut.c_str());
		assert(ofs);
		
		//Write ppm header
		cout << "PPM file opened for writing: " << filenameOut << endl;
		ofs << "P6" << endl
			<< "#." << endl
			<< width_ << " " << height_ << endl
			<< ((bytesPerSample == 1) ? 255 : 65535) << endl;
		
		//Write the data in interlaced RGB order
		if (bytesPerSample == 1)
		{
			cout << "Writing 8-bit PPM file..." << endl;
			
			char pixelC;
			for (int i = 0; i < totalPixels_; i++)
				for (int j = 0; j < 3; j++)
				{
					pixelC = roundToNearestInt((*img_[j])[i] * 255.0, 0, 255);
					ofs.put(pixelC);
				}
		}
		else    //bytesPerSample == 2
		{
			if (flipEndian)
				cout << "Writing Endian-flipped 16-bit PPM file..." << endl;
			else cout << "Writing 16-bit PPM file..." << endl;
			
			char *charPtr;
			unsigned short pixelS;
			for (int i = 0; i < totalPixels_; i++)
				for (int j = 0; j < 3; j++)
				{
					pixelS = roundToNearestInt((*img_[j])[i] * 65535.0, 0, 65535);
					
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
		
		cout << endl;
		ofs.close();
	}
	
	return string("");//filenameOut;
}

#include "main.h"
#include "Random.h"
#include "Bayer.h"
#include "SteerablePyramid.h"
#include "FourierTransform.h"
#include "ImagePairStats.h"
#include "ImageProcs.h"
#include "DebayerProcs.h"
#include "ImageReal.h"
#include "Image3ChnlReal.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;
struct ImageStats
{
	Image3ChnlReal rgbImg;
	Image3ChnlReal labImg;
	
	Float_t rPsnr, gPsnr, bPsnr;
	Float_t rgbPsnr, labPsnr, labDeltaE;
	
	bool minRPsnr, minGPsnr, minBPsnr, minRGBPsnr, minLabPsnr, minLabDeltaE;
	bool maxRPsnr, maxGPsnr, maxBPsnr, maxRGBPsnr, maxLabPsnr, maxLabDeltaE;
	
	char labels[4];
	
	ImageStats() :
		rPsnr(0), gPsnr(0), bPsnr(0), rgbPsnr(0), labPsnr(0), labDeltaE(0),
		minRPsnr(false), minGPsnr(false), minBPsnr(false), minRGBPsnr(false), minLabPsnr(false), minLabDeltaE(false),
		maxRPsnr(false), maxGPsnr(false), maxBPsnr(false), maxRGBPsnr(false), maxLabPsnr(false), maxLabDeltaE(false)
	{
		for (int i = 0; i < 4; i++)
			labels[i] = '\0';
	}
	
	struct RPsnrSorter
	{
		bool operator()(const ImageStats& is1, const ImageStats& is2)
		{
			return (is1.rPsnr > is2.rPsnr);
		}
	};
	struct GPsnrSorter
	{
		bool operator()(const ImageStats& is1, const ImageStats& is2)
		{
			return (is1.gPsnr > is2.gPsnr);
		}
	};
	struct BPsnrSorter
	{
		bool operator()(const ImageStats& is1, const ImageStats& is2)
		{
			return (is1.bPsnr > is2.bPsnr);
		}
	};
	struct RGBPsnrSorter
	{
		bool operator()(const ImageStats& is1, const ImageStats& is2)
		{
			return (is1.rgbPsnr > is2.rgbPsnr);
		}
	};
	struct LabPsnrSorter
	{
		bool operator()(const ImageStats& is1, const ImageStats& is2)
		{
			return (is1.labPsnr < is2.labPsnr);
		}
	};
	struct LabDeltaESorter
	{
		bool operator()(const ImageStats& is1, const ImageStats& is2)
		{
			return (is1.labDeltaE < is2.labDeltaE);
		}
	};
};

//******************************************************************************
//Extern Globals
extern int gBitDepth;

//******************************************************************************
//Global Declarations
enum ProgramBehavior {
	PB_NONE,
	PB_TEST_FOURIER, PB_TEST_FOURIER_CONVOLUTION, PB_TEST_RGB_XYZ_LAB,
	PB_NOISIFY, PB_DENOISE,
	PB_GRAY_IMG_PAIR_STATS, PB_RGB_IMG_PAIR_STATS,
	PB_NUM_PROGRAM_BEHAVIORS };

enum IdealShrinkageMethod {
	INDEPENDENT, ASCEND, DESCEND, NUM_IDEAL_SHRINKAGE_METHODS };

static const string skProgramBehaviorLabels[PB_NUM_PROGRAM_BEHAVIORS + 1] = {
	"PB_NONE",
	"PB_TEST_FOURIER", "PB_TEST_FOURIER_CONVOLUTION", "PB_TEST_RGB_XYZ_LAB",
	"PB_NOISIFY", "PB_DENOISE",
	"PB_GRAY_IMG_PAIR_STATS", "PB_RGB_IMG_PAIR_STATS",
	"PB_NUM_PROGRAM_BEHAVIORS" };

static const string skIdealShrinkageMethodLabels[NUM_IDEAL_SHRINKAGE_METHODS + 1] = {
	"INDEPENDENT", "ASCEND", "DESCEND", "NUM_IDEAL_SHRINKAGE_METHODS" };

bool gVerbose = false;
bool gGenerateIntermediateImages = false;
static ProgramBehavior sProgramBehavior = PB_NONE;
static bool sBayerize = false;
static string sHotpixelMapFilename = "";
static string sNoiseParamsFilename = "";
static vector<string> sInputFilenames;

//One vector per RGB color, each vector a set of levels starting with H, then each descending B after that.
vector<Float_t> gManualShrinkageThresholds[3];
string sManualShrinkageThresholdStr = "";	//Used to generate filename

//Noise parameter combinations, as indicated in the filenames
//First parameter is a label for the filename
static NoiseParam sNoiseParams;

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

#pragma mark -

void die(const char *msg)
{
	cerr << msg << endl;
	exit(1);
}

void warning(const char *msg)
{
	cout << endl;
	cout << "****************************************" << endl;
	cout << "  !!! WARNING !!!" << endl;
	cout << "  " << msg << endl;
	cout << "****************************************" << endl;
	cout << endl;
}

void myAssert(bool expression)
{
	assert(expression);
	if (!expression)
	{
		die("myAssert failure");
	}
}

void assertFloatEqual(const Float_t &f1, const Float_t &f2)
{
	Float_t diff = f1 - f2;
	myAssert(fabs(diff) < 1E-10);
}

//Pass the same value for minV and maxV to indicate no clipping 
int roundToNearestInt(Float_t valueF, int minV, int maxV)
{
	/*
	int value = 0;
	if (((int)(valueF * 10) % 10) < 5)
		value = valueF;
	value = valueF + 1;
	*/
	int value = (int)(valueF + .5);
	
	if (minV != maxV)
	{
		if (value < minV)
			value = minV;
		else if (value > maxV)
			value = maxV;
	}
	
	return value;
}

bool floatsEqual(Float_t f1, Float_t f2)
{
	return (fabs(f1 - f2) < .0001);
}

unsigned long getCurrentTime_ms()
{
	static timeval tp;
	gettimeofday(&tp, NULL);
	
	static unsigned long currentTime_ms;
	currentTime_ms = tp.tv_sec * CURRENT_TIME_PER_SEC;	//Convert seconds to milliseconds
	currentTime_ms += tp.tv_usec / CURRENT_TIME_PER_SEC;	//Convert microseconds to milliseconds
	
	return currentTime_ms;
}

void displayHelp()
{
	cout << endl << endl << "Flags:" << endl;
	cout << setw(4) << "" << setw(20) << left << "-h" << "Display program help" << endl;
	cout << setw(4) << "" << setw(20) << left << "-v" << "Verbose output" << endl;
	cout << setw(4) << "" << setw(20) << left << "-i" << "Produce all intermediate images" << endl;
	
	//========================================
	
	cout << endl;
	cout << setw(4) << "" << setw(20) << left << "--testFT" << "Test Fourier Transform Output and Speed" << endl;
	cout << setw(4) << "" << setw(20) << left << "--testFTconv" << "Test Fourier Convolution" << endl;
	cout << setw(4) << "" << setw(20) << left << "--testColor" << "Test RGB XYZ LAB Color Routines" << endl;
	
	//========================================
	
	cout << endl;
	
	cout << setw(4) << "" << setw(20) << left << "-N" << "Add simulated sensor noise" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "In:  1+ pgm/ppm" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "Out: 1+ pgm/ppm" << endl;
	
	cout << setw(4) << "" << setw(20) << left << "-b" << "Convert RGB input to Bayer gray before processing" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "Only valid when following -N in the arg list" << endl;
	
	cout << setw(4) << "" << setw(20) << left << "-n" << "Denoise images" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "In:  1+ pgm/ppm" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "Out: 1+ pgm/ppm" << endl;
	
	cout << setw(4) << "" << setw(20) << left << "-s [threshold]" << "Specify a shrinkage threshold (see README for syntax)" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "Only valid when following -n in the arg list" << endl;
	
	cout << setw(4) << "" << setw(20) << left << "--grayStat" << "Calculate grayscale image statistics" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "In:  1  pgm, original" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "In:  1+ pgm, result" << endl;
	
	cout << setw(4) << "" << setw(20) << left << "--rgbStat" << "Calculate RGB image statistics" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "In:  1  ppm, original" << endl;
	cout << setw(4) << "" << setw(20) << "" << setw(4) << "" << "In:  1+ ppm, result" << endl;
	
	cout << endl;
	
	exit(0);
}

void processArgs(int argc, char **argv)
{
	if (argc < 2)
		die("Enter -h flag to display help.");
	
	for (int argIdx = 1; argIdx < argc; argIdx++)
	{
		string prevArg(argv[argIdx - 1]);
		string argStr(argv[argIdx]);
		
		if (argStr[0] == '-')// && isalpha(argStr[1]))	//Flag
		{
			if (argStr.length() < 2)
				die("Flag is too short");
			if (argStr[1] != '-')	//A short flag
			{
				for (int charIdx = 1; charIdx < argStr.length(); charIdx++)
				{
					switch (argStr[charIdx])
					{
						case 'h':
							cout << "-" << argStr[charIdx] << " flag detected" << endl;
							displayHelp();
							exit(0);	//displayHelp() should exit the program
							break;
						case 'v':
							gVerbose = true;
							cout << "-" << argStr[charIdx] << " flag detected" << endl;
							break;
						case 'i':
							gGenerateIntermediateImages = true;
							cout << "-" << argStr[charIdx] << " flag detected" << endl;
							break;
						case 'N':	//Add noise to a pgm image
							cout << "-" << argStr[charIdx] << " flag detected" << endl;
							if (sProgramBehavior != PB_NONE)
								die("Please specify only one program behavior flag");
							sProgramBehavior = PB_NOISIFY;
							break;
						case 'n':	//Denoise a set of images
							cout << "-" << argStr[charIdx] << " flag detected" << endl;
							if (sProgramBehavior != PB_NONE && sProgramBehavior != PB_DENOISE)
								die("Please specify only one program behavior flag");
							sProgramBehavior = PB_DENOISE;
							break;
						case 'b':
							cout << "-" << argStr[charIdx] << " flag detected" << endl;
							if (sProgramBehavior != PB_NOISIFY)
								die("-b only works with -N");
							sBayerize = true;
							break;
						case 's':
							cout << "-" << argStr[charIdx] << " flag detected" << endl;
							if (sProgramBehavior != PB_DENOISE)
								die("-s only works with -n");
							break;
						default:
							die(string("Unknown flag: " + argStr).c_str());
					}	//switch argStr[charIdx]
				}	//Loop over flag chars
			}
			else	//argStr[1] == '-'
			{
				assert(argStr.length() >= 2);
				assert(argStr[0] == '-' && argStr[1] == '-');
				if (argStr.length() < 3)
					die("Flag is too short");
				
				if (argStr == "--testFT")
				{
					cout << argStr << " flag detected" << endl;
					if (sProgramBehavior != PB_NONE)
						die("Please specify only one program behavior flag");
					sProgramBehavior = PB_TEST_FOURIER;
				}
				else if (argStr == "--testFTconv")
				{
					cout << argStr << " flag detected" << endl;
					if (sProgramBehavior != PB_NONE)
						die("Please specify only one program behavior flag");
					sProgramBehavior = PB_TEST_FOURIER_CONVOLUTION;
				}
				else if (argStr == "--testColor")
				{
					cout << argStr << " flag detected" << endl;
					if (sProgramBehavior != PB_NONE)
						die("Please specify only one program behavior flag");
					sProgramBehavior = PB_TEST_RGB_XYZ_LAB;
				}
				else if (argStr == "--grayStat")
				{
					cout << argStr << " flag detected" << endl;
					if (sProgramBehavior != PB_NONE)
						die("Please specify only one program behavior flag");
					sProgramBehavior = PB_GRAY_IMG_PAIR_STATS;
				}
				else if (argStr == "--rgbStat")
				{
					cout << argStr << " flag detected" << endl;
					if (sProgramBehavior != PB_NONE)
						die("Please specify only one program behavior flag");
					sProgramBehavior = PB_RGB_IMG_PAIR_STATS;
				}
				else die(string("Unknown flag: " + argStr).c_str());
			}
		}
		else	//Nonflag
		{
			if (prevArg == "-N")
			{
				//Noise params file
				sNoiseParamsFilename = argStr;
				cout << "Noise params filename arg read: " << sNoiseParamsFilename << endl;
			}
			else if (prevArg == "-s")
			{
				//Denoise threshold
				
				string color = argStr.substr(0, 1);
				if (color[0] < '0' || color[0] > '3')
					die("Invalid shrinkage color argument");
				int colorIdx = color[0] - '1';
				
				string levelStr = argStr.substr(1, 2);
				
				char levelType = levelStr[0];
				if (levelType != 'A' && levelType != 'H' && levelType != 'B')
					die("Invalid shrinkage level type ('A' or 'H' or 'B')");
				
				int level = 0;
				if (levelType == 'B')
					level = levelStr[1] - '0';
				if (levelType == 'B' && (level < 0 || level > 9))
					die("Invalid shrinkage level (only levels 0-9 can be used with level type 'B')");
				
				sManualShrinkageThresholdStr += "-";
				sManualShrinkageThresholdStr += argStr;
				
				int thresholdStartPos = 2;
				if (levelType == 'B')
					thresholdStartPos = 3;
				stringstream ss(argStr.substr(thresholdStartPos, argStr.length() - thresholdStartPos));
				Float_t threshold;
				ss >> threshold;
				
				if (levelType == 'A')	//All levels
				{
					if (colorIdx > -1)
					{
						if (gManualShrinkageThresholds[colorIdx].size() < 11)
							gManualShrinkageThresholds[colorIdx].resize(11);
						for (int i = 0; i <= 9; i++)
							gManualShrinkageThresholds[colorIdx][i] = threshold;
					}
					else
					{
						for (int i = 0; i < 3; i++)
						{
							if (gManualShrinkageThresholds[i].size() < 11)
								gManualShrinkageThresholds[i].resize(11);
							for (int j = 0; j <= 9; j++)
								gManualShrinkageThresholds[i][j] = threshold;
						}
					}
					string colorName = "RGB";
					if (colorIdx == 0)
						colorName = "red";
					else if (colorIdx == 1)
						colorName = "green";
					else if (colorIdx == 2)
						colorName = "blue";
					cout << "All levels threshold argument read for color " << colorName
						<< ": " << threshold << endl;
				}
				else if (levelType == 'H')
				{
					if (colorIdx > -1)
					{
						if (gManualShrinkageThresholds[colorIdx].size() < 1)
							gManualShrinkageThresholds[colorIdx].resize(1);
						gManualShrinkageThresholds[colorIdx][0] = threshold;
					}
					else
					{
						for (int i = 0; i < 3; i++)
						{
							if (gManualShrinkageThresholds[i].size() < 1)
								gManualShrinkageThresholds[i].resize(1);
							gManualShrinkageThresholds[i][0] = threshold;
						}
					}
					string colorName = "RGB";
					if (colorIdx == 0)
						colorName = "red";
					else if (colorIdx == 1)
						colorName = "green";
					else if (colorIdx == 2)
						colorName = "blue";
					cout << "H threshold argument read for color " << colorName
						<< ": " << threshold << endl;
				}
				else
				{
					assert(levelType == 'B');
					
					if (colorIdx > -1)
					{
						if (gManualShrinkageThresholds[colorIdx].size() < level + 2)
							gManualShrinkageThresholds[colorIdx].resize(level + 2);
						gManualShrinkageThresholds[colorIdx][level + 1] = threshold;
					}
					else
					{
						for (int i = 0; i < 3; i++)
						{
							if (gManualShrinkageThresholds[i].size() < level + 2)
								gManualShrinkageThresholds[i].resize(level + 2);
							gManualShrinkageThresholds[i][level + 1] = threshold;
						}
					}
					string colorName = "RGB";
					if (colorIdx == 0)
						colorName = "red";
					else if (colorIdx == 1)
						colorName = "green";
					else if (colorIdx == 2)
						colorName = "blue";
					cout << "B threshold argument read for color " << colorName
						<< ", level " << level << ": " << threshold << endl;
				}
			}
			else
			{
				string filename(argv[argIdx]);
				sInputFilenames.push_back(filename);
				cout << "Input filename arg read: " << filename << endl;
			}
		}
	}
}

bool processArgsFromFile(int argc, char **argv)
{
	if (argc < 2)
		return false;
	
	string argListFile("argList.txt");
	
	ifstream ifs;
	ifs.open(argListFile.c_str());
	if (!ifs)
		return false;
	
	cout << "Building argument list from file " << argListFile << endl;
	
	vector<string> args;
	args.push_back(string(argv[0]));
	while (ifs)
	{
		string s;
		while (ifs && !isspace(ifs.peek()))
		{
			char c = ifs.get();
			if (ifs)
				s += c;
		}
		
		//cout << "    Arg read: \"" << s << "\"" << endl;
		
		args.push_back(s);
		
		//Skip spaces/endlines/tabs/etc. to the next argument
		while (ifs && isspace(ifs.peek()))
			ifs.get();
	}
	ifs.close();
	
	//cout << "Num args in file: " << args.size() << endl;
	//for (int i = 0; i < args.size(); i++)
	//	cout << "    Arg: \"" << args[i] << "\"" << endl;
	
	int argc2 = args.size();
	char **argv2 = new char*[argc2];
	for (int i = 0; i < argc2; i++)
	{
		argv2[i] = new char[256];
		strcpy(argv2[i], args[i].c_str());
		//cout << "    Arg " << i << " length: " << strlen(argv2[i]) << endl;
	}
	
	//cout << "Argument list built" << endl;
	//for (int i = 0; i < argc2; i++)
	//	cout << "    Arg: \"" << argv2[i] << "\"" << endl;
	
	processArgs(argc2, argv2);
	
	return true;
}

void writeSeparator(char sepChar, int rowLength, int numRows)
{
	cout << endl << endl;
	for (int r = 0; r < numRows; r++)
	{
		for (int c = 0; c < rowLength; c++)
			cout << sepChar;
		cout << endl;
	}
	cout << endl;
}

void readNoiseParamsFile()
{
	ifstream ifs;
	ifs.open(sNoiseParamsFilename.c_str());
	if (!ifs)
		die(string("Noise params file not found: " + sNoiseParamsFilename).c_str());
	
	for (int paramIdx = 0; paramIdx < 19; paramIdx++)
	{
		while (true)
		{
			while (ifs && isspace(ifs.peek()))	//Skip any whitespace (blank lines, leading spaces, etc).
				ifs.get();
			if (ifs && ifs.peek() != '#')	//Next line isn't a comment, so proceed with reading data from it
				break;
			//Next line is a comment, so skip the entire line
			while (ifs && ifs.get() != '\n');
		}
		
		if (ifs)
		{
			switch (paramIdx)
			{
				case 0:	ifs >> sNoiseParams.label_;	break;	//Note, this label string cannot contain whitespace
				case 1:	ifs >> sNoiseParams.lux_;	break;
				case 2:	ifs >> sNoiseParams.minReflectance_;	break;
				case 3:	ifs >> sNoiseParams.maxReflectance_;	break;
				case 4:	ifs >> sNoiseParams.temperatureC_;	break;
				case 5:	ifs >> sNoiseParams.pixelPitchMicrons_;	break;
				case 6:	ifs >> sNoiseParams.fullWellCapacity_;	break;
				case 7:	ifs >> sNoiseParams.readNoise_;	break;
				case 8:	ifs >> sNoiseParams.ampVarFactor_;	break;
				case 9:	ifs >> sNoiseParams.darkCurrentRoomTempElectronsPerPixPerSec_;	break;
				case 10:	ifs >> sNoiseParams.hotPixelFactor_;	break;
				case 11:	ifs >> sHotpixelMapFilename;	break;	//Note, this filename string cannot contain whitespace
				case 12:	ifs >> sNoiseParams.quantumEfficiency_[R];	break;
				case 13:	ifs >> sNoiseParams.quantumEfficiency_[G];	break;
				case 14:	ifs >> sNoiseParams.quantumEfficiency_[B];	break;
				case 15:	ifs >> sNoiseParams.ampGlowElectronsPerPixPerSec_;	break;
				case 16:	ifs >> sNoiseParams.exposureSecs_;	break;
				case 17:	ifs >> sNoiseParams.ampGain_;	break;
				case 18:	ifs >> sNoiseParams.focalRatio_;	break;
			}
		}
		
		if (!ifs)
		{
			stringstream ss;
			ss << "Invalid noise params file " << sNoiseParamsFilename << " at row " << paramIdx << " (0 indexed)";
			die(ss.str().c_str());
		}
		
		while (ifs && ifs.get() != '\n');	//The rest of the line can be a comment
	}
	
	ifs.close();
	
	sNoiseParams.selfCheck();
	
	//Read the hot pixel map
	if (sHotpixelMapFilename != "" && sHotpixelMapFilename != "NONE")
	{
		ImageReal *hotpixelMapImg = new ImageReal();
		hotpixelMapImg->readFile(sHotpixelMapFilename);
		sNoiseParams.hotPixelScalars_ = hotpixelMapImg;
	}
	
	//cout << "Noise params read:" << endl;
	//sNoiseParams.dump();
}

//======================================================================
#pragma mark -

void programBehaviorNoisify()
{
	readNoiseParamsFile();
	
	if (sInputFilenames.size() < 1)
		die("Supply at least one image as argument:\n  -- PGM - original gray image\n  -- PPM - original RGB image\n");
	
	for (int imageIdx = 0; imageIdx < sInputFilenames.size(); imageIdx++)
	{
		string fileExtension = sInputFilenames[imageIdx].substr(sInputFilenames[imageIdx].length() - 3, 3);
		myAssert(fileExtension == "pgm" || fileExtension == "ppm");
		
		if (fileExtension == "pgm")	//Noisify a grayscale image
		{
			ImageReal grayImg;
			grayImg.readFile(sInputFilenames[imageIdx]);
			
			if (sNoiseParams.ampGain_ >= 0)	//No auto-ISO
				addSensorNoise(grayImg, false, sNoiseParams);
			else addSensorNoiseWithAutoISO(grayImg, false, sNoiseParams);	//Auto-ISO
			
			cout << endl;
			grayImg.writeFile(1);
			//Write 16-bit file --- grayImg.writeFile(2);
		}
		else if (fileExtension == "ppm")	//Noisify a color image
		{
			Image3ChnlReal rgbImg;
			rgbImg.readFile(sInputFilenames[imageIdx]);
			
			if (!sBayerize)
			{
				if (sNoiseParams.ampGain_ >= 0)	//No auto-ISO
					addSensorNoise(rgbImg, false, sNoiseParams);
				else addSensorNoiseWithAutoISO(rgbImg, false, sNoiseParams);	//Auto-ISO
				
				cout << endl;
				rgbImg.writeFile(1);
				//Write 16-bit file --- rgbImg.writeFile(2);
			}
			else	//Bayerize
			{
				//Bayerize the input image (convert from rgb to gray)
				ImageReal grayImg;
				grayImg.resize(rgbImg.width(), rgbImg.height());
				bayer(rgbImg, grayImg);
				
				if (gGenerateIntermediateImages)
				{
					grayImg.writeFile(1);
					//Write 16-bit file --- grayImg.writeFile(2);
				}
				
				if (gGenerateIntermediateImages)
				{
					//Immediately deBayer
					//Strictly speaking, this isn't an "intermediate" image (ala gGenerateIntermediateImages),
					//but it might be interesting to see.
					deBayer(grayImg, rgbImg);
					rgbImg.writeFile(1);
				}
				
				//Denoise the gray Bayer image
				if (sNoiseParams.ampGain_ >= 0)	//No auto-ISO
					addSensorNoise(grayImg, true, sNoiseParams);
				else addSensorNoiseWithAutoISO(grayImg, true, sNoiseParams);	//Auto-ISO
				
				if (gGenerateIntermediateImages)
				{
					//Strictly speaking, this isn't an "intermediate" image (ala gGenerateIntermediateImages),
					//but it might be interesting to see.
					cout << endl;
					grayImg.writeFile(1);
					//Write 16-bit file --- grayImg.writeFile(2);
				}
				
				//Debayerize back to color
				deBayer(grayImg, rgbImg);
				rgbImg.writeFile(1);
				//Write 16-bit file --- rgbImg.writeFile(2);
			}
		}
	}
}

void programBehaviorDenoise()
{
	if (sInputFilenames.size() < 1)
		die("Supply at least one image as argument:\n  -- PGM - original gray image\n  -- PPM - original RGB image\n");
	
	cout << "Denoising with " << sManualShrinkageThresholdStr << endl;
	for (int imageIdx = 0; imageIdx < sInputFilenames.size(); imageIdx++)
	{
		string fileExtension = sInputFilenames[imageIdx].substr(sInputFilenames[imageIdx].length() - 3, 3);
		myAssert(fileExtension == "pgm" || fileExtension == "ppm");
		
		if (fileExtension == "pgm")
			denoiseGrayImage(sInputFilenames[imageIdx], sManualShrinkageThresholdStr);
		else if (fileExtension == "ppm")
			denoiseRGBImage(sInputFilenames[imageIdx], sManualShrinkageThresholdStr);
	}
}

void programBehavior1ChnlImagePairStats()
{
	ImageReal grayImg1;
	grayImg1.readFile(sInputFilenames[0]);
	int width = grayImg1.width(), height = grayImg1.height();
	
	//The bit depth should be 12 for RAW files from cameras,
	//but 16 for images created (by scaling by 256) from 8 bit images.
	ImageReal grayImg2(width, height);
	vector<vector<Float_t> > stats(sInputFilenames.size());
	Float_t mse, psnr;
	for (int imageIdx = 1; imageIdx < sInputFilenames.size(); imageIdx++)
	{
		stats[imageIdx].resize(3);
		grayImg2.readFile(sInputFilenames[imageIdx]);
		
		cout << endl;
		
		calc1ChnlMSE_PSNR(grayImg1, grayImg2, mse, psnr);
		stats[imageIdx][0] = mse;
		stats[imageIdx][1] = sqrt(mse);
		stats[imageIdx][2] = psnr;
	}
	
	for (int i = 0; i < 3; i++)
	{
		switch (i)
		{
			case 0:  cout << "MSE:              "; break;
			case 1:  cout << "RMSE:             "; break;
			case 2:  cout << "PSNR (dB):        "; break;
		}
		
		for (int j = 1; j < sInputFilenames.size(); j++)
			cout << setw(12) << stats[j][i];
		cout << endl;
	}
}

void programBehaviorRGBImagePairStats()
{
	//Keep track of the first image, it is the baseline for stat calculations
	Image3ChnlReal rgbImg1;
	rgbImg1.readFile(sInputFilenames[0]);
	int width = rgbImg1.width(), height = rgbImg1.height();
	int totalPixels = width * height;
	
	//RGB to Lab the first image
	Image3ChnlReal labImg1(width, height);
	rgbImgToLabImg(1, rgbImg1, labImg1);
	
	//Normalize Lab, this is for Lab PSNR.
	//Instead of normalizing and getting Lab PSNR,
	//it would be better to get Lab MSE instead.
	Image3ChnlReal labImg1Norm = labImg1;
	Float_t minVal = (*labImg1Norm[0])[0];
	for (int c = 0; c < 3; c++)
		for (int i = 0; i < totalPixels; i++)
			if ((*labImg1Norm[c])[i] < minVal)
				minVal = (*labImg1Norm[c])[i];
	if (minVal < 0)
		for (int c = 0; c < 3; c++)
			for (int i = 0; i < totalPixels; i++)
				(*labImg1Norm[c])[i] -= minVal;
	labImg1Norm.normalize(1.0, true, false);
	
	//Calculate the stats for each image
	Image3ChnlReal labImg2Norm(width, height);
	vector<ImageStats> imageStats(sInputFilenames.size() - 1);
	Float_t mse;
	for (int imageIdx = 1; imageIdx < sInputFilenames.size(); imageIdx++)
	{
		imageStats[imageIdx - 1].rgbImg.readFile(sInputFilenames[imageIdx]);
		
		calcRGB_MSE_PSNR(rgbImg1, imageStats[imageIdx - 1].rgbImg, mse, imageStats[imageIdx - 1].rgbPsnr);
		cout << "Image Index " << imageIdx
			<< "    mse: " << mse
			<< "    psnr: " << imageStats[imageIdx - 1].rgbPsnr;
		
		calc1ChnlMSE_PSNR(*rgbImg1[0], *imageStats[imageIdx - 1].rgbImg[0], mse, imageStats[imageIdx - 1].rPsnr);
		
		calc1ChnlMSE_PSNR(*rgbImg1[1], *imageStats[imageIdx - 1].rgbImg[1], mse, imageStats[imageIdx - 1].gPsnr);
		
		calc1ChnlMSE_PSNR(*rgbImg1[2], *imageStats[imageIdx - 1].rgbImg[2], mse, imageStats[imageIdx - 1].bPsnr);
		
		rgbImgToLabImg(1, imageStats[imageIdx - 1].rgbImg, imageStats[imageIdx - 1].labImg);
		labImg2Norm = imageStats[imageIdx - 1].labImg;
		minVal = (*labImg2Norm[0])[0];
		for (int c = 0; c < 3; c++)
			for (int i = 0; i < totalPixels; i++)
				if ((*labImg2Norm[c])[i] < minVal)
					minVal = (*labImg2Norm[c])[i];
		if (minVal < 0)
			for (int c = 0; c < 3; c++)
				for (int i = 0; i < totalPixels; i++)
					(*labImg2Norm[c])[i] -= minVal;
		labImg2Norm.normalize(1.0, true, false);
		
		//LAB PSNR isn't a justified stat in the literature.  I just made it up for my own curiosity.
		//calcRGB_MSE_PSNR(labImg1Norm, labImg2Norm, mse, imageStats[imageIdx - 1].labPsnr);
		//cout << "    labPsnr: " << imageStats[imageIdx - 1].labPsnr;
		
		imageStats[imageIdx - 1].labDeltaE = calcImgLabDeltaE(labImg1, imageStats[imageIdx - 1].labImg);
		cout << "    labDeltaE: " << imageStats[imageIdx - 1].labDeltaE;
		cout << endl;
	}
	
	//Mark the min and max stats
	Float_t minV, maxV;
	int minIdx, maxIdx;
	
	minV = maxV = imageStats[0].rPsnr;
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].rPsnr < minV)
			minV = imageStats[i].rPsnr;
		if (imageStats[i].rPsnr > maxV)
			maxV = imageStats[i].rPsnr;
	}
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].rPsnr == minV)
			imageStats[i].minRPsnr = true;
		if (imageStats[i].rPsnr == maxV)
			imageStats[i].maxRPsnr = true;
	}
	
	minV = maxV = imageStats[0].gPsnr;
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].gPsnr < minV)
			minV = imageStats[i].gPsnr;
		if (imageStats[i].gPsnr > maxV)
			maxV = imageStats[i].gPsnr;
	}
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].gPsnr == minV)
			imageStats[i].minGPsnr = true;
		if (imageStats[i].gPsnr == maxV)
			imageStats[i].maxGPsnr = true;
	}
	
	minV = maxV = imageStats[0].bPsnr;
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].bPsnr < minV)
			minV = imageStats[i].bPsnr;
		if (imageStats[i].bPsnr > maxV)
			maxV = imageStats[i].bPsnr;
	}
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].bPsnr == minV)
			imageStats[i].minBPsnr = true;
		if (imageStats[i].bPsnr == maxV)
			imageStats[i].maxBPsnr = true;
	}
	
	minV = maxV = imageStats[0].rgbPsnr;
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].rgbPsnr < minV)
			minV = imageStats[i].rgbPsnr;
		if (imageStats[i].rgbPsnr > maxV)
			maxV = imageStats[i].rgbPsnr;
	}
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].rgbPsnr == minV)
			imageStats[i].minRGBPsnr = true;
		if (imageStats[i].rgbPsnr == maxV)
			imageStats[i].maxRGBPsnr = true;
	}
	
	minV = maxV = imageStats[0].labPsnr;
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].labPsnr < minV)
			minV = imageStats[i].labPsnr;
		if (imageStats[i].labPsnr > maxV)
			maxV = imageStats[i].labPsnr;
	}
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].labPsnr == minV)
			imageStats[i].minLabPsnr = true;
		if (imageStats[i].labPsnr == maxV)
			imageStats[i].maxLabPsnr = true;
	}
	
	minV = maxV = imageStats[0].labDeltaE;
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].labDeltaE < minV)
			minV = imageStats[i].labDeltaE;
		if (imageStats[i].labDeltaE > maxV)
			maxV = imageStats[i].labDeltaE;
	}
	for (int i = 0; i < imageStats.size(); i++)
	{
		if (imageStats[i].labDeltaE == minV)
			imageStats[i].minLabDeltaE = true;
		if (imageStats[i].labDeltaE == maxV)
			imageStats[i].maxLabDeltaE = true;
	}
	
	//Write the file names
	cout << endl << "Filename order:" << endl;
	for (int imageIdx = 1; imageIdx < sInputFilenames.size(); imageIdx++)
		cout << (char)('A' + (imageIdx - 1)) << ": " << sInputFilenames[imageIdx] << endl;
	
	//Write the results
	cout << endl << "Results:" << endl;
	cout << "                      ";
	for (int imageIdx = 1; imageIdx < sInputFilenames.size(); imageIdx++)
		cout << setw(12) << (char)('A' + (imageIdx - 1));
	cout << endl;
	
	for (int i = 0; i < 6; i++)
	{
		//LAB PSNR isn't a justified stat in the literature.  I just made it up for my own curiosity.
		if (i == 4)
			continue;
		
		switch (i)
		{
			case 0: cout << "Red PSNR (dB):        "; break;
			case 1: cout << "Green PSNR (dB):      "; break;
			case 2: cout << "Blue PSNR (dB):       "; break;
			case 3: cout << "RGB PSNR (dB):        "; break;
			case 4: cout << "Lab PSNR (dB):        "; break;
			case 5: cout << "Lab Delta E:          "; break;
		}
		
		//Write the stats
		for (int j = 0; j < imageStats.size(); j++)
		{
			switch (i)
			{
				case 0:
					cout << setw(12) << imageStats[j].rPsnr;
					break;
				case 1:
					cout << setw(12) << imageStats[j].gPsnr;
					break;
				case 2:
					cout << setw(12) << imageStats[j].bPsnr;
					break;
				case 3:
					cout << setw(12) << imageStats[j].rgbPsnr;
					break;
				case 4:
					cout << setw(12) << imageStats[j].labPsnr;
					break;
				case 5:
					cout << setw(12) << imageStats[j].labDeltaE;
					break;
			}
		}
		cout << endl;
		
		//Mark the best and worst
		cout << "                      ";
		for (int j = 0; j < imageStats.size(); j++)
		{
			bool markWorst = false, markBest = false;
			switch (i)
			{
				case 0:
					if (imageStats[j].minRPsnr)
						markWorst = true;
					if (imageStats[j].maxRPsnr)
						markBest = true;
					break;
				case 1:
					if (imageStats[j].minGPsnr)
						markWorst = true;
					if (imageStats[j].maxGPsnr)
						markBest = true;
					break;
				case 2:
					if (imageStats[j].minBPsnr)
						markWorst = true;
					if (imageStats[j].maxBPsnr)
						markBest = true;
					break;
				case 3:
					if (imageStats[j].minRGBPsnr)
						markWorst = true;
					if (imageStats[j].maxRGBPsnr)
						markBest = true;
					break;
				case 4:
					if (imageStats[j].minLabPsnr)
						markWorst = true;
					if (imageStats[j].maxLabPsnr)
						markBest = true;
					break;
				case 5:
					if (imageStats[j].maxLabDeltaE)
						markWorst = true;
					if (imageStats[j].minLabDeltaE)
						markBest = true;
					break;
			}
			if (markWorst && markBest)
				markWorst = false;
			
			if (markWorst)
				cout << setw(12) << "-- WORST --";
			else if (markBest)
				cout << setw(12) << "=== BEST ===";
			else cout << setw(12) << "";
		}
		cout << endl;
	}
	
	//Write the results sorted
	cout << endl << "Name sorted:" << endl;
	cout << setw(15) << "RGB PSNR" << setw(15) << "Lab Delta E" << endl;
	for (int i = 0; i < imageStats.size(); i++)
		cout << setw(15) << imageStats[i].rgbPsnr << setw(15) << imageStats[i].labDeltaE << "    " << imageStats[i].rgbImg.filename() << endl;
	
	cout << endl << "Name sorted, tab delimited (for easy copy/paste to spread sheet):" << endl;
	cout << "RGB PSNR" << '\t' << "Lab Delta E" << endl;
	for (int i = 0; i < imageStats.size(); i++)
		cout << imageStats[i].rgbPsnr << '\t' << imageStats[i].labDeltaE << endl;
	
	cout << endl << "RGB Psnr sorted:" << endl;
	cout << setw(15) << "RGB PSNR" << setw(15) << "Lab Delta E" << endl;
	sort(imageStats.begin(), imageStats.end(), ImageStats::RGBPsnrSorter());
	for (int i = 0; i < imageStats.size(); i++)
		cout << setw(15) << imageStats[i].rgbPsnr << setw(15) << imageStats[i].labDeltaE << "    " << imageStats[i].rgbImg.filename() << endl;
	
	cout << endl << "RGB Psnr sorted, tab delimited (for easy copy/paste to spread sheet):" << endl;
	cout << "RGB PSNR" << '\t' << "Lab Delta E" << endl;
	for (int i = 0; i < imageStats.size(); i++)
		cout << imageStats[i].rgbPsnr << '\t' << imageStats[i].labDeltaE << endl;
	
	cout << endl << "Lab Delta E sorted:" << endl;
	cout << setw(15) << "RGB PSNR" << setw(15) << "Lab Delta E" << endl;
	sort(imageStats.begin(), imageStats.end(), ImageStats::LabDeltaESorter());
	for (int i = 0; i < imageStats.size(); i++)
		cout << setw(15) << imageStats[i].rgbPsnr << setw(15) << imageStats[i].labDeltaE << "    " << imageStats[i].rgbImg.filename() << endl;
	
	cout << endl << "Lab Delta E sorted, tab delimited (for easy copy/paste to spread sheet):" << endl;
	cout << "RGB PSNR" << '\t' << "Lab Delta E" << endl;
	for (int i = 0; i < imageStats.size(); i++)
		cout << imageStats[i].rgbPsnr << '\t' << imageStats[i].labDeltaE << endl;
}

//======================================================================
#pragma mark -

void testRGBtoLABtoRGB()
{
	//PPM
	for (int imageIdx = 0; imageIdx < sInputFilenames.size(); imageIdx++)
	{
		//Convert RGB to LAB
		Image3ChnlReal rgbImg;
		rgbImg.readFile(sInputFilenames[imageIdx]);
		int width = rgbImg.width(), height = rgbImg.height();
		int halfWidth = width / 2, halfHeight = height / 2;
		int totalPixels = rgbImg.totalPixels();
		
		Image3ChnlReal labImg;
		rgbImgToLabImg(1, rgbImg, labImg);
		
		Float_t minValL, maxValL, minValA, maxValA, minValB, maxValB;
		
		minValL = maxValL = (*labImg[0])[0];
		for (int i = 0; i < totalPixels; i++)
		{
			if ((*labImg[0])[i] > maxValL)
				maxValL = (*labImg[0])[i];
			if ((*labImg[0])[i] < minValL)
				minValL = (*labImg[0])[i];
		}
		labImg[0]->normalize(1, true);
		
		minValA = maxValA = (*labImg[1])[0];
		for (int i = 0; i < totalPixels; i++)
		{
			if ((*labImg[1])[i] > maxValA)
				maxValA = (*labImg[1])[i];
			if ((*labImg[1])[i] < minValA)
				minValA = (*labImg[1])[i];
		}
		if (minValA < 0)
			for (int i = 0; i < totalPixels; i++)
				(*labImg[1])[i] -= minValA;
		labImg[1]->normalize(1, true);
		
		minValB = maxValB = (*labImg[2])[0];
		for (int i = 0; i < totalPixels; i++)
		{
			if ((*labImg[2])[i] > maxValB)
				maxValB = (*labImg[2])[i];
			if ((*labImg[2])[i] < minValB)
				minValB = (*labImg[2])[i];
		}
		if (minValB < 0)
			for (int i = 0; i < totalPixels; i++)
				(*labImg[2])[i] -= minValB;
		labImg[2]->normalize(1, true);
		
		labImg.setFilename(rgbImg.filename());
		labImg.appendFilename("_lab");
		labImg.writeFile(1);
		
		//========================================
		/*
		//Denoise rgb and lab images
		Image3ChnlReal rgbImgDenoised, labImgDenoised;
		
		//еее
		//Need to reuse optimal denoising routines here
		*/
		//========================================
		
		//Convert LAB to RGB
		
		labImg[0]->scale(maxValL);
		labImg[1]->scale(maxValA - minValA);
		if (minValA < 0)
			for (int i = 0; i < totalPixels; i++)
				(*labImg[1])[i] += minValA;
		labImg[2]->scale(maxValB - minValB);
		if (minValB < 0)
			for (int i = 0; i < totalPixels; i++)
				(*labImg[2])[i] += minValB;
		
		Image3ChnlReal rgbImgSynth;
		labImgToRgbImg(1, labImg, rgbImgSynth);
		
		rgbImgSynth.setFilename(rgbImg.filename());
		rgbImgSynth.appendFilename("_rgb");
		rgbImgSynth.writeFile(1);
		
		Float_t mse, psnr;
		calcRGB_MSE_PSNR(rgbImg, rgbImgSynth, mse, psnr);
		cout << "MSE/PSNR: " << mse << setw(14) << psnr << endl;
		calcRGB_MSE_PSNR(rgbImg, rgbImgSynth, mse, psnr);
		cout << "MSE/PSNR: " << mse << setw(14) << psnr << endl;
	}
}

void testPGMFourierAnalysis()
{	
	ImageReal firstImagePowerSpectrum;
	
	//PGM
	for (int imageIdx = 0; imageIdx < sInputFilenames.size(); imageIdx++)
	{
		//Open pgm
		ImageReal grayImg;
		grayImg.readFile(sInputFilenames[imageIdx]);
		int width = grayImg.width(), height = grayImg.height();
		int halfWidth = width / 2, halfHeight = height / 2;
		
		ImageReal fourierTransform;
		bool err = doFourierAnalysis(grayImg, fourierTransform);
		
		//Make an image of the power spectrum
		ImageReal powerSpectrum;
		generatePowerSpectrum(fourierTransform, powerSpectrum, GT_LOG10, true);
		
		powerSpectrum.setFilename(grayImg.filename());
		powerSpectrum.appendFilename("_ps");
		powerSpectrum.writeFile(1);
		
		bool shiftCenter = false;
		/*
		//Make a 1D histogram of the power spectrum, binned by frequency
		myAssert(!shiftCenter);
		generatePowerSpectrum(fourierTransform, powerSpectrum, GT_NONE, shiftCenter);
		if (imageIdx == 0)
			firstImagePowerSpectrum = powerSpectrum;
		
		int numBins = sqrt(halfWidth * halfWidth + halfHeight * halfHeight);
		cout << "Num bins: " << numBins << endl;
		int numBinsScalar = 4;
		numBins /= numBinsScalar;
		numBins++;
		cout << "Num bins scaled: " << numBins << endl;
		
		Float_t *powerSpectrumHistrogram = new Float_t[numBins];
		memset(powerSpectrumHistrogram, 0, numBins * sizeof(Float_t));
		
		int *binCounts = new int[numBins];
		memset(binCounts, 0, numBins * sizeof(int));
		
		//This code only works if the power spectrum has NOT been center-shifted.
		for (int y = 0; y < height; y++)
		{
			int yCenter = (y <= height / 2) ? 0 : height;
			Float_t yDistSqr = pow((Float_t)y - (Float_t)yCenter, 2);
			for (int x = 0; x <= halfWidth; x++)	//Step one column past the center (Nyquist)
			{
				if (y == 0 && x == 0)
					continue;	//Skip the DC term
				
				Float_t xDistSqr = pow((Float_t)x, 2);
				int dcTermDistance = roundToNearestInt(sqrt(yDistSqr + xDistSqr), 0, 0);
				dcTermDistance--;	//Since we skipped the DC term
				dcTermDistance /= numBinsScalar;
				myAssert(dcTermDistance >= 0 && dcTermDistance < numBins);
				
				powerSpectrumHistrogram[dcTermDistance] += powerSpectrum(x, y);
				binCounts[dcTermDistance]++;
			}
		}
		
		//Divide each bin by its number of contributing Fourier coefficients.
		//The result will be a histogram in which the value in each bin
		//represents the mean power at the frequency corresponding to that bin,
		//i.e., the mean power in a very thin annulus surrounding the DC term
		//in the power spectrum.
		int maxBinCount = 0;
		for (int i = 0; i < numBins; i++)
			if (binCounts[i] > maxBinCount)
				maxBinCount = binCounts[i];
		for (int i = 0; i < numBins; i++)
			powerSpectrumHistrogram[i] /= (Float_t)maxBinCount;
		
		string filename = sInputFilenames[imageIdx] + "_hist.txt";
		ofstream ofs;
		ofs.open(filename.c_str());
		myAssert(ofs);
		
		for (int i = 0; i < numBins; i++)
			ofs << i << '\t' << powerSpectrumHistrogram[i] << endl;
		
		ofs.close();
		
		delete [] powerSpectrumHistrogram;
		*/
		//Diff the power spectra
		shiftCenter = true;
		generatePowerSpectrum(fourierTransform, powerSpectrum, GT_NONE, shiftCenter);
		if (imageIdx == 0)
			firstImagePowerSpectrum = powerSpectrum;
		else
		{
			ImageReal powerSpectrumDiff;
			diffImages(firstImagePowerSpectrum, powerSpectrum, powerSpectrumDiff);
			powerSpectrumDiff.writeFile(1);
		}
	}
}

void testFourierConvolution(int method)
{
	if (method == 0)
	{
		if (sInputFilenames.size() < 2)
			die("Please enter two PGM files to convolve");
		
		ImageReal img1, img2;
		img1.readFile(sInputFilenames[0]);
		img2.readFile(sInputFilenames[1]);
		
		ImageReal convolvedImgReal;
		convolveImages(img1, img2, convolvedImgReal);
		
		convolvedImgReal.normalize(1, false);
		
		convolvedImgReal.setFilename(img1.filename());
		convolvedImgReal.appendFilename("_conv");
		convolvedImgReal.writeFile(1);
	}
	else if (method == 1)
	{
		if (sInputFilenames.size() < 1)
			die("Please enter one PGM file to convolve with the steerable pyramid filters");
		
		ImageReal img1;
		img1.readFile(sInputFilenames[0]);
		
		SteerablePyramid steerPyr;
		steerPyr.DoFullWaveletAnalysis(img1.img(), img1.width(), img1.height(), false);
		
		ImageReal img1SynthReal(img1.width(), img1.height());
		steerPyr.DoFullWaveletSynthesis(img1SynthReal.img());
		
		Float_t mse, psnr;
		calc1ChnlMSE_PSNR(img1, img1SynthReal, mse, psnr);
		cout << "MSE/PSNR: " << mse << setw(14) << psnr << endl;
		
		//---------------------
		
		steerPyr.DoFullWaveletAnalysis(img1.img(), img1.width(), img1.height(), true);
		
		steerPyr.DoFullWaveletSynthesis(img1SynthReal.img());
		
		calc1ChnlMSE_PSNR(img1, img1SynthReal, mse, psnr);
		cout << "MSE/PSNR: " << mse << setw(14) << psnr << endl;
	}
}

//======================================================================
#pragma mark -

int main(int argc, char* argv[])
{
	//Sorta kinda verify the enum labels
	assert(skProgramBehaviorLabels[PB_NUM_PROGRAM_BEHAVIORS] == "PB_NUM_PROGRAM_BEHAVIORS");
	assert(skIdealShrinkageMethodLabels[NUM_IDEAL_SHRINKAGE_METHODS] == "NUM_IDEAL_SHRINKAGE_METHODS");
	
	//setiosflags(ios::showpoint | ios::fixed);
	cout.setf(ios::fixed);
	cout.precision(6);	//std out default I believe, therefore a pointless, but illustrative, command
	
	//Uncomment this line to get a random seed.
	//Pass a seed into SeedRandom to set the seed explicitly.
	//SeedRandom();
	
	if (!processArgsFromFile(argc, argv))
		processArgs(argc, argv);
	
	//Clearly mark a new run in the run log or console window
	writeSeparator('#', 100, 4);
	
	switch (sProgramBehavior)
	{
		case PB_TEST_FOURIER:
			fourierTest();  //FourierTransform.cpp
			break;
		case PB_TEST_FOURIER_CONVOLUTION:
			testFourierConvolution(0);
			break;
		case PB_TEST_RGB_XYZ_LAB:
			testRGBXYZLABRoutines();	//ImageProcs.cpp
			break;
		//========================================
		case PB_NOISIFY:
			programBehaviorNoisify();
			break;
		case PB_DENOISE:
			programBehaviorDenoise();
			break;
		case PB_GRAY_IMG_PAIR_STATS:
			programBehavior1ChnlImagePairStats();
			break;
		case PB_RGB_IMG_PAIR_STATS:
			programBehaviorRGBImagePairStats();
			break;
		//========================================
		case PB_NONE:
		default:
			die("Please specify a program behavior flag.");
	}
	
	cout << endl << "Program is complete" << endl;
	
	return 0;
}

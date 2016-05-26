#include "Bayer.h"
#include <assert.h>

//******************************************************************************
//Extern Globals

//******************************************************************************
//Global Declarations
DebayerMethod gDeBayerMethod = DB_ALTERNATING_PROJECTIONS;

Color kRGGB[4] = { R, G, G, B };
Color kBGGR[4] = { B, G, G, R };
Color kGRBG[4] = { G, R, B, G };
Color kGBRG[4] = { G, B, R, G };

Color *gBayerPattern = kRGGB;

//******************************************************************************
//Function Prototypes
	//	In this file
	//	In external files

//Determine if a pixel is red, green, or blue in the Bayer mosaic
Color determineBayerColor(int x, int y)
{
	if (x < 0)	//A flag for nonBayer images
		return R;
	
	int x2 = x % 2, y2 = y % 2;
	
	Color color = gBayerPattern[y2 * 2 + x2];
	
	return color;
}

//Determine which position a pixel is in within its 2x2 block
BlockCorner determineBlockCorner(int x, int y)
{
	int x2 = x % 2, y2 = y % 2;
	
	if (x2 == 0 && y2 == 0)
		return UL;
	if (x2 == 1 && y2 == 0)
		return UR;
	if (x2 == 0 && y2 == 1)
		return LL;
	//if (x2 == 1 && y2 == 1)
		return LR;
}

//Determine which position a color is in within a 2x2 block
void determineBlockCoord(Color color, int &x, int &y)
{
	assert(color != G);
	if (color == R)
	{
		if (gBayerPattern == kRGGB || gBayerPattern == kGRBG)
			y = 0;
		else y = 1;
		
		if (gBayerPattern == kRGGB || gBayerPattern == kGBRG)
			x = 0;
		else x = 1;
	}
	else //if (color == B)
	{
		if (gBayerPattern == kBGGR || gBayerPattern == kGBRG)
			y = 0;
		else y = 1;
		
		if (gBayerPattern == kBGGR || gBayerPattern == kGRBG)
			x = 0;
		else x = 1;
	}
}

//Determine the positions of all three colors (RGB) in the Bayer mosaic.
//In the case of green, the first position will be returned (in the order UL, UR, LL, LR);
void determineBlockColorOrder(int &rx, int &ry, int &g1x, int &g1y, int &bx, int &by,
							int &redArrayPos, int &greenArrayPos, int &blueArrayPos,
							int &secondGreenUnusedPos)
{
	determineBlockCoord(R, rx, ry);
	determineBlockCoord(B, bx, by);
	assert(rx != bx && ry != by);
	
	g1x = 0;
	g1y = 0;
	if ((rx == 0 && ry == 0) || (bx == 0 && by == 0))
		g1x = 1;
	
	if (rx == 0 && ry == 0)
	{
		redArrayPos = 0;
		greenArrayPos = 1;
		blueArrayPos = 3;
		secondGreenUnusedPos = 2;
	}
	else if (rx == 1 && ry == 0)
	{
		redArrayPos = 1;
		greenArrayPos = 0;
		blueArrayPos = 2;
		secondGreenUnusedPos = 3;
	}
	else if (bx == 0 && by == 0)
	{
		redArrayPos = 3;
		greenArrayPos = 1;
		blueArrayPos = 0;
		secondGreenUnusedPos = 2;
	}
	else if (bx == 1 && by == 0)
	{
		redArrayPos = 2;
		greenArrayPos = 0;
		blueArrayPos = 1;
		secondGreenUnusedPos = 3;
	}
	else assert(false);
}

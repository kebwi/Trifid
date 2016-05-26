#ifndef __BAYER__
#define __BAYER__

enum Color { R, G, B };

enum BlockCorner { UL, UR, LL, LR };	//What position is a particular pixel in, relative to its 2x2 block

enum DebayerMethod { DB_NEAREST_NEIGHBOR, DB_BILINEAR, DB_SMOOTH_HUE, DB_ALTERNATING_PROJECTIONS };

extern Color kRGGB[4];
extern Color kBGGR[4];
extern Color kGRBG[4];
extern Color kGBRG[4];

extern DebayerMethod gDeBayerMethod;
extern Color *gBayerPattern;

//Determine the color of a pixel in the Bayer mosaic
Color determineBayerColor(int x, int y);

//Determine which position a pixel is in within its 2x2 block
BlockCorner determineBlockCorner(int x, int y);

//Determine which position a color is in within a 2x2 block
void determineBlockCoord(Color color, int &x, int &y);

//Determine the positions of all three colors (RGB) in the Bayer mosaic.
//In the case of green, the first position will be returned (in the order UL, UR, LL, LR);
void determineBlockColorOrder(int &rx, int &ry, int &g1x, int &g1y, int &bx, int &by,
							int &redArrayPos, int &greenArrayPos, int &blueArrayPos,
							int &secondGreenUnusedPos);

#endif

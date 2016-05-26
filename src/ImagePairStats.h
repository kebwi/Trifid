#ifndef __IMAGE_PAIR_STATS__
#define __IMAGE_PAIR_STATS__

#include "FloatType.h"
#include "ImageReal.h"
#include "Image3ChnlReal.h"

Float_t calcLabDeltaE(const Float_t l1, const Float_t a1, const Float_t b1,
						const Float_t l2, const Float_t a2, const Float_t b2);
Float_t calcImgLabDeltaE(const Image3ChnlReal &img1, const Image3ChnlReal &img2);

Float_t calc1ChnlME(const ImageReal &img1, const ImageReal &img2);
Float_t calc1ChnlMSE(const ImageReal &img1, const ImageReal &img2);

Float_t calcRGB_ME(const Image3ChnlReal &img1, const Image3ChnlReal &img2);
Float_t calcRGB_MSE(const Image3ChnlReal &img1, const Image3ChnlReal &img2);

void calc1ChnlMSE_PSNR(const ImageReal &img1, const ImageReal &img2,
						Float_t &mse, Float_t &psnr);
void calcRGB_MSE_PSNR(const Image3ChnlReal &img1, const Image3ChnlReal &img2,
						Float_t &mse, Float_t &psnr);

Float_t calculate1ChnlPsnrScore(const ImageReal &imgOri, const ImageReal &img);
Float_t calculateRgbPsnrScore(const Image3ChnlReal &imgOri, const Image3ChnlReal &img);
Float_t calculateLabScore(const Image3ChnlReal &imgOri, const Image3ChnlReal &img);

#endif

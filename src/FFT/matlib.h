/* a few routines from a vector/matrix library */

#include "FloatType.h"

void xpose(Float_t *indata, long iRsiz, Float_t *outdata, long oRsiz, long Nrows, long Ncols);
/* not in-place matrix transpose	*/
/* INPUTS */
/* *indata = input data array	*/
/* iRsiz = offset to between rows of input data array	*/
/* oRsiz = offset to between rows of output data array	*/
/* Nrows = number of rows in input data array	*/
/* Ncols = number of columns in input data array	*/
/* OUTPUTS */
/* *outdata = output data array	*/

void cxpose(Float_t *indata, long iRsiz, Float_t *outdata, long oRsiz, long Nrows, long Ncols);
/* not in-place complex matrix transpose	*/
/* INPUTS */
/* *indata = input data array	*/
/* iRsiz = offset to between rows of input data array	*/
/* oRsiz = offset to between rows of output data array	*/
/* Nrows = number of rows in input data array	*/
/* Ncols = number of columns in input data array	*/
/* OUTPUTS */
/* *outdata = output data array	*/

void cvprod(Float_t *a, Float_t *b, Float_t *out, long N);
/* complex vector product, can be in-place */
/* product of complex vector *a times complex vector *b */
/* INPUTS */
/* N vector length */
/* *a complex vector length N complex numbers */
/* *b complex vector length N complex numbers */
/* OUTPUTS */
/* *out complex vector length N */

#ifndef FP_H
#define FP_H

#include <stdbool.h>

#include "params.h"

void fp_set(fp *x, uint64_t y);
void fp_cswap(fp *x, fp *y, bool c);

void fp_enc(fp *x, uint_c const *y); /* encode to Montgomery representation */
void fp_dec(uint_c *x, fp const *y); /* decode from Montgomery representation */

void fp_add2(fp *x, fp const *y);
void fp_sub2(fp *x, fp const *y);
void fp_mul2(fp *x, fp const *y);

void fp_add3(fp *x, fp const *y, fp const *z);
void fp_sub3(fp *x, fp const *y, fp const *z);
void fp_mul3(fp *x, fp const *y, fp const *z);

void fp_sq1(fp *x);
void fp_sq2(fp *x, fp const *y);
void fp_inv(fp *x);
bool fp_issquare(fp *x); /* destroys input! */

void fp_random(fp *x);


//------print version-----//

void fp_mul3_p(fp *x, fp const *y, fp const *z);

void fp_sq1_p(fp *x);
void fp_sq2_p(fp *x, fp const *y);
void fp_inv_p(fp *x);
void fp_mul2_p(fp *x, fp const *y);
static void fp_pow_p(fp *x, uint_c const *e);



#endif

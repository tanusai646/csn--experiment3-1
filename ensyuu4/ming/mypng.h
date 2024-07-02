#ifndef MYPNG_H
#define MYPNG_H

#include <stdio.h>
#include <stdlib.h>
#include "xydata.h"

extern int data_ans_to_png_img(FILE *fp,
			       const double rx, const double ry,
			       const struct xy *g, size_t gsz,
			       const struct xy *b, size_t bsz,
			       const struct xy *data, size_t sz);

extern int data_to_png_img(FILE *fp, const struct xy *data, size_t sz);

#endif

#ifndef GENPNG_H
#define GENPNG_H
#include <stdio.h>
#include <png.h>

#define PNG_WIDTH 400
#define PNG_HIGHT 400

/* 
   イメージデータは png_byte image[PNG_HIGHT][3*PNG_WIDTH] 決め打ちとする．
*/

/* png_byte image[PNG_HIGHT][3*PNG_WIDTH]; */
typedef png_byte image_t[PNG_HIGHT][3*PNG_WIDTH];

int gen_png(FILE *fp, image_t image);

#endif

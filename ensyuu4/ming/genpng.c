#include <stdio.h>
#include <stdlib.h>
#include "genpng.h"

/* 
   #define PNG_WIDTH 400
   #define PNG_HIGHT 400
   イメージデータは png_byte image[PNG_HIGHT][3*PNG_WIDTH] 決め打ちとする．
*/


/* image 画像を，fp に PNG データとして書き込む */

int
gen_png(FILE *fp, image_t image)
{
  png_uint_32 j, width = PNG_WIDTH, height = PNG_HIGHT;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  volatile png_bytepp rows = NULL;
  if (image == NULL) return -1;
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) goto error;
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) goto error;
  if (setjmp(png_jmpbuf(png_ptr))) goto error;
  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
	       PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
	       PNG_FILTER_TYPE_DEFAULT);
  rows = (png_bytepp)png_malloc(png_ptr, sizeof(png_bytep) * height);
  if (rows == NULL) goto error;
  for (j = 0; j < height; j++)
    rows[j] = &(image[j][0]);
  png_set_rows(png_ptr, info_ptr, rows);
  png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  return 0;
 error:
  if (rows)
    {
      png_bytepp tmp = rows;
      rows = NULL;
      png_free(png_ptr, tmp); // 万一 longjmpしても次は rows == NULL
    }
  png_destroy_write_struct(&png_ptr, &info_ptr);
  return -1;
}

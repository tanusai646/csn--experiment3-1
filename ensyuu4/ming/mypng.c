#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "genpng.h"
#include "mypng.h"

/* メインデータ
   RGBそれぞれ渡す Rのみ1点  */
int
data_ans_to_png_img(FILE *fp,
		    const double rx, const double ry,
		    const struct xy *g, size_t gsz,
		    const struct xy *b, size_t bsz,
		    const struct xy *data, size_t sz)
{
  size_t i;
  image_t *img = (image_t *)calloc(sizeof(image_t), sizeof(png_byte));
  if (img == 0)
    return -1;
  // 枠 (円周)
  for (i = 0; i < 1600; i++)
    {
      int ix = (int)(( cos(i*M_PI/800) + 1.0) * 200.0);
      int iy = (int)((-sin(i*M_PI/800) + 1.0) * 200.0);
      ix = (ix < 0) ? 0 : (ix < 400) ? ix : 399;
      iy = (iy < 0) ? 0 : (iy < 400) ? iy : 399;
      (*img)[iy][ix*3+0] = 255;
      (*img)[iy][ix*3+1] = 255;
    }
  // 星雲
  for (i = 0; i < sz; i++)
    {
      int ix = (int)(( data[i].x + 1.0) * 200.0);
      int iy = (int)((-data[i].y + 1.0) * 200.0);
      ix = (ix < 0) ? 0 : (ix < 400) ? ix : 399;
      iy = (iy < 0) ? 0 : (iy < 400) ? iy : 399;
      (*img)[iy][ix*3+2] = 255;
      (*img)[iy][ix*3+1] = 255;
      if (--((*img)[iy][ix*3+0]) == 0)
	(*img)[iy][ix*3+0] = 1;
    }
  // blue
  for (i = 0; i < bsz; i++)
    {
      int ix = (int)(( b[i].x + 1.0) * 200.0);
      int iy = (int)((-b[i].y + 1.0) * 200.0);
      ix = (ix < 0) ? 0 : (ix < 400) ? ix : 399;
      iy = (iy < 0) ? 0 : (iy < 400) ? iy : 399;
      (*img)[iy][ix*3+0] = 0;
      (*img)[iy][ix*3+1] = 0;
      (*img)[iy][ix*3+2] = 255;
    }
  // green
  for (i = 0; i < gsz; i++)
    {
      int ix = (int)(( g[i].x + 1.0) * 200.0);
      int iy = (int)((-g[i].y + 1.0) * 200.0);
      ix = (ix < 0) ? 0 : (ix < 400) ? ix : 399;
      iy = (iy < 0) ? 0 : (iy < 400) ? iy : 399;
      (*img)[iy][ix*3+0] = 0;
      (*img)[iy][ix*3+1] = 255;
      (*img)[iy][ix*3+2] = 0;
    }
  // red
  {
    int ix = (int)(( rx + 1.0) * 200.0);
    int iy = (int)((-ry + 1.0) * 200.0);
    ix = (ix < 0) ? 0 : (ix < 400) ? ix : 399;
    iy = (iy < 0) ? 0 : (iy < 400) ? iy : 399;
    (*img)[iy][ix*3+0] = 255;
    (*img)[iy][ix*3+1] = 0;
    (*img)[iy][ix*3+2] = 0;
  }
  // とりあえず．このまま
  return gen_png(fp, *img);
}

/* データのみ  */
int
data_to_png_img(FILE *fp, const struct xy *data, size_t sz)
{
  size_t i;
  image_t *img = (image_t *)calloc(sizeof(image_t), sizeof(png_byte));
  if (img == 0)
    return -1;
  // 枠 (円周)
  for (i = 0; i < 1600; i++)
    {
      int ix = (int)(( cos(i*M_PI/800) + 1.0) * 200.0);
      int iy = (int)((-sin(i*M_PI/800) + 1.0) * 200.0);
      ix = (ix < 0) ? 0 : (ix < 400) ? ix : 399;
      iy = (iy < 0) ? 0 : (iy < 400) ? iy : 399;
      (*img)[iy][ix*3+0] = 255;
      (*img)[iy][ix*3+1] = 255;
    }
  // 星雲
  for (i = 0; i < sz; i++)
    {
      int ix = (int)(( data[i].x + 1.0) * 200.0);
      int iy = (int)((-data[i].y + 1.0) * 200.0);
      ix = (ix < 0) ? 0 : (ix < 400) ? ix : 399;
      iy = (iy < 0) ? 0 : (iy < 400) ? iy : 399;
      (*img)[iy][ix*3+2] = 255;
      (*img)[iy][ix*3+1] = 255;
      if (--((*img)[iy][ix*3+0]) == 0)
	(*img)[iy][ix*3+0] = 1;
    }
  // とりあえず．このまま
  return gen_png(fp, *img);
}

//runtime 実行時間指定
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// 時間測定のための関数
static double elapsed_time(struct timeval tv[2]) {
  return tv[1].tv_sec - tv[0].tv_sec +
    (1e-6) * (tv[1].tv_usec - tv[0].tv_usec);
}

// （四分）円の面積の近似計算
double getArea(int n)
{
  int i, j, s1;
  double x, y, s = 0.0, d2 = 0.5/n;
#pragma omp parallel for private(j,x,y,s1) reduction(+:s) schedule(runtime)
  for (i = 0; i < n; i++) {
    y = (2 * i + 1) * d2;     // 正方形の中心 y座標
    s1 = 0;
    for (j = 0; j < n; j++) {
      x = (2 * j + 1) * d2;   // 正方形の中心 x座標
      if (x * x + y * y < 1.0) s1++;   // カウント
    }
    s += s1;
  }
  return s / n / n;  
}

int main(int argc, char *argv[])
{
  struct timeval tv[2];
  gettimeofday(tv + 0, 0);    // 開始時刻記録
  double a = getArea(100000); // 計算本体
  gettimeofday(tv + 1, 0);    // 終了時刻記録
  printf("area * 4.0 = %.13f\n", a * 4);
  printf("time: %f\n", elapsed_time(tv));
}


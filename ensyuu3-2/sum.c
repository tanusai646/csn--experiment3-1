#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

// 時間測定のための関数
static double elapsed_time(struct timeval tv[2]) {
  return tv[1].tv_sec - tv[0].tv_sec +
    (1e-6) * (tv[1].tv_usec - tv[0].tv_usec);
}

// n までの総和
long long sumUp(long long n)
{
  long long i, s = 0;
#pragma omp parallel for reduction(+:s) schedule(runtime)
  for (i = 1; i <= n; i++) {
    s += i;
  }
  return s;  
}

int main(int argc, char *argv[])
{
  struct timeval tv[2];
  gettimeofday(tv + 0, 0);    // 開始時刻記録
  long long s = sumUp(100000000LL); // 計算本体
  gettimeofday(tv + 1, 0);    // 終了時刻記録
  printf("sumup(100000000) = %lld\n", s);
  printf("time: %f\n", elapsed_time(tv));
}


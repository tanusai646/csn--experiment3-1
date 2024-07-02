#ifndef JUDGE_H
#define JUDGE_H

#include <stdlib.h>
#include "xydata.h"

/* この時間制限以下で，別途，競技用の時間制限を */
#define TIME_LIMIT 200.0

/* 以下のAPI
   全体として，マルチスレッドでの並行呼び出しは不可とする．
   single や crititalで保護すること */

/* ユーザ側のメモリに data 生成．審判内にも同じデータを持つ．
   配列 data のサイズ(要素数)は sz
   sz は 2^d と仮定．d は偶数と仮定 */
int judge_gendata(struct xy *data, size_t sz, int d, double key1, double key2);

/* judge_gendataしてからの時間(秒) */
double judge_gettime();

/* 報告．中間報告許す．評価値は審判として再計算 */
int judge_report(double x0, double y0);

/* 終了処理 */
int judge_finalize();

#endif

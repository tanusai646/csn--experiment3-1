#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "judge.h"
#include "mypng.h"

/* 近すぎるときの発散防止 */
#define EPARAM 1e-40

////////////////////////////////////////
// 評価値を計算 (どちらも正となるように符号を考える)
//   評価値1: 重力ポテンシャル(の総和)/総質量
//   評価値2: 重力加速度(の総和)の絶対値/総質量 × 系の半径

static int get_vals(double *val1, double *val2, double x0, double y0, struct xy *data, size_t sz){
	double poten = 0.0;
	double accelx = 0.0, accely = 0.0;

	//並列化
    #pragma omp parallel
    {
        double local_poten = 0.0;
        double local_accelx = 0.0;
        double local_accely = 0.0;

        #pragma omp for nowait
        for (size_t i = 0; i < sz; i++) {
            double x = data[i].x - x0;
            double y = data[i].y - y0;
            double x2 = x * x;
            double y2 = y * y;
            double r2 = x2 + y2 + EPARAM;
            double r_2 = 1.0 / r2;
            double r_1 = sqrt(r_2);
            double r_3 = r_1 * r_2;

            // 重力ポテンシャルの総和
            local_poten += r_1;
            // 加速度ベクトルの総和
            local_accelx += r_3 * x;
            local_accely += r_3 * y;
        }

        // スレッドごとに計算した結果を集約
        #pragma omp atomic
        poten += local_poten;
        #pragma omp atomic
        accelx += local_accelx;
        #pragma omp atomic
        accely += local_accely;
    }

	*val1 = poten / sz;
	*val2 = sqrt(accelx * accelx + accely * accely) / sz;
	return 0;
}


////////////////////////////////////////
/*
  候補位置 (*xp, *yp) をできるだけ改良．改良した回数を返す．最大 maxi回改良
  ポテンシャルは使わず，加速度が0になるところを目指す
*/

// 勾配降下法を使って候補位置 (*xp, *yp) をできるだけ改良．改良した回数を返す．最大 maxi回改良
int improve(double *xp, double *yp, int maxi, struct xy *data, size_t sz) {
    // 候補位置
    double x0 = *xp, y0 = *yp;
    // 勾配ベクトル
    double gradx = 0.0, grady = 0.0;
    // ステップサイズ
    double step_size = 0.01;
    double step_size_min = 1e-5; // 最小ステップサイズ
    double step_size_max = 0.1;  // 最大ステップサイズ
    // 更新量の大きさ^2，更新後候補位置の中心からの距離^2
    double u2, nr2;
    // 改良用反復の番号
    int ii = 0;
    // 前回の評価値
    double prev_val1 = 1e50, prev_val2 = 1e50;

    while (ii < maxi) {
        // 勾配を総和
        #pragma omp parallel
        {
            double local_gradx = 0.0, local_grady = 0.0;

            #pragma omp for
            for (size_t i = 0; i < sz; i++) {
                double x = data[i].x - x0;
                double y = data[i].y - y0;
                double x2 = x * x;
                double y2 = y * y;
                double r2 = x2 + y2 + EPARAM;
                double r_2 = 1.0 / r2;
                double r_1 = sqrt(r_2);
                double r_3 = r_1 * r_2;

                // 勾配の総和
                local_gradx += r_3 * x;
                local_grady += r_3 * y;
            }

            #pragma omp critical
            {
                gradx += local_gradx;
                grady += local_grady;
            }
        }

        // 勾配のノルムを計算
        double grad_norm = sqrt(gradx * gradx + grady * grady);

        // 勾配が非常に小さい場合は収束とみなす
        if (grad_norm < 1e-6) {
            break;
        }

        // ステップサイズの調整
        step_size = fmin(fmax(step_size / grad_norm, step_size_min), step_size_max);

        // 勾配に基づいて更新
        x0 -= step_size * gradx;
        y0 -= step_size * grady;

        // 更新量の大きさ^2
        u2 = gradx * gradx + grady * grady;

        // 更新後候補位置の中心からの距離^2
        nr2 = x0 * x0 + y0 * y0;

        // 系の枠外 (中心からの距離が1を超える) なら，枠内ぎりぎりに
        if (nr2 > 1.0) {
            double nr_1 = sqrt(1.0 / nr2);
            x0 *= nr_1;
            y0 *= nr_1;
        }

        // 次の反復へ
        ii++;

        // 次の反復で勾配を求める準備
        gradx = 0.0;
        grady = 0.0;

        // 早期に終了させる処理
        double val1, val2;
        get_vals(&val1, &val2, x0, y0, data, sz);
        if (fabs(val1 - prev_val1) < 1e-6 && fabs(val2 - prev_val2) < 1e-6) {
            break;  // Early stopping if improvement is very small
        }
        prev_val1 = val1;
        prev_val2 = val2;
    }

    *xp = x0;
    *yp = y0;
    return ii;
}


// 時間制限を満たすために処理をスキップするかの判断のため
#define MY_TIME_LIMIT 28.2

// 系内での「適当な」候補位置についてその改善をして探す．
int scanarea(int s, double key1, double key2, struct xy *data, size_t sz){
	double ssa = fabs(cos(key2)); // 適当な差分
	double ss0 = fabs(sin(key2)); // 適当なベース
	size_t ssz = (((size_t)1) << s); // スキャンする候補位置の数
	double minval = 1e+50;     // これまで最良の評価値
	double minval1 = 1e+50;    // これまで最良の評価値1
	double x1 = 0.0, y1 = 0.0; // これまで最良の位置
	struct xy *b = (struct xy *)calloc(ssz, sizeof(struct xy)); // 改良前
	if (b == 0) return -2;
	struct xy *c = (struct xy *)calloc(ssz, sizeof(struct xy)); // 改良後
	if (c == 0) return -2;
	int i0;

	int trial_count = 0;

	//並列化処理
	#pragma omp parallel
	{
		double local_minval = minval;
		double local_minval1 = minval1;
		double local_x1=0, local_y1=0;
		int local_trial_count = 0;
		double density_factor = 0.03;

		#pragma omp for
		for (i0 = 0; i0 < ssz; i0++){ // 各候補位置について
			if (judge_gettime() < MY_TIME_LIMIT){ // 時間が残っていれば処理
				double ss = ss0 + ssa * i0;
				ss = ss - (long long)ss;
				ss = ss - (int)ss;    // 小数部分を残す (0から1)
				double rr = sqrt(ss); // 中心からの距離，密度(面積あたり)を一定に
				// 中心から方向を表す単位ベクトル (xkey1, ykey1)
				double xkey1 = cos(key1 * i0), ykey1 = sin(key1 * i0);
				double x0 = rr * xkey1, y0 = rr * ykey1; // 候補位置 (距離，方向)
				b[i0].x = x0; 
				b[i0].y = y0; // 表示用
				double val1, val2;
				get_vals(&val1, &val2, x0, y0, data, sz); // 評価値を得る
				double val = val1 + val2 * 1.0;
				fprintf(stderr, "*%d: (%f, %f) with %f %f %f\n", i0, x0, y0, val, val1, val2); // 初期評価値を表示
				
				// できるだけ改良
				int ii = improve(&x0, &y0, 100, data, sz);
				c[i0].x = x0; 
				c[i0].y = y0;
				get_vals(&val1, &val2, x0, y0, data, sz);
				val = val1 + val2 * 1.0;
				fprintf(stderr, ":%d(%d): (%f, %f) with %f %f %f\n", i0, ii, x0, y0, val, val1, val2); // 改良後評価値を表示
				
				local_trial_count++;

                // これまでより良い評価なら，最良値の情報を更新
                if (val < local_minval) {
                    local_x1 = x0;
                    local_y1 = y0;
                    local_minval = val;
                    fprintf(stderr, "!:%d\n", i0);
                }
                if (val1 < local_minval1) { // 評価値1 (のみ) についても参考までに
                    local_minval1 = val1;
                    fprintf(stderr, "?:%d\n", i0);
                }
			}
		}

		#pragma omp critical
		{
			// 各スレッドのローカル値を全体での最良値と比較して更新
			if (local_minval < minval) {
				minval = local_minval;
				x1 = local_x1; 
				y1 = local_y1;
			}
			if (local_minval1 < minval1) {
				minval1 = local_minval1;
			}
		}
	}
	// 仕上げ
	improve(&x1, &y1, 3, data, sz);
	// まず審判に報告
	judge_report(x1, y1);
	// 以下は，情報表示
	// 評価値表示
	double val1, val2;
	get_vals(&val1, &val2, x1, y1, data, sz);
	double val = val1 + val2 * 1.0;
	fprintf(stderr, "min: (%f, %f) with %f %f %f\n", x1, y1, val, val1, val2);
	fprintf(stderr, "Total trials: %d\n", trial_count); // 試行回数を出力

	// PNG 生成 とりあえず標準出力に (ファイルをfopenする手もあるが)
	int r = data_ans_to_png_img(stdout, x1, y1, c, ssz, b, ssz, data, sz);
	// 近傍の点を表示
	size_t i;
	//並列処理
	#pragma omp parallel for
	for (i = 0; i < sz; i++){
		double r2 = (data[i].x - x1) * (data[i].x - x1) + (data[i].y - y1) * (data[i].y - y1);
		if (r2 < 0.0001){
			fprintf(stderr, "neib: %f, %f dist: %f\n", data[i].x, data[i].y, sqrt(r2));
		}
	}
	
	free(c);
	free(b);
	return r;
}

int main(int argc, char *argv[]) {
	size_t sz = 0;
	struct xy *data = 0;
	if (argc > 4){
		// 問題を得る
		{
			int d = atoi(argv[1]); // 問題サイズは 2^d (dは偶数)
			if (d < 0 || d > 50 || (d & 1)) return -1;
			double key1 = 0.0, key2 = 0.0;
			int r;
			r = sscanf(argv[2], "%lf", &key1);
			if (r != 1) return -1;
			r = sscanf(argv[3], "%lf", &key2);
			if (r != 1) return -1;
			sz = (((size_t)1) << d);
			data = (struct xy *)malloc(sizeof(struct xy) << d);
			if (data == 0) return -1;
			// ログ情報は標準エラー出力に
			fprintf(stderr, "data keys: %ld %d %f %f\n", (long)sz, d, key1, key2);
			// 審判から問題のデータを得る
			r = judge_gendata(data, sz, d, key1, key2);
			if (r != 0) return -1;
		}
		// 問題を解いて画像出力，または問題の画像出力
		{
			int s = atoi(argv[4]);
			if (s < 0 || s > 50) return -1;
			// スキャン用のデータ．実験して良さそうなものを選んでいる
			double key1 = 1.1, key2 = 0.5;
			int r;
			if (s)
			// s が 0 でないなら，問題を解き，judge_reportし，画像出力
			r = scanarea(s, key1, key2, data, sz);
			else
			// PNGデータは標準出力に (ファイルをfopenする手もあるが)
			r = data_to_png_img(stdout, data, sz);
			// 審判の終了処理
			r |= judge_finalize();
			return r;
		}
		}
	else
		{
		fprintf(stderr,
			"Usage: %s LOG2SIZE KEY1 KEY2 LOG2SIZE_SCAN\n"
			"LOG2SIZE: 整数d．問題サイズが 2^d．dは偶数 (例: 20)\n"
			"KEY1: 倍精度浮動小数点数．問題生成の鍵1    (例: 1.07)\n"
			"KEY2: 倍精度浮動小数点数．問題生成の鍵2    (例: 0.7)\n"
			"LOG2SIZE_SCAN: 整数s．スキャン回数 2^s     (例: 12)\n",
			argv[0]);
		return -1;
		}
}

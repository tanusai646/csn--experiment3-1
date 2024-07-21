#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "judge.h"
#include "mypng.h"

/* 近すぎるときの発散防止 */
#define EPARAM 1e-40

//Barnes-Hutアルゴリズム用のデータ構造を追加
typedef struct QuadNode {
    double x, y;  // 中心座標
    double mass;  // ノードの質量
    struct QuadNode *children[4];  // 4つの子ノード
    struct QuadNode *parent;  // 親ノード
    struct xy *points;  // このノードに含まれる点
    size_t num_points;  // 点の数
	double x_min, x_max;  // ノードの領域の境界
    double y_min, y_max;
} QuadNode;

QuadNode* create_quad_node(double x_min, double x_max, double y_min, double y_max){
	QuadNode *node = (QuadNode *)malloc(sizeof(QuadNode));
    node->x = 0.0;
    node->y = 0.0;
    node->mass = 0.0;
    node->points = NULL;
    node->num_points = 0;
    node->x_min = x_min;
    node->x_max = x_max;
    node->y_min = y_min;
    node->y_max = y_max;
    for (int i = 0; i < 4; i++) {
        node->children[i] = NULL;
    }
    return node;
};

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

int improve(double *xp, double *yp, int maxi, struct xy *data, size_t sz){
	// 候補位置
	double x0 = *xp, y0 = *yp;
	// 加速度ベクトル
	double accelx = 0.0, accely = 0.0;
	// 加速度ベクトルの勾配，2x2のヘッセ行列．2要素分は同じなので3要素で
	double derivxx = 0.0, derivxy = 0.0, derivyy = 0.0;
	// 更新量の大きさ^2，更新後候補位置の中心からの距離^2
	double u2, nr2;
	// 改良用反復の番号
	int ii = 0;
	//
	double prev_val1 = 1e50, prev_val2 = 1e50;

	while (ii < maxi){
      // 加速度ベクトルとその勾配を総和
	#pragma omp parallel
        {
            // 各スレッドのローカル変数を定義
            double local_accelx = 0.0, local_accely = 0.0;
            double local_derivxx = 0.0, local_derivxy = 0.0, local_derivyy = 0.0;

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

                // 加速度ベクトルの総和
                local_accelx += r_3 * x;
                local_accely += r_3 * y;

                // 加速度ベクトルの勾配の総和
                local_derivxx += 3.0 * r_3 * (x2 * r_2 - 1.0 / 3.0);
                local_derivxy += 3.0 * r_3 * x * y * r_2;
                local_derivyy += 3.0 * r_3 * (y2 * r_2 - 1.0 / 3.0);
            }

            #pragma omp critical
            {
                accelx += local_accelx;
                accely += local_accely;
                derivxx += local_derivxx;
                derivxy += local_derivxy;
                derivyy += local_derivyy;
            }
        }
      // x0, y0 を更新．u2, nr2 によって打ち切るか決められるようにする
    	{
#if 0
		fprintf(stderr, "pos=(%f,%f) on %d\n", x0, y0, ii);
		fprintf(stderr, "accel=(%f,%f) on %d\n", accelx, accely, ii);
		fprintf(stderr, "deriv=(%f,%f,%f) on %d\n", derivxx, derivxy, derivyy, ii);
#endif
		// 2x2 ヘッセ行列: a, b, c, d
		double a = derivxx, b = derivxy, c = derivxy, d = derivyy;
		// 行列式を用いて，逆行列を扱う
		double det = a * d - b * c;

		if(fabs(det) < 1e-20){
			det = (det < 0 ? -1 : 1) * 1e-20;
		}

		// 加速度の総和を0にするための更新量を求める (ニュートン法)
		double ux = - 1.0 / det * (  d * accelx - b * accely);
		double uy = - 1.0 / det * (- c * accelx + a * accely);
#if 0
		fprintf(stderr, "update %f,%f\n", ux, uy);
#endif
		// 更新量の大きさ^2
		u2 = ux * ux + uy * uy;

        // ステップサイズのスケーリングファクタを計算
        double scale = 1.0 / sqrt(u2);
        double max_step = 0.1; // 最大ステップサイズを設定

        if (scale > max_step) {
            scale = max_step;
        }

		double nx = x0 + ux, ny = y0 + uy;
		// 更新後候補位置の中心からの距離^2
		nr2 = nx * nx + ny * ny;
#if 0
		double pax = accelx + a * (nx - x0) + b * (ny - y0);
		double pay = accely + c * (nx - x0) + d * (ny - y0);
		fprintf(stderr, "accelp=(%f,%f)\n", pax, pay);
#endif
		x0 = nx; 
		y0 = ny;
		if (nr2 > 1.0) { 
			// 系の枠外 (中心からの距離が1を超える) なら，枠内ぎりぎりに
			double nr_1 = sqrt(1.0 / nr2);
			x0 = nx * nr_1; 
			y0 = ny * nr_1;
		}
		// 次の反復へ． 
		ii++;
		// 次の反復で加速度ベクトルとその勾配を求める準備
		accelx = 0.0; 
		accely = 0.0;
		derivxx = 0.0; 
		derivxy = 0.0; 
		derivyy = 0.0;
		}
		// u2, nr2 によって次の反復は実行せずに打ち切るか決める
		if (nr2 > 1.0 || u2 < 1e-20){
			break;
		}

		//早期に終了させる処理
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
				//fprintf(stderr, "*%d: (%f, %f) with %f %f %f\n", i0, x0, y0, val, val1, val2); // 初期評価値を表示
				// できるだけ改良
				int ii = improve(&x0, &y0, 100, data, sz);
				c[i0].x = x0; 
				c[i0].y = y0;
				get_vals(&val1, &val2, x0, y0, data, sz);
				val = val1 + val2 * 1.0;
				//fprintf(stderr, ":%d(%d): (%f, %f) with %f %f %f\n", i0, ii, x0, y0, val, val1, val2); // 改良後評価値を表示
				
				local_trial_count++;

				// まず初期評価値で一定以下の場合にのみ改良を行う
                if (val < minval * 1.5) { // 例として最良値の1.5倍以下とする
                    // できるだけ改良
                    int ii = improve(&x0, &y0, 100, data, sz);
                    c[i0].x = x0;
                    c[i0].y = y0;
                    get_vals(&val1, &val2, x0, y0, data, sz);
                    val = val1 + val2 * 1.0;
                    //fprintf(stderr, ":%d(%d): (%f, %f) with %f %f %f\n", i0, ii, x0, y0, val, val1, val2); // 改良後評価値を表示

                    // これまでより良い評価なら，最良値の情報を更新
                    if (val < local_minval) {
                        local_x1 = x0;
                        local_y1 = y0;
                        local_minval = val;
                        //fprintf(stderr, "!:%d\n", i0);
                    }
                    if (val1 < local_minval1) { // 評価値1 (のみ) についても参考までに
                        local_minval1 = val1;
                        //fprintf(stderr, "?:%d\n", i0);
                    }
					trial_count += local_trial_count;
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

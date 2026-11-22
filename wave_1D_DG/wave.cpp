#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 1.0 // Length
#ifndef Ne
#define Ne 50
#endif
#define NP 2 // nodes per element;
#define C 340.0 // sound speed
#define T 0.1
#define x0 (0.5*L)
#define sigma 0.05
#define CFL 0.2

double sigma_inv = 1 / sigma;
double L_inv = 1 / L;

// initial gaussian wave
double gaussian(double x){
    return exp(-(x-x0)*(x-x0) * sigma_inv * sigma_inv);
}

// sx 
double dgaussian(double x){
    return -2.0 * (x-x0) * sigma_inv * sigma_inv *  gaussian(x);
}

// analytical solotion (D’Alembert) for fixed-end boundary 
// u = 0.5 * [gaussian]
double fixed_end(double s){
    s -= 2.0 * L * floor((s+L) * 0.5 * L_inv);
    return copysign(gaussian(fabs(s)), s);
}

double analytical_sol(double x, double t){
    return 0.5 * (fixed_end(x - C * t) + fixed_end(x + C * t));
}

void make_f(double *v, double *s, double *f_v, double *f_s){
    for (int i = 0; i < 2; i++){
        f_v[i] = -C*C*s[i];
        f_s[i] = -v[i];
    }
}

void upwind_flux(double vl, double sl, double vr, double sr, double *Fv, double *Fs){
    *Fv = 0.5 * (-C*C * sl - C*C* sr) - 0.5 * C * (vr - vl);
    *Fs = 0.5 * (-vl - vr) - 0.5 * C * (sr - sl);
}

void calc_flux(double v[][2], double s[][2], double *flux_v, double *flux_s){
    // j=0のupwindをここで呼べばいい、for のなかの分岐は減らした方がいい
    upwind_flux(-v[0][0], s[0][0],v[0][0], s[0][0], &flux_v[0], &flux_s[0]);
    for (int j = 1; j < Ne; j++){
        upwind_flux(v[j-1][1], s[j-1][1], v[j][0], s[j][0], &flux_v[j], &flux_s[j]);
    }
    // メモリアクセス的には順番通りにした方がいい
    upwind_flux(v[Ne-1][1], s[Ne-1][1],-v[Ne-1][1], s[Ne-1][1], &flux_v[Ne], &flux_s[Ne]);
}

void make_flux_vec(double *flux, double *flux_mat, int e){
        flux_mat[0] = flux[e];
        flux_mat[1] = -flux[e+1];   
}
// matrix * vec multipy (2*2 * 2*1)
void matvec(double x[2][2], double y[2], double z[2]){
    for (int i = 0; i < 2; i++){
            z[i] = x[i][0] * y[0] + x[i][1] * y[1];
    }
}

// GEMM X*Y + Z
void gemv(double x[2][2], double y[2], double z[2], double w[2]){
    for (int i = 0; i < 2; i++){
            w[i] = x[i][0] * y[0] + x[i][1] * y[1] + z[i]; 
    }
}

// M_inv * (Se * F + flux)
void rhs_mat_calc(double Me_inv[2][2], double Se[2][2], double f[2], double flux_mat[2], double dq[2]){
    double tmp[2];
    gemv(Se, f, flux_mat, tmp);
    matvec(Me_inv, tmp, dq);
}

void rhs(double v[][2], double dv[][2], double s[][2], double ds[][2], double *flux_v, double *flux_s, double Se[2][2], double Me_inv[2][2]){
    calc_flux(v,s,flux_v,flux_s);
    double f_v[2], f_s[2], flux_v_mat[2], flux_s_mat[2];
    for (int e = 0; e < Ne; e++){
        make_f(v[e], s[e], f_v, f_s);
        make_flux_vec(flux_v, flux_v_mat, e);    
        rhs_mat_calc(Me_inv, Se, f_v, flux_v_mat, dv[e]); // calc dv/dt

        make_flux_vec(flux_s, flux_s_mat, e);  
        rhs_mat_calc(Me_inv, Se, f_s, flux_s_mat, ds[e]); // calc ds/dt
    }
}

// 工夫が必要　アスキーにしないで、バイナリファイルに書き出すとか（そうすると早くなる、桁落ちしなくなる）
// hdf5というライブラリを使うことが多い。pythonも読み方を知っている。配列要素、型とかを書きだしてくれるやつがある
// 時間ごとに別のファイルに書き出すするといい（壊れたときの対策となる）
// hdf5だと一つの時間ごとに、別々の物理量を出すことができる ただのバイナリだとやりにくいからこっちの方が遥かに便利になるという感じ
// メッシュは多次元配列っぽく書くか一次元配列っぽく書くか　
void snapshot(FILE *fp, double u[Ne][2], double h, double t){
    fprintf(fp, "# t = %.6e\n", t);
    for (int e = 0; e < Ne; e++){
        for (int i = 0; i < 2; i++){
            double x = (e + i) * h;
            fprintf(fp, "%.6e %.6e %.6e\n", x, u[e][i], analytical_sol(x, t));
        } 
    }
    fprintf(fp, "\n\n");
}

int main(){
    static double v[Ne][2], v1[Ne][2], v2[Ne][2];
    static double s[Ne][2], s1[Ne][2], s2[Ne][2];
    static double f_v[Ne][2], f_s[Ne][2]; // flux F
    static double flux_v[Ne+1], flux_s[Ne+1]; // numeric flux 
    static double vk1[Ne][2], vk2[Ne][2], vk3[Ne][2]; // slope of v at each stage in runge-kutta
    static double sk1[Ne][2], sk2[Ne][2], sk3[Ne][2]; // slope of s at each stage in runge-kutta
    static double u[Ne][2];

    double h = L / Ne; // the wise of each element
    double Me[2][2] = {{h/3, h/6}, {h/6, h/3}}; // mass-matrix of each elememt
    double Me_inv[2][2] = {{4/h, -2/h}, {-2/h, 4/h}};
    double Se[2][2] = {{-0.5, -0.5}, {0.5, 0.5}}; // stiffness-matrix of each element

    double dt = (h / C) * CFL; // CFL condition もう少し精度よくする

    double t = 0;
    int step = 0;
    int interval = (int)Ne/50;

    // initial condition u(0)
    for (int i = 0; i < Ne ; i++){
        double xl = i*h;
        double xr = (i+1)*h;
        u[i][0] = gaussian(xl);
        u[i][1] = gaussian(xr);
        s[i][0] = dgaussian(xl); // s_left
        s[i][1] = dgaussian(xr); // s_right
    }
    FILE *fp = fopen("wave.dat", "w");
    snapshot(fp, u, h, t);

    // runge-kutta stage3 
    // 時間方向は四次精度、ＤＧ法は理論的に空間は無限精度になるけど倍精度じゃ足りなくなる可能性がある。時間と空間精度はそろえるようにしてね
    // 時間精度はsemi-ラグランジアン法とかがあったりする、DG法と組み合わせると積分精度を改善できる可能性がある。
    while (t < T){
        // stage 1 (calc k1)
        rhs(v,vk1,s,sk1,flux_v, flux_s, Se, Me_inv);
        for (int e = 0; e < Ne; e++)
            for (int i = 0; i < 2; i++){
                    v1[e][i] = v[e][i] + dt * vk1[e][i];
                    s1[e][i] = s[e][i] + dt * sk1[e][i];
            }
        
        // stage 2 (calc k2)
        rhs(v1,vk2,s1,sk2,flux_v, flux_s,Se, Me_inv);
        for (int e = 0; e < Ne; e++)
            for (int i = 0; i < 2; i++){
                    v2[e][i] = v[e][i] + 0.25 * dt * (vk1[e][i] + vk2[e][i]);
                    s2[e][i] = s[e][i] + 0.25 * dt * (sk1[e][i] + sk2[e][i]);
            }
        // stage 3 (calc k3)
        rhs(v2,vk3,s2,sk3,flux_v, flux_s,Se, Me_inv);

        // update u(du/dt = v)
        for (int e = 0; e < Ne; e++)
            for (int i = 0; i < 2; i++)
                u[e][i] += dt * (v[e][i] / 6.0 + v1[e][i] / 6.0 + 2.0 * v2[e][i] / 3.0);
        
        // update  v, s
        for (int e = 0; e < Ne; e++)
            for (int i = 0; i < 2; i++){
                v[e][i] += dt * (vk1[e][i] / 6.0 + vk2[e][i]/ 6.0 + 2.0 * vk3[e][i] / 3.0);
                s[e][i] += dt * (sk1[e][i] / 6.0 + sk2[e][i]/ 6.0 + 2.0 * sk3[e][i] / 3.0);

            }
        t += dt;
        step += 1;

        if (step % interval == 0 && t < 0.02 ) snapshot(fp, u, h, t);
    }
    double err_max = 0.0;
    for (int e = 0; e < Ne; e++){
        for (int k = 0; k < 2; k++){
            double err = u[e][k] - analytical_sol((e+k)*h, t);
            if (fabs(err) > err_max) err_max = fabs(err);
        }
    }
    printf("Elements = %d, error = %.6e\n", Ne, err_max);
    fclose(fp);
}
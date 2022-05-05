#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <time.h>

const int NX = 200;
const int NY = 200;
const int NZ = 320;
const int nx = 126;
const int ny = 126;
const int nz = 126;

const double D = 30;
const int dir = 19;

int fineindex(int i, int j, int k, int a)
{
    int c;
    c = ((i * (ny + 2) + j) * (nz + 2) + k) * dir + a;
    return c;
}

int fineisnindex1(int i, int j, int k)
{
    return (i * (ny + 2) + j) * (nz + 2) + k;
}

int topbottom(int i, int k, int a)
{
    return (i * (nz + 2) + k) * dir + a;
}

int leftright(int j, int k, int a)
{
    return (j * (nz + 2) + k) * dir + a;
}

int updown(int i, int j, int a)
{
    return (i * (ny + 2) + j) * dir + a;
}

int coarseindex(int i, int j, int k, int a)
{
    int c;
    c = ((i * NY + j) * NZ + k) * dir + a;
    return c;
}

int coarseisnindex1(int i, int j, int k)
{
    return (i * NY + j) * NZ + k;
}

void interpolation_fine_boundary(double *ff, double *feqf, double *fneqf_b, double *fneqf_r, double *fneqf_t, double *fneqf_l, double *fneqf_u, double *fneqf_d, int side, int line, int ex[19], int ey[19], int ez[19], double w[19], double *rhof, double *jxf, double *jyf, double *jzf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F);
void dist_transfer_cf(double *ff, double *fc, int Start1, int End1, int Start2, int End2, int face, int side, int ex[19], int ey[19], int ez[19], double *rhoc, double *jxc, double *jyc, double *jzc, double *feqc, double w[19], double tauc, double tauf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F);
void dist_transfer_fc(double *ff, double *fc, int Start1, int End1, int Start2, int End2, int face, int side, int ex[19], int ey[19], int ez[19], double *rhof, double *jxf, double *jyf, double *jzf, double *feqf, double w[19], double tauc, double tauf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F);
void transfer_mid_boundary(double *fc, double *ff, int Start1, int End1, int Start2, int End2, int side, int start, double *feqc, double *feqf, double *jxc, double *jyc, double *jzc, double *jxf, double *jyf, double *jzf, double *rhoc, double *rhof, double tauc, double tauf, int ex[19], int ey[19], int ez[19], double w[19], int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F);
void interpolation_mid_boundary(double *ff, double *feqf, double *fneqf_b, double *fneqf_r, double *fneqf_t, double *fneqf_l, double *fneqf_u, double *fneqf_d, int side, int face, int ex[19], int ey[19], int ez[19], double w[19], double *rhof, double *jxf, double *jyf, double *jzf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F);
void bicubic_interpolation(double *ff, double *feqf, double *fneqf_b, double *fneqf_r, double *fneqf_t, double *fneqf_l, double *fneqf_u, double *fneqf_d, int face, int side, double *jxf, double *jyf, double *jzf, double *rhof, int ex[19], int ey[19], int ez[19], double w[19], int alpha);

FILE *in, *out;

int main()
{
    int a, i, j, k, v, ia, ja, ka, IA, JA, KA, I, J, K, count, neb;
    int ex[19], ey[19], ez[19], kb[19];
    double dist_l[19], dist_eq_l[19], fneqf[19];
    double sc[19], sf[19];
    double *ff, *ftf, *feqf, *fc, *ftc, *feqc;
    int *isnf, *isnf1, *isnc, *wb;
    double *jxc, *jyc, *jzc, *rhoc, *jxf, *jyf, *jzf, *rhof, *fneqf_b, *fneqf_l, *fneqf_r, *fneqf_t, *fneqf_u, *fneqf_d;
    double rho0 = 1.0;
    double CX[5], CY[5], CZ[5];
    double cx[5], cy[5], cz[5], ucx[5], ucy[5], ucz[5], omega_x[5], omega_y[5], omega_z[5];
    double ubx, uby, ubz;
    int alpha = 2;
    double tauc = 0.8, tauf = 0.5 + alpha * (tauc - 0.5); // Re=20
    double visc = (tauc - 0.5) / 3.0;
    double w[19];
    double H = NY;
    int ts = 1;
    double tb;
    double dist, r1, thita;
    double pi = acos(-1);
    double Fx, Fy, Fz, Tx, Ty, Tz;
    double fx, fy, fz, rx, ry, rz;
    double Pxx, Pxy, Pxz, Pyy, Pyz, Pzz, ux, uy, uz;
    int q;
    double rho_av;
    int num;
    double rhos = 1.164241164;
    double g = 0.000888734;  // gravity
    double Ct = 0.000212831; // time step calculation
    double Cl = 0.0005;

    double mass, Inertia, Fg;
    mass = rhos * (pi / 6.0) * D * D * D;            // mass
    Inertia = mass * D * D / 10.0;                   // polar moment of inertia
    Fg = (pi / 6.0) * D * D * D * (rho0 - rhos) * g; // Gravity Force

    double rho_l, jx_l, jy_l, rho, rho1;
    double dist_neq_l[4][19];

    double delxf = 0.5, delxc = 1.0, delyc = 1.0, delzc = 1.0; // lattice units in coarse and fine block
    double deltf = 0.5, deltc = 1.0;                           // Time units in coarse and fine block
    double Df = D * alpha;                                     // Diameter in fine block
    int CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F;

    // clock_t start, end;
    // double cpu_time_used;

    double M[19][19] = {
        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        {-30.0, -11.0, -11.0, -11.0, -11.0, -11.0, -11.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0},
        {12.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
        {0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0},
        {0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, -4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0},
        {0.0, 2.0, 2.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, -2.0, -2.0},
        {0.0, -4.0, -4.0, 2.0, 2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -2.0, -2.0, -2.0, -2.0},
        {0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, -2.0, -2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0}};

    double Minv[19][19] = {
        {1.0 / 19.0, -5.0 / 399.0, 1.0 / 21.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0, 1.0 / 10.0, -1.0 / 10.0, 0, 0, 0, 0, 1.0 / 18.0, -1.0 / 18.0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0, -1.0 / 10.0, 1.0 / 10.0, 0, 0, 0, 0, 1.0 / 18.0, -1.0 / 18.0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0, 0, 0, 1.0 / 10.0, -1.0 / 10.0, 0, 0, -1.0 / 36.0, 1.0 / 36.0, 1.0 / 12.0, -1.0 / 12.0, 0, 0, 0, 0, 0, 0},
        {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0, 0, 0, -1.0 / 10.0, 1.0 / 10.0, 0, 0, -1.0 / 36.0, 1.0 / 36.0, 1.0 / 12.0, -1.0 / 12.0, 0, 0, 0, 0, 0, 0},
        {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0, 0, 0, 0, 0, 1.0 / 10.0, -1.0 / 10.0, -1.0 / 36.0, 1.0 / 36.0, -1.0 / 12.0, 1.0 / 12.0, 0, 0, 0, 0, 0, 0},
        {1.0 / 19.0, -11.0 / 2394.0, -1.0 / 63.0, 0, 0, 0, 0, -1.0 / 10.0, 1.0 / 10.0, -1.0 / 36.0, 1.0 / 36.0, -1.0 / 12.0, 1.0 / 12.0, 0, 0, 0, 0, 0, 0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 1.0 / 10.0, 1.0 / 40.0, 1.0 / 10.0, 1.0 / 40.0, 0, 0, 1.0 / 36.0, 1.0 / 72.0, 1.0 / 12.0, 1.0 / 24.0, 1.0 / 4.0, 0, 0, 1.0 / 8.0, -1.0 / 8.0, 0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0, 1.0 / 10.0, 1.0 / 40.0, 0, 0, 1.0 / 36.0, 1.0 / 72.0, 1.0 / 12.0, 1.0 / 24.0, -1.0 / 4.0, 0, 0, -1.0 / 8.0, -1.0 / 8.0, 0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 1.0 / 10.0, 1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0, 0, 0, 1.0 / 36.0, 1.0 / 72.0, 1.0 / 12.0, 1.0 / 24.0, -1.0 / 4.0, 0, 0, 1.0 / 8.0, 1.0 / 8.0, 0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0, 0, 0, 1.0 / 36.0, 1.0 / 72.0, 1.0 / 12.0, 1.0 / 24.0, 1.0 / 4.0, 0, 0, -1.0 / 8.0, 1.0 / 8.0, 0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 1.0 / 10.0, 1.0 / 40.0, 0, 0, 1.0 / 10.0, 1.0 / 40.0, 1.0 / 36.0, 1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0, 0, 1.0 / 4.0, -1.0 / 8.0, 0, 1.0 / 8.0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0, 0, 0, 1.0 / 10.0, 1.0 / 40.0, 1.0 / 36.0, 1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0, 0, -1.0 / 4.0, 1.0 / 8.0, 0, 1.0 / 8.0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 1.0 / 10.0, 1.0 / 40.0, 0, 0, -1.0 / 10.0, -1.0 / 40.0, 1.0 / 36.0, 1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0, 0, -1.0 / 4.0, -1.0 / 8.0, 0, -1.0 / 8.0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, -1.0 / 10.0, -1.0 / 40.0, 0, 0, -1.0 / 10.0, -1.0 / 40.0, 1.0 / 36.0, 1.0 / 72.0, -1.0 / 12.0, -1.0 / 24.0, 0, 0, 1.0 / 4.0, 1.0 / 8.0, 0, -1.0 / 8.0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 0, 0, 1.0 / 10.0, 1.0 / 40.0, 1.0 / 10.0, 1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0, 0, 0, 1.0 / 4.0, 0, 0, 1.0 / 8.0, -1.0 / 8.0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 0, 0, -1.0 / 10.0, -1.0 / 40.0, 1.0 / 10.0, 1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0, 0, 0, -1.0 / 4.0, 0, 0, -1.0 / 8.0, -1.0 / 8.0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 0, 0, 1.0 / 10.0, 1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0, 0, 0, -1.0 / 4.0, 0, 0, 1.0 / 8.0, 1.0 / 8.0},
        {1.0 / 19.0, 4.0 / 1197.0, 1.0 / 252.0, 0, 0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 10.0, -1.0 / 40.0, -1.0 / 18.0, -1.0 / 36.0, 0, 0, 0, 1.0 / 4.0, 0, 0, -1.0 / 8.0, 1.0 / 8.0}};

    double rc[19], reqc[19], rf[19], reqf[19];

    sf[0] = 0.0;
    sf[1] = 1.0;
    sf[2] = 1.0;
    sf[3] = 0.0;
    sf[4] = 1.2;
    sf[5] = 0.0;
    sf[6] = 1.2;
    sf[7] = 0.0;
    sf[8] = 1.2;
    sf[9] = 1.0 / tauf;
    sf[10] = 1.2;
    sf[11] = 1.0 / tauf;
    sf[12] = 1.2;
    sf[13] = 1.0 / tauf;
    sf[14] = 1.0 / tauf;
    sf[15] = 1.0 / tauf;
    sf[16] = 1.2;
    sf[17] = 1.2;
    sf[18] = 1.2;

    sc[0] = 0.0;
    sc[3] = 0.0;
    sc[5] = 0.0;
    sc[7] = 0.0;
    sc[9] = 1.0 / tauc;
    sc[11] = 1.0 / tauc;
    sc[13] = 1.0 / tauc;
    sc[14] = 1.0 / tauc;
    sc[15] = 1.0 / tauc;

    sc[1] = pow((1.0 / sf[1] - 0.5) / alpha + 0.5, -1);
    sc[2] = pow((1.0 / sf[2] - 0.5) / alpha + 0.5, -1);
    sc[4] = pow((1.0 / sf[4] - 0.5) / alpha + 0.5, -1);
    sc[6] = pow((1.0 / sf[6] - 0.5) / alpha + 0.5, -1);
    sc[8] = pow((1.0 / sf[8] - 0.5) / alpha + 0.5, -1);
    sc[10] = pow((1.0 / sf[10] - 0.5) / alpha + 0.5, -1);
    sc[12] = pow((1.0 / sf[12] - 0.5) / alpha + 0.5, -1);
    sc[16] = pow((1.0 / sf[16] - 0.5) / alpha + 0.5, -1);
    sc[17] = pow((1.0 / sf[17] - 0.5) / alpha + 0.5, -1);
    sc[18] = pow((1.0 / sf[18] - 0.5) / alpha + 0.5, -1);

    //******* discrete velocity vector *********

    ex[0] = 0;
    ey[0] = 0;
    ez[0] = 0;
    ex[1] = 1;
    ey[1] = 0;
    ez[1] = 0;
    ex[2] = -1;
    ey[2] = 0;
    ez[2] = 0;
    ex[3] = 0;
    ey[3] = 1;
    ez[3] = 0;
    ex[4] = 0;
    ey[4] = -1;
    ez[4] = 0;
    ex[5] = 0;
    ey[5] = 0;
    ez[5] = 1;
    ex[6] = 0;
    ey[6] = 0;
    ez[6] = -1;
    ex[7] = 1;
    ey[7] = 1;
    ez[7] = 0;
    ex[8] = -1;
    ey[8] = 1;
    ez[8] = 0;
    ex[9] = 1;
    ey[9] = -1;
    ez[9] = 0;
    ex[10] = -1;
    ey[10] = -1;
    ez[10] = 0;
    ex[11] = 1;
    ey[11] = 0;
    ez[11] = 1;
    ex[12] = -1;
    ey[12] = 0;
    ez[12] = 1;
    ex[13] = 1;
    ey[13] = 0;
    ez[13] = -1;
    ex[14] = -1;
    ey[14] = 0;
    ez[14] = -1;
    ex[15] = 0;
    ey[15] = 1;
    ez[15] = 1;
    ex[16] = 0;
    ey[16] = -1;
    ez[16] = 1;
    ex[17] = 0;
    ey[17] = 1;
    ez[17] = -1;
    ex[18] = 0;
    ey[18] = -1;
    ez[18] = -1;

    kb[0] = 0;
    kb[1] = 2;
    kb[2] = 1;
    kb[3] = 4;
    kb[4] = 3;
    kb[5] = 6;
    kb[6] = 5;
    kb[7] = 10;
    kb[8] = 9;
    kb[9] = 8;
    kb[10] = 7;
    kb[11] = 14;
    kb[12] = 13;
    kb[13] = 12;
    kb[14] = 11;
    kb[15] = 18;
    kb[16] = 17;
    kb[17] = 16;
    kb[18] = 15;

    fc = (double *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * dir * sizeof(double));
    ftc = (double *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * dir * sizeof(double));
    feqc = (double *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * dir * sizeof(double));

    ff = (double *)malloc((nx + 5) * (ny + 5) * (nz + 5) * dir * sizeof(double));
    ftf = (double *)malloc((nx + 5) * (ny + 5) * (nz + 5) * dir * sizeof(double));
    feqf = (double *)malloc((nx + 5) * (ny + 5) * (nz + 5) * dir * sizeof(double));

    isnf = (int *)malloc((nx + 5) * (ny + 5) * (nz + 5) * sizeof(int));
    isnf1 = (int *)malloc((nx + 5) * (ny + 5) * (nz + 5) * sizeof(int));
    wb = (int *)malloc((nx + 5) * (ny + 5) * (nz + 5) * sizeof(int));
    rhof = (double *)malloc((nx + 5) * (ny + 5) * (nz + 5) * sizeof(double));
    jxf = (double *)malloc((nx + 5) * (ny + 5) * (nz + 5) * sizeof(double));
    jyf = (double *)malloc((nx + 5) * (ny + 5) * (nz + 5) * sizeof(double));
    jzf = (double *)malloc((nx + 5) * (ny + 5) * (nz + 5) * sizeof(double));

    isnc = (int *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * sizeof(int));
    rhoc = (double *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * sizeof(double));
    jxc = (double *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * sizeof(double));
    jyc = (double *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * sizeof(double));
    jzc = (double *)malloc((NX + 5) * (NY + 5) * (NZ + 5) * sizeof(double));

    fneqf_b = (double *)malloc((nx + 5) * (ny + 5) * dir * sizeof(double));
    fneqf_t = (double *)malloc((nx + 5) * (ny + 5) * dir * sizeof(double));
    fneqf_l = (double *)malloc((nx + 5) * (ny + 5) * dir * sizeof(double));
    fneqf_r = (double *)malloc((nx + 5) * (ny + 5) * dir * sizeof(double));
    fneqf_u = (double *)malloc((nx + 5) * (ny + 5) * dir * sizeof(double));
    fneqf_d = (double *)malloc((nx + 5) * (ny + 5) * dir * sizeof(double));

    //********* weight funnctions ********

    for (a = 0; a < 19; a++)
    {
        if (a == 0)
            w[a] = 1.0 / 3.0;

        else if ((a >= 1) && (a <= 6))
            w[a] = 1.0 / 18.0;

        else
            w[a] = 1.0 / 36.0;

        printf("wt[%i]=%f\n", a, w[a]);
        fflush(stdout);
    }

    CX[0] = 100.5;
    CX[1] = 100.5;
    CY[0] = 100.5;
    CY[1] = 100.5;
    CZ[0] = 255.5;
    CZ[1] = 255.5;

    CX_I = 70;
    CX_F = 131;
    CY_I = 70;
    CY_F = 131;
    CZ_I = 225;
    CZ_F = 286;

    cx[0] = 0.5 * nx;
    cy[0] = 0.5 * ny;
    cz[0] = 0.5 * nz;
    cx[1] = 0.5 * nx;
    cy[1] = 0.5 * ny;
    cz[1] = 0.5 * nz;

    double x_cen0 = CX[1], x_cen1, y_cen0 = CY[1], y_cen1, z_cen0 = CZ[1], z_cen1;
    int x_cen2 = 0, y_cen2 = 0, z_cen2 = 0;

    ucx[0] = 0.0, ucy[0] = 0.0, ucz[0] = 0.0;
    ucx[1] = 0.0, ucy[1] = 0.0, ucz[1] = 0.0;
    omega_x[0] = 0.0;
    omega_x[1] = 0.0;
    omega_y[0] = 0.0;
    omega_y[1] = 0.0;
    omega_z[0] = 0.0;
    omega_z[1] = 0.0;

    /********** setting the values of isn in coarse and fine blocks ********/
    for (I = 0; I <= NX + 1; I++)
    {
        for (J = 0; J <= NY + 1; J++)
        {
            for (K = 0; K <= NZ + 1; K++)
            {
                isnc[coarseisnindex1(I, J, K)] = 0;
            }
        }
    }

    // ******* Setting initial values of ftf[a][n][] in fine block ********
    for (i = 0; i <= nx; i++)
    {
        for (j = 0; j <= ny; j++)
        {
            for (k = 0; k <= nz; k++)
            {
                rf[0] = rho0;
                rf[1] = -11.0 * rho0;
                rf[2] = 3.0 * rho0;
                rf[3] = 0.0;
                rf[4] = 0.0;
                rf[5] = 0.0;
                rf[6] = 0.0;
                rf[7] = 0.0;
                rf[8] = 0.0;
                rf[9] = 0.0;
                rf[10] = 0.0;
                rf[11] = 0.0;
                rf[12] = 0.0;
                rf[13] = 0.0;
                rf[14] = 0.0;
                rf[15] = 0.0;
                rf[16] = 0.0;
                rf[17] = 0.0;
                rf[18] = 0.0;

                for (a = 0; a < 19; a++)
                {
                    ftf[fineindex(i, j, k, a)] = 0.0;

                    for (v = 0; v < 19; v++)
                    {
                        ftf[fineindex(i, j, k, a)] += Minv[a][v] * rf[v];
                    }
                }
            }
        }
    }

    // ******* Setting initial values of ftc[a][n][] in coarse block ********
    for (I = 0; I <= NX + 1; I++)
    {
        for (J = 0; J <= NY + 1; J++)
        {
            for (K = 0; K <= NZ + 1; K++)
            {
                rc[0] = rho0;
                rc[1] = -11.0 * rho0;
                rc[2] = 3.0 * rho0;
                rc[3] = 0.0;
                rc[4] = 0.0;
                rc[5] = 0.0;
                rc[6] = 0.0;
                rc[7] = 0.0;
                rc[8] = 0.0;
                rc[9] = 0.0;
                rc[10] = 0.0;
                rc[11] = 0.0;
                rc[12] = 0.0;
                rc[13] = 0.0;
                rc[14] = 0.0;
                rc[15] = 0.0;
                rc[16] = 0.0;
                rc[17] = 0.0;
                rc[18] = 0.0;

                for (a = 0; a < 19; a++)
                {
                    ftc[coarseindex(I, J, K, a)] = 0.0;

                    for (v = 0; v < 19; v++)
                    {
                        ftc[coarseindex(I, J, K, a)] += Minv[a][v] * rc[v];
                    }
                }
            }
        }
    }

    FILE *out1;
    out1 = fopen("Velocity_vs_time.dat", "w");

    FILE *out3;
    char sol[30] = "shearflow_", str1[] = ".dat", str[30], sol1[5];
    int solnumber = 0;
    double distx, disty, distz;

    FILE *out4;
    out4 = fopen("LD_vs_time.dat", "w");
    FILE *out5;
    out5 = fopen("Fz_vs_time.dat", "w");
    // FILE *time;
    // time = fopen("Time_CPU", "w");

    // start = clock();
    do
    {
        // ******** Streaming in the coarse block to get fc[a][n+1][] *********

        for (I = 1; I <= NX; I++)
        {
            for (J = 1; J <= NY; J++)
            {
                for (K = 1; K <= NZ; K++)
                {
                    for (a = 0; a < 19; a++)
                    {
                        IA = I + ex[a];
                        JA = J + ey[a];
                        KA = K + ez[a];
                        fc[coarseindex(IA, JA, KA, a)] = ftc[coarseindex(I, J, K, a)]; // streaming step
                    }
                }
            }
        }

        //********* boundary conditions in coarse block **********
        for (I = 1, J = 1; J <= NY; J++)
        {
            for (K = 1; K <= NZ; K++)
            {
                for (a = 0; a < 19; a++)
                {
                    IA = I - ex[a];
                    JA = J - ey[a];
                    KA = K - ez[a];

                    if (IA == 0)
                    {
                        fc[coarseindex(I, J, K, a)] = ftc[coarseindex(I, J, K, kb[a])];
                    }
                }
            }
        }

        for (I = NX, J = 1; J <= NY; J++)
        {
            for (K = 1; K <= NZ; K++)
            {
                for (a = 0; a < 19; a++)
                {
                    IA = I - ex[a];
                    JA = J - ey[a];
                    KA = K - ez[a];

                    if (IA == NX + 1)
                    {
                        fc[coarseindex(I, J, K, a)] = ftc[coarseindex(I, J, K, kb[a])];
                    }
                }
            }
        }

        for (J = NY, I = 1; I <= NX; I++)
        {
            for (K = 1; K <= NZ; K++)
            {
                for (a = 0; a < 19; a++)
                {
                    IA = I - ex[a];
                    JA = J - ey[a];
                    KA = K - ez[a];

                    if (JA == NY + 1)
                    {
                        fc[coarseindex(I, J, K, a)] = ftc[coarseindex(I, J, K, kb[a])];
                    }
                }
            }
        }

        for (J = 1, I = 1; I <= NX; I++)
        {
            for (K = 1; K <= NZ; K++)
            {
                for (a = 0; a < 19; a++)
                {
                    IA = I - ex[a];
                    JA = J - ey[a];
                    KA = K - ez[a];

                    if (JA == 0)
                    {
                        fc[coarseindex(I, J, K, a)] = ftc[coarseindex(I, J, K, kb[a])];
                    }
                }
            }
        }

        for (K = 1, I = 1; I <= NX; I++)
        {
            for (J = 1; J <= NY; J++)
            {
                for (a = 0; a < 19; a++)
                {
                    IA = I - ex[a];
                    JA = J - ey[a];
                    KA = K - ez[a];

                    if (KA == 0)
                    {
                        fc[coarseindex(I, J, K, a)] = ftc[coarseindex(I, J, K, kb[a])];
                    }
                }
            }
        }

        for (K = NZ, I = 1; I <= NX; I++)
        {
            for (J = 1; J <= NY; J++)
            {
                for (a = 0; a < 19; a++)
                {
                    IA = I - ex[a];
                    JA = J - ey[a];
                    KA = K - ez[a];

                    if (KA == NZ + 1)
                    {
                        fc[coarseindex(I, J, K, a)] = ftc[coarseindex(I, J, K, kb[a])];
                    }
                }
            }
        }

        // ****** streaming in the fine block to get ff[a][n+1/2][] *******

        for (i = 0; i <= nx; i++)
        {
            for (j = 0; j <= ny; j++)
            {
                for (k = 0; k <= nz; k++)
                {
                    for (a = 0; a < 19; a++)
                    {
                        ia = i + ex[a];
                        ja = j + ey[a];
                        ka = k + ez[a];
                        ff[fineindex(ia, ja, ka, a)] = ftf[fineindex(i, j, k, a)]; // streaming step
                    }
                }
            }
        }

        //******** updating the centre positions and isn values in fine block *********

        cx[1] = (CX[1] - CX_I + 1) * alpha;
        cy[1] = (CY[1] - CY_I + 1) * alpha;
        cz[1] = (CZ[1] - CZ_I + 1) * alpha;

        for (i = 1; i <= nx - 1; i++)
        {
            for (j = 1; j <= ny - 1; j++)
            {
                for (k = 1; k <= nz - 1; k++)
                {
                    wb[fineisnindex1(i, j, k)] = 0;
                    dist = sqrt(pow((i - cx[1]), 2) + pow((j - cy[1]), 2) + pow((k - cz[1]), 2));
                    if (dist <= (Df / 2.0))
                    {
                        isnf[fineisnindex1(i, j, k)] = 1;
                    }
                    else
                    {
                        isnf[fineisnindex1(i, j, k)] = 0;
                    }
                }
            }
        }

        for (i = 1; i <= nx - 1; i++)
        {
            for (j = 1; j <= ny - 1; j++)
            {
                for (k = 1; k <= nz - 1; k++)
                {
                    if (isnf[fineisnindex1(i, j, k)] == 0)
                    {
                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1)
                            {
                                wb[fineisnindex1(i, j, k)] = 1;
                                wb[fineisnindex1(ia, ja, ka)] = 2; // solid
                            }
                        }
                    }
                }
            }
        }

        //***** Force Calculation at time n+1/2 ******
        Fx = 0.0, Fy = 0.0, Fz = 0.0, Tx = 0.0, Ty = 0.0, Tz = 0.0;

        for (i = 1; i <= nx - 1; i++)
        {
            for (j = 1; j <= ny - 1; j++)
            {
                for (k = 1; k <= nz - 1; k++)
                {
                    if (wb[fineisnindex1(i, j, k)] == 1) // fluid boundary node
                    {
                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1)
                            { // solid boundary

                                rx = (i - cx[1]) * delxf;
                                ry = (j - cy[1]) * delxf;
                                rz = (k - cz[1]) * delxf;

                                ubx = ucx[1] + omega_y[1] * rz - omega_z[1] * ry;
                                uby = ucy[1] - omega_x[1] * rz + omega_z[1] * rx;
                                ubz = ucz[1] + omega_x[1] * ry - omega_y[1] * rx;

                                ff[fineindex(i, j, k, kb[a])] = ftf[fineindex(i, j, k, a)] + 6 * w[a] * rho0 * (ubx * ex[kb[a]] + uby * ey[kb[a]] + ubz * ez[kb[a]]); // bounceback distribution with momentum transfer
                            }
                        }

                        rho = 0.0;
                        Pxx = 0.0;
                        Pyy = 0.0;
                        Pxy = 0.0;
                        Pxz = 0.0;
                        Pyz = 0.0;
                        Pzz = 0.0;

                        for (a = 0; a < 19; a++)
                        {
                            rho += ff[fineindex(i, j, k, a)];
                            Pxx += ex[a] * ex[a] * ff[fineindex(i, j, k, a)];
                            Pxy += ex[a] * ey[a] * ff[fineindex(i, j, k, a)];
                            Pxz += ex[a] * ez[a] * ff[fineindex(i, j, k, a)];
                            Pyy += ey[a] * ey[a] * ff[fineindex(i, j, k, a)];
                            Pyz += ey[a] * ez[a] * ff[fineindex(i, j, k, a)];
                            Pzz += ez[a] * ez[a] * ff[fineindex(i, j, k, a)];
                        }

                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1) // solid
                            {
                                ff[fineindex(i, j, k, kb[a])] = w[a] * (rho + 3.0 * rho * (ubx * ex[kb[a]] + uby * ey[kb[a]] + ubz * ez[kb[a]]) + 4.5 * ((Pxx - rho / 3) * (ex[kb[a]] * ex[kb[a]] - 1.0 / 3.0) + 2 * Pxy * ex[kb[a]] * ey[kb[a]] + 2 * Pxz * ex[kb[a]] * ez[kb[a]] + 2 * Pyz * ey[kb[a]] * ez[kb[a]] + (Pyy - rho / 3) * (ey[kb[a]] * ey[kb[a]] - 1.0 / 3.0) + (Pzz - rho / 3) * (ez[kb[a]] * ez[kb[a]] - 1.0 / 3.0)));
                            }
                        }

                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1) // solid
                            {

                                Fx += (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ex[a] * delxf * delxf;
                                Fy += (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ey[a] * delxf * delxf;
                                Fz += (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ez[a] * delxf * delxf;

                                fx = (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ex[a] * delxf * delxf;
                                fy = (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ey[a] * delxf * delxf;
                                fz = (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ez[a] * delxf * delxf;

                                Tx += ry * fz - rz * fy;
                                Ty += rz * fx - rx * fz;
                                Tz += rx * fy - ry * fx;
                            }
                        }
                    }
                }
            }
        }

        Fz = Fz + Fg;

        ucx[2] = ucx[1] + deltf * Fx / mass;
        ucy[2] = ucy[1] + deltf * Fy / mass;
        ucz[2] = ucz[1] + deltf * Fz / mass;

        omega_x[2] = omega_x[1] + deltf * Tx / Inertia;
        omega_y[2] = omega_y[1] + deltf * Ty / Inertia;
        omega_z[2] = omega_z[1] + deltf * Tz / Inertia;

        CX[2] = CX[1] + deltf * ucx[2];
        CY[2] = CY[1] + deltf * ucy[2];
        CZ[2] = CZ[1] + deltf * ucz[2];

        fprintf(out1, "%f\t%f\n", ts * Ct, ucz[2] * Cl / Ct);
        fflush(out1);
        fprintf(out4, "%f\t%f\n", ts * Ct, (CZ[2] - 15) / D);
        fflush(out4);
        fprintf(out5, "%f\t%f\n", ts * Ct, Fz);
        fflush(out5);

        cx[1] = (CX[2] - CX_I + 1) * alpha;
        cy[1] = (CY[2] - CY_I + 1) * alpha;
        cz[1] = (CZ[2] - CZ_I + 1) * alpha;

        printf("%d Fx = %.8f cx = %.8f ucx = %.8f Fy = %.8f cy = %.8f ucy = %.8f Fz = %.8f cz = %.8f ucz = %.8f\n", ts, Fx, cx[1], ucx[1], Fy, cy[1], ucy[1], Fz, cz[1], ucz[1]);
        fflush(stdout);

        //******* Locating any uncovered node by storing isn parameter in a new variable (isn1) *****

        for (i = 1; i <= nx - 1; i++)
        {
            for (j = 1; j <= ny - 1; j++)
            {
                for (k = 1; k <= nz - 1; k++)
                {
                    dist = sqrt(pow((i - cx[1]), 2) + pow((j - cy[1]), 2) + pow((k - cz[1]), 2));
                    if (dist <= (Df / 2.0))
                    {
                        isnf1[fineisnindex1(i, j, k)] = 1;
                    }

                    else
                    {
                        isnf1[fineisnindex1(i, j, k)] = 0;
                    }
                }
            }
        }

        for (i = 1; i <= nx - 1; i++)
        {
            for (j = 1; j <= ny - 1; j++)
            {
                for (k = 1; k <= nz - 1; k++)
                {
                    if (isnf1[fineisnindex1(i, j, k)] == 0 && isnf[fineisnindex1(i, j, k)] == 1) // new uncovered node from solid to fluid
                    {

                        neb = 0; // count of neigbouring fluid nodes
                        rho1 = 0.0;
                        ux = 0.0;
                        uy = 0.0;
                        uz = 0.0;
                        Pxx = 0.0;
                        Pyy = 0.0;
                        Pxy = 0.0;
                        Pxz = 0.0;
                        Pyz = 0.0;
                        Pzz = 0.0;

                        // printf("i=%d,j=%d,k=%d\n",i,j,k);
                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 0) // neighbouring old fluid nodes
                            {
                                neb++;
                                for (a = 0; a < 19; a++)
                                {
                                    rho1 += ff[fineindex(ia, ja, ka, a)];
                                    ux += ex[a] * ff[fineindex(ia, ja, ka, a)];
                                    uy += ey[a] * ff[fineindex(ia, ja, ka, a)];
                                    uz += ez[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pxx += ex[a] * ex[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pxy += ex[a] * ey[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pxz += ex[a] * ez[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pyy += ey[a] * ey[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pyz += ey[a] * ez[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pzz += ez[a] * ez[a] * ff[fineindex(ia, ja, ka, a)];
                                }
                            }
                        }

                        ux = ux / rho1;
                        uy = uy / rho1;
                        uz = uz / rho1;
                        rho1 = rho1 / neb;
                        Pxx = Pxx / neb;
                        Pxy = Pxy / neb;
                        Pxz = Pxz / neb;
                        Pyy = Pyy / neb;
                        Pyz = Pyz / neb;
                        Pzz = Pzz / neb;

                        for (a = 0; a < 19; a++)
                        {
                            ff[fineindex(i, j, k, a)] = w[a] * (rho1 + 3.0 * rho1 * (ux * ex[a] + uy * ey[a] + uz * ez[a]) + 4.5 * ((Pxx - rho1 / 3) * (ex[a] * ex[a] - 1.0 / 3.0) + 2 * Pxy * ex[a] * ey[a] + 2 * Pxz * ex[a] * ez[a] + 2 * Pyz * ey[a] * ez[a] + (Pyy - rho1 / 3) * (ey[a] * ey[a] - 1.0 / 3.0) + (Pzz - rho1 / 3) * (ez[a] * ez[a] - 1.0 / 3.0))); // uncovered node
                                                                                                                                                                                                                                                                                                                                                                      //    printf("ff[%d][%d][%d][%d]=%f\n",a,i,j,k,ff[fineindex(i, j, k,a)]);
                        }
                    }
                }
            }
        }

        //    printf("neb=%d\n",neb);

        CX[1] = CX[2];
        CY[1] = CY[2];
        CZ[1] = CZ[2];

        ucx[1] = ucx[2];
        ucy[1] = ucy[2];
        ucz[1] = ucz[2];

        omega_x[1] = omega_x[2];
        omega_y[1] = omega_y[2];
        omega_z[1] = omega_z[2];

        // ****** calculating physical variables ux[n+1/2],uy[n+1/2],rho[n+1/2] and coLlision to get ftf[a][n+1/2][] in fine block ******

        // den_av=0.0;
        count = 0;
        rho_av = 0.0;
        for (i = 1; i <= nx - 1; i++)
        {
            for (j = 1; j <= ny - 1; j++)
            {
                for (k = 1; k <= nz - 1; k++)
                {
                    if (isnf1[fineisnindex1(i, j, k)] == 0)
                    {
                        for (a = 0; a < 19; a++)
                        {
                            rf[a] = 0.0;
                            for (v = 0; v < 19; v++)
                            {
                                rf[a] += M[a][v] * ff[fineindex(i, j, k, v)];
                            }
                        }

                        rhof[fineisnindex1(i, j, k)] = rf[0];
                        jxf[fineisnindex1(i, j, k)] = rf[3];
                        jyf[fineisnindex1(i, j, k)] = rf[5];
                        jzf[fineisnindex1(i, j, k)] = rf[7];

                        reqf[0] = 0.0;
                        reqf[1] = -11.0 * rhof[fineisnindex1(i, j, k)] + 19.0 * (jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[2] = 3.0 * rhof[fineisnindex1(i, j, k)] - 11.0 * (jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / (2.0 * rho0);
                        reqf[3] = 0.0;
                        reqf[4] = -(2.0 / 3.0) * jxf[fineisnindex1(i, j, k)];
                        reqf[5] = 0.0;
                        reqf[6] = -(2.0 / 3.0) * jyf[fineisnindex1(i, j, k)];
                        reqf[7] = 0.0;
                        reqf[8] = -(2.0 / 3.0) * jzf[fineisnindex1(i, j, k)];
                        reqf[9] = (3.0 * jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] - (jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)])) / rho0;
                        reqf[10] = -0.5 * reqf[9]; //-(3.0*jx[i][j][k]*jx[i][j][k] - (jx[i][j][k]*jx[i][j][k] + jy[i][j][k]*jy[i][j][k] + jz[i][j][k]*jz[i][j][k]))/(2.0*rho0);
                        reqf[11] = (jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] - jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[12] = -0.5 * reqf[11]; //-(jy[i][j][k]*jy[i][j][k] - jz[i][j][k]*jz[i][j][k])/(2.0*rho0);
                        reqf[13] = (jxf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[14] = (jyf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[15] = (jxf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[16] = 0.0;
                        reqf[17] = 0.0;
                        reqf[18] = 0.0;

                        count++;
                        rho_av += rhof[fineisnindex1(i, j, k)];

                        for (a = 0; a < 19; a++)
                        {
                            rf[a] = rf[a] - sf[a] * (rf[a] - reqf[a]); // collision
                        }

                        for (a = 0; a < 19; a++)
                        {
                            ftf[fineindex(i, j, k, a)] = 0.0;

                            for (v = 0; v < 19; v++)
                            {
                                ftf[fineindex(i, j, k, a)] += Minv[a][v] * rf[v];
                            }
                        }
                    }
                }
            }
        }

        // printf("Average density fine block = %f\tCount=%d\n",rho_av/count,count);
        // ****** streaming in fine block to get ff[a][n+1][] ******
        for (i = 1; i <= nx - 1; i++)
        {
            for (j = 1; j <= ny - 1; j++)
            {
                for (k = 1; k <= nz - 1; k++)
                {
                    for (a = 0; a < 19; a++)
                    {
                        ia = i + ex[a];
                        ja = j + ey[a];
                        ka = k + ez[a];
                        ff[fineindex(ia, ja, ka, a)] = ftf[fineindex(i, j, k, a)]; // streaming step
                    }
                }
            }
        }

        // printf("CZ=%f\n",CZ[1]);

        //******** updating the centre positions and isn values in fine block *********

        cx[1] = (CX[1] - CX_I + 1) * alpha;
        cy[1] = (CY[1] - CY_I + 1) * alpha;
        cz[1] = (CZ[1] - CZ_I + 1) * alpha;

        for (i = 2; i <= nx - 2; i++)
        {
            for (j = 2; j <= ny - 2; j++)
            {
                for (k = 2; k <= nz - 2; k++)
                {
                    wb[fineisnindex1(i, j, k)] = 0;
                    dist = sqrt(pow((i - cx[1]), 2) + pow((j - cy[1]), 2) + pow((k - cz[1]), 2));
                    if (dist <= (Df / 2.0))
                    {
                        isnf[fineisnindex1(i, j, k)] = 1;
                    }
                    else
                    {
                        isnf[fineisnindex1(i, j, k)] = 0;
                    }
                }
            }
        }

        for (i = 2; i <= nx - 2; i++)
        {
            for (j = 2; j <= ny - 2; j++)
            {
                for (k = 2; k <= nz - 2; k++)
                {
                    if (isnf[fineisnindex1(i, j, k)] == 0)
                    {
                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1)
                            {
                                wb[fineisnindex1(i, j, k)] = 1;    // fluid boundary node
                                wb[fineisnindex1(ia, ja, ka)] = 2; // solid
                            }
                        }
                    }
                }
            }
        }

        //***** Force Calculation at time n+1 ******
        Fx = 0.0, Fy = 0.0, Fz = 0.0, Tx = 0.0, Ty = 0.0, Tz = 0.0;

        for (i = 2; i <= nx - 2; i++)
        {
            for (j = 2; j <= ny - 2; j++)
            {
                for (k = 2; k <= nz - 2; k++)
                {
                    if (wb[fineisnindex1(i, j, k)] == 1) // fluid boundary node
                    {
                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1)
                            { // solid boundary

                                rx = (i - cx[1]) * delxf;
                                ry = (j - cy[1]) * delxf;
                                rz = (k - cz[1]) * delxf;

                                ubx = ucx[1] + omega_y[1] * rz - omega_z[1] * ry;
                                uby = ucy[1] - omega_x[1] * rz + omega_z[1] * rx;
                                ubz = ucz[1] + omega_x[1] * ry - omega_y[1] * rx;

                                ff[fineindex(i, j, k, kb[a])] = ftf[fineindex(i, j, k, a)] + 6 * w[a] * rho0 * (ubx * ex[kb[a]] + uby * ey[kb[a]] + ubz * ez[kb[a]]); // bounceback distribution with momentum transfer
                            }
                        }

                        rho = 0.0;
                        Pxx = 0.0;
                        Pyy = 0.0;
                        Pxy = 0.0;
                        Pxz = 0.0;
                        Pyz = 0.0;
                        Pzz = 0.0;

                        for (a = 0; a < 19; a++)
                        {
                            rho += ff[fineindex(i, j, k, a)];
                            Pxx += ex[a] * ex[a] * ff[fineindex(i, j, k, a)];
                            Pxy += ex[a] * ey[a] * ff[fineindex(i, j, k, a)];
                            Pxz += ex[a] * ez[a] * ff[fineindex(i, j, k, a)];
                            Pyy += ey[a] * ey[a] * ff[fineindex(i, j, k, a)];
                            Pyz += ey[a] * ez[a] * ff[fineindex(i, j, k, a)];
                            Pzz += ez[a] * ez[a] * ff[fineindex(i, j, k, a)];
                        }

                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1) // solid
                            {
                                ff[fineindex(i, j, k, kb[a])] = w[a] * (rho + 3.0 * rho * (ubx * ex[kb[a]] + uby * ey[kb[a]] + ubz * ez[kb[a]]) + 4.5 * ((Pxx - rho / 3) * (ex[kb[a]] * ex[kb[a]] - 1.0 / 3.0) + 2 * Pxy * ex[kb[a]] * ey[kb[a]] + 2 * Pxz * ex[kb[a]] * ez[kb[a]] + 2 * Pyz * ey[kb[a]] * ez[kb[a]] + (Pyy - rho / 3) * (ey[kb[a]] * ey[kb[a]] - 1.0 / 3.0) + (Pzz - rho / 3) * (ez[kb[a]] * ez[kb[a]] - 1.0 / 3.0)));
                            }
                        }

                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 1) // solid
                            {

                                Fx += (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ex[a] * delxf * delxf;
                                Fy += (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ey[a] * delxf * delxf;
                                Fz += (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ez[a] * delxf * delxf;

                                fx = (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ex[a] * delxf * delxf;
                                fy = (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ey[a] * delxf * delxf;
                                fz = (ff[fineindex(i, j, k, kb[a])] + ftf[fineindex(i, j, k, a)]) * ez[a] * delxf * delxf;

                                Tx += ry * fz - rz * fy;
                                Ty += rz * fx - rx * fz;
                                Tz += rx * fy - ry * fx;
                            }
                        }
                    }
                }
            }
        }

        Fz = Fz + Fg;

        ucx[2] = ucx[1] + deltf * Fx / mass;
        ucy[2] = ucy[1] + deltf * Fy / mass;
        ucz[2] = ucz[1] + deltf * Fz / mass;

        omega_x[2] = omega_x[1] + deltf * Tx / Inertia;
        omega_y[2] = omega_y[1] + deltf * Ty / Inertia;
        omega_z[2] = omega_z[1] + deltf * Tz / Inertia;

        CX[2] = CX[1] + deltf * ucx[2];
        CY[2] = CY[1] + deltf * ucy[2];
        CZ[2] = CZ[1] + deltf * ucz[2];

        fprintf(out1, "%f\t%f\n", (ts + 0.5) * Ct, ucz[2] * Cl / Ct);
        fflush(out1);
        fprintf(out4, "%f\t%f\n", (ts + 0.5) * Ct, (CZ[2] - 15) / D);
        fflush(out4);
        fprintf(out5, "%f\t%f\n", (ts + 0.5) * Ct, Fz);
        fflush(out5);

        cx[1] = (CX[2] - CX_I + 1) * alpha;
        cy[1] = (CY[2] - CY_I + 1) * alpha;
        cz[1] = (CZ[2] - CZ_I + 1) * alpha;

        printf("%f Fx = %.8f cx = %.8f ucx = %.8f Fy = %.8f cy = %.8f ucy = %.8f Fz = %.8f cz = %.8f ucz = %.8f\n", ts + 0.5, Fx, cx[1], ucx[1], Fy, cy[1], ucy[1], Fz, cz[1], ucz[1]);
        fflush(stdout);

        //******* Locating any uncovered node by storing isn parameter in a new variable (isn1) *****

        for (i = 2; i <= nx - 2; i++)
        {
            for (j = 2; j <= ny - 2; j++)
            {
                for (k = 2; k <= nz - 2; k++)
                {
                    dist = sqrt(pow((i - cx[1]), 2) + pow((j - cy[1]), 2) + pow((k - cz[1]), 2));
                    if (dist <= (Df / 2.0))
                    {
                        isnf1[fineisnindex1(i, j, k)] = 1;
                    }

                    else
                    {
                        isnf1[fineisnindex1(i, j, k)] = 0;
                    }
                }
            }
        }

        for (i = 2; i <= nx - 2; i++)
        {
            for (j = 2; j <= ny - 2; j++)
            {
                for (k = 2; k <= nz - 2; k++)
                {
                    if (isnf1[fineisnindex1(i, j, k)] == 0 && isnf[fineisnindex1(i, j, k)] == 1) // new uncovered node from solid to fluid
                    {

                        neb = 0; // count of neigbouring fluid nodes
                        rho1 = 0.0;
                        ux = 0.0;
                        uy = 0.0;
                        uz = 0.0;
                        Pxx = 0.0;
                        Pyy = 0.0;
                        Pxy = 0.0;
                        Pxz = 0.0;
                        Pyz = 0.0;
                        Pzz = 0.0;

                        for (a = 0; a < 19; a++)
                        {
                            ia = i + ex[a];
                            ja = j + ey[a];
                            ka = k + ez[a];

                            if (isnf[fineisnindex1(ia, ja, ka)] == 0) // neighbouring old fluid nodes
                            {
                                neb++;
                                for (a = 0; a < 19; a++)
                                {
                                    rho1 += ff[fineindex(ia, ja, ka, a)];
                                    ux += ex[a] * ff[fineindex(ia, ja, ka, a)];
                                    uy += ey[a] * ff[fineindex(ia, ja, ka, a)];
                                    uz += ez[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pxx += ex[a] * ex[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pxy += ex[a] * ey[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pxz += ex[a] * ez[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pyy += ey[a] * ey[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pyz += ey[a] * ez[a] * ff[fineindex(ia, ja, ka, a)];
                                    Pzz += ez[a] * ez[a] * ff[fineindex(ia, ja, ka, a)];
                                }
                            }
                        }

                        ux = ux / rho1;
                        uy = uy / rho1;
                        uz = uz / rho1;
                        rho1 = rho1 / neb;
                        Pxx = Pxx / neb;
                        Pxy = Pxy / neb;
                        Pxz = Pxz / neb;
                        Pyy = Pyy / neb;
                        Pyz = Pyz / neb;
                        Pzz = Pzz / neb;

                        for (a = 0; a < 19; a++)
                        {
                            ff[fineindex(i, j, k, a)] = w[a] * (rho1 + 3.0 * rho1 * (ux * ex[a] + uy * ey[a] + uz * ez[a]) + 4.5 * ((Pxx - rho1 / 3) * (ex[a] * ex[a] - 1.0 / 3.0) + 2 * Pxy * ex[a] * ey[a] + 2 * Pxz * ex[a] * ez[a] + 2 * Pyz * ey[a] * ez[a] + (Pyy - rho1 / 3) * (ey[a] * ey[a] - 1.0 / 3.0) + (Pzz - rho1 / 3) * (ez[a] * ez[a] - 1.0 / 3.0))); // uncovered node
                        }
                    }
                }
            }
        }

        CX[1] = CX[2];
        CY[1] = CY[2];
        CZ[1] = CZ[2];

        ucx[1] = ucx[2];
        ucy[1] = ucy[2];
        ucz[1] = ucz[2];

        omega_x[1] = omega_x[2];
        omega_y[1] = omega_y[2];
        omega_z[1] = omega_z[2];

        //******* transfer of ff[a][n+1][] to fc[a][n+1][] on coarse block boundary ******
        dist_transfer_fc(ff, fc, CX_I, CX_F, CZ_I, CZ_F, CY_I, 1, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_fc(ff, fc, CY_I, CY_F, CZ_I, CZ_F, CX_F, 2, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_fc(ff, fc, CX_I, CX_F, CZ_I, CZ_F, CY_F, 3, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_fc(ff, fc, CY_I, CY_F, CZ_I, CZ_F, CX_I, 4, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_fc(ff, fc, CX_I, CX_F, CY_I, CY_F, CZ_I, 5, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_fc(ff, fc, CX_I, CX_F, CY_I, CY_F, CZ_F, 6, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

        //******* transfer of fc[a][n+1][] to ff[a][n+1][] on fine block boundary *******
        dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CZ_I - 1, CZ_F + 1, CY_I - 1, 1, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_cf(ff, fc, CY_I - 1, CY_F + 1, CZ_I - 1, CZ_F + 1, CX_F + 1, 2, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CZ_I - 1, CZ_F + 1, CY_F + 1, 3, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_cf(ff, fc, CY_I - 1, CY_F + 1, CZ_I - 1, CZ_F + 1, CX_I - 1, 4, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CY_I - 1, CY_F + 1, CZ_I - 1, 5, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CY_I - 1, CY_F + 1, CZ_F + 1, 6, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        //*********Spatial interpolation to find the fine block points at fine block boundary **********
        interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 0, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 0, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 0, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 0, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 0, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 0, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);

        //******* transfer of fc[a][n+1][] to ff[a][n+1][] on coarse block at (j=-2, i=0, i=nx+2 & j=ny+2) for 4-point langranginan interpolation *******
        transfer_mid_boundary(fc, ff, CX_I, CX_F, CZ_I, CZ_F, 1, 2, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        transfer_mid_boundary(fc, ff, CY_I, CY_F, CZ_I, CZ_F, 2, 2, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        transfer_mid_boundary(fc, ff, CX_I, CX_F, CZ_I, CZ_F, 3, 2, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        transfer_mid_boundary(fc, ff, CY_I, CY_F, CZ_I, CZ_F, 4, 2, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        transfer_mid_boundary(fc, ff, CX_I, CX_F, CY_I, CY_F, 5, 2, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        transfer_mid_boundary(fc, ff, CX_I, CX_F, CY_I, CY_F, 6, 2, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

        //*********Spatial interpolation to find the fine block points at fine block mid boundary **********
        interpolation_mid_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_mid_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx - 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_mid_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny - 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_mid_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_mid_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);
        interpolation_mid_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz - 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx - 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny - 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz - 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
        // printf("Check7\n");
        x_cen1 = CX[1] - x_cen0;
        y_cen1 = CY[1] - y_cen0;
        z_cen1 = CZ[1] - z_cen0;

        if (CZ_I >= 3.0)
        {
            //********* shifting of finer grid in X direction **********
            if (fabs(x_cen1 - x_cen2) >= delxc)
            {
                if (x_cen1 > x_cen2)
                {
                    CX_I++;
                    CX_F++;
                    x_cen2++;

                    for (i = 0; i <= nx - alpha; i++)
                    {
                        for (j = 0; j <= ny; j++)
                        {
                            for (k = 0; k <= nz; k++)
                            {
                                for (a = 0; a < 19; a++)
                                {
                                    ia = i + alpha;
                                    ff[fineindex(i, j, k, a)] = ff[fineindex(ia, j, k, a)]; // ff[a][i + alpha][j][k];
                                }
                            }
                        }
                    }
                    //********** transfer of ff[a][n+1][] to fc[a][n+1][] on a newly created coarse block boundary ***********
                    dist_transfer_fc(ff, fc, CY_I, CY_F, CZ_I, CZ_F, CX_I, 4, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //********* transfer of fc[a][n+1][] on the newly added fine block boundary ********
                    dist_transfer_cf(ff, fc, CY_I - 1, CY_F + 1, CZ_I - 1, CZ_F + 1, CX_F + 1, 2, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*******Spatial Interpolation to find the fine block points at the find block boundary*********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);

                    //********* transfer of fc[a][n+1][] on the newly added i=nx-1 ********
                    transfer_mid_boundary(fc, ff, CY_I - 1, CY_F + 1, CZ_I - 1, CZ_F + 1, 2, 0, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*******Spatial Interpolation to find the fine block points at the find block mid boundary*********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx - 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 2, nx - 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
                }

                if (x_cen1 < x_cen2)
                {
                    CX_I--;
                    CX_F--;
                    x_cen2--;

                    for (i = nx; i >= alpha; i--)
                    {
                        for (j = 0; j <= ny; j++)
                        {
                            for (k = 0; k <= nz; k++)
                            {
                                for (a = 0; a < 19; a++)
                                {
                                    ia = i - alpha;
                                    ff[fineindex(i, j, k, a)] = ff[fineindex(ia, j, k, a)];
                                }
                            }
                        }
                    }

                    //********** transfer of ff[a][n+1][] to fc[a][n+1][] on a newly created coarse block boundary ***********
                    dist_transfer_fc(ff, fc, CY_I, CY_F, CZ_I, CZ_F, CX_F, 2, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //********** transfer of fc[a][n+1][] to ff[a][n+1][] on a newly created fine block boundary ***********
                    dist_transfer_cf(ff, fc, CY_I - 1, CY_F + 1, CZ_I - 1, CZ_F + 1, CX_I - 1, 4, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*******Calculation of microscopic quantities and non-equilibrium distribution for interpolation*********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 0, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 0, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);

                    //*****Interpolation on fine grids***
                    transfer_mid_boundary(fc, ff, CY_I - 1, CY_F + 1, CZ_I - 1, CZ_F + 1, 4, 0, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //********** transfer of fc[a][n+1][] to ff[a][n+1][] on a newly created i=1 ***********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 4, 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
                }
            }

            //********* shifting of finer grid in Y direction **********
            if (fabs(y_cen1 - y_cen2) >= delyc)
            {
                if (y_cen1 > y_cen2)
                {
                    CY_I++;
                    CY_F++;
                    y_cen2++;

                    for (j = 0; j <= ny - alpha; j++)
                    {
                        for (i = 0; i <= nx; i++)
                        {
                            for (k = 0; k <= nz; k++)
                            {
                                for (a = 0; a < 19; a++)
                                {
                                    ja = j + alpha;
                                    ff[fineindex(i, j, k, a)] = ff[fineindex(i, ja, k, a)];
                                }
                            }
                        }
                    }

                    //******transfer of ff[a][n+1][] to fc[a][n+1][] on coarse block boundary ******
                    dist_transfer_fc(ff, fc, CX_I, CX_F, CZ_I, CZ_F, CY_I, 1, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //******transfer of fc[a][n+1][] to ff[a][n+1][] on fine block boundary ******
                    dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CZ_I - 1, CZ_F + 1, CY_F + 1, 3, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //***Calculation of macroscopic variables at fine points on j=ny ********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);

                    //*****Interpolation to find coarse point distribution at j=ny-1 *********
                    transfer_mid_boundary(fc, ff, CX_I - 1, CX_F + 1, CZ_I - 1, CZ_F + 1, 3, 0, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*******Interpolation of fine points at j=ny-1 *********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny - 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 3, ny - 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
                }

                if (y_cen1 < y_cen2)
                {
                    CY_I--;
                    CY_F--;
                    y_cen2--;

                    for (j = ny; j >= alpha; j--)
                    {
                        for (i = 0; i <= nx; i++)
                        {
                            for (k = 0; k <= nz; k++)
                            {
                                for (a = 0; a < 19; a++)
                                {
                                    ja = j - alpha;
                                    ff[fineindex(i, j, k, a)] = ff[fineindex(i, ja, k, a)];
                                }
                            }
                        }
                    }

                    //******* transfer of ff[a][n+1][] to fc[a][n+1][] on newly formed coarse block boundary ******
                    dist_transfer_fc(ff, fc, CX_I, CX_F, CZ_I, CZ_F, CY_F, 3, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    // ******* transfer of fc[a][n+1][] to ff[a][n+1][] on fine block boundary *******
                    dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CZ_I - 1, CZ_F + 1, CY_I - 1, 1, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*********Spatial interpolation to find the fine block points at fine block boundary **********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 0, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 0, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);

                    //******** Finding ff[a][i][j] at coarse points at j=1 ********
                    transfer_mid_boundary(fc, ff, CX_I - 1, CX_F + 1, CZ_I - 1, CZ_F + 1, 1, 0, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*********Spatial interpolation to find the fine block points at j=1 **********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 1, 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
                }
            }

            //********* shifting of finer grid in Y direction **********
            if (fabs(z_cen1 - z_cen2) >= delzc)
            {
                if (z_cen1 > z_cen2)
                {
                    CZ_I++;
                    CZ_F++;
                    z_cen2++;

                    for (k = 0; k <= nz - alpha; k++)
                    {
                        for (i = 0; i <= nx; i++)
                        {
                            for (j = 0; j <= ny; j++)
                            {
                                for (a = 0; a < 19; a++)
                                {
                                    ka = k + alpha;
                                    ff[fineindex(i, j, k, a)] = ff[fineindex(i, j, ka, a)];
                                }
                            }
                        }
                    }

                    //******* transfer of ff[a][n+1][] to fc[a][n+1][] on newly formed coarse block boundary ******
                    dist_transfer_fc(ff, fc, CX_I, CX_F, CY_I, CY_F, CZ_I, 5, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    // ******* transfer of fc[a][n+1][] to ff[a][n+1][] on fine block boundary *******
                    dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CY_I - 1, CY_F + 1, CZ_F + 1, 6, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*********Spatial interpolation to find the fine block points at fine block boundary **********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);

                    //******** Finding ff[a][i][j] at coarse points at j=1 ********
                    transfer_mid_boundary(fc, ff, CX_I - 1, CX_F + 1, CY_I - 1, CY_F + 1, 6, 0, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*********Spatial interpolation to find the fine block points at j=1 **********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz - 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 6, nz - 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
                }

                if (z_cen1 < z_cen2)
                {
                    CZ_I--;
                    CZ_F--;
                    z_cen2--;

                    for (k = nz; k >= alpha; k--)
                    {
                        for (i = 0; i <= nx; i++)
                        {
                            for (j = 0; j <= ny; j++)
                            {
                                for (a = 0; a < 19; a++)
                                {
                                    ka = k - alpha;
                                    ff[fineindex(i, j, k, a)] = ff[fineindex(i, j, ka, a)];
                                }
                            }
                        }
                    }

                    printf("CZ_I=%d\n", CZ_I);
                    //******* transfer of ff[a][n+1][] to fc[a][n+1][] on newly formed coarse block boundary ******
                    dist_transfer_fc(ff, fc, CX_I, CX_F, CY_I, CY_F, CZ_F, 6, ex, ey, ez, rhof, jxf, jyf, jzf, feqf, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    // ******* transfer of fc[a][n+1][] to ff[a][n+1][] on fine block boundary *******
                    dist_transfer_cf(ff, fc, CX_I - 1, CX_F + 1, CY_I - 1, CY_F + 1, CZ_I - 1, 5, ex, ey, ez, rhoc, jxc, jyc, jzc, feqc, w, tauc, tauf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*********Spatial interpolation to find the fine block points at fine block boundary **********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 0, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 0, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);

                    //******** Finding ff[a][i][j] at coarse points at j=1 ********
                    transfer_mid_boundary(fc, ff, CX_I - 1, CX_F + 1, CY_I - 1, CY_F + 1, 5, 0, feqc, feqf, jxc, jyc, jzc, jxf, jyf, jzf, rhoc, rhof, tauc, tauf, ex, ey, ez, w, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    //*********Spatial interpolation to find the fine block points at j=1 **********
                    interpolation_fine_boundary(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 1, ex, ey, ez, w, rhof, jxf, jyf, jzf, alpha, CX_I, CX_F, CY_I, CY_F, CZ_I, CZ_F);

                    bicubic_interpolation(ff, feqf, fneqf_b, fneqf_r, fneqf_t, fneqf_l, fneqf_u, fneqf_d, 5, 1, jxf, jyf, jzf, rhof, ex, ey, ez, w, alpha);
                }
            }
        }
        printf("CX[1]=%f\tCY[1]=%f\tCZ[1]=%f\n", CX[1], CY[1], CZ[1]);
        cx[1] = (CX[1] - CX_I + 1) * alpha;
        cy[1] = (CY[1] - CY_I + 1) * alpha;
        cz[1] = (CZ[1] - CZ_I + 1) * alpha;

        for (i = 0; i <= nx; i++)
        {
            for (j = 0; j <= ny; j++)
            {
                for (k = 0; k <= nz; k++)
                {
                    dist = sqrt(pow((i - cx[1]), 2) + pow((j - cy[1]), 2) + pow((k - cz[1]), 2));
                    if (dist <= (Df / 2.0))
                    {
                        isnf[fineisnindex1(i, j, k)] = 1;
                    }
                    else
                    {
                        isnf[fineisnindex1(i, j, k)] = 0;
                    }
                }
            }
        }

        // ****** calculating physical variables ux[n+1],uy[n+1],rho[n+1] and collision to get ftf[a][n+1][] in fine block ******
        count = 0;
        rho_av = 0.0;
        for (i = 0; i <= nx; i++)
        {
            for (j = 0; j <= ny; j++)
            {
                for (k = 0; k <= nz; k++)
                {
                    if (isnf[fineisnindex1(i, j, k)] == 0)
                    {
                        for (a = 0; a < 19; a++)
                        {
                            rf[a] = 0.0;
                            for (v = 0; v < 19; v++)
                            {
                                rf[a] += M[a][v] * ff[fineindex(i, j, k, v)];
                            }
                        }

                        rhof[fineisnindex1(i, j, k)] = rf[0];
                        jxf[fineisnindex1(i, j, k)] = rf[3];
                        jyf[fineisnindex1(i, j, k)] = rf[5];
                        jzf[fineisnindex1(i, j, k)] = rf[7];

                        reqf[0] = 0.0;
                        reqf[1] = -11.0 * rhof[fineisnindex1(i, j, k)] + 19.0 * (jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[2] = 3.0 * rhof[fineisnindex1(i, j, k)] - 11.0 * (jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / (2.0 * rho0);
                        reqf[3] = 0.0;
                        reqf[4] = -(2.0 / 3.0) * jxf[fineisnindex1(i, j, k)];
                        reqf[5] = 0.0;
                        reqf[6] = -(2.0 / 3.0) * jyf[fineisnindex1(i, j, k)];
                        reqf[7] = 0.0;
                        reqf[8] = -(2.0 / 3.0) * jzf[fineisnindex1(i, j, k)];
                        reqf[9] = (3.0 * jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] - (jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)])) / rho0;
                        reqf[10] = -0.5 * reqf[9]; //-(3.0*jx[i][j][k]*jx[i][j][k] - (jx[i][j][k]*jx[i][j][k] + jy[i][j][k]*jy[i][j][k] + jz[i][j][k]*jz[i][j][k]))/(2.0*rho0);
                        reqf[11] = (jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] - jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[12] = -0.5 * reqf[11]; //-(jy[i][j][k]*jy[i][j][k] - jz[i][j][k]*jz[i][j][k])/(2.0*rho0);
                        reqf[13] = (jxf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[14] = (jyf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[15] = (jxf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)]) / rho0;
                        reqf[16] = 0.0;
                        reqf[17] = 0.0;
                        reqf[18] = 0.0;

                        rho_av += rhof[fineisnindex1(i, j, k)];
                        count++;

                        for (a = 0; a < 19; a++)
                        {
                            rf[a] = rf[a] - sf[a] * (rf[a] - reqf[a]); // collision
                        }

                        for (a = 0; a < 19; a++)
                        {
                            ftf[fineindex(i, j, k, a)] = 0.0;

                            for (v = 0; v < 19; v++)
                            {
                                ftf[fineindex(i, j, k, a)] += Minv[a][v] * rf[v];
                            }
                        }
                    }
                }
            }
        }

        //     printf("Average density fine block = %f\tCount=%d\n",rho_av/count,count);

        count = 0;
        rho_av = 0.0;
        // ******** Calculating physical variables ux[n],uy[n],rho[n] and collision to get ftc[a][n+1][] in the coarse block *******
        for (I = 1; I <= NX; I++)
        {
            for (J = 1; J <= NY; J++)
            {
                for (K = 1; K <= NZ; K++)
                {
                    for (a = 0; a < 19; a++)
                    {
                        rc[a] = 0.0;
                        for (v = 0; v < 19; v++)
                        {
                            rc[a] += M[a][v] * fc[coarseindex(I, J, K, v)];
                        }
                    }

                    rhoc[coarseisnindex1(I, J, K)] = rc[0];
                    jxc[coarseisnindex1(I, J, K)] = rc[3];
                    jyc[coarseisnindex1(I, J, K)] = rc[5];
                    jzc[coarseisnindex1(I, J, K)] = rc[7];

                    reqc[0] = 0.0;
                    reqc[1] = -11.0 * rhoc[coarseisnindex1(I, J, K)] + 19.0 * (jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)]) / rho0;
                    reqc[2] = 3.0 * rhoc[coarseisnindex1(I, J, K)] - 11.0 * (jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)]) / (2.0 * rho0);
                    reqc[3] = 0.0;
                    reqc[4] = -(2.0 / 3.0) * jxc[coarseisnindex1(I, J, K)];
                    reqc[5] = 0.0;
                    reqc[6] = -(2.0 / 3.0) * jyc[coarseisnindex1(I, J, K)];
                    reqc[7] = 0.0;
                    reqc[8] = -(2.0 / 3.0) * jzc[coarseisnindex1(I, J, K)];
                    reqc[9] = (3.0 * jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] - (jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)])) / rho0;
                    reqc[10] = -0.5 * reqc[9]; //-(3.0*jx[i][j][k]*jx[i][j][k] - (jx[i][j][k]*jx[i][j][k] + jy[i][j][k]*jy[i][j][k] + jz[i][j][k]*jz[i][j][k]))/(2.0*rho0);
                    reqc[11] = (jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] - jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)]) / rho0;
                    reqc[12] = -0.5 * reqc[11]; //-(jy[i][j][k]*jy[i][j][k] - jz[i][j][k]*jz[i][j][k])/(2.0*rho0);
                    reqc[13] = (jxc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)]) / rho0;
                    reqc[14] = (jyc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)]) / rho0;
                    reqc[15] = (jxc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)]) / rho0;
                    reqc[16] = 0.0;
                    reqc[17] = 0.0;
                    reqc[18] = 0.0;

                    rho_av += rhoc[coarseisnindex1(I, J, K)];
                    count++;

                    for (a = 0; a < 19; a++)
                    {
                        rc[a] = rc[a] - sc[a] * (rc[a] - reqc[a]); // collision
                    }

                    for (a = 0; a < 19; a++)
                    {
                        ftc[coarseindex(I, J, K, a)] = 0.0;

                        for (v = 0; v < 19; v++)
                        {
                            ftc[coarseindex(I, J, K, a)] += Minv[a][v] * rc[v];
                        }
                    }
                }
            }
        }
        // printf("Average density coarse block = %f\tCount=%d\n",rho_av/count,count);
        /*
if(solnumber>=1)
 {
   solnumber++;
   strcpy(str,sol);
   sprintf(sol1,"%d",solnumber);
   strcat(str,sol1);
   strcat(str,str1);
   out3 = fopen(str,"w");

    fprintf(out3,"VARIABLES = X, Y, Z, ISN, RHO, UX, UY, UZ\n");

    fprintf(out3,"ZONE I = %d, J= %d, K= %d\n",NX,NY,NZ);
    for(K=1;K<=NZ;K++){
      for(J=1;J<=NY;J++){
         for(I=1;I<=NX;I++)
               {
                rhoc[coarseisnindex1(I ,J, K)]=0.0;
                jxc[coarseisnindex1(I ,J, K)]=0.0;
                jyc[coarseisnindex1(I ,J, K)]=0.0;
                jzc[coarseisnindex1(I ,J, K)]=0.0;

            for(a=0;a<19;a++)
            {
                rhoc[coarseisnindex1(I ,J, K)]+=fc[coarseindex(I, J, K,a)];
                jxc[coarseisnindex1(I ,J, K)]+=ex[a]*fc[coarseindex(I, J, K,a)];
                jyc[coarseisnindex1(I ,J, K)]+=ey[a]*fc[coarseindex(I, J, K,a)];
                jzc[coarseisnindex1(I ,J, K)]+=ez[a]*fc[coarseindex(I, J, K,a)];
            }
                  fprintf(out3,"%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n",I,J,K,isnc[I][J][K],rhoc[coarseisnindex1(I ,J, K)],jxc[coarseisnindex1(I ,J, K)],jyc[coarseisnindex1(I ,J, K)],jzc[coarseisnindex1(I ,J, K)]);
            }
           }
       }

    fprintf(out3,"ZONE I =%d, J=%d, K=%d\n",nx+1,ny+1,nz+1);
    for(k=0,distz=CZ_I-1;k<=nz;k++,distz+=0.5){
      for (j=0,disty=CY_I-1;j<=ny;j++,disty+=0.5){
        for(i=0,distx=CX_I-1;i<=nx;i++,distx+=0.5)
         {
          rhof[fineisnindex1(i, j, k)]=0.0;
          jxf[fineisnindex1(i, j, k)]=0.0;
          jyf[fineisnindex1(i, j, k)]=0.0;
          jzf[fineisnindex1(i, j, k)]=0.0;
          if(isnf1[i][j][k]==0)
          {
            for(a=0;a<19;a++)
                {
                rhof[fineisnindex1(i, j, k)]+=ff[fineindex(i, j, k,a)];
                jxf[fineisnindex1(i, j, k)]+=ex[a]*ff[fineindex(i, j, k,a)];
                jyf[fineisnindex1(i, j, k)]+=ey[a]*ff[fineindex(i, j, k,a)];
                jzf[fineisnindex1(i, j, k)]+=ez[a]*ff[fineindex(i, j, k,a)];
            }
          }
            fprintf(out3,"%f\t%f\t%f\t%d\t%f\t%f\t%f\t%f\n",distx,disty,distz,isnf1[i][j][k],rhof[fineisnindex1(i, j, k)],jxf[fineisnindex1(i, j, k)],jyf[fineisnindex1(i, j, k)],jzf[fineisnindex1(i, j, k)]);
         }
        }
    }

fclose(out3);
} */

        ts++;
    } while (CZ[1] >= 18.0);

    /*   end = clock();
       cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
       fprintf(time, "CPU_time_used=%f\n", cpu_time_used);
       fflush(time);
   */
    //   fclose(time);
    fclose(out1);
    fclose(out4);
    fclose(out5);
    return 0;
}

void dist_transfer_fc(double *ff, double *fc, int Start1, int End1, int Start2, int End2, int face, int side, int ex[19], int ey[19], int ez[19], double *rhof, double *jxf, double *jyf, double *jzf, double *feqf, double w[19], double tauc, double tauf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F)
{
    int I, J, K, i, j, k, a;
    double l, m, n, p, u2, t;

    if (side == 1)
    {
        for (J = face, j = 2, I = Start1, i = 2; I <= End1; I++, i += alpha)
        {
            for (K = Start2, k = 2; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    fc[coarseindex(I, J, K, a)] = feqf[fineindex(i, j, k, a)] + alpha * (tauc / tauf) * (ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]);
                }
            }
        }
    }
    if (side == 2)
    {
        for (J = Start1, j = 2, I = face, i = nx - 2; J <= End1; J++, j += alpha)
        {
            for (K = Start2, k = 2; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    fc[coarseindex(I, J, K, a)] = feqf[fineindex(i, j, k, a)] + alpha * (tauc / tauf) * (ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]);
                }
            }
        }
    }
    if (side == 3)
    {
        for (J = face, j = ny - 2, I = Start1, i = 2; I <= End1; I++, i += alpha)
        {
            for (K = Start2, k = 2; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    fc[coarseindex(I, J, K, a)] = feqf[fineindex(i, j, k, a)] + alpha * (tauc / tauf) * (ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]);
                }
            }
        }
    }
    if (side == 4)
    {
        for (J = Start1, j = 2, I = face, i = 2; J <= End1; J++, j += alpha)
        {
            for (K = Start2, k = 2; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    fc[coarseindex(I, J, K, a)] = feqf[fineindex(i, j, k, a)] + alpha * (tauc / tauf) * (ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]);
                }
            }
        }
    }
    if (side == 5)
    {
        for (K = face, k = 2, I = Start1, i = 2; I <= End1; I++, i += alpha)
        {
            for (J = Start2, j = 2; J <= End2; J++, j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    fc[coarseindex(I, J, K, a)] = feqf[fineindex(i, j, k, a)] + alpha * (tauc / tauf) * (ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]);
                }
            }
        }
    }
    if (side == 6)
    {
        for (K = face, k = nz - 2, I = Start1, i = 2; I <= End1; I++, i += alpha)
        {
            for (J = Start2, j = 2; J <= End2; J++, j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    fc[coarseindex(I, J, K, a)] = feqf[fineindex(i, j, k, a)] + alpha * (tauc / tauf) * (ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]);
                }
            }
        }
    }
}

void dist_transfer_cf(double *ff, double *fc, int Start1, int End1, int Start2, int End2, int face, int side, int ex[19], int ey[19], int ez[19], double *rhoc, double *jxc, double *jyc, double *jzc, double *feqc, double w[19], double tauc, double tauf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F)
{
    int I, J, K, i, j, k, a;
    double l, m, n, p, u2, t;

    if (side == 1)
    {
        for (J = face, j = 0, I = Start1, i = 0; I <= End1; I++, i += alpha)
        {
            for (K = Start2, k = 0; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    ff[fineindex(i, j, k, a)] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }
            }
        }
    }
    if (side == 4)
    {
        for (I = face, i = 0, J = Start1, j = 0; J <= End1; J++, j += alpha)
        {
            for (K = Start2, k = 0; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    ff[fineindex(i, j, k, a)] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }
            }
        }
    }
    if (side == 3)
    {
        for (J = face, j = ny, I = Start1, i = 0; I <= End1; I++, i += alpha)
        {
            for (K = Start2, k = 0; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    ff[fineindex(i, j, k, a)] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }
            }
        }
    }
    if (side == 2)
    {
        for (I = face, i = nx, J = Start1, j = 0; J <= End1; J++, j += alpha)
        {
            for (K = Start2, k = 0; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    ff[fineindex(i, j, k, a)] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }
            }
        }
    }
    if (side == 5)
    {
        for (I = Start1, i = 0, K = face, k = 0; I <= End1; I++, i += alpha)
        {
            for (J = Start2, j = 0; J <= End2; J++, j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    ff[fineindex(i, j, k, a)] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }
            }
        }
    }
    if (side == 6)
    {
        for (I = Start1, i = 0, K = face, k = nz; I <= End1; I++, i += alpha)
        {
            for (J = Start2, j = 0; J <= End2; J++, j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);

                    ff[fineindex(i, j, k, a)] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }
            }
        }
    }
}

void interpolation_fine_boundary(double *ff, double *feqf, double *fneqf_b, double *fneqf_r, double *fneqf_t, double *fneqf_l, double *fneqf_u, double *fneqf_d, int side, int face, int ex[19], int ey[19], int ez[19], double w[19], double *rhof, double *jxf, double *jyf, double *jzf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F)
{
    int I, J, K, i, j, k, a;
    double l, m, n, p, u2, t;

    if (side == 1)
    {
        // for (i = 0, j = face; i <= nx; i += alpha)
        // {
        //     for (k = 0; k <= nz; k += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_b[topbottom(i, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (j = face, k = 0; k <= nz; k += alpha)
        {
            for (i = 3; i <= nx - 3; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_b[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_b[topbottom(i + 1, k, a)] + fneqf_b[topbottom(i - 1, k, a)]) - (1.0 / 16.0) * (fneqf_b[topbottom(i + 3, k, a)] + fneqf_b[topbottom(i - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (j = face, i = 0; i <= nx; i += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_b[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_b[topbottom(i, k + 1, a)] + fneqf_b[topbottom(i, k - 1, a)]) - (1.0 / 16.0) * (fneqf_b[topbottom(i, k + 3, a)] + fneqf_b[topbottom(i, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }

        //****Corner points topbottom side
        for (k = 0; k <= nz; k += alpha)
        {
            i = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i - 1, k, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i + 1, k, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            i = nx - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i + 1, k, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i - 1, k, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }

        for (i = 0; i <= nx; i += alpha)
        {
            k = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i, k - 1, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i, k + 1, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            k = nz - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i, k + 1, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i, k - 1, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 2)
    {
        // for (k = 0, i = face; k <= nz; k += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_r[leftright(j, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (i = face, k = 0; k <= nz; k += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_r[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_r[leftright(j + 1, k, a)] + fneqf_r[leftright(j - 1, k, a)]) - (1.0 / 16.0) * (fneqf_r[leftright(j + 3, k, a)] + fneqf_r[leftright(j - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }

        for (i = face, j = 0; j <= ny; j += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_r[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_r[leftright(j, k + 1, a)] + fneqf_r[leftright(j, k - 1, a)]) - (1.0 / 16.0) * (fneqf_r[leftright(j, k + 3, a)] + fneqf_r[leftright(j, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }

        //******Corner points leftright Side
        for (k = 0; k <= nz; k += alpha)
        {
            i = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j - 1, k, a)] + (3.0 / 4.0) * fneqf_r[leftright(j + 1, k, a)] - (1.0 / 8.0) * fneqf_r[leftright(j + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            i = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j + 1, k, a)] + (3.0 / 4.0) * fneqf_r[leftright(j - 1, k, a)] - (1.0 / 8.0) * fneqf_r[leftright(j - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }
        for (j = 0; j <= ny; j += alpha)
        {
            i = face;
            k = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j, k - 1, a)] + (3.0 / 4.0) * fneqf_r[leftright(j, k + 1, a)] - (1.0 / 8.0) * fneqf_r[leftright(j, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            i = face;
            k = nz - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j, k + 1, a)] + (3.0 / 4.0) * fneqf_r[leftright(j, k - 1, a)] - (1.0 / 8.0) * fneqf_r[leftright(j, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];
                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 3)
    {
        // for (i = 0, j = face; i <= nx; i += alpha)
        // {
        //     for (k = 0; k <= nz; k += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_t[topbottom(i, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (j = face, k = 0; k <= nz; k += alpha)
        {
            for (i = 3; i <= nx - 3; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_t[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_t[topbottom(i + 1, k, a)] + fneqf_t[topbottom(i - 1, k, a)]) - (1.0 / 16.0) * (fneqf_t[topbottom(i + 3, k, a)] + fneqf_t[topbottom(i - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (j = face, i = 0; i <= nx; i += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_t[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_t[topbottom(i, k + 1, a)] + fneqf_t[topbottom(i, k - 1, a)]) - (1.0 / 16.0) * (fneqf_t[topbottom(i, k + 3, a)] + fneqf_t[topbottom(i, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }

        //********Top side
        for (k = 0; k <= nz; k += alpha)
        {
            i = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i - 1, k, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i + 1, k, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            i = nx - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i + 1, k, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i - 1, k, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }

        for (i = 0; i <= nx; i += alpha)
        {
            k = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i, k - 1, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i, k + 1, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            k = nz - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i, k + 1, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i, k - 1, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 4)
    {
        // for (k = 0, i = face; k <= nz; k += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_l[leftright(j, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (i = face, k = 0; k <= nz; k += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_l[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_l[leftright(j + 1, k, a)] + fneqf_l[leftright(j - 1, k, a)]) - (1.0 / 16.0) * (fneqf_l[leftright(j + 3, k, a)] + fneqf_l[leftright(j - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }

        for (i = face, j = 0; j <= ny; j += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_l[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_l[leftright(j, k + 1, a)] + fneqf_l[leftright(j, k - 1, a)]) - (1.0 / 16.0) * (fneqf_l[leftright(j, k + 3, a)] + fneqf_l[leftright(j, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }
        //********Corner points Left side
        for (k = 0; k <= nz; k += alpha)
        {
            i = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j - 1, k, a)] + (3.0 / 4.0) * fneqf_l[leftright(j + 1, k, a)] - (1.0 / 8.0) * fneqf_l[leftright(j + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            i = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j + 1, k, a)] + (3.0 / 4.0) * fneqf_l[leftright(j - 1, k, a)] - (1.0 / 8.0) * fneqf_l[leftright(j - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }

        for (j = 0; j <= ny; j += alpha)
        {
            i = face;
            k = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j, k - 1, a)] + (3.0 / 4.0) * fneqf_l[leftright(j, k + 1, a)] - (1.0 / 8.0) * fneqf_l[leftright(j, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            i = face;
            k = nz - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j, k + 1, a)] + (3.0 / 4.0) * fneqf_l[leftright(j, k - 1, a)] - (1.0 / 8.0) * fneqf_l[leftright(j, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 5)
    {
        // for (i = 0, k = face; i <= nx; i += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_d[updown(i, j, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (k = face, j = 0; j <= ny; j += alpha)
        {
            for (i = 3; i <= nx - 3; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_d[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_d[updown(i + 1, j, a)] + fneqf_d[updown(i - 1, j, a)]) - (1.0 / 16.0) * (fneqf_d[updown(i + 3, j, a)] + fneqf_d[updown(i - 3, j, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (k = face, i = 0; i <= nx; i += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_d[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_d[updown(i, j + 1, a)] + fneqf_d[updown(i, j - 1, a)]) - (1.0 / 16.0) * (fneqf_d[updown(i, j + 3, a)] + fneqf_d[updown(i, j - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }
        //********Corner points Left side
        for (i = 0; i <= nx; i += alpha)
        {
            k = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i, j - 1, a)] + (3.0 / 4.0) * fneqf_d[updown(i, j + 1, a)] - (1.0 / 8.0) * fneqf_d[updown(i, j + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            k = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i, j + 1, a)] + (3.0 / 4.0) * fneqf_d[updown(i, j - 1, a)] - (1.0 / 8.0) * fneqf_d[updown(i, j - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }

        for (j = 0; j <= ny; j += alpha)
        {
            k = face;
            i = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i - 1, j, a)] + (3.0 / 4.0) * fneqf_d[updown(i + 1, j, a)] - (1.0 / 8.0) * fneqf_d[updown(i + 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            k = face;
            i = nx - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i + 1, j, a)] + (3.0 / 4.0) * fneqf_d[updown(i - 1, j, a)] - (1.0 / 8.0) * fneqf_d[updown(i - 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }
    }

    if (side == 6)
    {
        // for (i = 0, k = face; i <= nx; i += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_u[updown(i, j, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (k = face, j = 0; j <= ny; j += alpha)
        {
            for (i = 3; i <= nx - 3; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_u[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_u[updown(i + 1, j, a)] + fneqf_u[updown(i - 1, j, a)]) - (1.0 / 16.0) * (fneqf_u[updown(i + 3, j, a)] + fneqf_u[updown(i - 3, j, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (k = face, i = 0; i <= nx; i += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_u[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_u[updown(i, j + 1, a)] + fneqf_u[updown(i, j - 1, a)]) - (1.0 / 16.0) * (fneqf_u[updown(i, j + 3, a)] + fneqf_u[updown(i, j - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }
        //********Corner points Left side
        for (i = 0; i <= nx; i += alpha)
        {
            k = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i, j - 1, a)] + (3.0 / 4.0) * fneqf_u[updown(i, j + 1, a)] - (1.0 / 8.0) * fneqf_u[updown(i, j + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            k = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i, j + 1, a)] + (3.0 / 4.0) * fneqf_u[updown(i, j - 1, a)] - (1.0 / 8.0) * fneqf_u[updown(i, j - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }

        for (j = 0; j <= ny; j += alpha)
        {
            k = face;
            i = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i - 1, j, a)] + (3.0 / 4.0) * fneqf_u[updown(i + 1, j, a)] - (1.0 / 8.0) * fneqf_u[updown(i + 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            k = face;
            i = nx - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i + 1, j, a)] + (3.0 / 4.0) * fneqf_u[updown(i - 1, j, a)] - (1.0 / 8.0) * fneqf_u[updown(i - 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }
    }
}

void transfer_mid_boundary(double *fc, double *ff, int Start1, int End1, int Start2, int End2, int side, int start, double *feqc, double *feqf, double *jxc, double *jyc, double *jzc, double *jxf, double *jyf, double *jzf, double *rhoc, double *rhof, double tauc, double tauf, int ex[19], int ey[19], int ez[19], double w[19], int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F)
{
    int I, J, K, i, j, k, a, q;
    double l, m, n, p, u2, t, rho_l, jx_l, jy_l, jz_l;
    double dist_l[19], fneqf[19];
    double dist_neq_l[4][19];

    if (side == 1)
    {
        for (J = CY_I - 2, I = Start1, i = start; I <= End1; I++, i += alpha)
        {
            for (K = Start2, k = start; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    dist_l[a] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                    // dist_neq_l[0][a] = (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }

                // l = 0.0;
                // m = 0.0;
                // n = 0.0;
                // p = 0.0;

                // for (a = 0; a < 19; a++)
                // {
                //     l += dist_l[a];
                //     m += ex[a] * dist_l[a];
                //     n += ey[a] * dist_l[a];
                //     p += ez[a] * dist_l[a];
                // }
                // rho_l = l;
                // jx_l = m / rho_l;
                // jy_l = n / rho_l;
                // jz_l = p / rho_l;

                // //******Find the non equilibrium points for j=0,j=2,j=4
                // for (q = 1, j = 0; j <= 4; j += 2, q++)
                // {
                //     l = 0.0;
                //     m = 0.0;
                //     n = 0.0;
                //     p = 0.0;

                //     for (a = 0; a < 19; a++)
                //     {
                //         l += ff[fineindex(i, j, k, a)];
                //         m += ex[a] * ff[fineindex(i, j, k, a)];
                //         n += ey[a] * ff[fineindex(i, j, k, a)];
                //         p += ez[a] * ff[fineindex(i, j, k, a)];
                //     }
                //     rhof[fineisnindex1(i, j, k)] = l;
                //     jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                //     jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                //     jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                //     u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                //     for (a = 0; a < 19; a++)
                //     {
                //         t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                //         feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                //         dist_neq_l[q][a] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]; // Non equilibrium distribution for j=-2 point
                //     }
                // }

                // //*****4 point langrangian interpolation of the conserved quantities ****

                // jxf[fineisnindex1(i, 1, k)] = (-jx_l + 9.0 * jxf[fineisnindex1(i, 0, k)] + 9.0 * jxf[fineisnindex1(i, 2, k)] - jxf[fineisnindex1(i, 4, k)]) / 16.0;
                // jyf[fineisnindex1(i, 1, k)] = (-jy_l + 9.0 * jyf[fineisnindex1(i, 0, k)] + 9.0 * jyf[fineisnindex1(i, 2, k)] - jyf[fineisnindex1(i, 4, k)]) / 16.0;
                // jzf[fineisnindex1(i, 1, k)] = (-jz_l + 9.0 * jzf[fineisnindex1(i, 0, k)] + 9.0 * jzf[fineisnindex1(i, 2, k)] - jzf[fineisnindex1(i, 4, k)]) / 16.0;
                // rhof[fineisnindex1(i, 1, k)] = (-rho_l + 9.0 * rhof[fineisnindex1(i, 0, k)] + 9.0 * rhof[fineisnindex1(i, 2, k)] - rhof[fineisnindex1(i, 4, k)]) / 16.0;

                // u2 = jxf[fineisnindex1(i, 1, k)] * jxf[fineisnindex1(i, 1, k)] + jyf[fineisnindex1(i, 1, k)] * jyf[fineisnindex1(i, 1, k)] + jzf[fineisnindex1(i, 1, k)] * jzf[fineisnindex1(i, 1, k)];

                // //*****Equilibrium distribution at j=1 ********
                // for (a = 0; a < 19; a++)
                // {
                //     t = ex[a] * jxf[fineisnindex1(i, 1, k)] + ey[a] * jyf[fineisnindex1(i, 1, k)] + ez[a] * jzf[fineisnindex1(i, 1, k)];
                //     feqf[fineindex(i, 1, k, a)] = w[a] * rhof[fineisnindex1(i, 1, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // }

                // //*****Find Non-equilibrium distribution at j=1 using 4 point langrangian interpolation*******
                // for (a = 0; a < 19; a++)
                // {
                //     fneqf[a] = (-dist_neq_l[0][a] + 9.0 * dist_neq_l[1][a] + 9.0 * dist_neq_l[2][a] - dist_neq_l[3][a]) / 16.0;
                // }

                //*****Find f = feq+fneq******
                for (j = 1, a = 0; a < 19; a++)
                {
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf[a];

                    ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, 0, k, a)] + 9.0 * ff[fineindex(i, 2, k, a)] - ff[fineindex(i, 4, k, a)]) / 16.0;
                }
            }
        }
    }
    if (side == 2)
    {
        for (I = CX_F + 2, J = Start1, j = start; J <= End1; J++, j += alpha) // store for i=nx+2
        {
            for (K = Start2, k = start; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    dist_l[a] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                    // dist_neq_l[0][a] = (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }

                // l = 0.0;
                // m = 0.0;
                // n = 0.0;
                // p = 0.0;

                // for (a = 0; a < 19; a++)
                // {
                //     l += dist_l[a];
                //     m += ex[a] * dist_l[a];
                //     n += ey[a] * dist_l[a];
                //     p += ez[a] * dist_l[a];
                // }
                // rho_l = l;
                // jx_l = m / rho_l;
                // jy_l = n / rho_l;
                // jz_l = p / rho_l;

                // //******Find the non equilibrium points for i=nx, i=nx-2, i=nx-4
                // for (q = 1, i = nx; i >= nx - 4; i -= 2, q++)
                // {
                //     l = 0.0;
                //     m = 0.0;
                //     n = 0.0;
                //     p = 0.0;

                //     for (a = 0; a < 19; a++)
                //     {
                //         l += ff[fineindex(i, j, k, a)];
                //         m += ex[a] * ff[fineindex(i, j, k, a)];
                //         n += ey[a] * ff[fineindex(i, j, k, a)];
                //         p += ez[a] * ff[fineindex(i, j, k, a)];
                //     }
                //     rhof[fineisnindex1(i, j, k)] = l;
                //     jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                //     jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                //     jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                //     u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                //     for (a = 0; a < 19; a++)
                //     {
                //         t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                //         feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                //         dist_neq_l[q][a] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]; // Non equilibrium distribution for j=-2 point
                //     }
                // }

                // //*****4 point langrangian interpolation of the conserved quantities****

                // jxf[fineisnindex1(nx - 1, j, k)] = (-jx_l + 9.0 * jxf[fineisnindex1(nx, j, k)] + 9.0 * jxf[fineisnindex1(nx - 2, j, k)] - jxf[fineisnindex1(nx - 4, j, k)]) / 16.0;
                // jyf[fineisnindex1(nx - 1, j, k)] = (-jy_l + 9.0 * jyf[fineisnindex1(nx, j, k)] + 9.0 * jyf[fineisnindex1(nx - 2, j, k)] - jyf[fineisnindex1(nx - 4, j, k)]) / 16.0;
                // jzf[fineisnindex1(nx - 1, j, k)] = (-jz_l + 9.0 * jzf[fineisnindex1(nx, j, k)] + 9.0 * jzf[fineisnindex1(nx - 2, j, k)] - jzf[fineisnindex1(nx - 4, j, k)]) / 16.0;
                // rhof[fineisnindex1(nx - 1, j, k)] = (-rho_l + 9.0 * rhof[fineisnindex1(nx, j, k)] + 9.0 * rhof[fineisnindex1(nx - 2, j, k)] - rhof[fineisnindex1(nx - 4, j, k)]) / 16.0;

                // u2 = jxf[fineisnindex1(nx - 1, j, k)] * jxf[fineisnindex1(nx - 1, j, k)] + jyf[fineisnindex1(nx - 1, j, k)] * jyf[fineisnindex1(nx - 1, j, k)] + jzf[fineisnindex1(nx - 1, j, k)] * jzf[fineisnindex1(nx - 1, j, k)];

                // //****** Equilibrium distribution at i=nx-1*******
                // for (a = 0; a < 19; a++)
                // {
                //     t = ex[a] * jxf[fineisnindex1(nx - 1, j, k)] + ey[a] * jyf[fineisnindex1(nx - 1, j, k)];
                //     feqf[fineindex(nx - 1, j, k, a)] = w[a] * rhof[fineisnindex1(nx - 1, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // }

                // //*****Find Non-equilibrium distribution at i=nx-1 using 4 point langrangian interpolation*******
                // for (a = 0; a < 19; a++)
                // {
                //     fneqf[a] = (-dist_neq_l[0][a] + 9.0 * dist_neq_l[1][a] + 9.0 * dist_neq_l[2][a] - dist_neq_l[3][a]) / 16.0;
                // }

                //*****Find f = feq+fneq******
                for (i = nx - 1, a = 0; a < 9; a++)
                {
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf[a];

                    ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(nx, j, k, a)] + 9.0 * ff[fineindex(nx - 2, j, k, a)] - ff[fineindex(nx - 4, j, k, a)]) / 16.0;
                }
            }
        }
    }
    if (side == 3)
    {
        for (J = CY_F + 2, I = Start1, i = start; I <= End1; I++, i += alpha)
        {
            for (K = Start2, k = start; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    dist_l[a] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                    // dist_neq_l[0][a] = (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }

                // l = 0.0;
                // m = 0.0;
                // n = 0.0;
                // p = 0.0;

                // for (a = 0; a < 19; a++)
                // {
                //     l += dist_l[a];
                //     m += ex[a] * dist_l[a];
                //     n += ey[a] * dist_l[a];
                //     p += ez[a] * dist_l[a];
                // }
                // rho_l = l;
                // jx_l = m / rho_l;
                // jy_l = n / rho_l;
                // jz_l = p / rho_l;

                // //******Find the non equilibrium points for j=ny, j=ny-2, j=ny-4
                // for (q = 1, j = ny; j >= ny - 4; j -= 2, q++)
                // {
                //     l = 0.0;
                //     m = 0.0;
                //     n = 0.0;
                //     p = 0.0;

                //     for (a = 0; a < 19; a++)
                //     {
                //         l += ff[fineindex(i, j, k, a)];
                //         m += ex[a] * ff[fineindex(i, j, k, a)];
                //         n += ey[a] * ff[fineindex(i, j, k, a)];
                //         p += ez[a] * ff[fineindex(i, j, k, a)];
                //     }
                //     rhof[fineisnindex1(i, j, k)] = l;
                //     jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                //     jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                //     jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                //     u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                //     for (a = 0; a < 19; a++)
                //     {
                //         t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                //         feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                //         dist_neq_l[q][a] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]; // Non equilibrium distribution for j=-2 point
                //     }
                // }

                // //*****4 point langrangian interpolation of the conserved quantities****

                // jxf[fineisnindex1(i, ny - 1, k)] = (-jx_l + 9.0 * jxf[fineisnindex1(i, ny, k)] + 9.0 * jxf[fineisnindex1(i, ny - 2, k)] - jxf[fineisnindex1(i, ny - 4, k)]) / 16.0;
                // jyf[fineisnindex1(i, ny - 1, k)] = (-jy_l + 9.0 * jyf[fineisnindex1(i, ny, k)] + 9.0 * jyf[fineisnindex1(i, ny - 2, k)] - jyf[fineisnindex1(i, ny - 4, k)]) / 16.0;
                // jzf[fineisnindex1(i, ny - 1, k)] = (-jz_l + 9.0 * jzf[fineisnindex1(i, ny, k)] + 9.0 * jzf[fineisnindex1(i, ny - 2, k)] - jzf[fineisnindex1(i, ny - 4, k)]) / 16.0;
                // rhof[fineisnindex1(i, ny - 1, k)] = (-rho_l + 9.0 * rhof[fineisnindex1(i, ny, k)] + 9.0 * rhof[fineisnindex1(i, ny - 2, k)] - rhof[fineisnindex1(i, ny - 4, k)]) / 16.0;

                // u2 = jxf[fineisnindex1(i, ny - 1, k)] * jxf[fineisnindex1(i, ny - 1, k)] + jyf[fineisnindex1(i, ny - 1, k)] * jyf[fineisnindex1(i, ny - 1, k)] + jzf[fineisnindex1(i, ny - 1, k)] * jzf[fineisnindex1(i, ny - 1, k)];

                // //*****Equilibrium distribution at j=ny-1 ********
                // for (a = 0; a < 19; a++)
                // {
                //     t = ex[a] * jxf[fineisnindex1(i, ny - 1, k)] + ey[a] * jyf[fineisnindex1(i, ny - 1, k)] + ez[a] * jzf[fineisnindex1(i, ny - 1, k)];
                //     feqf[fineindex(i, ny - 1, k, a)] = w[a] * rhof[fineisnindex1(i, ny - 1, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // }

                // //*****Find Non-equilibrium distribution at j=ny-1 using 4 point langrangian interpolation*******
                // for (a = 0; a < 19; a++)
                // {
                //     fneqf[a] = (-dist_neq_l[0][a] + 9.0 * dist_neq_l[1][a] + 9.0 * dist_neq_l[2][a] - dist_neq_l[3][a]) / 16.0;
                // }

                //*****Find f = feq+fneq******
                for (j = ny - 1, a = 0; a < 19; a++)
                {
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf[a];

                    ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, ny, k, a)] + 9.0 * ff[fineindex(i, ny - 2, k, a)] - ff[fineindex(i, ny - 4, k, a)]) / 16.0;
                }
            }
        }
    }
    if (side == 4)
    {
        for (J = Start1, j = start, I = CX_I - 2; J <= End1; J++, j += alpha)
        {
            for (K = Start2, k = start; K <= End2; K++, k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    dist_l[a] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                    // dist_neq_l[0][a] = (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }

                // l = 0.0;
                // m = 0.0;
                // n = 0.0;
                // p = 0.0;

                // for (a = 0; a < 19; a++)
                // {
                //     l += dist_l[a];
                //     m += ex[a] * dist_l[a];
                //     n += ey[a] * dist_l[a];
                //     p += ez[a] * dist_l[a];
                // }
                // rho_l = l;
                // jx_l = m / rho_l;
                // jy_l = n / rho_l;
                // jz_l = p / rho_l;

                // //******Find the non equilibrium points for i=0,i=2,i=4
                // for (q = 1, i = 0; i <= 4; i += 2, q++)
                // {
                //     l = 0.0;
                //     m = 0.0;
                //     n = 0.0;
                //     p = 0.0;

                //     for (a = 0; a < 19; a++)
                //     {
                //         l += ff[fineindex(i, j, k, a)];
                //         m += ex[a] * ff[fineindex(i, j, k, a)];
                //         n += ey[a] * ff[fineindex(i, j, k, a)];
                //         p += ez[a] * ff[fineindex(i, j, k, a)];
                //     }
                //     rhof[fineisnindex1(i, j, k)] = l;
                //     jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                //     jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                //     jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                //     u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                //     for (a = 0; a < 19; a++)
                //     {
                //         t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                //         feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                //         dist_neq_l[q][a] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]; // Non equilibrium distribution for j=-2 point
                //     }
                // }

                // //*****4 point langrangian interpolation of the conserved quantities****

                // jxf[fineisnindex1(1, j, k)] = (-jx_l + 9.0 * jxf[fineisnindex1(0, j, k)] + 9.0 * jxf[fineisnindex1(2, j, k)] - jxf[fineisnindex1(4, j, k)]) / 16.0;
                // jyf[fineisnindex1(1, j, k)] = (-jy_l + 9.0 * jyf[fineisnindex1(0, j, k)] + 9.0 * jyf[fineisnindex1(2, j, k)] - jyf[fineisnindex1(4, j, k)]) / 16.0;
                // jzf[fineisnindex1(1, j, k)] = (-jz_l + 9.0 * jzf[fineisnindex1(0, j, k)] + 9.0 * jzf[fineisnindex1(2, j, k)] - jzf[fineisnindex1(4, j, k)]) / 16.0;
                // rhof[fineisnindex1(1, j, k)] = (-rho_l + 9.0 * rhof[fineisnindex1(0, j, k)] + 9.0 * rhof[fineisnindex1(2, j, k)] - rhof[fineisnindex1(4, j, k)]) / 16.0;

                // u2 = jxf[fineisnindex1(1, j, k)] * jxf[fineisnindex1(1, j, k)] + jyf[fineisnindex1(1, j, k)] * jyf[fineisnindex1(1, j, k)] + jzf[fineisnindex1(1, j, k)] * jzf[fineisnindex1(1, j, k)];

                // //****** Equilibrium distribution at i=1*******
                // for (a = 0; a < 19; a++)
                // {
                //     t = ex[a] * jxf[fineisnindex1(1, j, k)] + ey[a] * jyf[fineisnindex1(1, j, k)] + ez[a] * jzf[fineisnindex1(1, j, k)];
                //     feqf[fineindex(1, j, k, a)] = w[a] * rhof[fineisnindex1(1, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // }

                // //*****Find Non-equilibrium distribution at i=1 using 4 point langrangian interpolation*******
                // for (a = 0; a < 19; a++)
                // {
                //     fneqf[a] = (-dist_neq_l[0][a] + 9.0 * dist_neq_l[1][a] + 9.0 * dist_neq_l[2][a] - dist_neq_l[3][a]) / 16.0;
                // }

                //*****Find f = feq+fneq******
                for (i = 1, a = 0; a < 19; a++)
                {
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf[a];

                    ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(0, j, k, a)] + 9.0 * ff[fineindex(2, j, k, a)] - ff[fineindex(4, j, k, a)]) / 16.0;
                }
            }
        }
    }
    if (side == 5)
    {
        for (I = Start1, i = start, K = CZ_I - 2; I <= End1; I++, i += alpha)
        {
            for (J = Start2, j = start; J <= End2; J++, j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    dist_l[a] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                    // dist_neq_l[0][a] = (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }

                // l = 0.0;
                // m = 0.0;
                // n = 0.0;
                // p = 0.0;

                // for (a = 0; a < 19; a++)
                // {
                //     l += dist_l[a];
                //     m += ex[a] * dist_l[a];
                //     n += ey[a] * dist_l[a];
                //     p += ez[a] * dist_l[a];
                // }
                // rho_l = l;
                // jx_l = m / rho_l;
                // jy_l = n / rho_l;
                // jz_l = p / rho_l;

                // //******Find the non equilibrium points for k=0,k=2,k=4
                // for (q = 1, k = 0; k <= 4; k += 2, q++)
                // {
                //     l = 0.0;
                //     m = 0.0;
                //     n = 0.0;
                //     p = 0.0;

                //     for (a = 0; a < 19; a++)
                //     {
                //         l += ff[fineindex(i, j, k, a)];
                //         m += ex[a] * ff[fineindex(i, j, k, a)];
                //         n += ey[a] * ff[fineindex(i, j, k, a)];
                //         p += ez[a] * ff[fineindex(i, j, k, a)];
                //     }
                //     rhof[fineisnindex1(i, j, k)] = l;
                //     jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                //     jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                //     jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                //     u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                //     for (a = 0; a < 19; a++)
                //     {
                //         t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                //         feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                //         dist_neq_l[q][a] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]; // Non equilibrium distribution for j=-2 point
                //     }
                // }

                // //*****4 point langrangian interpolation of the conserved quantities****

                // jxf[fineisnindex1(i, j, 1)] = (-jx_l + 9.0 * jxf[fineisnindex1(i, j, 0)] + 9.0 * jxf[fineisnindex1(i, j, 2)] - jxf[fineisnindex1(i, j, 4)]) / 16.0;
                // jyf[fineisnindex1(i, j, 1)] = (-jy_l + 9.0 * jyf[fineisnindex1(i, j, 0)] + 9.0 * jyf[fineisnindex1(i, j, 2)] - jyf[fineisnindex1(i, j, 4)]) / 16.0;
                // jzf[fineisnindex1(i, j, 1)] = (-jz_l + 9.0 * jzf[fineisnindex1(i, j, 0)] + 9.0 * jzf[fineisnindex1(i, j, 2)] - jzf[fineisnindex1(i, j, 4)]) / 16.0;
                // rhof[fineisnindex1(i, j, 1)] = (-rho_l + 9.0 * rhof[fineisnindex1(i, j, 0)] + 9.0 * rhof[fineisnindex1(i, j, 2)] - rhof[fineisnindex1(i, j, 4)]) / 16.0;

                // u2 = jxf[fineisnindex1(i, j, 1)] * jxf[fineisnindex1(i, j, 1)] + jyf[fineisnindex1(i, j, 1)] * jyf[fineisnindex1(i, j, 1)] + jzf[fineisnindex1(i, j, 1)] * jzf[fineisnindex1(i, j, 1)];

                // //****** Equilibrium distribution at i=1*******
                // for (a = 0; a < 19; a++)
                // {
                //     t = ex[a] * jxf[fineisnindex1(i, j, 1)] + ey[a] * jyf[fineisnindex1(i, j, 1)] + ez[a] * jzf[fineisnindex1(i, j, 1)];
                //     feqf[fineindex(i, j, 1, a)] = w[a] * rhof[fineisnindex1(i, j, 1)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // }

                // //*****Find Non-equilibrium distribution at i=1 using 4 point langrangian interpolation*******
                // for (a = 0; a < 19; a++)
                // {
                //     fneqf[a] = (-dist_neq_l[0][a] + 9.0 * dist_neq_l[1][a] + 9.0 * dist_neq_l[2][a] - dist_neq_l[3][a]) / 16.0;
                // }

                //*****Find f = feq+fneq******
                for (k = 1, a = 0; a < 19; a++)
                {
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf[a];

                    ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, j, 0, a)] + 9.0 * ff[fineindex(i, j, 2, a)] - ff[fineindex(i, j, 4, a)]) / 16.0;
                }
            }
        }
    }
    if (side == 6)
    {
        for (I = Start1, i = start, K = CZ_F + 2; I <= End1; I++, i += alpha)
        {
            for (J = Start2, j = start; J <= End2; J++, j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += fc[coarseindex(I, J, K, a)];
                    m += ex[a] * fc[coarseindex(I, J, K, a)];
                    n += ey[a] * fc[coarseindex(I, J, K, a)];
                    p += ez[a] * fc[coarseindex(I, J, K, a)];
                }
                rhoc[coarseisnindex1(I, J, K)] = l;
                jxc[coarseisnindex1(I, J, K)] = m / rhoc[coarseisnindex1(I, J, K)];
                jyc[coarseisnindex1(I, J, K)] = n / rhoc[coarseisnindex1(I, J, K)];
                jzc[coarseisnindex1(I, J, K)] = p / rhoc[coarseisnindex1(I, J, K)];

                u2 = jxc[coarseisnindex1(I, J, K)] * jxc[coarseisnindex1(I, J, K)] + jyc[coarseisnindex1(I, J, K)] * jyc[coarseisnindex1(I, J, K)] + jzc[coarseisnindex1(I, J, K)] * jzc[coarseisnindex1(I, J, K)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxc[coarseisnindex1(I, J, K)] + ey[a] * jyc[coarseisnindex1(I, J, K)] + ez[a] * jzc[coarseisnindex1(I, J, K)];
                    feqc[coarseindex(I, J, K, a)] = w[a] * rhoc[coarseisnindex1(I, J, K)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    dist_l[a] = feqc[coarseindex(I, J, K, a)] + (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                    // dist_neq_l[0][a] = (tauf / tauc) * (fc[coarseindex(I, J, K, a)] - feqc[coarseindex(I, J, K, a)]) / alpha;
                }

                // l = 0.0;
                // m = 0.0;
                // n = 0.0;
                // p = 0.0;

                // for (a = 0; a < 19; a++)
                // {
                //     l += dist_l[a];
                //     m += ex[a] * dist_l[a];
                //     n += ey[a] * dist_l[a];
                //     p += ez[a] * dist_l[a];
                // }
                // rho_l = l;
                // jx_l = m / rho_l;
                // jy_l = n / rho_l;
                // jz_l = p / rho_l;

                // //******Find the non equilibrium points for k=0,k=2,k=4
                // for (q = 1, k = nz; k >= nz - 4; k -= 2, q++)
                // {
                //     l = 0.0;
                //     m = 0.0;
                //     n = 0.0;
                //     p = 0.0;

                //     for (a = 0; a < 19; a++)
                //     {
                //         l += ff[fineindex(i, j, k, a)];
                //         m += ex[a] * ff[fineindex(i, j, k, a)];
                //         n += ey[a] * ff[fineindex(i, j, k, a)];
                //         p += ez[a] * ff[fineindex(i, j, k, a)];
                //     }
                //     rhof[fineisnindex1(i, j, k)] = l;
                //     jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                //     jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                //     jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                //     u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                //     for (a = 0; a < 19; a++)
                //     {
                //         t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                //         feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                //         dist_neq_l[q][a] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)]; // Non equilibrium distribution for j=-2 point
                //     }
                // }

                // //*****4 point langrangian interpolation of the conserved quantities****

                // jxf[fineisnindex1(i, j, nz - 1)] = (-jx_l + 9.0 * jxf[fineisnindex1(i, j, nz)] + 9.0 * jxf[fineisnindex1(i, j, nz - 2)] - jxf[fineisnindex1(i, j, nz - 4)]) / 16.0;
                // jyf[fineisnindex1(i, j, nz - 1)] = (-jy_l + 9.0 * jyf[fineisnindex1(i, j, nz)] + 9.0 * jyf[fineisnindex1(i, j, nz - 2)] - jyf[fineisnindex1(i, j, nz - 4)]) / 16.0;
                // jzf[fineisnindex1(i, j, nz - 1)] = (-jz_l + 9.0 * jzf[fineisnindex1(i, j, nz)] + 9.0 * jzf[fineisnindex1(i, j, nz - 2)] - jzf[fineisnindex1(i, j, nz - 4)]) / 16.0;
                // rhof[fineisnindex1(i, j, nz - 1)] = (-rho_l + 9.0 * rhof[fineisnindex1(i, j, nz)] + 9.0 * rhof[fineisnindex1(i, j, nz - 2)] - rhof[fineisnindex1(i, j, nz - 4)]) / 16.0;

                // u2 = jxf[fineisnindex1(i, j, nz - 1)] * jxf[fineisnindex1(i, j, nz - 1)] + jxf[fineisnindex1(i, j, nz - 1)] * jxf[fineisnindex1(i, j, nz - 1)] + jzf[fineisnindex1(i, j, nz - 1)] * jzf[fineisnindex1(i, j, nz - 1)];

                // //****** Equilibrium distribution at i=1*******
                // for (a = 0; a < 19; a++)
                // {
                //     t = ex[a] * jxf[fineisnindex1(i, j, nz - 1)] + ey[a] * jyf[fineisnindex1(i, j, nz - 1)] + ez[a] * jzf[fineisnindex1(i, j, nz - 1)];
                //     feqf[fineindex(i, j, nz - 1, a)] = w[a] * rhof[fineisnindex1(i, j, nz - 1)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // }

                // //*****Find Non-equilibrium distribution at i=1 using 4 point langrangian interpolation*******
                // for (a = 0; a < 19; a++)
                // {
                //     fneqf[a] = (-dist_neq_l[0][a] + 9.0 * dist_neq_l[1][a] + 9.0 * dist_neq_l[2][a] - dist_neq_l[3][a]) / 16.0;
                // }

                //*****Find f = feq+fneq******
                for (k = nz - 1, a = 0; a < 19; a++)
                {
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf[a];

                    ff[fineindex(i, j, k, a)] = (-dist_l[a] + 9.0 * ff[fineindex(i, j, nz, a)] + 9.0 * ff[fineindex(i, j, nz - 2, a)] - ff[fineindex(i, j, nz - 4, a)]) / 16.0;
                }
            }
        }
    }
}

void interpolation_mid_boundary(double *ff, double *feqf, double *fneqf_b, double *fneqf_r, double *fneqf_t, double *fneqf_l, double *fneqf_u, double *fneqf_d, int side, int face, int ex[19], int ey[19], int ez[19], double w[19], double *rhof, double *jxf, double *jyf, double *jzf, int alpha, int CX_I, int CX_F, int CY_I, int CY_F, int CZ_I, int CZ_F)
{
    int I, J, K, i, j, k, a;
    double l, m, n, p, u2, t;

    if (side == 1)
    {
        // for (i = 0, j = face; i <= nx; i += alpha)
        // {
        //     for (k = 0; k <= nz; k += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_b[topbottom(i, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (j = face, k = 2; k <= nz - 2; k += alpha)
        {
            for (i = 3; i <= nx - 3; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_b[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_b[topbottom(i + 1, k, a)] + fneqf_b[topbottom(i - 1, k, a)]) - (1.0 / 16.0) * (fneqf_b[topbottom(i + 3, k, a)] + fneqf_b[topbottom(i - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (j = face, i = 2; i <= nx - 2; i += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_b[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_b[topbottom(i, k + 1, a)] + fneqf_b[topbottom(i, k - 1, a)]) - (1.0 / 16.0) * (fneqf_b[topbottom(i, k + 3, a)] + fneqf_b[topbottom(i, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }

        //****Corner points topbottom side
        for (k = 2; k <= nz - 2; k += alpha)
        {
            i = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i - 1, k, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i + 1, k, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            i = nx - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i + 1, k, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i - 1, k, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }

        for (i = 2; i <= nx - 2; i += alpha)
        {
            k = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i, k - 1, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i, k + 1, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            k = nz - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i, k + 1, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i, k - 1, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 2)
    {
        // for (k = 0, i = face; k <= nz; k += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_r[leftright(j, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (i = face, k = 2; k <= nz - 2; k += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_r[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_r[leftright(j + 1, k, a)] + fneqf_r[leftright(j - 1, k, a)]) - (1.0 / 16.0) * (fneqf_r[leftright(j + 3, k, a)] + fneqf_r[leftright(j - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }

        for (i = face, j = 2; j <= ny - 2; j += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_r[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_r[leftright(j, k + 1, a)] + fneqf_r[leftright(j, k - 1, a)]) - (1.0 / 16.0) * (fneqf_r[leftright(j, k + 3, a)] + fneqf_r[leftright(j, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }

        //******Corner points leftright Side
        for (k = 2; k <= nz - 2; k += alpha)
        {
            i = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j - 1, k, a)] + (3.0 / 4.0) * fneqf_r[leftright(j + 1, k, a)] - (1.0 / 8.0) * fneqf_r[leftright(j + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            i = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j + 1, k, a)] + (3.0 / 4.0) * fneqf_r[leftright(j - 1, k, a)] - (1.0 / 8.0) * fneqf_r[leftright(j - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }
        for (j = 2; j <= ny - 2; j += alpha)
        {
            i = face;
            k = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j, k - 1, a)] + (3.0 / 4.0) * fneqf_r[leftright(j, k + 1, a)] - (1.0 / 8.0) * fneqf_r[leftright(j, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            i = face;
            k = nz - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j, k + 1, a)] + (3.0 / 4.0) * fneqf_r[leftright(j, k - 1, a)] - (1.0 / 8.0) * fneqf_r[leftright(j, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 3)
    {
        // for (i = 0, j = face; i <= nx; i += alpha)
        // {
        //     for (k = 0; k <= nz; k += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_t[topbottom(i, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (j = face, k = 2; k <= nz - 2; k += alpha)
        {
            for (i = 3; i <= nx - 3; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_t[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_t[topbottom(i + 1, k, a)] + fneqf_t[topbottom(i - 1, k, a)]) - (1.0 / 16.0) * (fneqf_t[topbottom(i + 3, k, a)] + fneqf_t[topbottom(i - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (j = face, i = 2; i <= nx - 2; i += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_t[topbottom(i, k, a)] = (9.0 / 16.0) * (fneqf_t[topbottom(i, k + 1, a)] + fneqf_t[topbottom(i, k - 1, a)]) - (1.0 / 16.0) * (fneqf_t[topbottom(i, k + 3, a)] + fneqf_t[topbottom(i, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }

        //********Corner points Top side
        for (k = 2; k <= nz - 2; k += alpha)
        {
            i = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i - 1, k, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i + 1, k, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            i = nx - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i + 1, k, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i - 1, k, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }

        for (i = 2; i <= nx - 2; i += alpha)
        {
            k = 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i, k - 1, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i, k + 1, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            k = nz - 1;
            j = face;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i, k + 1, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i, k - 1, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 4)
    {
        // for (k = 0, i = face; k <= nz; k += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_l[leftright(j, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (i = face, k = 2; k <= nz - 2; k += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_l[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_l[leftright(j + 1, k, a)] + fneqf_l[leftright(j - 1, k, a)]) - (1.0 / 16.0) * (fneqf_l[leftright(j + 3, k, a)] + fneqf_l[leftright(j - 3, k, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }

        for (i = face, j = 2; j <= ny - 2; j += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 1)] + jxf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j, k + 3)] + jxf[fineisnindex1(i, j, k - 3)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 1)] + jyf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j, k + 3)] + jyf[fineisnindex1(i, j, k - 3)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 1)] + jzf[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j, k + 3)] + jzf[fineisnindex1(i, j, k - 3)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 1)] + rhof[fineisnindex1(i, j, k - 1)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j, k + 3)] + rhof[fineisnindex1(i, j, k - 3)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_l[leftright(j, k, a)] = (9.0 / 16.0) * (fneqf_l[leftright(j, k + 1, a)] + fneqf_l[leftright(j, k - 1, a)]) - (1.0 / 16.0) * (fneqf_l[leftright(j, k + 3, a)] + fneqf_l[leftright(j, k - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j, k + 1, a)] + ff[fineindex(i, j, k - 1, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j, k + 3, a)] + ff[fineindex(i, j, k - 3, a)]);
                }
            }
        }
        //********Corner points Left side
        for (k = 2; k <= nz - 2; k += alpha)
        {
            i = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j - 1, k, a)] + (3.0 / 4.0) * fneqf_l[leftright(j + 1, k, a)] - (1.0 / 8.0) * fneqf_l[leftright(j + 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            i = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j + 1, k, a)] + (3.0 / 4.0) * fneqf_l[leftright(j - 1, k, a)] - (1.0 / 8.0) * fneqf_l[leftright(j - 3, k, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }

        for (j = 2; j <= ny - 2; j += alpha)
        {
            i = face;
            k = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j, k - 1, a)] + (3.0 / 4.0) * fneqf_l[leftright(j, k + 1, a)] - (1.0 / 8.0) * fneqf_l[leftright(j, k + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k - 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k + 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k + 3, a)];
            }

            i = face;
            k = nz - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j, k + 1, a)] + (3.0 / 4.0) * fneqf_l[leftright(j, k - 1, a)] - (1.0 / 8.0) * fneqf_l[leftright(j, k - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j, k + 1, a)] + (3.0 / 4.0) * ff[fineindex(i, j, k - 1, a)] - (1.0 / 8.0) * ff[fineindex(i, j, k - 3, a)];
            }
        }
    }
    if (side == 5)
    {
        // for (i = 0, k = face; i <= nx; i += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_d[updown(i, j, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (k = face, i = 3; i <= nx - 3; i += alpha)
        {
            for (j = 2; j <= ny - 2; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_d[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_d[updown(i + 1, j, a)] + fneqf_d[updown(i - 1, j, a)]) - (1.0 / 16.0) * (fneqf_d[updown(i + 3, j, a)] + fneqf_d[updown(i - 3, j, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (k = face, j = 3; j <= ny - 3; j += alpha)
        {
            for (i = 2; i <= nx - 2; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_d[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_d[updown(i, j + 1, a)] + fneqf_d[updown(i, j - 1, a)]) - (1.0 / 16.0) * (fneqf_d[updown(i, j + 3, a)] + fneqf_d[updown(i, j - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }
        //********Corner points down side
        for (i = 2; i <= nx - 2; i += alpha)
        {
            k = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i, j - 1, a)] + (3.0 / 4.0) * fneqf_d[updown(i, j + 1, a)] - (1.0 / 8.0) * fneqf_d[updown(i, j + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            k = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i, j + 1, a)] + (3.0 / 4.0) * fneqf_d[updown(i, j - 1, a)] - (1.0 / 8.0) * fneqf_d[updown(i, j - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }

        for (j = 2; j <= ny - 2; j += alpha)
        {
            k = face;
            i = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i - 1, j, a)] + (3.0 / 4.0) * fneqf_d[updown(i + 1, j, a)] - (1.0 / 8.0) * fneqf_d[updown(i + 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            k = face;
            i = nx - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i + 1, j, a)] + (3.0 / 4.0) * fneqf_d[updown(i - 1, j, a)] - (1.0 / 8.0) * fneqf_d[updown(i - 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }
    }

    if (side == 6)
    {
        // for (i = 0, k = face; i <= nx; i += alpha)
        // {
        //     for (j = 0; j <= ny; j += alpha)
        //     {
        //         l = 0.0;
        //         m = 0.0;
        //         n = 0.0;
        //         p = 0.0;

        //         for (a = 0; a < 19; a++)
        //         {
        //             l += ff[fineindex(i, j, k, a)];
        //             m += ex[a] * ff[fineindex(i, j, k, a)];
        //             n += ey[a] * ff[fineindex(i, j, k, a)];
        //             p += ez[a] * ff[fineindex(i, j, k, a)];
        //         }
        //         rhof[fineisnindex1(i, j, k)] = l;
        //         jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
        //         jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
        //         jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

        //         u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

        //         for (a = 0; a < 19; a++)
        //         {
        //             t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
        //             feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
        //             fneqf_u[updown(i, j, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
        //         }
        //     }
        // }

        for (k = face, i = 3; i <= nx - 3; i += alpha)
        {
            for (j = 2; j <= ny - 2; j += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i + 1, j, k)] + jxf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i + 3, j, k)] + jxf[fineisnindex1(i - 3, j, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i + 1, j, k)] + jyf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i + 3, j, k)] + jyf[fineisnindex1(i - 3, j, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i + 1, j, k)] + jzf[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i + 3, j, k)] + jzf[fineisnindex1(i - 3, j, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i + 1, j, k)] + rhof[fineisnindex1(i - 1, j, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i + 3, j, k)] + rhof[fineisnindex1(i - 3, j, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_u[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_u[updown(i + 1, j, a)] + fneqf_u[updown(i - 1, j, a)]) - (1.0 / 16.0) * (fneqf_u[updown(i + 3, j, a)] + fneqf_u[updown(i - 3, j, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i + 1, j, k, a)] + ff[fineindex(i - 1, j, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i + 3, j, k, a)] + ff[fineindex(i - 3, j, k, a)]);
                }
            }
        }

        for (k = face, j = 3; j <= ny - 3; j += alpha)
        {
            for (i = 2; i <= nx - 2; i += alpha)
            {
                // jxf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jxf[fineisnindex1(i, j + 1, k)] + jxf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jxf[fineisnindex1(i, j + 3, k)] + jxf[fineisnindex1(i, j - 3, k)]);
                // jyf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jyf[fineisnindex1(i, j + 1, k)] + jyf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jyf[fineisnindex1(i, j + 3, k)] + jyf[fineisnindex1(i, j - 3, k)]);
                // jzf[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (jzf[fineisnindex1(i, j + 1, k)] + jzf[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (jzf[fineisnindex1(i, j + 3, k)] + jzf[fineisnindex1(i, j - 3, k)]);
                // rhof[fineisnindex1(i, j, k)] = (9.0 / 16.0) * (rhof[fineisnindex1(i, j + 1, k)] + rhof[fineisnindex1(i, j - 1, k)]) - (1.0 / 16.0) * (rhof[fineisnindex1(i, j + 3, k)] + rhof[fineisnindex1(i, j - 3, k)]);

                // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    // fneqf_u[updown(i, j, a)] = (9.0 / 16.0) * (fneqf_u[updown(i, j + 1, a)] + fneqf_u[updown(i, j - 1, a)]) - (1.0 / 16.0) * (fneqf_u[updown(i, j + 3, a)] + fneqf_u[updown(i, j - 3, a)]);

                    // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                    ff[fineindex(i, j, k, a)] = (9.0 / 16.0) * (ff[fineindex(i, j + 1, k, a)] + ff[fineindex(i, j - 1, k, a)]) - (1.0 / 16.0) * (ff[fineindex(i, j + 3, k, a)] + ff[fineindex(i, j - 3, k, a)]);
                }
            }
        }
        //********Corner points up side
        for (i = 2; i <= nx - 2; i += alpha)
        {
            k = face;
            j = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i, j - 1, a)] + (3.0 / 4.0) * fneqf_u[updown(i, j + 1, a)] - (1.0 / 8.0) * fneqf_u[updown(i, j + 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j - 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j + 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j + 3, k, a)];
            }

            k = face;
            j = ny - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i, j + 1, a)] + (3.0 / 4.0) * fneqf_u[updown(i, j - 1, a)] - (1.0 / 8.0) * fneqf_u[updown(i, j - 3, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i, j + 1, k, a)] + (3.0 / 4.0) * ff[fineindex(i, j - 1, k, a)] - (1.0 / 8.0) * ff[fineindex(i, j - 3, k, a)];
            }
        }

        for (j = 2; j <= ny - 2; j += alpha)
        {
            k = face;
            i = 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i - 1, j, a)] + (3.0 / 4.0) * fneqf_u[updown(i + 1, j, a)] - (1.0 / 8.0) * fneqf_u[updown(i + 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i - 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i + 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i + 3, j, k, a)];
            }

            k = face;
            i = nx - 1;

            // jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            // jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            // jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            // rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            // u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                // fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i + 1, j, a)] + (3.0 / 4.0) * fneqf_u[updown(i - 1, j, a)] - (1.0 / 8.0) * fneqf_u[updown(i - 3, j, a)];

                // t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                // feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                // ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];

                ff[fineindex(i, j, k, a)] = (3.0 / 8.0) * ff[fineindex(i + 1, j, k, a)] + (3.0 / 4.0) * ff[fineindex(i - 1, j, k, a)] - (1.0 / 8.0) * ff[fineindex(i - 3, j, k, a)];
            }
        }
    }
}

void bicubic_interpolation(double *ff, double *feqf, double *fneqf_b, double *fneqf_r, double *fneqf_t, double *fneqf_l, double *fneqf_u, double *fneqf_d, int side, int face, double *jxf, double *jyf, double *jzf, double *rhof, int ex[19], int ey[19], int ez[19], double w[19], int alpha)
{
    int i, j, k, a;
    double l, m, n, p, u2, t;
    if (side == 1)
    {
        for (i = 0, j = face; i <= nx; i += alpha)
        {
            for (k = 0; k <= nz; k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    fneqf_b[topbottom(i, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
                }
            }
        }

        for (j = face, i = 3; i <= nx - 3; i += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                jxf[fineisnindex1(i, j, k)] = (jxf[fineisnindex1(i - 3, j, k + 3)] - 9.0 * jxf[fineisnindex1(i - 3, j, k + 1)] - 9.0 * jxf[fineisnindex1(i - 3, j, k - 1)] + jxf[fineisnindex1(i - 3, j, k - 3)] - 9.0 * jxf[fineisnindex1(i - 1, j, k + 3)] + 81.0 * jxf[fineisnindex1(i - 1, j, k + 1)] + 81.0 * jxf[fineisnindex1(i - 1, j, k - 1)] - 9.0 * jxf[fineisnindex1(i - 1, j, k - 3)] - 9.0 * jxf[fineisnindex1(i + 1, j, k + 3)] + 81.0 * jxf[fineisnindex1(i + 1, j, k + 1)] + 81.0 * jxf[fineisnindex1(i + 1, j, k - 1)] - 9.0 * jxf[fineisnindex1(i + 1, j, k - 3)] + jxf[fineisnindex1(i + 3, j, k + 3)] - 9.0 * jxf[fineisnindex1(i + 3, j, k + 1)] - 9.0 * jxf[fineisnindex1(i + 3, j, k - 1)] + jxf[fineisnindex1(i + 3, j, k - 3)]) / 256.0;
                jyf[fineisnindex1(i, j, k)] = (jyf[fineisnindex1(i - 3, j, k + 3)] - 9.0 * jyf[fineisnindex1(i - 3, j, k + 1)] - 9.0 * jyf[fineisnindex1(i - 3, j, k - 1)] + jyf[fineisnindex1(i - 3, j, k - 3)] - 9.0 * jyf[fineisnindex1(i - 1, j, k + 3)] + 81.0 * jyf[fineisnindex1(i - 1, j, k + 1)] + 81.0 * jyf[fineisnindex1(i - 1, j, k - 1)] - 9.0 * jyf[fineisnindex1(i - 1, j, k - 3)] - 9.0 * jyf[fineisnindex1(i + 1, j, k + 3)] + 81.0 * jyf[fineisnindex1(i + 1, j, k + 1)] + 81.0 * jyf[fineisnindex1(i + 1, j, k - 1)] - 9.0 * jyf[fineisnindex1(i + 1, j, k - 3)] + jyf[fineisnindex1(i + 3, j, k + 3)] - 9.0 * jyf[fineisnindex1(i + 3, j, k + 1)] - 9.0 * jyf[fineisnindex1(i + 3, j, k - 1)] + jyf[fineisnindex1(i + 3, j, k - 3)]) / 256.0;
                jzf[fineisnindex1(i, j, k)] = (jzf[fineisnindex1(i - 3, j, k + 3)] - 9.0 * jzf[fineisnindex1(i - 3, j, k + 1)] - 9.0 * jzf[fineisnindex1(i - 3, j, k - 1)] + jzf[fineisnindex1(i - 3, j, k - 3)] - 9.0 * jzf[fineisnindex1(i - 1, j, k + 3)] + 81.0 * jzf[fineisnindex1(i - 1, j, k + 1)] + 81.0 * jzf[fineisnindex1(i - 1, j, k - 1)] - 9.0 * jzf[fineisnindex1(i - 1, j, k - 3)] - 9.0 * jzf[fineisnindex1(i + 1, j, k + 3)] + 81.0 * jzf[fineisnindex1(i + 1, j, k + 1)] + 81.0 * jzf[fineisnindex1(i + 1, j, k - 1)] - 9.0 * jzf[fineisnindex1(i + 1, j, k - 3)] + jzf[fineisnindex1(i + 3, j, k + 3)] - 9.0 * jzf[fineisnindex1(i + 3, j, k + 1)] - 9.0 * jzf[fineisnindex1(i + 3, j, k - 1)] + jzf[fineisnindex1(i + 3, j, k - 3)]) / 256.0;
                rhof[fineisnindex1(i, j, k)] = (rhof[fineisnindex1(i - 3, j, k + 3)] - 9.0 * rhof[fineisnindex1(i - 3, j, k + 1)] - 9.0 * rhof[fineisnindex1(i - 3, j, k - 1)] + rhof[fineisnindex1(i - 3, j, k - 3)] - 9.0 * rhof[fineisnindex1(i - 1, j, k + 3)] + 81.0 * rhof[fineisnindex1(i - 1, j, k + 1)] + 81.0 * rhof[fineisnindex1(i - 1, j, k - 1)] - 9.0 * rhof[fineisnindex1(i - 1, j, k - 3)] - 9.0 * rhof[fineisnindex1(i + 1, j, k + 3)] + 81.0 * rhof[fineisnindex1(i + 1, j, k + 1)] + 81.0 * rhof[fineisnindex1(i + 1, j, k - 1)] - 9.0 * rhof[fineisnindex1(i + 1, j, k - 3)] + rhof[fineisnindex1(i + 3, j, k + 3)] - 9.0 * rhof[fineisnindex1(i + 3, j, k + 1)] - 9.0 * rhof[fineisnindex1(i + 3, j, k - 1)] + rhof[fineisnindex1(i + 3, j, k - 3)]) / 256.0;

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    fneqf_b[topbottom(i, k, a)] = (fneqf_b[topbottom(i - 3, k + 3, a)] - 9.0 * fneqf_b[topbottom(i - 3, k + 1, a)] - 9.0 * fneqf_b[topbottom(i - 3, k - 1, a)] + fneqf_b[topbottom(i - 3, k - 3, a)] - 9.0 * fneqf_b[topbottom(i - 1, k + 3, a)] + 81.0 * fneqf_b[topbottom(i - 1, k + 1, a)] + 81.0 * fneqf_b[topbottom(i - 1, k - 1, a)] - 9.0 * fneqf_b[topbottom(i - 1, k - 3, a)] - 9.0 * fneqf_b[topbottom(i + 1, k + 3, a)] + 81.0 * fneqf_b[topbottom(i + 1, k + 1, a)] + 81.0 * fneqf_b[topbottom(i + 1, k - 1, a)] - 9.0 * fneqf_b[topbottom(i + 1, k - 3, a)] + fneqf_b[topbottom(i + 3, k + 3, a)] - 9.0 * fneqf_b[topbottom(i + 3, k + 1, a)] - 9.0 * fneqf_b[topbottom(i + 3, k - 1, a)] + fneqf_b[topbottom(i + 3, k - 3, a)]) / 256.0;

                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];
                }
            }
        }

        //****Corner points topbottom side
        for (k = 1; k <= nz - 1; k += alpha)
        {
            i = 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i - 1, k, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i + 1, k, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i + 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];
            }

            i = nx - 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i + 1, k, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i - 1, k, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i - 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];
            }
        }

        for (i = 1; i <= nx - 1; i += alpha)
        {
            k = 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i, k - 1, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i, k + 1, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i, k + 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];
            }

            k = nz - 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_b[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_b[topbottom(i, k + 1, a)] + (3.0 / 4.0) * fneqf_b[topbottom(i, k - 1, a)] - (1.0 / 8.0) * fneqf_b[topbottom(i, k - 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_b[topbottom(i, k, a)];
            }
        }
    }
    if (side == 2)
    {

        for (k = 0, i = face; k <= nz; k += alpha)
        {
            for (j = 0; j <= ny; j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    fneqf_r[leftright(j, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
                }
            }
        }

        for (i = face, j = 3; j <= ny - 3; j += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                jxf[fineisnindex1(i, j, k)] = (jxf[fineisnindex1(i, j - 3, k + 3)] - 9.0 * jxf[fineisnindex1(i, j - 3, k + 1)] - 9.0 * jxf[fineisnindex1(i, j - 3, k - 1)] + jxf[fineisnindex1(i, j - 3, k - 3)] - 9.0 * jxf[fineisnindex1(i, j - 1, k + 3)] + 81.0 * jxf[fineisnindex1(i, j - 1, k + 1)] + 81.0 * jxf[fineisnindex1(i, j - 1, k - 1)] - 9.0 * jxf[fineisnindex1(i, j - 1, k - 3)] - 9.0 * jxf[fineisnindex1(i, j + 1, k + 3)] + 81.0 * jxf[fineisnindex1(i, j + 1, k + 1)] + 81.0 * jxf[fineisnindex1(i, j + 1, k - 1)] - 9.0 * jxf[fineisnindex1(i, j + 1, k - 3)] + jxf[fineisnindex1(i, j + 3, k + 3)] - 9.0 * jxf[fineisnindex1(i, j + 3, k + 1)] - 9.0 * jxf[fineisnindex1(i, j + 3, k - 1)] + jxf[fineisnindex1(i, j + 3, k - 3)]) / 256.0;
                jyf[fineisnindex1(i, j, k)] = (jyf[fineisnindex1(i, j - 3, k + 3)] - 9.0 * jyf[fineisnindex1(i, j - 3, k + 1)] - 9.0 * jyf[fineisnindex1(i, j - 3, k - 1)] + jyf[fineisnindex1(i, j - 3, k - 3)] - 9.0 * jyf[fineisnindex1(i, j - 1, k + 3)] + 81.0 * jyf[fineisnindex1(i, j - 1, k + 1)] + 81.0 * jyf[fineisnindex1(i, j - 1, k - 1)] - 9.0 * jyf[fineisnindex1(i, j - 1, k - 3)] - 9.0 * jyf[fineisnindex1(i, j + 1, k + 3)] + 81.0 * jyf[fineisnindex1(i, j + 1, k + 1)] + 81.0 * jyf[fineisnindex1(i, j + 1, k - 1)] - 9.0 * jyf[fineisnindex1(i, j + 1, k - 3)] + jyf[fineisnindex1(i, j + 3, k + 3)] - 9.0 * jyf[fineisnindex1(i, j + 3, k + 1)] - 9.0 * jyf[fineisnindex1(i, j + 3, k - 1)] + jyf[fineisnindex1(i, j + 3, k - 3)]) / 256.0;
                jzf[fineisnindex1(i, j, k)] = (jzf[fineisnindex1(i, j - 3, k + 3)] - 9.0 * jzf[fineisnindex1(i, j - 3, k + 1)] - 9.0 * jzf[fineisnindex1(i, j - 3, k - 1)] + jzf[fineisnindex1(i, j - 3, k - 3)] - 9.0 * jzf[fineisnindex1(i, j - 1, k + 3)] + 81.0 * jzf[fineisnindex1(i, j - 1, k + 1)] + 81.0 * jzf[fineisnindex1(i, j - 1, k - 1)] - 9.0 * jzf[fineisnindex1(i, j - 1, k - 3)] - 9.0 * jzf[fineisnindex1(i, j + 1, k + 3)] + 81.0 * jzf[fineisnindex1(i, j + 1, k + 1)] + 81.0 * jzf[fineisnindex1(i, j + 1, k - 1)] - 9.0 * jzf[fineisnindex1(i, j + 1, k - 3)] + jzf[fineisnindex1(i, j + 3, k + 3)] - 9.0 * jzf[fineisnindex1(i, j + 3, k + 1)] - 9.0 * jzf[fineisnindex1(i, j + 3, k - 1)] + jzf[fineisnindex1(i, j + 3, k - 3)]) / 256.0;
                rhof[fineisnindex1(i, j, k)] = (rhof[fineisnindex1(i, j - 3, k + 3)] - 9.0 * rhof[fineisnindex1(i, j - 3, k + 1)] - 9.0 * rhof[fineisnindex1(i, j - 3, k - 1)] + rhof[fineisnindex1(i, j - 3, k - 3)] - 9.0 * rhof[fineisnindex1(i, j - 1, k + 3)] + 81.0 * rhof[fineisnindex1(i, j - 1, k + 1)] + 81.0 * rhof[fineisnindex1(i, j - 1, k - 1)] - 9.0 * rhof[fineisnindex1(i, j - 1, k - 3)] - 9.0 * rhof[fineisnindex1(i, j + 1, k + 3)] + 81.0 * rhof[fineisnindex1(i, j + 1, k + 1)] + 81.0 * rhof[fineisnindex1(i, j + 1, k - 1)] - 9.0 * rhof[fineisnindex1(i, j + 1, k - 3)] + rhof[fineisnindex1(i, j + 3, k + 3)] - 9.0 * rhof[fineisnindex1(i, j + 3, k + 1)] - 9.0 * rhof[fineisnindex1(i, j + 3, k - 1)] + rhof[fineisnindex1(i, j + 3, k - 3)]) / 256.0;

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    fneqf_r[leftright(j, k, a)] = (fneqf_r[leftright(j - 3, k + 3, a)] - 9.0 * fneqf_r[leftright(j - 3, k + 1, a)] - 9.0 * fneqf_r[leftright(j - 3, k - 1, a)] + fneqf_r[leftright(j - 3, k - 3, a)] - 9.0 * fneqf_r[leftright(j - 1, k + 3, a)] + 81.0 * fneqf_r[leftright(j - 1, k + 1, a)] + 81.0 * fneqf_r[leftright(j - 1, k - 1, a)] - 9.0 * fneqf_r[leftright(j - 1, k - 3, a)] - 9.0 * fneqf_r[leftright(j + 1, k + 3, a)] + 81.0 * fneqf_r[leftright(j + 1, k + 1, a)] + 81.0 * fneqf_r[leftright(j + 1, k - 1, a)] - 9.0 * fneqf_r[leftright(j + 1, k - 3, a)] + fneqf_r[leftright(j + 3, k + 3, a)] - 9.0 * fneqf_r[leftright(j + 3, k + 1, a)] - 9.0 * fneqf_r[leftright(j + 3, k - 1, a)] + fneqf_r[leftright(j + 3, k - 3, a)]) / 256.0;

                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];
                }
            }
        }

        //******Corner points leftright Side
        for (k = 1; k <= nz - 1; k += alpha)
        {
            i = face;
            j = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j - 1, k, a)] + (3.0 / 4.0) * fneqf_r[leftright(j + 1, k, a)] - (1.0 / 8.0) * fneqf_r[leftright(j + 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];
            }

            i = face;
            j = ny - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j + 1, k, a)] + (3.0 / 4.0) * fneqf_r[leftright(j - 1, k, a)] - (1.0 / 8.0) * fneqf_r[leftright(j - 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];
            }
        }
        for (j = 1; j <= ny - 1; j += alpha)
        {
            i = face;
            k = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j, k - 1, a)] + (3.0 / 4.0) * fneqf_r[leftright(j, k + 1, a)] - (1.0 / 8.0) * fneqf_r[leftright(j, k + 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];
            }

            i = face;
            k = nz - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_r[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_r[leftright(j, k + 1, a)] + (3.0 / 4.0) * fneqf_r[leftright(j, k - 1, a)] - (1.0 / 8.0) * fneqf_r[leftright(j, k - 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_r[leftright(j, k, a)];
            }
        }
    }
    if (side == 3)
    {
        for (i = 0, j = face; i <= nx; i += alpha)
        {
            for (k = 0; k <= nz; k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    fneqf_t[topbottom(i, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
                }
            }
        }

        for (j = face, i = 3; i <= nx - 3; i += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                jxf[fineisnindex1(i, j, k)] = (jxf[fineisnindex1(i - 3, j, k + 3)] - 9.0 * jxf[fineisnindex1(i - 3, j, k + 1)] - 9.0 * jxf[fineisnindex1(i - 3, j, k - 1)] + jxf[fineisnindex1(i - 3, j, k - 3)] - 9.0 * jxf[fineisnindex1(i - 1, j, k + 3)] + 81.0 * jxf[fineisnindex1(i - 1, j, k + 1)] + 81.0 * jxf[fineisnindex1(i - 1, j, k - 1)] - 9.0 * jxf[fineisnindex1(i - 1, j, k - 3)] - 9.0 * jxf[fineisnindex1(i + 1, j, k + 3)] + 81.0 * jxf[fineisnindex1(i + 1, j, k + 1)] + 81.0 * jxf[fineisnindex1(i + 1, j, k - 1)] - 9.0 * jxf[fineisnindex1(i + 1, j, k - 3)] + jxf[fineisnindex1(i + 3, j, k + 3)] - 9.0 * jxf[fineisnindex1(i + 3, j, k + 1)] - 9.0 * jxf[fineisnindex1(i + 3, j, k - 1)] + jxf[fineisnindex1(i + 3, j, k - 3)]) / 256.0;
                jyf[fineisnindex1(i, j, k)] = (jyf[fineisnindex1(i - 3, j, k + 3)] - 9.0 * jyf[fineisnindex1(i - 3, j, k + 1)] - 9.0 * jyf[fineisnindex1(i - 3, j, k - 1)] + jyf[fineisnindex1(i - 3, j, k - 3)] - 9.0 * jyf[fineisnindex1(i - 1, j, k + 3)] + 81.0 * jyf[fineisnindex1(i - 1, j, k + 1)] + 81.0 * jyf[fineisnindex1(i - 1, j, k - 1)] - 9.0 * jyf[fineisnindex1(i - 1, j, k - 3)] - 9.0 * jyf[fineisnindex1(i + 1, j, k + 3)] + 81.0 * jyf[fineisnindex1(i + 1, j, k + 1)] + 81.0 * jyf[fineisnindex1(i + 1, j, k - 1)] - 9.0 * jyf[fineisnindex1(i + 1, j, k - 3)] + jyf[fineisnindex1(i + 3, j, k + 3)] - 9.0 * jyf[fineisnindex1(i + 3, j, k + 1)] - 9.0 * jyf[fineisnindex1(i + 3, j, k - 1)] + jyf[fineisnindex1(i + 3, j, k - 3)]) / 256.0;
                jzf[fineisnindex1(i, j, k)] = (jzf[fineisnindex1(i - 3, j, k + 3)] - 9.0 * jzf[fineisnindex1(i - 3, j, k + 1)] - 9.0 * jzf[fineisnindex1(i - 3, j, k - 1)] + jzf[fineisnindex1(i - 3, j, k - 3)] - 9.0 * jzf[fineisnindex1(i - 1, j, k + 3)] + 81.0 * jzf[fineisnindex1(i - 1, j, k + 1)] + 81.0 * jzf[fineisnindex1(i - 1, j, k - 1)] - 9.0 * jzf[fineisnindex1(i - 1, j, k - 3)] - 9.0 * jzf[fineisnindex1(i + 1, j, k + 3)] + 81.0 * jzf[fineisnindex1(i + 1, j, k + 1)] + 81.0 * jzf[fineisnindex1(i + 1, j, k - 1)] - 9.0 * jzf[fineisnindex1(i + 1, j, k - 3)] + jzf[fineisnindex1(i + 3, j, k + 3)] - 9.0 * jzf[fineisnindex1(i + 3, j, k + 1)] - 9.0 * jzf[fineisnindex1(i + 3, j, k - 1)] + jzf[fineisnindex1(i + 3, j, k - 3)]) / 256.0;
                rhof[fineisnindex1(i, j, k)] = (rhof[fineisnindex1(i - 3, j, k + 3)] - 9.0 * rhof[fineisnindex1(i - 3, j, k + 1)] - 9.0 * rhof[fineisnindex1(i - 3, j, k - 1)] + rhof[fineisnindex1(i - 3, j, k - 3)] - 9.0 * rhof[fineisnindex1(i - 1, j, k + 3)] + 81.0 * rhof[fineisnindex1(i - 1, j, k + 1)] + 81.0 * rhof[fineisnindex1(i - 1, j, k - 1)] - 9.0 * rhof[fineisnindex1(i - 1, j, k - 3)] - 9.0 * rhof[fineisnindex1(i + 1, j, k + 3)] + 81.0 * rhof[fineisnindex1(i + 1, j, k + 1)] + 81.0 * rhof[fineisnindex1(i + 1, j, k - 1)] - 9.0 * rhof[fineisnindex1(i + 1, j, k - 3)] + rhof[fineisnindex1(i + 3, j, k + 3)] - 9.0 * rhof[fineisnindex1(i + 3, j, k + 1)] - 9.0 * rhof[fineisnindex1(i + 3, j, k - 1)] + rhof[fineisnindex1(i + 3, j, k - 3)]) / 256.0;

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    fneqf_t[topbottom(i, k, a)] = (fneqf_t[topbottom(i - 3, k + 3, a)] - 9.0 * fneqf_t[topbottom(i - 3, k + 1, a)] - 9.0 * fneqf_t[topbottom(i - 3, k - 1, a)] + fneqf_t[topbottom(i - 3, k - 3, a)] - 9.0 * fneqf_t[topbottom(i - 1, k + 3, a)] + 81.0 * fneqf_t[topbottom(i - 1, k + 1, a)] + 81.0 * fneqf_t[topbottom(i - 1, k - 1, a)] - 9.0 * fneqf_t[topbottom(i - 1, k - 3, a)] - 9.0 * fneqf_t[topbottom(i + 1, k + 3, a)] + 81.0 * fneqf_t[topbottom(i + 1, k + 1, a)] + 81.0 * fneqf_t[topbottom(i + 1, k - 1, a)] - 9.0 * fneqf_t[topbottom(i + 1, k - 3, a)] + fneqf_t[topbottom(i + 3, k + 3, a)] - 9.0 * fneqf_t[topbottom(i + 3, k + 1, a)] - 9.0 * fneqf_t[topbottom(i + 3, k - 1, a)] + fneqf_t[topbottom(i + 3, k - 3, a)]) / 256.0;

                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];
                }
            }
        }

        //********Corner points Top side
        for (k = 1; k <= nz - 1; k += alpha)
        {
            i = 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i - 1, k, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i + 1, k, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i + 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];
            }

            i = nx - 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i + 1, k, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i - 1, k, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i - 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];
            }
        }

        for (i = 1; i <= nx - 1; i += alpha)
        {
            k = 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i, k - 1, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i, k + 1, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i, k + 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];
            }

            k = nz - 1;
            j = face;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_t[topbottom(i, k, a)] = (3.0 / 8.0) * fneqf_t[topbottom(i, k + 1, a)] + (3.0 / 4.0) * fneqf_t[topbottom(i, k - 1, a)] - (1.0 / 8.0) * fneqf_t[topbottom(i, k - 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_t[topbottom(i, k, a)];
            }
        }
    }
    if (side == 4)
    {
        for (j = 0, i = face; j <= ny; j += alpha)
        {
            for (k = 0; k <= nz; k += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    fneqf_l[leftright(j, k, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
                }
            }
        }

        for (i = face, j = 3; j <= ny - 3; j += alpha)
        {
            for (k = 3; k <= nz - 3; k += alpha)
            {
                jxf[fineisnindex1(i, j, k)] = (jxf[fineisnindex1(i, j - 3, k + 3)] - 9.0 * jxf[fineisnindex1(i, j - 3, k + 1)] - 9.0 * jxf[fineisnindex1(i, j - 3, k - 1)] + jxf[fineisnindex1(i, j - 3, k - 3)] - 9.0 * jxf[fineisnindex1(i, j - 1, k + 3)] + 81.0 * jxf[fineisnindex1(i, j - 1, k + 1)] + 81.0 * jxf[fineisnindex1(i, j - 1, k - 1)] - 9.0 * jxf[fineisnindex1(i, j - 1, k - 3)] - 9.0 * jxf[fineisnindex1(i, j + 1, k + 3)] + 81.0 * jxf[fineisnindex1(i, j + 1, k + 1)] + 81.0 * jxf[fineisnindex1(i, j + 1, k - 1)] - 9.0 * jxf[fineisnindex1(i, j + 1, k - 3)] + jxf[fineisnindex1(i, j + 3, k + 3)] - 9.0 * jxf[fineisnindex1(i, j + 3, k + 1)] - 9.0 * jxf[fineisnindex1(i, j + 3, k - 1)] + jxf[fineisnindex1(i, j + 3, k - 3)]) / 256.0;
                jyf[fineisnindex1(i, j, k)] = (jyf[fineisnindex1(i, j - 3, k + 3)] - 9.0 * jyf[fineisnindex1(i, j - 3, k + 1)] - 9.0 * jyf[fineisnindex1(i, j - 3, k - 1)] + jyf[fineisnindex1(i, j - 3, k - 3)] - 9.0 * jyf[fineisnindex1(i, j - 1, k + 3)] + 81.0 * jyf[fineisnindex1(i, j - 1, k + 1)] + 81.0 * jyf[fineisnindex1(i, j - 1, k - 1)] - 9.0 * jyf[fineisnindex1(i, j - 1, k - 3)] - 9.0 * jyf[fineisnindex1(i, j + 1, k + 3)] + 81.0 * jyf[fineisnindex1(i, j + 1, k + 1)] + 81.0 * jyf[fineisnindex1(i, j + 1, k - 1)] - 9.0 * jyf[fineisnindex1(i, j + 1, k - 3)] + jyf[fineisnindex1(i, j + 3, k + 3)] - 9.0 * jyf[fineisnindex1(i, j + 3, k + 1)] - 9.0 * jyf[fineisnindex1(i, j + 3, k - 1)] + jyf[fineisnindex1(i, j + 3, k - 3)]) / 256.0;
                jzf[fineisnindex1(i, j, k)] = (jzf[fineisnindex1(i, j - 3, k + 3)] - 9.0 * jzf[fineisnindex1(i, j - 3, k + 1)] - 9.0 * jzf[fineisnindex1(i, j - 3, k - 1)] + jzf[fineisnindex1(i, j - 3, k - 3)] - 9.0 * jzf[fineisnindex1(i, j - 1, k + 3)] + 81.0 * jzf[fineisnindex1(i, j - 1, k + 1)] + 81.0 * jzf[fineisnindex1(i, j - 1, k - 1)] - 9.0 * jzf[fineisnindex1(i, j - 1, k - 3)] - 9.0 * jzf[fineisnindex1(i, j + 1, k + 3)] + 81.0 * jzf[fineisnindex1(i, j + 1, k + 1)] + 81.0 * jzf[fineisnindex1(i, j + 1, k - 1)] - 9.0 * jzf[fineisnindex1(i, j + 1, k - 3)] + jzf[fineisnindex1(i, j + 3, k + 3)] - 9.0 * jzf[fineisnindex1(i, j + 3, k + 1)] - 9.0 * jzf[fineisnindex1(i, j + 3, k - 1)] + jzf[fineisnindex1(i, j + 3, k - 3)]) / 256.0;
                rhof[fineisnindex1(i, j, k)] = (rhof[fineisnindex1(i, j - 3, k + 3)] - 9.0 * rhof[fineisnindex1(i, j - 3, k + 1)] - 9.0 * rhof[fineisnindex1(i, j - 3, k - 1)] + rhof[fineisnindex1(i, j - 3, k - 3)] - 9.0 * rhof[fineisnindex1(i, j - 1, k + 3)] + 81.0 * rhof[fineisnindex1(i, j - 1, k + 1)] + 81.0 * rhof[fineisnindex1(i, j - 1, k - 1)] - 9.0 * rhof[fineisnindex1(i, j - 1, k - 3)] - 9.0 * rhof[fineisnindex1(i, j + 1, k + 3)] + 81.0 * rhof[fineisnindex1(i, j + 1, k + 1)] + 81.0 * rhof[fineisnindex1(i, j + 1, k - 1)] - 9.0 * rhof[fineisnindex1(i, j + 1, k - 3)] + rhof[fineisnindex1(i, j + 3, k + 3)] - 9.0 * rhof[fineisnindex1(i, j + 3, k + 1)] - 9.0 * rhof[fineisnindex1(i, j + 3, k - 1)] + rhof[fineisnindex1(i, j + 3, k - 3)]) / 256.0;

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    fneqf_l[leftright(j, k, a)] = (fneqf_l[leftright(j - 3, k + 3, a)] - 9.0 * fneqf_l[leftright(j - 3, k + 1, a)] - 9.0 * fneqf_l[leftright(j - 3, k - 1, a)] + fneqf_l[leftright(j - 3, k - 3, a)] - 9.0 * fneqf_l[leftright(j - 1, k + 3, a)] + 81.0 * fneqf_l[leftright(j - 1, k + 1, a)] + 81.0 * fneqf_l[leftright(j - 1, k - 1, a)] - 9.0 * fneqf_l[leftright(j - 1, k - 3, a)] - 9.0 * fneqf_l[leftright(j + 1, k + 3, a)] + 81.0 * fneqf_l[leftright(j + 1, k + 1, a)] + 81.0 * fneqf_l[leftright(j + 1, k - 1, a)] - 9.0 * fneqf_l[leftright(j + 1, k - 3, a)] + fneqf_l[leftright(j + 3, k + 3, a)] - 9.0 * fneqf_l[leftright(j + 3, k + 1, a)] - 9.0 * fneqf_l[leftright(j + 3, k - 1, a)] + fneqf_l[leftright(j + 3, k - 3, a)]) / 256.0;

                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];
                }
            }
        }

        //********Corner points Left side
        for (k = 1; k <= nz - 1; k += alpha)
        {
            i = face;
            j = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j - 1, k, a)] + (3.0 / 4.0) * fneqf_l[leftright(j + 1, k, a)] - (1.0 / 8.0) * fneqf_l[leftright(j + 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];
            }

            i = face;
            j = ny - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j + 1, k, a)] + (3.0 / 4.0) * fneqf_l[leftright(j - 1, k, a)] - (1.0 / 8.0) * fneqf_l[leftright(j - 3, k, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];
            }
        }

        for (j = 1; j <= ny - 1; j += alpha)
        {
            i = face;
            k = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k + 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k + 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k + 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k - 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k + 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k + 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j, k - 1, a)] + (3.0 / 4.0) * fneqf_l[leftright(j, k + 1, a)] - (1.0 / 8.0) * fneqf_l[leftright(j, k + 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];
            }

            i = face;
            k = nz - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j, k - 3)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j, k - 3)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j, k - 3)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j, k + 1)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j, k - 1)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j, k - 3)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_l[leftright(j, k, a)] = (3.0 / 8.0) * fneqf_l[leftright(j, k + 1, a)] + (3.0 / 4.0) * fneqf_l[leftright(j, k - 1, a)] - (1.0 / 8.0) * fneqf_l[leftright(j, k - 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_l[leftright(j, k, a)];
            }
        }
    }
    if (side == 5)
    {
        for (i = 0, k = face; i <= nx; i += alpha)
        {
            for (j = 0; j <= ny; j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    fneqf_d[updown(i, j, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
                }
            }
        }

        for (k = face, i = 3; i <= nx - 3; i += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                jxf[fineisnindex1(i, j, k)] = (jxf[fineisnindex1(i - 3, j + 3, k)] - 9.0 * jxf[fineisnindex1(i - 3, j + 1, k)] - 9.0 * jxf[fineisnindex1(i - 3, j - 1, k)] + jxf[fineisnindex1(i - 3, j - 3, k)] - 9.0 * jxf[fineisnindex1(i - 1, j + 3, k)] + 81.0 * jxf[fineisnindex1(i - 1, j + 1, k)] + 81.0 * jxf[fineisnindex1(i - 1, j - 1, k)] - 9.0 * jxf[fineisnindex1(i - 1, j - 3, k)] - 9.0 * jxf[fineisnindex1(i + 1, j + 3, k)] + 81.0 * jxf[fineisnindex1(i + 1, j + 1, k)] + 81.0 * jxf[fineisnindex1(i + 1, j - 1, k)] - 9.0 * jxf[fineisnindex1(i + 1, j - 3, k)] + jxf[fineisnindex1(i + 3, j + 3, k)] - 9.0 * jxf[fineisnindex1(i + 3, j + 1, k)] - 9.0 * jxf[fineisnindex1(i + 3, j - 1, k)] + jxf[fineisnindex1(i + 3, j - 3, k)]) / 256.0;
                jyf[fineisnindex1(i, j, k)] = (jyf[fineisnindex1(i - 3, j + 3, k)] - 9.0 * jyf[fineisnindex1(i - 3, j + 1, k)] - 9.0 * jyf[fineisnindex1(i - 3, j - 1, k)] + jyf[fineisnindex1(i - 3, j - 3, k)] - 9.0 * jyf[fineisnindex1(i - 1, j + 3, k)] + 81.0 * jyf[fineisnindex1(i - 1, j + 1, k)] + 81.0 * jyf[fineisnindex1(i - 1, j - 1, k)] - 9.0 * jyf[fineisnindex1(i - 1, j - 3, k)] - 9.0 * jyf[fineisnindex1(i + 1, j + 3, k)] + 81.0 * jyf[fineisnindex1(i + 1, j + 1, k)] + 81.0 * jyf[fineisnindex1(i + 1, j - 1, k)] - 9.0 * jyf[fineisnindex1(i + 1, j - 3, k)] + jyf[fineisnindex1(i + 3, j + 3, k)] - 9.0 * jyf[fineisnindex1(i + 3, j + 1, k)] - 9.0 * jyf[fineisnindex1(i + 3, j - 1, k)] + jyf[fineisnindex1(i + 3, j - 3, k)]) / 256.0;
                jzf[fineisnindex1(i, j, k)] = (jzf[fineisnindex1(i - 3, j + 3, k)] - 9.0 * jzf[fineisnindex1(i - 3, j + 1, k)] - 9.0 * jzf[fineisnindex1(i - 3, j - 1, k)] + jzf[fineisnindex1(i - 3, j - 3, k)] - 9.0 * jzf[fineisnindex1(i - 1, j + 3, k)] + 81.0 * jzf[fineisnindex1(i - 1, j + 1, k)] + 81.0 * jzf[fineisnindex1(i - 1, j - 1, k)] - 9.0 * jzf[fineisnindex1(i - 1, j - 3, k)] - 9.0 * jzf[fineisnindex1(i + 1, j + 3, k)] + 81.0 * jzf[fineisnindex1(i + 1, j + 1, k)] + 81.0 * jzf[fineisnindex1(i + 1, j - 1, k)] - 9.0 * jzf[fineisnindex1(i + 1, j - 3, k)] + jzf[fineisnindex1(i + 3, j + 3, k)] - 9.0 * jzf[fineisnindex1(i + 3, j + 1, k)] - 9.0 * jzf[fineisnindex1(i + 3, j - 1, k)] + jzf[fineisnindex1(i + 3, j - 3, k)]) / 256.0;
                rhof[fineisnindex1(i, j, k)] = (rhof[fineisnindex1(i - 3, j + 3, k)] - 9.0 * rhof[fineisnindex1(i - 3, j + 1, k)] - 9.0 * rhof[fineisnindex1(i - 3, j - 1, k)] + rhof[fineisnindex1(i - 3, j - 3, k)] - 9.0 * rhof[fineisnindex1(i - 1, j + 3, k)] + 81.0 * rhof[fineisnindex1(i - 1, j + 1, k)] + 81.0 * rhof[fineisnindex1(i - 1, j - 1, k)] - 9.0 * rhof[fineisnindex1(i - 1, j - 3, k)] - 9.0 * rhof[fineisnindex1(i + 1, j + 3, k)] + 81.0 * rhof[fineisnindex1(i + 1, j + 1, k)] + 81.0 * rhof[fineisnindex1(i + 1, j - 1, k)] - 9.0 * rhof[fineisnindex1(i + 1, j - 3, k)] + rhof[fineisnindex1(i + 3, j + 3, k)] - 9.0 * rhof[fineisnindex1(i + 3, j + 1, k)] - 9.0 * rhof[fineisnindex1(i + 3, j - 1, k)] + rhof[fineisnindex1(i + 3, j - 3, k)]) / 256.0;

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    fneqf_d[updown(i, j, a)] = (fneqf_d[updown(i - 3, j + 3, a)] - 9.0 * fneqf_d[updown(i - 3, j + 1, a)] - 9.0 * fneqf_d[updown(i - 3, j - 1, a)] + fneqf_d[updown(i - 3, j - 3, a)] - 9.0 * fneqf_d[updown(i - 1, j + 3, a)] + 81.0 * fneqf_d[updown(i - 1, j + 1, a)] + 81.0 * fneqf_d[updown(i - 1, j - 1, a)] - 9.0 * fneqf_d[updown(i - 1, j - 3, a)] - 9.0 * fneqf_d[updown(i + 1, j + 3, a)] + 81.0 * fneqf_d[updown(i + 1, j + 1, a)] + 81.0 * fneqf_d[updown(i + 1, j - 1, a)] - 9.0 * fneqf_d[updown(i + 1, j - 3, a)] + fneqf_d[updown(i + 3, j + 3, a)] - 9.0 * fneqf_d[updown(i + 3, j + 1, a)] - 9.0 * fneqf_d[updown(i + 3, j - 1, a)] + fneqf_d[updown(i + 3, j - 3, a)]) / 256.0;

                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];
                }
            }
        }

        //********Corner points down side
        for (i = 1; i <= nx - 1; i += alpha)
        {
            k = face;
            j = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i, j - 1, a)] + (3.0 / 4.0) * fneqf_d[updown(i, j + 1, a)] - (1.0 / 8.0) * fneqf_d[updown(i, j + 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];
            }

            k = face;
            j = ny - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i, j + 1, a)] + (3.0 / 4.0) * fneqf_d[updown(i, j - 1, a)] - (1.0 / 8.0) * fneqf_d[updown(i, j - 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];
            }
        }

        for (j = 1; j <= ny - 1; j += alpha)
        {
            k = face;
            i = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i - 1, j, a)] + (3.0 / 4.0) * fneqf_d[updown(i + 1, j, a)] - (1.0 / 8.0) * fneqf_d[updown(i + 3, j, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];
            }

            k = face;
            i = nx - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_d[updown(i, j, a)] = (3.0 / 8.0) * fneqf_d[updown(i + 1, j, a)] + (3.0 / 4.0) * fneqf_d[updown(i - 1, j, a)] - (1.0 / 8.0) * fneqf_d[updown(i - 3, j, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_d[updown(i, j, a)];
            }
        }
    }
    if (side == 6)
    {
        for (i = 0, k = face; i <= nx; i += alpha)
        {
            for (j = 0; j <= ny; j += alpha)
            {
                l = 0.0;
                m = 0.0;
                n = 0.0;
                p = 0.0;

                for (a = 0; a < 19; a++)
                {
                    l += ff[fineindex(i, j, k, a)];
                    m += ex[a] * ff[fineindex(i, j, k, a)];
                    n += ey[a] * ff[fineindex(i, j, k, a)];
                    p += ez[a] * ff[fineindex(i, j, k, a)];
                }
                rhof[fineisnindex1(i, j, k)] = l;
                jxf[fineisnindex1(i, j, k)] = m / rhof[fineisnindex1(i, j, k)];
                jyf[fineisnindex1(i, j, k)] = n / rhof[fineisnindex1(i, j, k)];
                jzf[fineisnindex1(i, j, k)] = p / rhof[fineisnindex1(i, j, k)];

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    fneqf_u[updown(i, j, a)] = ff[fineindex(i, j, k, a)] - feqf[fineindex(i, j, k, a)];
                }
            }
        }

        for (k = face, i = 3; i <= nx - 3; i += alpha)
        {
            for (j = 3; j <= ny - 3; j += alpha)
            {
                jxf[fineisnindex1(i, j, k)] = (jxf[fineisnindex1(i - 3, j + 3, k)] - 9.0 * jxf[fineisnindex1(i - 3, j + 1, k)] - 9.0 * jxf[fineisnindex1(i - 3, j - 1, k)] + jxf[fineisnindex1(i - 3, j - 3, k)] - 9.0 * jxf[fineisnindex1(i - 1, j + 3, k)] + 81.0 * jxf[fineisnindex1(i - 1, j + 1, k)] + 81.0 * jxf[fineisnindex1(i - 1, j - 1, k)] - 9.0 * jxf[fineisnindex1(i - 1, j - 3, k)] - 9.0 * jxf[fineisnindex1(i + 1, j + 3, k)] + 81.0 * jxf[fineisnindex1(i + 1, j + 1, k)] + 81.0 * jxf[fineisnindex1(i + 1, j - 1, k)] - 9.0 * jxf[fineisnindex1(i + 1, j - 3, k)] + jxf[fineisnindex1(i + 3, j + 3, k)] - 9.0 * jxf[fineisnindex1(i + 3, j + 1, k)] - 9.0 * jxf[fineisnindex1(i + 3, j - 1, k)] + jxf[fineisnindex1(i + 3, j - 3, k)]) / 256.0;
                jyf[fineisnindex1(i, j, k)] = (jyf[fineisnindex1(i - 3, j + 3, k)] - 9.0 * jyf[fineisnindex1(i - 3, j + 1, k)] - 9.0 * jyf[fineisnindex1(i - 3, j - 1, k)] + jyf[fineisnindex1(i - 3, j - 3, k)] - 9.0 * jyf[fineisnindex1(i - 1, j + 3, k)] + 81.0 * jyf[fineisnindex1(i - 1, j + 1, k)] + 81.0 * jyf[fineisnindex1(i - 1, j - 1, k)] - 9.0 * jyf[fineisnindex1(i - 1, j - 3, k)] - 9.0 * jyf[fineisnindex1(i + 1, j + 3, k)] + 81.0 * jyf[fineisnindex1(i + 1, j + 1, k)] + 81.0 * jyf[fineisnindex1(i + 1, j - 1, k)] - 9.0 * jyf[fineisnindex1(i + 1, j - 3, k)] + jyf[fineisnindex1(i + 3, j + 3, k)] - 9.0 * jyf[fineisnindex1(i + 3, j + 1, k)] - 9.0 * jyf[fineisnindex1(i + 3, j - 1, k)] + jyf[fineisnindex1(i + 3, j - 3, k)]) / 256.0;
                jzf[fineisnindex1(i, j, k)] = (jzf[fineisnindex1(i - 3, j + 3, k)] - 9.0 * jzf[fineisnindex1(i - 3, j + 1, k)] - 9.0 * jzf[fineisnindex1(i - 3, j - 1, k)] + jzf[fineisnindex1(i - 3, j - 3, k)] - 9.0 * jzf[fineisnindex1(i - 1, j + 3, k)] + 81.0 * jzf[fineisnindex1(i - 1, j + 1, k)] + 81.0 * jzf[fineisnindex1(i - 1, j - 1, k)] - 9.0 * jzf[fineisnindex1(i - 1, j - 3, k)] - 9.0 * jzf[fineisnindex1(i + 1, j + 3, k)] + 81.0 * jzf[fineisnindex1(i + 1, j + 1, k)] + 81.0 * jzf[fineisnindex1(i + 1, j - 1, k)] - 9.0 * jzf[fineisnindex1(i + 1, j - 3, k)] + jzf[fineisnindex1(i + 3, j + 3, k)] - 9.0 * jzf[fineisnindex1(i + 3, j + 1, k)] - 9.0 * jzf[fineisnindex1(i + 3, j - 1, k)] + jzf[fineisnindex1(i + 3, j - 3, k)]) / 256.0;
                rhof[fineisnindex1(i, j, k)] = (rhof[fineisnindex1(i - 3, j + 3, k)] - 9.0 * rhof[fineisnindex1(i - 3, j + 1, k)] - 9.0 * rhof[fineisnindex1(i - 3, j - 1, k)] + rhof[fineisnindex1(i - 3, j - 3, k)] - 9.0 * rhof[fineisnindex1(i - 1, j + 3, k)] + 81.0 * rhof[fineisnindex1(i - 1, j + 1, k)] + 81.0 * rhof[fineisnindex1(i - 1, j - 1, k)] - 9.0 * rhof[fineisnindex1(i - 1, j - 3, k)] - 9.0 * rhof[fineisnindex1(i + 1, j + 3, k)] + 81.0 * rhof[fineisnindex1(i + 1, j + 1, k)] + 81.0 * rhof[fineisnindex1(i + 1, j - 1, k)] - 9.0 * rhof[fineisnindex1(i + 1, j - 3, k)] + rhof[fineisnindex1(i + 3, j + 3, k)] - 9.0 * rhof[fineisnindex1(i + 3, j + 1, k)] - 9.0 * rhof[fineisnindex1(i + 3, j - 1, k)] + rhof[fineisnindex1(i + 3, j - 3, k)]) / 256.0;

                u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

                for (a = 0; a < 19; a++)
                {
                    fneqf_u[updown(i, j, a)] = (fneqf_u[updown(i - 3, j + 3, a)] - 9.0 * fneqf_u[updown(i - 3, j + 1, a)] - 9.0 * fneqf_u[updown(i - 3, j - 1, a)] + fneqf_u[updown(i - 3, j - 3, a)] - 9.0 * fneqf_u[updown(i - 1, j + 3, a)] + 81.0 * fneqf_u[updown(i - 1, j + 1, a)] + 81.0 * fneqf_u[updown(i - 1, j - 1, a)] - 9.0 * fneqf_u[updown(i - 1, j - 3, a)] - 9.0 * fneqf_u[updown(i + 1, j + 3, a)] + 81.0 * fneqf_u[updown(i + 1, j + 1, a)] + 81.0 * fneqf_u[updown(i + 1, j - 1, a)] - 9.0 * fneqf_u[updown(i + 1, j - 3, a)] + fneqf_u[updown(i + 3, j + 3, a)] - 9.0 * fneqf_u[updown(i + 3, j + 1, a)] - 9.0 * fneqf_u[updown(i + 3, j - 1, a)] + fneqf_u[updown(i + 3, j - 3, a)]) / 256.0;

                    t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                    feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                    ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];
                }
            }
        }

        //********Corner points up side
        for (i = 1; i <= nx - 1; i += alpha)
        {
            k = face;
            j = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j + 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j + 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j + 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j - 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j + 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j + 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i, j - 1, a)] + (3.0 / 4.0) * fneqf_u[updown(i, j + 1, a)] - (1.0 / 8.0) * fneqf_u[updown(i, j + 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];
            }

            k = face;
            j = ny - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i, j - 3, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i, j - 3, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i, j - 3, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i, j + 1, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i, j - 1, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i, j - 3, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i, j + 1, a)] + (3.0 / 4.0) * fneqf_u[updown(i, j - 1, a)] - (1.0 / 8.0) * fneqf_u[updown(i, j - 3, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];
            }
        }

        for (j = 1; j <= ny - 1; j += alpha)
        {
            k = face;
            i = 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i + 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i + 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i + 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i - 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i + 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i + 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i - 1, j, a)] + (3.0 / 4.0) * fneqf_u[updown(i + 1, j, a)] - (1.0 / 8.0) * fneqf_u[updown(i + 3, j, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];
            }

            k = face;
            i = nx - 1;

            jxf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jxf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jxf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jxf[fineisnindex1(i - 3, j, k)];
            jyf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jyf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jyf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jyf[fineisnindex1(i - 3, j, k)];
            jzf[fineisnindex1(i, j, k)] = (3.0 / 8.0) * jzf[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * jzf[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * jzf[fineisnindex1(i - 3, j, k)];
            rhof[fineisnindex1(i, j, k)] = (3.0 / 8.0) * rhof[fineisnindex1(i + 1, j, k)] + (3.0 / 4.0) * rhof[fineisnindex1(i - 1, j, k)] - (1.0 / 8.0) * rhof[fineisnindex1(i - 3, j, k)];

            u2 = jxf[fineisnindex1(i, j, k)] * jxf[fineisnindex1(i, j, k)] + jyf[fineisnindex1(i, j, k)] * jyf[fineisnindex1(i, j, k)] + jzf[fineisnindex1(i, j, k)] * jzf[fineisnindex1(i, j, k)];

            for (a = 0; a < 19; a++)
            {
                fneqf_u[updown(i, j, a)] = (3.0 / 8.0) * fneqf_u[updown(i + 1, j, a)] + (3.0 / 4.0) * fneqf_u[updown(i - 1, j, a)] - (1.0 / 8.0) * fneqf_u[updown(i - 3, j, a)];

                t = ex[a] * jxf[fineisnindex1(i, j, k)] + ey[a] * jyf[fineisnindex1(i, j, k)] + ez[a] * jzf[fineisnindex1(i, j, k)];
                feqf[fineindex(i, j, k, a)] = w[a] * rhof[fineisnindex1(i, j, k)] * (1 + 3 * t + 4.5 * t * t - 1.5 * u2);
                ff[fineindex(i, j, k, a)] = feqf[fineindex(i, j, k, a)] + fneqf_u[updown(i, j, a)];
            }
        }
    }
}
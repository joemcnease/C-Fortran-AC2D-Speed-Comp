#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>


void fd2d(float *p, const int nx, const int nz, const float dx, const float dz,
          const int nt, const float dt, const float *stf, const int sx, const int sz,
          const float *c, bool save);


void print_array(const float *array, const int size);


void save_pressure(const float *p, const int nx, const int nz, const char *filename);


void fill(float *array, const int size, const float val);


void fd2d(float *p, const int nx, const int nz, const float dx, const float dz,
          const int nt, const float dt, const float *stf, const int sx, const int sz,
          const float *c, bool save) {

    float pOld[nx*nz];
    float pNew[nx*nz];
    float d2px[nx*nz];
    float d2pz[nx*nz];

    fill(pOld, nx*nz, 0.);
    fill(pNew, nx*nz, 0.);
    fill(d2px, nx*nz, 0.);
    fill(d2pz, nx*nz, 0.);

    for (int it=0; it<nt-1; it++) {
        for (int i=0; i<nz-1; i++) {
            for (int j=1; j<nx-2; j++) {
                int idx = i*nx + j;
                d2px[idx] = (p[idx-1] - 2*p[idx] + p[idx+1]) / (dx*dx);
            }
        }
        for (int j=0; j<nx-1; j++) {
            for (int i=1; i<nz-2; i++) {
                int idx = i*nx + j;
                d2pz[idx] = (p[idx-nx] - 2*p[idx] + p[idx+nx]) / (dz*dz);
            }
        }
        for (int i=0; i<(nx*nz)-1; i++) {
            pNew[i] = 2.*p[i] - pOld[i] + (dt*dt) * (c[i]*c[i]) * (d2px[i] + d2pz[i]);
        }

        pNew[sz*nx + sx] = pNew[sz*nx + sx] + stf[it];
        for (int i=0; i<nz-1; i++) {
            for (int j=0; j<nx-1; j++) {
                int idx = i*nx + j;
                pOld[idx] = p[idx];
                p[idx] = pNew[idx];
            }
        }

        if (save && (it%10==0)) {
            char fn[50];
            sprintf(fn, "pressure/pressure%08d", it);
            save_pressure(p, nx, nz, fn);
        }

        printf("Time step: %d", it);
        printf("\n");
    }
}


void print_array(const float *array, const int size) {
    for (int i=0; i<size; i++) {
        printf("array[%i] = %10.5f \n", i, array[i]);
    }
}


void save_pressure(const float *p, const int nx, const int nz, const char *filename) {
    FILE *fptr;
    fptr = fopen(filename, "w");

    for (int i=0; i<nz-1; ++i) {
        for (int j=0; j<nx-1; ++j) {
            fprintf(fptr, "%6.3f,", p[i*nx + j]);
        }
        fputs("\n", fptr);
    }
    fclose(fptr);
}


void fill(float *array, const int size, const float val) {
    for (int i=0; i<size-1; i++) {
        array[i] = val;
    }
}


int main() {

#define nx 301
#define nz 201
#define dx 5.0
#define dz 5.0
#define sx floor(nx/2)
#define sz floor(nz/2)
#define nt 1000
#define dt 0.0005
#define t0 0.01
#define f0 100.0
#define c0 3000.0

    int res = mkdir("pressure", 0777);

    float c[nx*nz];
    float p[nx*nz];
    float stf[nt];

    fill(c, nx*nz, c0);
    fill(p, nx*nz, 0.);
    fill(stf, nt, 0.);

    for (int i=0; i<nt-1; i++) {
        stf[i] = exp(-(f0*f0)*(i*dt-t0)*(i*dt-t0));
    }

    clock_t start, end;
    double elapsed;

    start = clock();
    fd2d(p, nx, nz, dx, dz, nt, dt, stf, sx, sz, c, true);
    end = clock();

    elapsed = (double)(end-start) / CLOCKS_PER_SEC;
    printf("Time to compute: %f \n", elapsed);

    return 0;
}

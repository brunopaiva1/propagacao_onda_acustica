#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void initializeSource(float *s, float f, float dt, int nt) {
    float t;
    float pi = 3.14;

    for (int i = 0; i < nt; i++){
        t = i * dt;
        s[i] = (1 - 2 * pi * pi * f * f * t * t) * exp(-pi * pi * f * f * t * t);
    }
    
}

void propagateWave(float *s, float c, float dx, float dy, float dz, float dt,
                    int nx, int ny, int nz, int nt, int xs, int ys, int zs) {
    
    float dEx, dEy, dEz;
    float *uAnterior = (float*) malloc(nx * ny * nz * sizeof(float));
    float *uProximo = (float*) malloc(nx * ny * nz * sizeof(float));
    float *u = (float*) malloc(nx * ny * nz * sizeof(float));

    memset(u, 0, nx * ny * nz * sizeof(float));
    memset(uAnterior, 0, nx * ny * nz * sizeof(float));
    memset(uProximo, 0, nx * ny * nz * sizeof(float));

    for (int t = 0; t < nt; t++) {
        for (int idx = 0; idx < (nx - 4) * (ny - 4) * (nz - 4); idx++) {
            int x = 2 + idx / ((ny - 4) * (nz - 4));
            int y = 2 + (idx / (nz - 4)) % (ny - 4);
            int z = 2 + idx % (nz - 4);

            dEx = ((-1.0/12.0)*uAnterior[(x - 2) * ny * nz + y * nz + z] +
                    (4.0/3.0)*uAnterior[(x - 1) * ny * nz + y * nz + z] -
                    (5.0/2.0)*uAnterior[x * ny * nz + y * nz + z] +
                    (4.0/3.0)*uAnterior[(x + 1) * ny * nz + y * nz + z] -
                    (1.0/12.0)*uAnterior[(x + 2) * ny * nz + y * nz + z]) / (dx * dx);
            
            dEy = ((-1.0/12.0)*uAnterior[x * ny * nz + (y - 2) * nz + z] +
                    (4.0/3.0)*uAnterior[x * ny * nz + (y - 1) * nz + z] -
                    (5.0/2.0)*uAnterior[x * ny * nz + y * nz + z] +
                    (4.0/3.0)*uAnterior[x * ny * nz + (y + 1) * nz + z] -
                    (1.0/12.0)*uAnterior[x * ny * nz + (y + 2) * nz + z]) / (dy * dy);

            dEz = ((-1.0/12.0)*uAnterior[x * ny * nz + y * nz + (z - 2)] +
                    (4.0/3.0)*uAnterior[x * ny * nz + y * nz + (z - 1)] -
                    (5.0/2.0)*uAnterior[x * ny * nz + y * nz + z] +
                    (4.0/3.0)*uAnterior[x * ny * nz + y * nz + (z + 1)] -
                    (1.0/12.0)*uAnterior[x * ny * nz + y * nz + (z + 2)]) / (dz * dz);

            uProximo[x * ny * nz + y * nz + z] = c * c * dt * dt * (dEx + dEy + dEz) - uAnterior[x * ny * nz + y * nz + z] + 2 * u[x * ny * nz + y * nz + z];
                    
        }


        uProximo[xs * ny * nz + ys * nz + zs] -= c * c * dt * dt * s[t];

        float *temp = u;
        u = uProximo;
        uProximo = uAnterior;
        uAnterior = temp;

        // if (t % 50 == 0)
        // {
        //     char filename[50];
        //     sprintf(filename, "samples/sample_t%d.bin", t); // Cria um nome de arquivo único para cada tempo
        //     FILE *file = fopen(filename, "wb");
        //     if (file != NULL) {
        //         // Escreva os dados de uProximo no arquivo binário
        //         fwrite(uProximo, sizeof(float), nx * ny * nz, file);
        //         fclose(file);
        //     } else {
        //         printf("Erro ao abrir o arquivo para escrita.\n");
        //     }
        // }
        

    }
    
    free(uAnterior);
    free(uProximo);
                    }

int main(int argc, char* argv[]) {
    int xs = 15, ys = 15, zs = 15;
    float dx = 10, dy = 10, dz = 10;
    float dt = 0.001;
    int nx = 50, ny = 50, nz = 50;
    int nt = 10000;
    float f = 10;
    float c = 1500.0;
    int;

    float *s = (float *)malloc(nt * sizeof(float));

    initializeSource(s, f, dt, nt);
    propagateWave(s, c, dx, dy, dz, dt, nx, ny, nz, nt, xs, ys, zs);

    free(s);

    return 0;
}
/*
 * wav2.c
 *
 * Copyright 2015 M. Uemoto <uemoto@mp.es.osaka-u.ac.jp>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.141592653589793
#define FACT1 3.81 /* hbar^2/2m */
#define NGRID 64
#define EMAX 2.5

typedef union {
  struct {
    double x;
    double y;
    double z;
  };
  double v[3];
} vec3d;

typedef union {
  struct {
    vec3d a;
    vec3d b;
    vec3d c;
  };
  vec3d u[3];
  double v[9];
} lattice;

typedef struct {
  double n_record;
  double n_spin;
  double n_prec;
} wavecar_header;

typedef struct {
  double n_kpoint;
  double n_band;
  double e_cutoff;
  lattice a;
  double e_fermi;
} wavecar_info;

typedef struct {
  double n_coeff;
  vec3d k;
} kpoint_info;

typedef struct {
  double energy;
  double dummy;  // The meaning of this part is not sure.... (:P
  double occup;
} state_info;

static int etbl[3][2] = {{1, 2}, {2, 0}, {0, 1}};

int ireduce(int i, int cycle) {
  int half = cycle / 2;
  return (i + half) % cycle - half;
}

void scalar3d(double a, vec3d *p, vec3d *q) {
  for (int i = 0; i < 3; i++) {
    q->v[i] = a * p->v[i];
  }
}

void add3d(vec3d *p, vec3d *q, vec3d *r) {
  for (int i = 0; i < 3; i++) {
    r->v[i] = p->v[i] + q->v[i];
  }
}

double inner3d(vec3d *p, vec3d *q) {
  double r = 0;
  for (int i = 0; i < 3; i++) {
    r += p->v[i] * q->v[i];
  }
  return r;
}

void outer3d(vec3d *a, vec3d *b, vec3d *c) {
  for (int i = 0; i < 3; i++) {
    int j = etbl[i][0], k = etbl[i][1];
    c->v[i] = (a->v[j] * b->v[k]) - (a->v[k] * b->v[j]);
  }
}

void convert3d(vec3d *p, lattice *a, vec3d *q) {
  for (int i = 0; i < 3; i++) {
    q->v[i] = 0.0;
    for (int j = 0; j < 3; j++) {
      q->v[i] += a->u[j].v[i] * p->v[j];
    }
  }
}

double recipl(lattice *a, lattice *b) {
  vec3d aa[3];
  for (int i = 0; i < 3; i++) {
    int j = etbl[i][0], k = etbl[i][1];
    outer3d(&a->u[i], &a->u[j], &aa[k]);
  }
  double v = inner3d(&a->u[0], &aa[0]);
  double f = 2.0 * PI / v;
  for (int i = 0; i < 3; i++) {
    scalar3d(f, &aa[i], &b->u[i]);
  }
  return v;
}

double _Complex prod(float _Complex *a, float _Complex *b, int n) {
  double _Complex s = 0;
  for (int i = 0; i < n; i++) {
    s += conj(a[i]) * b[i];
  }
  return s;
}

int create_gtable(lattice *b, double e_cut, int n_plain, vec3d *gtbl) {
  vec3d gi;
  int n_count = 0;
  for (int i0 = 0; i0 < NGRID; i0++) {
    gi.v[0] = ireduce(i0, NGRID);
    for (int i1 = 0; i1 < NGRID; i1++) {
      gi.v[1] = ireduce(i1, NGRID);
      for (int i2 = 0; i2 < NGRID; i2++) {
        gi.v[2] = ireduce(i2, NGRID);
        // Calculate plane wave energy
        vec3d gk;
        convert3d(&gi, b, &gk);
        double e = FACT1 * inner3d(&gk, &gk);
        // g must be inside the cutoff sphere:
        if (e < e_cut) {
          if (n_count < n_plain) {
            gtbl[n_count] = gk;
            n_count++;
          } else {
            printf("# Broken NGRID is too small\n");
            return -1;
          }
        }
      }
    }
  }
  return n_count;
}

double read_wave(FILE *fh, int n_record, int n_coeff, int i,
                 float _Complex *wf) {
  fseek(fh, (3 + i) * n_record, SEEK_SET);
  fread((void *)wf, sizeof(float _Complex), n_coeff, fh);
  return 1.0 / sqrt(creal(prod(wf, wf, n_coeff)));  // Normalize const.
}

// Transition probablity by linear polarized light
void matrix1(int n_plain, vec3d *gtbl, float _Complex *wf1, float _Complex *wf2,
             double _Complex *mx, double _Complex *my) {
  // Initialize
  *mx = 0;
  *my = 0;
  // Summension over all diagonal element:
  for (int i = 0; i < n_plain; i++) {
    vec3d g = gtbl[i];
    double _Complex wf1u = wf1[i];
    double _Complex wf2u = wf2[i];
    double _Complex wf1d = wf1[i + n_plain];
    double _Complex wf2d = wf2[i + n_plain];
    // First term
    *mx += conj(wf1u) * g.x * wf2u;
    *mx += conj(wf1d) * g.x * wf2d;
    *my += conj(wf1u) * g.y * wf2u;
    *my += conj(wf1d) * g.y * wf2d;
  }
}

int main(int argc, char **argv) {
  FILE *fh = fopen("WAVECAR", "r");

  // Read Wavecar "HEADER" section
  fseek(fh, 0, SEEK_SET);
  wavecar_header header;
  fread(&header, sizeof(wavecar_header), 1, fh);
  int n_record = (int)header.n_record;
  // Error check
  if ((int)header.n_spin != 1) {
    printf("# Error: not spinor!\n");
    return -1;
  }
  if ((int)header.n_prec != 45200) {
    printf("# Error: not double precision!\n");
    return -1;
  }

  // Read Wavecar "INFOMATION" section
  fseek(fh, n_record, SEEK_SET);
  wavecar_info info;
  fread(&info, sizeof(wavecar_info), 1, fh);
  int n_band = info.n_band;
  // Calculate reciplocal lattice vector: b
  lattice b;
  recipl(&info.a, &b);

  // Read "GAMMA POINT" section
  fseek(fh, n_record * 2, SEEK_SET);
  kpoint_info kpoint;
  state_info *state = malloc(n_band * sizeof(state_info));
  fread(&kpoint, sizeof(kpoint_info), 1, fh);    // infomation
  fread(state, sizeof(state_info), n_band, fh);  // eigenenergies
  int n_coeff = (int)kpoint.n_coeff;
  int n_plain = n_coeff / 2;  // 2 is spin degeneracy

  // Create table info
  vec3d *gtbl = malloc(n_plain * sizeof(vec3d));
  int n_count = create_gtable(&b, info.e_cutoff, n_plain, gtbl);
  if (n_count != n_plain) {
    printf("# BROKEN WAVECAR!\n");
    return -1;
  }

  // Write out all band energies
  printf("# eigen energies on k=0\n");
  for (int i = 0; i < n_band; i++) {
    printf("eigenstate i=%d e=%f n=%f\n", i, state[i].energy, state[i].occup);
  }

  // Write out all considerable transittion
  printf("# transition at k=0\n");
  // Wavefunction: wfi, wfj
  float _Complex *wfi = malloc(sizeof(float _Complex) * n_coeff);
  float _Complex *wfj = malloc(sizeof(float _Complex) * n_coeff);

  for (int i = 0; i < n_band - 1; i++) {
    if (state[i].energy < info.e_fermi) {
      // Read wavefunction: i
      double ci = read_wave(fh, n_record, n_coeff, i, wfi);

      for (int j = i + 1; j < n_band; j++) {
        if (info.e_fermi < state[j].energy) {
          double ediff = state[j].energy - state[i].energy;
          double ndiff = state[i].occup - state[j].occup;
          // Upper limit of e-MAX energy
          if (ediff < EMAX) {
            // Read wavefunction: j
            double cj = read_wave(fh, n_record, n_coeff, j, wfj);

            // Calculate
            double _Complex mx, my;
            matrix1(n_plain, gtbl, wfi, wfj, &mx, &my);
            mx = mx * ci * cj;  // Normalize
            my = my * ci * cj;  // Normalize

            // Output Transition
            printf("i=%d j=%d e=%f n=%f x2=%e y2=%e \n", i, j, ediff, ndiff,
                   cabs(mx) * cabs(mx), cabs(my) * cabs(my));
          }
        }
      }
    }
  }
  return 0;
}

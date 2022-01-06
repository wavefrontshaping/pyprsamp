/*
 * Copyright 2016 Boshra Rajaei <b.rajaei@sadjad.ac.ir>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "prSAMP_Calibration.h"
//#include <time.h>


/* (Sequential) AMP for sparse matrices */

double cgamp (
        int n, int m, int strIndY, int strInda, double *y, int *ir, int *jc,
        double delta, double *p_prms,
        int t_max, double eps, double damp, double vnf, int disp, complex double *a
    ) {
    complex double *a_proj, *o, *g, *g_old, *r;
    complex double a_old, g_proj;
    double *v, *dg, *sig,*c_proj;
    //complex double *a_init;
    double c_old, dg_proj, v_old;
    double diff;//, diff_y=0;
    double *c;


    //time_t t1, t2;

    int i, mu, idx, t;
    int *seq, key;

    /* Alloc. structures */

    o = malloc(m*sizeof(*o));
    g = malloc(m*sizeof(*g));
    g_old = malloc(m*sizeof(*g_old));
    a_proj = malloc(m*sizeof(*a_proj));
    c_proj = malloc(m*sizeof(*c_proj));
    v = malloc(m*sizeof(*v));
    dg = malloc(m*sizeof(*dg));
    seq = malloc(n*sizeof(*seq));
    r = malloc(n*sizeof(*r));
    sig = malloc(n*sizeof(*sig));
    //a_init = malloc(n*sizeof(*a_init));
    c = malloc(n*sizeof(*c));

    for (i=0;i<n;i++)
        c[i] = 0.5*sqrt(0.1);


    //for(i=0;i<n;i++)
    //    a_init[i] = a[i];


    srand(time(NULL));

    //if (!a_proj || !c_proj || !o || !v || !g || !g_old || !dg || !seq)
    //    printf("Failure in allocating memory.");

    for (mu = 0; mu < m; mu++) g[mu] = 0;

    diff = 0.;
    //t1 = time(NULL);
    for (t = 0; t < t_max; t++) {
        // Generate random permutation
        for (key = 0; key < n; key++) seq[key] = key;
        sort_rand(n, seq);

        // Update a_proj and c_proj
        for (mu = 0; mu < m; mu++)
            a_proj[mu] = c_proj[mu] = 0;
        for (i = 0; i < n; i++) {
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += a[strInda+i];
                c_proj[ ir[idx] ] += c[i];
                //v[ ir[idx] ] += c[i];
            }
        }

        // Update w and v
        for (mu = 0; mu < m; mu++) {
            o[mu] = a_proj[mu] - c_proj[mu] * g[mu];
            v[mu] = c_proj[mu];
        }

        for (mu = 0; mu < m; mu++)
            channel(y[strIndY+mu], o[mu], v[mu], delta, &g[mu], &dg[mu], vnf);

        for (mu = 0; mu < m; mu++) g_old[mu] = g[mu];


        // Sweep over all n variables, in random order
        diff = 0.;

        for (key = 0; key < n; key++) {
            i = seq[key];

            // Update {sig, r}, {a, c}
            a_old = a[strInda+i], c_old = c[i];

            g_proj = dg_proj = 0.; // Dot products: F * g and -F^2 * dg
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                g_proj += g[ ir[idx] ];
                dg_proj -= dg[ ir[idx] ];
            }
            if(t>0)
                sig[i] = damp * sig[i] + (1 - damp) * (1. / dg_proj);
            else
                sig[i] = 1. / dg_proj;

            if(sig[i] < 0) sig[i] = 0.1*vnf;

            if(t>0)
                r[i] = damp * r[i] + (1 - damp) * (a[strInda+i] + sig[i] * g_proj);
            else
                r[i] = a[strInda+i] + sig[i] * g_proj;

            prior(r[i], sig[i], p_prms, &a[strInda+i], &c[i], vnf);


            // Update {w, v}, {g, dg}
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                mu = ir[idx];
                v_old = v[mu];

                v[mu] += (c[i] - c_old);
                o[mu] += (a[strInda+i] - a_old)
                        - (v[mu] - v_old) * g_old[mu];

                channel(y[strIndY+mu], o[mu], v[mu], delta, &g[mu], &dg[mu], vnf);
            }

            diff += cabs(a[strInda+i] - a_old);

        }
        diff /= n;
/*
        for (mu = 0; mu < m; mu++)
            a_proj[mu] = 0;
        for (i = 0; i < n; i++) {
            for (idx = jc[i]; idx < jc[i + 1]; idx++) {
                a_proj[ ir[idx] ] += a[i];
            }
        }
        diff_y = 0;
        for (mu = 0; mu < m; mu++)
            diff_y += abs(y[strIndY+mu] - cabs(a_proj[mu]));
        diff_y /= m;*/

        // Print some info.
        if (disp)
            PySys_WriteStdout("t: %3d; diff: %.4e\n", t, diff);
           // printf("t: %3d; diff: %.4e; diff_y: %.4e\n", t, diff, diff_y);
        //mexEvalString("drawnow");

        // Check for convergence
        if (diff < eps || isnan(diff)) break;

        //if(diff>1 && t>100)
        //  break;

    }

    if(isnan(diff))//diff>1e-1 ||
      printf("fail\n");


    //t2 = time(NULL);
    //printf("execution time: %d\n", t2-t1);

    free(seq);
    free(dg);
    free(g_old);
    free(g);
    free(v);
    free(o);
    free(c_proj);
    free(a_proj);
    free(r);
    free(sig);
    free(c);

    return diff;

}

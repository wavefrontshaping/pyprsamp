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

/* Random order */
void sort_rand( int n, int *seq ) {
    int i, key, exchg;

    for (key = 0; key < n - 1; key++) {
        exchg = key + rand_int(n - key);
        i = seq[key]; seq[key] = seq[exchg]; seq[exchg] = i;
    }
}

double drand(double dMin, double dMax)
{
    return dMin + ( dMax - dMin) * (rand() % RAND_MAX) / RAND_MAX;
}

double gaussrand()
{
//#pragma omp critical
//{
	static double V1, V2, S;
	static int phase = 0;
	double X, U1, U2;

	if(phase == 0) {
		do {
			U1 = (double)rand() / RAND_MAX;
			U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
//}
	return X;
}

double logsig(double x)
{
    return 1/(1+exp(-x));
}

//b=A*x  or b=A'*x
//depending on m, n (size of A) and p (size of x)
void project(double *A, double *x, double *b, int m, int n, int p)
{
    int i,j,ind;

   // for(i=0;i<m*n;i++)
   //     printf("index: %d, %.4e\n", i, A[i]);

    if(p==n)
    {
        for (i = 0; i < m; i++) b[i] = 0;
        for (i = 0; i < m; i++)//row
            for (j = 0; j < n; j++) {//column
                ind = i + j * m;
                //printf("index: %d, %.4e\n", ind, A[ind]);
                b[i] += A[ind] * x[j];
            }
    }
    else
    {
        for (i = 0; i < n; i++) b[i] = 0;
        for (i = 0; i < n; i++)//column
            for (j = 0; j < m; j++) {//row
                ind = j + i * m;
                b[i] += A[ind] * x[j];
            }
    }

}

void readFromText(FILE *fp, int n, int nnz, int *ir, int *jc, double *p_prms, int *t_max, double *prec, double *damp, double *vnf, double *delta, int *disp)
{
    int i;

    fscanf(fp, "%lf", delta);
    fscanf(fp, "%lf", &p_prms[0]);
    fscanf(fp, "%lf", &p_prms[1]);
    fscanf(fp, "%lf", &p_prms[2]);
    fscanf(fp, "%d", t_max);
    fscanf(fp, "%lf", prec);
    fscanf(fp, "%d", disp);
    fscanf(fp, "%lf", damp);
    fscanf(fp, "%lf", vnf);

    //printf("delta: %lf\n rho=%lf\n mean=%lf\n var=%lf\n iter=%d\n prec:%lf\n disp:%d\n damp:%lf\n vnf:%lf\n",*delta,p_prms[0],p_prms[1],p_prms[2], *t_max, *prec, *disp, *damp,*vnf);


    for(i=0;i<n+1;i++)
        fscanf(fp, "%d", &jc[i]);

    //printf("%d\n%d\n%d\n",jc[n-3],jc[n-2],jc[n-1]);

    for(i=0;i<nnz;i++)
        fscanf(fp, "%d", &ir[i]);

    //printf("%d\n%d\n%d\n",ir[nnz-3],ir[nnz-2],ir[nnz-1]);
}

void readFromTextY(FILE *fp, int m, int p, double *y)
{
    int i;

    for(i=0;i<m*p;i++)
        fscanf(fp, "%lf", &y[i]);
}

void readDimensions(FILE *fp, int *m, int *n, int *p, int *nnz)
{
    fscanf(fp, "%d", m);
    fscanf(fp, "%d", n);
    fscanf(fp, "%d", p);
    fscanf(fp, "%d", nnz);

    //printf("M=%d\n N=%d\n P=%d\n NNZ=%d\n",*m,*n,*p,*nnz);

}

void writeRow(FILE *fp, int n, complex double *a)
{
    int i;
    for(i=0;i<n;i++)
        if(isnan((float)a[i]))
            fprintf(fp, "0\n0\n");
        else
            fprintf(fp, "%0.16f\n%0.16f\n", creal(a[i]), cimag(a[i]));
}



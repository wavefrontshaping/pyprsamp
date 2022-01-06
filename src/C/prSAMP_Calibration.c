/*
 * Copyright 2016 Boshra Rajaei <b.rajaei@sadjad.ac.ir>
 * Modified for Python by Sebastien M. Popoff
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
#include <complex.h>

static PyObject *PrsampError;

static PyMethodDef Methods[] = {
	{"calibration",  (PyCFunction)calibration, METH_VARARGS, "Run the prSAMP calibration alogrithm."},
     {"calibrationFromFile",  (PyCFunction)calibrationFromFile, METH_VARARGS, "Run the prSAMP calibration alogrithm with data sent with a file."},
    	{NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef prsamp =
{
    PyModuleDef_HEAD_INIT,
    "prsamp", /* name of module */
    "usage: \n", /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    Methods
};

PyMODINIT_FUNC PyInit_prsamp(void)
{
    PyObject *m;
    m = PyModule_Create(&prsamp);
    if (m == NULL)
        return;
    PrsampError = PyErr_NewException("prsamp.error", NULL, NULL);
    Py_INCREF(PrsampError);
    PyModule_AddObject(m, "error", PrsampError);
    return m;
}

PyObject* calibration(PyObject* self, PyObject* args)
{
      int i, j;
      double p_prms[3], eps, damp, vnf, delta;
      double *y;
      int *ir, *jc;
      int m, n, p, nnz; //nnz: number of non zero cells in X
      int t_max;
    
    int disp = 0;
     
    PyObject *py_jc, *py_ir, *py_y;
    PyObject *py_H = PyList_New(0);
    PyObject *py_diff = PyList_New(0);
    PyObject *py_ret = PyList_New(0);

    if (!PyArg_ParseTuple(args, "O!O!O!d(ddd)iddd|i", 
        &PyList_Type, &py_jc, &PyList_Type, &py_ir,&PyList_Type, &py_y,
        &delta, &p_prms[0],&p_prms[1],&p_prms[2],&t_max, &eps, &damp, &vnf, &disp))
        {
            PyErr_SetString(PrsampError, "Invalid arguments.");
            return 0;
        }

 
    n = PyList_Size(py_jc)-1;
    nnz = PyList_Size(py_ir);
    m = PyList_Size(py_y);
    p = PyList_Size(PyList_GetItem(py_y,0));

    ir = malloc(nnz*sizeof(*ir));
    jc = malloc((n+1)*sizeof(*jc));
    y = malloc(p*m*sizeof(*y));

  
    // Get data from Python lits
    for (i = 0; i < m; i++)
    {
        PyObject *y_line = PyList_GetItem(py_y,i);
        for (j=0; j < p; j++)
            y[p*i+j] = PyFloat_AsDouble(PyList_GetItem(y_line,j));
    }
    for (i = 0; i < nnz; i++) 
        ir[i] = (int)PyLong_AsLong(PyList_GetItem(py_ir,i));
    for (i = 0; i < n+1; i++) 
        jc[i] = (int)PyLong_AsLong(PyList_GetItem(py_jc,i));

    
    
    srand(time(NULL));
    gaussrand();
    complex double *h = malloc(m*n*sizeof(*h));
    double *diff = malloc(m*sizeof(*diff));

    for(j=0;j<m*n;j++)
        h[j] = gaussrand()*sqrt(1.0/n)+gaussrand()*sqrt(1.0/n)*I;
    
    
    for(i = 0; i < m; i++) 
    {
        diff[i] = cgamp(n, p, i*p, i*n, y, ir, jc, delta, p_prms,
        t_max, eps, damp, vnf, disp, h);
    }
     
    
    

    for (i = 0; i < m; i++)
        for(j = 0; j < n; j++)
            PyList_Append(py_H,PyComplex_FromDoubles(crealf(h[i*n+j]),cimagf(h[i*n+j])));

    for (i = 0; i < m; i++)
         PyList_Append(py_diff,PyFloat_FromDouble(diff[i]));

    PyList_Append(py_ret, py_H);
    PyList_Append(py_ret, py_diff);

    return py_ret;


}

PyObject* calibrationFromFile(PyObject* self, PyObject* args)
{
      int i, j;
      double *p_prms, eps;
      double damp, vnf;
      double delta, *y;
      int *ir, *jc;
      int m, n, p, nnz; //nnz: number of non zero cells in X
      int disp, t_max;
      char * part;
     
    

      // Get the parametrs
    if (!PyArg_ParseTuple(args, "s", &part)) 
        {
            PyErr_SetString(PrsampError, "Invalid arguments!");
            return 0;
        }
     
      
        // CALIBRATION CODE START HERE
  
    FILE *fp,*fp_w;
    char name[30] = "input";
    strcat(name, part);
    strcat(name, ".txt");
    fp = fopen(name, "r");
    PySys_WriteStdout("Reading file: %s.\n",name);
    
    readDimensions(fp, &m, &n, &p, &nnz);

    
    p_prms = malloc(3 * sizeof *p_prms);
    ir = malloc(nnz*sizeof(*ir));
    jc = malloc((n+1)*sizeof(*jc));
    y = malloc(p*m*sizeof(*y));
      
    readFromText(fp, n, nnz, ir, jc, p_prms, &t_max, &eps, &damp, &vnf, &delta, &disp);
    readFromTextY(fp, m, p, y);
    fclose(fp);
    
    
    
    srand(time(NULL));
    gaussrand();
    complex double *h = malloc(m*n*sizeof(*h));
    double *diff = malloc(m*sizeof(*diff));

    for(j=0;j<m*n;j++)
        h[j] = gaussrand()*sqrt(1.0/n)+gaussrand()*sqrt(1.0/n)*I;
    
  
    for(i = 0; i < m; i++) {
        diff[i] = cgamp(n, p, i*p, i*n, y, ir, jc, delta, p_prms,
        t_max, eps, damp, vnf, disp, h);
     }
    

    char nameo[30] = "output";
    strcat(nameo, part);
    strcat(nameo, ".txt");
    fp_w = fopen(nameo, "w");
    writeRow(fp_w, m*n, h);
    //Write the diff values   
    for(i=0;i<n;i++)
        if(isnan((float)diff[i]))
            fprintf(fp_w, "-1.0\n");
        else
            fprintf(fp_w, "%0.16f\n", diff[i]);
    fclose(fp_w);
    PySys_WriteStdout("Writing file: %s.\n",nameo);

  Py_RETURN_NONE;
}



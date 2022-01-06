# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 17:43:13 2016

@author: Sebastien
"""

import prsamp
import numpy as np
import os
import math



__all__=['save_input_file','load_output_file','delete_files','Run','RunFromFile','prsampPrms']

#from joblib import Parallel, delayed
#import multiprocessing
class prsampPrms():
    def __init__(self,shape):
        self.m = shape[0]
        self.n = shape[1]
        self.delta = 1e-1 # noise variation parameter
        self.priorPrms = (1., 0., 1./np.sqrt(self.n)) # gaussian prior parameters 
                           # 1. sparsity, 0<rho<=1 (maximum number of nonzero entries in the
                           # resulting vector)
                           # 2. mean, mu
                           # 3. variance, v
        self.maxIter = 500*self.n/2; # maximum number of main AMP loop iteration
        self.prec = 1e-5; # convergence criterion. minimum difference between two consequent x
        self.display = 0; #display the temporary convergance results or not
        self.damp = 0.9; # how fast the algorithm converges to a local minima 0<=damp<1
        self.vnf = 0.5; # variance normalization factor
        
        
        
def save_input_file(X, Y, prms = None):
    '''' save2file_calib
    Generates input<blk_code>.txt file necessary to call ./prSAMP_Calibration in Unix or
    prSAMP_Calib.exe in Windows.
    '''

    # Generate a random id
    prblmid = int(np.round(np.random.rand()*1000))
    [n, p]  = X.shape
    m = Y.shape[0]

    if prms == None:
        prms = prsampPrms([m,n])
    
    nnz = np.sum(X) # Number of non-zero values in X
    jc = np.sum(X.transpose(),axis = 0)
    for i in range(1,len(jc)):
        jc[i] = jc[i]+jc[i-1]
    
    jc = np.insert(jc,0,0)
    jc = jc.tolist()
    
    ir = [i for j in range(n) for i in range(p) if X[j,i]]

    with open(r''.join(['input',str(prblmid),'.txt']),'w') as fileID:

        fileID.write('%.0f\n' % m)
        fileID.write('%.0f\n'%n)
        fileID.write('%.0f\n'%p)
        fileID.write('%.0f\n'%nnz)
        fileID.write('%f\n'%prms.delta)
    #    fileID.write('%f\n'%opt_priorPrms)
        for value in  prms.priorPrms:
            fileID.write('%f\n'%value)
    #    fileID.write('\n')
        fileID.write('%.0f\n'%prms.maxIter)
        fileID.write('%f\n'%prms.prec)
        fileID.write('%.0f\n'%prms.display)
        fileID.write('%.2f\n'%prms.damp)
        fileID.write('%.2f\n'%prms.vnf)
    #    fileID.write('%.0f\n'%jc)
        for value in  jc:
            fileID.write('%.0f\n'%value)
    #    fileID.write('%.0f\n'%ir)
        for value in ir:
            fileID.write('%.0f\n'%value)
    #    fileID.write('%.16f\n'%Y)
        for vec in  Y:
            for value in vec:
                fileID.write('%.16f\n'%value)   
                    
        fileID.close()
    return prblmid 



def load_output_file(filesuffix, m, n):
    '''
    loadfile_calib
    Loads output<filesuffix>.txt file which is the output of ./prSAMP_Calibration in Unix or
    prSAMP_Calib.exe in Windows.
    '''
    with open(r''.join(['output',str(filesuffix),'.txt']),'r') as f:
        H_hat = np.zeros([m,n],dtype = 'complex')
        for i in range(m):
            for j in range(n):
                H_hat[i,j] = float(f.next().replace(',','.'))+ complex(0,1)*float(f.next().replace(',','.'))

        f.close()       
    return H_hat

    
def delete_files(filesuffix):
#    pass
    os.remove(''.join(['input',str(filesuffix),'.txt']))
    os.remove(''.join(['output',str(filesuffix),'.txt']))

    

def RunFromFile(X,Y, prms = None, delete=True):
    
    n = X.shape[0]
    m = Y.shape[0]
    if prms == None:
        prms = prsampPrms([m,n])
        
    calibid = save_input_file(X, Y, prms)
    prsamp.calibrationFromFile(str(calibid))
    ret = load_output_file(calibid,m,n)
    
    if delete == True:
        delete_files(calibid)
        
    return ret
    
    
def RunUntilConv(X,Y,  prms, repeatMax = 10, maxIterIncrease = 1.):  
    n,p = X.shape
    m = Y.shape[0]
    flag = True
    counter = 0
    H_ret = np.zeros([m,n],dtype='complex')
    eps_ret = [1e8]*m
    indVec = np.arange(m)
    while flag:

        print('Proceeding iteration %d:\t %d rows left.' % (counter+1,len(indVec)))
        H_temp,eps_temp = Run(X,Y[indVec],prms)

               
#        indCopy = [ind for i,ind in enumerate(indVec) if eps_temp[i] < eps_ret[ind]]

        for i,ind in enumerate(indVec): 
            # We keep the values only if the results are better than the previous iterations 
            if( eps_temp[i] < eps_ret[ind] ):
                eps_ret[ind] = eps_temp[i]
                H_ret[ind] = H_temp[i]

#        H_ret[indCopy] = H_temp[:]
        counter += 1
        
        indVec = [x for x in indVec if (math.isnan(eps_ret[x]) or eps_ret[x] > prms.prec)]
        if not indVec:
            flag = False
            print("Convergence criteria met.")

        if not counter < repeatMax:
            flag = False
            print("Convergence criteria not met.")
            

            
        prms.maxIter = int(maxIterIncrease*prms.maxIter)
        
#        print(H_temp.shape)

            
        
    return H_ret,eps_ret
    
def Run(X,Y,  prms = None):
    
    n,p = X.shape
    m = Y.shape[0]

    ## Get the parameter of the X matrix as a sparse matrix (csc)
    ## c.f. csc_matrix of scipy.sparse, nnx, indptr and indices parameters
    #    Xsparse = scipy.sparse.csc_matrix(X.transpose())
    #    nnz = Xsparse.nnz
    #    jc = Xsparse.indptr
    #    ir = Xsparse.indices
    jc = np.sum(X.transpose(),axis = 0)
    for i in range(1,len(jc)):
        jc[i] = jc[i]+jc[i-1]
    
    jc = np.insert(jc,0,0)
    jc = jc.tolist()
    
    ir = [i for j in range(n) for i in range(p) if X[j,i]]
        
    if prms == None:
        prms = prsampPrms([m,n])
        
    H,eps = prsamp.calibration(jc, ir, Y.tolist(), prms.delta, prms.priorPrms, prms.maxIter, prms.prec, \
                       prms.damp, prms.vnf,prms.display)
    return np.array(H).reshape(m,n),eps

   


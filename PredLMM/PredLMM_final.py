import pandas as pd
import numpy as np
import sys
import os
import time
import h5py
from copy import copy
from scipy.sparse import csr_matrix, rand
from struct import unpack, calcsize
from numpy.linalg import inv
from scipy.optimize import newton
from random import sample
from scipy.linalg.blas import dgemm,sgemm,sgemv
from bed_reader import open_bed, sample_file

def Identity(p):
 Identity = np.zeros((p, p), int)
 np.fill_diagonal(Identity, 1)
 return(Identity)

def Compute_blocks(Genotype, sub_sample, weights):
    X = Genotype; N = Genotype.shape[0]; p = Genotype.shape[1]
    p_allele = np.mean(X, axis=0)/2
    X = (X - 2*p_allele) / np.sqrt(2*p_allele*(1-p_allele))
    weights_star = weights*p/np.sum(weights)
    X =  np.multiply(np.sqrt(weights_star)[np.newaxis,:], X)
    subsample_size = len(sub_sample)
    non_subsample = np.setdiff1d(range(0,N),sub_sample)
    indices = np.hstack((sub_sample,non_subsample))
    X = X[indices, :]
    X_1 =X[range(0,subsample_size),]
    X_2 =X[range(subsample_size, N),]
    G_1 = sgemm(alpha=1,a=X_1,b=X, trans_b = 1)/p
    diag_22 = np.sum(np.square(X_2), axis = 1)/p
    GRM_array = csr_matrix((N, N)); GRM_array = GRM_array.todense().A1.reshape(N, N)
    GRM_array[np.ix_(range(0,subsample_size), range(subsample_size, N))] = G_1[:, range(subsample_size, N)]; 
    GRM_array = GRM_array + GRM_array.T
    GRM_array[np.ix_(range(0,subsample_size), range(0, subsample_size))] = G_1[:, range(0,subsample_size)]
    np.fill_diagonal(GRM_array, np.hstack((np.diag(GRM_array[np.ix_(range(0,subsample_size), range(0, subsample_size))]), diag_22)))
    return GRM_array



def derivative_minim_sub(y_sub, X_sub, X_subT, G_selected, A_selc, subsample_size):
 def smaller_predproc_exponential(param):
  h = param
  C_inv = inv(h*G_selected+(1-h)*Identity(subsample_size))
  C_invX = sgemm(alpha=1,a=C_inv,b=X_sub)
  beta = sgemm(alpha=1,a=inv(sgemm(alpha=1,a=X_subT.reshape(1,subsample_size),b=C_invX)),b=sgemm(alpha=1,a=C_invX,b=y_sub,trans_a=1))
  residual = (y_sub.reshape(subsample_size,1) - np.matmul(X_sub.reshape(subsample_size,1),beta))
  C_invResid = sgemm(alpha=1,a=C_inv,b=residual,trans_b=0)
  qf = sgemm(alpha=1,a=residual,b=C_invResid,trans_a=1)
  diff1 = np.sum(np.multiply(C_inv, A_selc))-subsample_size/qf*sgemm(alpha=1,a=C_invResid.T,b=sgemm(alpha=1,a=A_selc,b=C_invResid))
  #print(h)
  return(diff1)

 start_time = time.time()
 try:
  pc_minimizer_easy = newton(smaller_predproc_exponential,0.5,tol=0.0000001)
 except:
  pc_minimizer_easy=np.array([0])
 if pc_minimizer_easy>1:
    pc_minimizer_easy = np.array([1])   
 if pc_minimizer_easy<0:
    pc_minimizer_easy = np.array([0])
 h = pc_minimizer_easy
 C_inv = inv(h*G_selected+(1-h)*Identity(subsample_size))
 C_invX = sgemm(alpha=1,a=C_inv,b=X_sub)
 beta = sgemm(alpha=1,a=inv(sgemm(alpha=1,a=X_subT.reshape(1,subsample_size),b=C_invX)),b=sgemm(alpha=1,a=C_invX,b=y_sub,trans_a=1))
 residual = (y_sub.reshape(subsample_size,1) - np.matmul(X_sub.reshape(subsample_size,1),beta))
 C_invResid = sgemm(alpha=1,a=C_inv,b=residual,trans_b=0)
 sigma = sgemm(alpha=1,a=residual,b=C_invResid,trans_a=1)/subsample_size
 GRM_array_sub = sgemm(alpha=1,a=C_inv,b=A_selc) #V_pp^-1 A_ppc
 W = np.maximum(GRM_array_sub, GRM_array_sub.transpose() )
 a = np.sum(np.multiply(W,W))
 del C_inv; del W;
 sd_sub = np.sqrt(2/a)
 t1 = (time.time() - start_time)
 #result = np.hstack((np.asscalar(pc_minimizer_easy),np.asscalar(sd_sub),np.asscalar(sigma),t1))
 result = {'Heritability estimate': h, 'SD of heritability estimate': sd_sub.round(4), 'Variance estimate':    np.asscalar(sigma.round(4)), 'Time taken': (np.array(t1).round(4)).item()}
 return(result)


def derivative_minim_full(y, X, X_T, Ct, id_diag, add, G_selected, GRM_array, N):
 def Identity(p):
  Identity = np.zeros((p, p), int)
  np.fill_diagonal(Identity, 1)
  return(Identity)

 def der_predproc_exponential(param):
  h = param
  addedId = np.reshape((1-h)+ h*add,N)
  addedId_invU = np.multiply((1/addedId)[:,np.newaxis], Ct.T)
  CTadded_Id_invC = sgemm(alpha=1,a=Ct,b=addedId_invU)
  C_inv = (-sgemm(alpha=1,a=h*addedId_invU, b=sgemm(alpha=1,a=inv(G_selected+h*CTadded_Id_invC),b=addedId_invU.T)))
  np.fill_diagonal(C_inv,(1/addedId + C_inv[id_diag]))
  C_invX = sgemm(alpha=1,a=C_inv,b=X)
  beta = sgemm(alpha=1,a=inv(sgemm(alpha=1,a=X_T,b=C_invX)),b=sgemm(alpha=1,a=C_invX,b=y,trans_a=1))
  residual = (y.reshape(N,1) - np.matmul(X,beta)).T
  C_invResid = sgemm(alpha=1,a=C_inv,b=residual,trans_b=1)
  qf = sgemm(alpha=1,a=residual,b=C_invResid,trans_a=0)
  diff1 = np.sum(np.multiply(C_inv, GRM_array))-N/qf*sgemm(alpha=1,a=C_invResid.T,b=sgemm(alpha=1,a=GRM_array,b=C_invResid))
  del C_inv,addedId,addedId_invU,CTadded_Id_invC
  #print(h)
  return(diff1)

 start_time = time.time()
 pc_minimizer_f = newton(der_predproc_exponential,0.5,tol=0.000005)
 if pc_minimizer_f>1:
    pc_minimizer_f = 1   
 if pc_minimizer_f<0:
    pc_minimizer_f = 0
 h = pc_minimizer_f
 addedId = np.reshape((1-h)+ h*add,N)
 addedId_invU = np.multiply((1/addedId)[:,np.newaxis], Ct.T)
 CTadded_Id_invC = sgemm(alpha=1,a=Ct,b=addedId_invU)
 C_inv = (-sgemm(alpha=1,a=h*addedId_invU, b=sgemm(alpha=1,a=inv(G_selected+h*CTadded_Id_invC),b=addedId_invU.T)))
 np.fill_diagonal(C_inv,(1/addedId + C_inv[id_diag]))
 C_invX = sgemm(alpha=1,a=C_inv,b=X)
 beta = sgemm(alpha=1,a=inv(sgemm(alpha=1,a=X_T,b=C_invX)),b=sgemm(alpha=1,a=C_invX,b=y,trans_a=1))
 residual = (y.reshape(N,1) - np.matmul(X,beta)).T
 C_invResid = sgemm(alpha=1,a=C_inv,b=residual,trans_b=1)
 sigma = sgemm(alpha=1,a=residual,b=C_invResid,trans_a=0)/N
 GRM_array= sgemm(alpha=1,a=C_inv,b=GRM_array) #V_pp^-1 A_ppc
 W = np.maximum(GRM_array, GRM_array.transpose())
 a = np.sum(np.multiply(W,W))
 #print(a)
 del C_inv;
 sd = np.sqrt(2/a)
 t1 = (time.time() - start_time)
 #result = np.hstack((np.asscalar(pc_minimizer_f),np.asscalar(sd),np.asscalar(sigma),t1))
 result = {'Heritability estimate': (np.array(h)), 'SD of heritability estimate': sd.round(4), 'Variance estimate':    np.asscalar(sigma.round(4)), 'Time taken': (np.array(t1).round(4)).item()}
 return(result)

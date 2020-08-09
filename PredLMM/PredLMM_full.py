import numpy as np
import time as time
from numpy.linalg import inv
from scipy.optimize import newton
from scipy.linalg.blas import dgemm,sgemm,sgemv


def derivative_minim(y, X, X_T, Ct, id_diag, add, G_selected, GRM_array, N):
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
  print(h)
  return(diff1)

 start_time = time.time()
 pc_minimizer_f = newton(der_predproc_exponential,0.5,tol=0.000005)
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
 t1 = (time.time() - start_time)
 result = np.hstack((np.asscalar(pc_minimizer_f),np.asscalar(sigma),t1))
 return(result)

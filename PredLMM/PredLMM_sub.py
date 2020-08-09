import numpy as np
import time as time
from numpy.linalg import inv
from scipy.optimize import newton, root
from scipy.linalg.blas import dgemm,sgemm,sgemv


def derivative_minim(y_sub, X_sub, X_subT, G_selected, A_selc, subsample_size):
 def Identity(p):
  Identity = np.zeros((p, p), int)
  np.fill_diagonal(Identity, 1)
  return(Identity)

 def smaller_predproc_exponential(param):
  h = param
  C_inv = inv(h*G_selected+(1-h)*Identity(subsample_size))
  C_invX = sgemm(alpha=1,a=C_inv,b=X_sub)
  beta = sgemm(alpha=1,a=inv(sgemm(alpha=1,a=X_subT.reshape(1,subsample_size),b=C_invX)),b=sgemm(alpha=1,a=C_invX,b=y_sub,trans_a=1))
  residual = (y_sub.reshape(subsample_size,1) - np.matmul(X_sub.reshape(subsample_size,1),beta))
  C_invResid = sgemm(alpha=1,a=C_inv,b=residual,trans_b=0)
  qf = sgemm(alpha=1,a=residual,b=C_invResid,trans_a=1)
  diff1 = np.sum(np.multiply(C_inv, A_selc))-subsample_size/qf*sgemm(alpha=1,a=C_invResid.T,b=sgemm(alpha=1,a=A_selc,b=C_invResid))
  print(h)
  return(diff1)

 start_time = time.time()
 try:
  pc_minimizer_easy = newton(smaller_predproc_exponential,0.5,tol=0.0000001)
 except:
  pc_minimizer_easy=0
 h = pc_minimizer_easy
 C_inv = inv(h*G_selected+(1-h)*Identity(subsample_size))
 C_invX = sgemm(alpha=1,a=C_inv,b=X_sub)
 beta = sgemm(alpha=1,a=inv(sgemm(alpha=1,a=X_subT.reshape(1,subsample_size),b=C_invX)),b=sgemm(alpha=1,a=C_invX,b=y_sub,trans_a=1))
 residual = (y_sub.reshape(subsample_size,1) - np.matmul(X_sub.reshape(subsample_size,1),beta))
 C_invResid = sgemm(alpha=1,a=C_inv,b=residual,trans_b=0)
 sigma = sgemm(alpha=1,a=residual,b=C_invResid,trans_a=1)/subsample_size
 t1 = (time.time() - start_time)
 result = np.hstack((np.asscalar(pc_minimizer_easy),np.asscalar(sigma),t1))
 return(result)

import os
import numpy as np
import ctypes

def clinreg(libpath, x, y):
    _path = os.path.dirname(__file__)
    clib = np.ctypeslib.load_library(libpath, _path)
    cfstat = clib.fit
    cfstat.restype = ctypes.c_int
    cfstat.argtypes = [np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED'),
                       ctypes.c_int,
                       ctypes.c_int,
                       ctypes.c_int,
                       np.ctypeslib.ndpointer(ctypes.c_double, ndim=1, flags='C_CONTIGUOUS, ALIGNED')
                      ]
    score = np.zeros(1)
    success = cfstat(x, y, 1, 1, 450, score)
    return score

def fstat(x, y):
    libfile = "/usr/users/sbanerj/trans-eQTL/codebase/tejaas/lib/linear_regression.so"
    return clinreg(libfile, x, y)

def zstat(x, y):
    libfile = "/usr/users/sbanerj/trans-eQTL/codebase/tejaas/lib/linear_regression_zstat.so"
    return clinreg(libfile, x, y)


datax = np.loadtxt("datax.txt")
datay = np.loadtxt("datay.txt")

print( "F-stat:", fstat(datax, datay))
print( "Z-stat:", zstat(datax, datay))

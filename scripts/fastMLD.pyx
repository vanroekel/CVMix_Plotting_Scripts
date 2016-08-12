#cython: boundscheck=False
from numpy import empty
cimport numpy as np
ctypedef np.float64_t DOUBLE

cdef interpBWlevels(y0,y1,x0,x1,xT):
    coeff2=(y1-y0)/(x1-x0)
    coeff1=y0 - coeff2*x0

    return coeff2*xT + coeff1

def computeThreshMLD(np.ndarray[DOUBLE,ndim=2] inField,double threshVal,\
                     np.ndarray[DOUBLE,ndim=2] zvals,double refDepth):
    """
    computes threshold based mixed layer depth

    This routine computes the mixed layer depth via a threshold method.  It
    is similar to that laid out in de Boyer Montegut et al. 2004.  It allows the
    user to specify the threshold (temperature or salinity or density).

    Parameters
    =======
    inField -- input data must be a 2d array of doubles with time dimension 
                then the vertical
    threshVal -- threshold to compute MLD on
    zvals -- depths (in meters or pressure)
    refDepth -- reference depth for threshold computation (must be in same unit 
                    zvals)

    NOTE
    =====
    The code does NOT verify units of different fields
    It also assumes that the row dimension of the input arrays are time and
    the column dimension is vertical levels
    """
    cdef int nz, nt
    cdef int k, refIndex, j
    cdef double dV, dVp1, mldT

    nz=inField.shape[1]
    nt=inField.shape[0]
    
    cdef np.ndarray[DOUBLE] mld
    mld = empty(shape = nt)
    for j in range(nt):
        for k in range(nz-1):
            if zvals[j, k+1] < refDepth:
                refVal = interpBWlevels(inField[j, k], inField[j, k+1], zvals[j, k],  \
                                      zvals[j, k+1], refDepth)
                refIndex = k
                break

        for k in range(refIndex, nz-1):
            if abs(inField[j, k+1]-refVal) >= threshVal:
                dV = abs(inField[j, k]-refVal)
                dVp1 = abs(inField[j, k+1]-refVal)

                mldT = interpBWlevels(zvals[j, k], zvals[j, k+1], dV, dVp1, threshVal)
                mldT = max(mldT, zvals[j, k+1])
                mld[j] = abs(min(mldT, zvals[j, k]))
            
#        print type(zvals[j, nz-1])
        mld[j] = abs(zvals[j, nz-1])
    return mld

def computeGradMLD(np.ndarray[DOUBLE, ndim = 2] inField, threshVal,   \
                   np.ndarray[DOUBLE, ndim = 2] zvals):
    """
    computes gradient based mixed layer depth

    This routine computes the mixed layer depth via a gradient method.  It allows the
    user to specify the gradient threshold (temperature or salinity or density).

    Parameters
    =======
    inField -- input data must be a 2d array of doubles with time dimension 
                then the vertical
    threshVal -- threshold to compute MLD on, gradient must exceed this value
    zvals -- depths (in meters)

    NOTE
    =====
    The code does NOT verify units of different fields
    It also assumes that the row dimension of the input arrays are time and
    the column dimension is vertical levels
    """
    cdef int j, k
    cdef int nz, nx, found
    
    nz = inField.shape[1]
    nt = inField.shape[0]
    
    cdef np.ndarray[DOUBLE] gradient, mld
    gradient = empty(shape = (nz))
    mld = empty(shape = (nt))
    if threshVal > 0:
        for j in range(nt):
            found=0
            for k in range(1, nz):
                dz = (zvals[j, k-1] - zvals[j, k])
                gradient[k] = (inField[j, k-1] - inField[j, k])/dz

            for k in range(1, nz-1):
                gradient[k]=(gradient[k-1] + gradient[k] + gradient[k+1])/3.0

            for k in range(1,nz-1):
                if gradient[k+1] >= threshVal:
                    mldT = interpBWlevels(zvals[j, k], zvals[j, k+1], gradient[k], gradient[k+1], threshVal)
                    mldT = max(mldT, zvals[j, k+1])
                    mld[j] = abs(min(mldT, zvals[j, k]))
                    found=1
                    break

            if found == 0:
                spot = gradient.argmax()
                mld[j] = abs(zvals[j, spot])
    else:
        for j in range(nt):
            found = 0
            for k in range(1, nz):
                dz = (zvals[j, k-1] - zvals[j, k])
                gradient[k] = (inField[j, k-1] - inField[j, k])/dz

            for k in range(1, nz-1):
                gradient[k]=(gradient[k-1] + gradient[k] + gradient[k+1])/3.0

            for k in range(1,nz-1):
                if gradient[k+1] <= threshVal:
                    mldT = interpBWlevels(zvals[j, k], zvals[j, k+1], gradient[k], gradient[k+1], threshVal)
                    mldT = max(mldT, zvals[j, k+1])
                    mld[j] = abs(min(mldT, zvals[j, k]))
                    found = 1
                    break
                    
            if found == 0:
                spot = gradient.argmin()
#            print 'end = ;',threshVal, gradient.max(), gradient.min(),gradient[0:3]
                mld[j] = abs(zvals[j, spot])
           
    return mld
            

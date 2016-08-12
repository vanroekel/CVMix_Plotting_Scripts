import numpy
cimport numpy as np
ctypedef np.float64_t DOUBLE
#cython: boundscheck=False

cdef double compMean(np.ndarray[DOUBLE, ndim=1] a):
    cdef int i
    cdef int n = a.shape[0]
    cdef double m = 0.0
    for i in range(n):
        m += a[i]
    m /= n
    return m

cdef double interpBWlevels(double y0,double y1,double x0,double x1,double xT):
    cdef int k,k2
#    cdef double Minv[3][3]
#    cdef double coeffs[3]
#    cdef double rhs[3]
    cdef np.ndarray Minv = numpy.zeros([3,3], dtype=numpy.float64)
    cdef np.ndarray coeffs = numpy.zeros(3, dtype=numpy.float64)
    cdef np.ndarray rhs = numpy.zeros(3, dtype=numpy.float64)
    cdef double det

    det=-((x1-x0)**2)
    rhs[0]=y1
    rhs[1]=y0
    rhs[2]=0.0
    Minv[0,0] = -1.0/det
    Minv[0,1] = 1.0/det
    Minv[0,2] = -1.0/(x1-x0)
    Minv[1,0] = 2.0*x0/det
    Minv[1,1] = -2.0*x0/det
    Minv[1,2] = (x1+x0)/(x1-x0)
    Minv[2,0] = -x0**2/det
    Minv[2,1] = x1*(2.0*x0-x1)/det
    Minv[2,2] = -x1*x0/(x1-x0)
    
    for k in range(3):
        for k2 in range(3):
            coeffs[k2]=coeffs[k2]+Minv[2-k2,k]*rhs[k]

    return coeffs[0]+coeffs[1]*xT+coeffs[2]*xT**2

cdef double compute_unresolved_vT(double crit_RI, double zmid, \
                                  double Cs, double N_cnt, \
                                  double ws, double sig):

    cdef double Vtc, Cv, kap
    cdef VT
    
    kap = 0.4
    Vtc = numpy.sqrt(0.2 / (Cs * sig))
        
    if N_cnt < 0.002:
        Cv = 2.1 - 200.0*N_cnt
    else:
        Cv = 1.7
    
    VT = -Cv * Vtc * abs(zmid) * N_cnt * ws / (kap**2 * crit_RI)
        
    return VT

cdef compute_phi(double zeta):
    
    cdef double a_s = numpy.sqrt(1.0 + 16.0)*(1.0 - 8.0)
    cdef double c_s = numpy.sqrt(1.0 + 16.0)*24.0
    
    cdef double phis = 0.0
    if zeta > -100:
        if zeta > 0:
            return c_s, 1.0 + 5.0*zeta
        elif zeta >= -1.0:
            return c_s, (1.0 - 16.0*zeta)**(-1.0/2.0)
        else:
            return c_s, (a_s - c_s*zeta)**(-1.0 / 3.0)
    else:
        return c_s, float(0.0)

cdef compute_w(double buoy_sfc, double ustar, double sig, double sigma, double BLD):
    cdef int k
    cdef double ws, zeta
    cdef double c_s, phis
    
    zeta = min(sigma, sig) * BLD * buoy_sfc * 0.4 / (ustar**3 + 1.E-15)
    
    c_s, phis = compute_phi(zeta)
    
    if ustar == 0:
        if buoy_sfc >= 0:
            ws = 0.0
        else:
            ws = -c_s*min(sig, sigma)*BLD*0.4*buoy_sfc
            ws = 0.4*(ws)**(1.0/3.0)
    else:
        ws = 0.4 * ustar / phis
    
    return c_s, ws
         
def computeBLD(int nscl, double crit_RI, np.ndarray[DOUBLE, ndim = 2] zmid, np.ndarray[DOUBLE,ndim = 2] ztop, \
               double ref_level, bint average, np.ndarray[DOUBLE, ndim = 2] temp, \
               np.ndarray[DOUBLE, ndim = 2] salt, np.ndarray[DOUBLE, ndim = 2] u,  \
               np.ndarray[DOUBLE, ndim = 2] v, np.ndarray[DOUBLE, ndim = 2] N2, double randval, \
               double buoy_sfc, double ustar, double refValT, double refValS, bint compVt):        

    cdef int k, t, j, spot
    cdef int nz, nt
    cdef double c_s, ws, sigVal, fact
    
    nt = zmid.shape[0]
    nz = zmid.shape[1]
#    cdef double ri[nz]
#    cdef double ztemp[nz]
#    cdef double bld[nt]
    
#    cdef np.ndarray[DOUBLE] buoy = numpy.empty(shape = nz)
    cdef np.ndarray[DOUBLE] ri = numpy.empty(shape = nz)
    cdef np.ndarray[DOUBLE] ztemp = numpy.empty(shape = nz)
    cdef np.ndarray[DOUBLE] ztempT = numpy.empty(shape = nz)
    cdef np.ndarray[DOUBLE] bld = numpy.empty(shape = nt)
#    cdef np.ndarray[DOUBLE] N_cnt = numpy.empty(shape = nz)
    cdef double refT, refS, refU, refV, Hval, deltab, deltaU
    
    refT = 0.
    refS = 0.
    refV = 0.
    refU = 0.
    
    if compVt:
        fact = 1.0
    else:
        fact = 0.0
                #find BLD via iterative method,  assume BLD is the depth of the layer
    for j in range(nt):    
        ztemp = abs(zmid[j, :])
        ztempT = abs(ztop[j,:])
         
        scaling = 1.0-ref_level*0.5
        for k in range(0,nz - 3):
            sigVal = ztemp[k] / ztempT[k+1]
            c_s, ws = compute_w(buoy_sfc, ustar, ref_level, sigVal, ztempT[k+1])
            VT = compute_unresolved_vT(crit_RI, ztempT[k+1], c_s, max(N2[j,k],0), ws, ref_level)

            VT = max(1E-10, abs(VT))
            
            f = numpy.where(ztemp < ref_level*ztempT[k+1])[0]
            if len(f) == 0:
                refT = temp[j,0]
                refS = salt[j,0]
                refU = u[j,0]
                refV = v[j,0]
            else:
                if average:
                    refT = compMean(temp[j, f])
                    refS = compMean(salt[j, f])
                    refU = compMean(u[j, f])
                    refV = compMean(v[j, f])
                else:
                    if randval is None:
                        spot=len(f)/2
                    else:
                        spot = randval*f[len(f)-1]
                        refT = temp[j, spot]
                        refS = salt[j, spot]
                        refU = u[j, spot]
                        refV = v[j, spot]

            
            deltaB = -9.81*2E-4*(temp[j,k]-refT) + 9.81*8E-4*(salt[j,k]-refS)
            deltaU = ((u[j,k]-refU)**2+(v[j,k]-refV)**2)
            Hval = ztemp[k]*scaling
            ri[k] = Hval*deltaB/(fact*VT + deltaU + 1E-15)
       
            if ri[k] > 10:
                break
        t = 0
        while ri[t] < crit_RI and t < len(ri) - 1:
            t += 1
  
        bld[j] = interpBWlevels(ztemp[t-1],ztemp[t],ri[t-1],ri[t],crit_RI)
    return bld     
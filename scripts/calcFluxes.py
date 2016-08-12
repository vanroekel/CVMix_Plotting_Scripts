def getFluxes(ntimes,var_array,diff_array,refZMid,non_local=0.0):
    import numpy as np
    diffsize=diff_array.shape
    if len(diffsize) != 2:
        print "Error you must choose one cell only at a time"
        return None
   
    nzm=refZMid.shape[1]
    nzt=nzm+1

    vert_grad=np.zeros((ntimes,nzt))

    for j in range(1,nzm):
        vert_grad[:,j] = diff_array[:,j]*(var_array[:,j-1] - var_array[:,j]) / (refZMid[:,j-1] - refZMid[:,j])

    localflux=vert_grad + non_local
#    return vert_grad    
#    return -np.multiply(diff_array,vert_grad) + non_local
    return vert_grad - non_local




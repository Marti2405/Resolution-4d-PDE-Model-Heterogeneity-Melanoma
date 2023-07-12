"""
    # divx2


    Purpose:
    Module takes as inputs the gradient of a 2-dimensional distribution
    and returns the derivative of this distribution through its
    spatial (1st & 2nd) dimensions, as a divergence          


    Input:
       3D vector of two coordinates ([.,.,.] [.,.,.])


    Output:
        3D divergence of that vector

"""
function divx2(in)   

    global dx
    

    dvx1 = in[[2:end...,end],:,:,1] - in[[1,1:end-1...],:,:,1] # central difference for unknown points
    dvx2 = in[:,[2:end...,end],:,2] - in[:,[1,1:end-1...],:,2] # and forward difference for boundaries
    
    # Boundaries 
    dvx1[1,:,:] *=2
    dvx1[end,:,:] *= 2
    dvx2[:,1,:] *= 2
    dvx2[:,end,:] *=2


    out = (1/(2*dx))*(dvx1 + dvx2)

    # corner boundaries
    out[1  ,1  ,:] /=2
    out[1  ,end,:] /=2
    out[end,1  ,:] /=2
    out[end,end,:] /=2
    
    return out
end
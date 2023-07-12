"""
    # divx


    Purpose:
    Module takes as inputs the gradient of a 4-dimensional distribution 
    and returns the derivative of this distribution through its structural 
    (1st & 2nd) dimensions, as a divergence             


    Input:
       4D vector of two coordinates ([.,.,.,.] [.,.,.,.])


    Output:
        4D divergence of that vector

"""

function divx(in)
    
    global dx
    
    dvx1 = in[[2:end...,end], :, :, :, 1] - in[[1,1:end-1...], :, :, :, 1] # central difference for all points except the boundaries where we do forward difference
    dvx2 = in[:, [2:end...,end], :, :, 2] - in[:, [1,1:end-1...], :, :, 2]

    ## Boundaries 
    # sides
    dvx1[1,:,:,:] *=2
    dvx1[end,:,:,:] *= 2
    dvx2[:,1,:,:] *= 2
    dvx2[:,end,:,:] *=2
    

    
    # Result
    out = (1 / (2 * dx))* (dvx1 + dvx2) 

    # corner boundaries
    out[1,1,:,:] /=2
    out[1,end,:,:] /=2
    out[end,1,:,:] /=2
    out[end,end,:,:] /=2

    
    
    return out
end


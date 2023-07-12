"""
    # divy


    Purpose:
    Module takes as inputs the gradient of a 4-dimensional distribution
    and returns the derivative of this distribution through its
    structural (3rd & 4th) dimensions, as a divergence

    Input:
       4D vector of two coordinates ([.,.,.,.] [.,.,.,.])


    Output:
        4D divergence of that vector

"""
function  divy(in)

    global dy

    dvy1 = in[:,:,[2:end...,end],:,1] - in[:,:,[1,1:end-1...],:,1] # Unkwnown points -> Central difference
    dvy2 = in[:,:,:,[2:end...,end],2] - in[:,:,:,[1,1:end-1...],2] # Boundaries -> Forward difference

    # Boundaries 
    dvy1[:,:,1,:] *= 2
    dvy1[:,:,end,:] *= 2
    dvy2[:,:,:,1] *= 2
    dvy2[:,:,:,end] *= 2


    out = (1/(2*dy))*(dvy1 + dvy2) 

    # corner boundaries
    out[:,:,1  ,1  ] /=2
    out[:,:,1  ,end] /=2
    out[:,:,end,1  ] /=2
    out[:,:,end,end] /=2
    

    return out
end
"""
    # grady


    Purpose:
    Module takes as inputs a 4-dimensional distribution and returns the
    derivative of this distribution through its structural (3rd & 4th)
    dimensions, as a gradient vector

    Input:
        4D array


    Output:
        4D gradient vector of two coordinates ([.,.,.,.] [.,.,.,.])

"""
function grady(in)

    global dy

    out = zeros(size(in)...,2)

    out[:,:,:,:,1] = in[:,:,[2:end...,end],:] - in[:,:,[1,1:end-1...],:] # Central diff
    out[:,:,:,:,2] = in[:,:,:,[2:end...,end]] - in[:,:,:,[1,1:end-1...]] # Forward diff for boundaries


    # Boundaries 
    out[:,:,1,:,1] *=2
    out[:,:,end,:,1] *= 2
    out[:,:,:,1,2] *= 2
    out[:,:,:,end,2] *=2

    
    out = (1/(2*dy))*out

    # corner boundaries
    out[:,:,1  ,1  ,:] /=2
    out[:,:,1  ,end,:] /=2
    out[:,:,end,1  ,:] /=2
    out[:,:,end,end,:] /=2

    return out

end
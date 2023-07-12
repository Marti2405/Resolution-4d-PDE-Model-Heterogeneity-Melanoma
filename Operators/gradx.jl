"""
    # gradx


    Purpose:
    Module takes as inputs a 4-dimensional distribution and returns the
    derivative of this distribution through its spatial (1st & 2nd)
    dimensions, as a gradient vector

    Input:
       4D array 


    Output:
        4D gradient vector of two coordinates of that array ([.,.,.,.] [.,.,.,.])

"""
function gradx(in)

    global dx

    out = zeros(size(in)...,2)

    out[:,:,:,:,1] = in[[2:end...,end],:,:,:] - in[[1,1:end-1...],:,:,:] # Central dif for all points
    out[:,:,:,:,2] = in[:,[2:end...,end],:,:] - in[:,[1,1:end-1...],:,:] # Forward diff for boundaries

    # Boundaries 
    out[1,:,:,:,1] *=2
    out[end,:,:,:,1] *= 2
    out[:,1,:,:,2] *= 2
    out[:,end,:,:,2] *=2


    out = (1/(2*dx))*out 

    # corner boundaries
    out[1  ,1  ,:,:] /=2
    out[1  ,end,:,:] /=2
    out[end,1  ,:,:] /=2
    out[end,end,:,:] /=2

end
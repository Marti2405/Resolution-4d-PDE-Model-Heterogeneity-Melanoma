"""
    # gradx2


    Purpose:
    Module takes as inputs a 2-dimensional distribution 
    and returns the derivative of this distribution through 
    its spatial (1st & 2nd) dimensions, as a gradient vector.

    Input:
       3D array[x1,x2,s]


    Output:
        4D array[x1-1,x2-1,s,1:2]  
        The last index represents the column index in the matrix. 
        Matrix given by the product of the gradient with 'in'

"""
function  gradx2(in)

    global dx

    out = zeros(tuple(size(in)..., 2))
    
    
    out[:,:,:,1] = in[[2:end..., end],:,:] - in[[1,1:end-1...],:,:] # Central diff
    out[:,:,:,2] = in[:,[2:end..., end],:] - in[:,[1,1:end-1...],:] # Forward diff for boundaries
    

    # Boundaries 
    out[1,:,:,1] *=2
    out[end,:,:,1] *= 2
    out[:,1,:,2] *= 2
    out[:,end,:,2] *=2


    out = (1/(2*dx))*out 

    # corner boundaries
    out[1  ,1  ,:,:] /=2
    out[1  ,end,:,:] /=2
    out[end,1  ,:,:] /=2
    out[end,end,:,:] /=2

    return out
end
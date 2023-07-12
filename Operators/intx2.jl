"""
    # intx2


    Purpose:
        Module takes as inputs a 2-dimensional distribution and returns the
        integral of this distribution through its spatial (1st & 2nd)
        dimensions

    Input:
        2D array


    Output:
        Float 

"""
function intx2(in)

    global dx

    out1 = dropdims((sum(in;dims=2) + sum(in[:,2:end-1];dims=2));dims=2) # midpoint integration method
    out = dropdims((sum(out1;dims=1) + sum(out1[2:end-1];dims=1));dims=1)

    out = ((0.5*dx)^2)*out # We divide by 2 for the midpoint method and multiply by dx
                           # to ensure correct normalization and account for the grid spacing
                           # ^2 to account for both dimensions

    return out[1] # to return a float and not an array containing a float
end
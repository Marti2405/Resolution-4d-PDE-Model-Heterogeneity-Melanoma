"""
    # intx


    Purpose:
        Module takes as inputs a 4-dimensional distribution and returns the
        integral of this distribution through its spatial (1st & 2nd)
        dimensions

    Input:
        4D array


    Output:
        2D integral of array

"""
function intx(in)

    global dx

    out1 = dropdims((sum(in;dims=2) + sum(in[:,2:end-1,:,:];dims=2));dims=2) # we use midpoint integration method
    out = dropdims((sum(out1;dims=1) + sum(out1[2:end-1,:,:];dims=1));dims=1)

    out = ((0.5*dx)^2)*out # We divide by 2 for the midpoint method and multiply by dx
                           # to ensure correct normalization and account for the grid spacing
                           # ^2 to account for both dimensions

    return out
end
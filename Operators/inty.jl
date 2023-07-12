"""
    # inty


    Purpose:
        Module takes as inputs a 4-dimensional distribution and returns the
        integral of this distribution through its structural (3rd & 4th)
        dimensions

    Input:
        4D array


    Output:
        2D integral result

"""
function inty(in)

    global dy

    out1 = sum(in;dims=4) + sum(in[:,:,:,2:end-1];dims=4)    # check if dropdims needed
    out = sum(out1;dims=3) + sum(out1[:,:,2:end-1];dims=3)
    
    out = ((0.5*dy)^2)*out # We divide by 2 for the midpoint method and multiply by dx
                           # to ensure correct normalization and account for the grid spacing
                           # ^2 to account for both dimensions

    return dropdims(dropdims(out;dims=4);dims=3)

end

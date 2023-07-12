"""
    # projx


    Purpose:
        Module takes as inputs a 2-dimensional spatial distribution and
        returns the numerical projection of this distribution through an
        additional 2 structural dimensions

    Input:
        2D array


    Output:
        4D projection of the array

"""
function projx(in)

    global y

    leny = length(y)

    out = repeat(in,1,1,leny,leny)

    return out

end
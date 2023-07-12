"""
    # projx


    Purpose:
        Module takes as inputs a 2-dimensional structural distribution and
        returns the numerical projection of this distribution through an
        additional 2 spatial dimensions

    Input:
        2D array


    Output:
        4D projection of the array

"""
function projy(in1)

    global x, y

    lenx = length(x)
    leny = length(y)

    iny = zeros(1,1,leny,leny)
    iny[1,1,:,:] = in1
    out = repeat(iny,lenx,lenx,1,1)

    return out
end
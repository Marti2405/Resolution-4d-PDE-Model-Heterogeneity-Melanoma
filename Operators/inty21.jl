"""
    # inty21


    Purpose:
        Module takes as inputs a 2-dimensional distribution and returns the
        integral of this distribution through its 2nd dimension

    Input:
        4D array


    Output:
        2D integral result

"""
function inty21(in)

    global dy

    out = sum((in[:,1:end-1] + in[:,2:end]);dims=2)  # check if dropdims needed

    out = (0.5*dy)*out

    return out
end
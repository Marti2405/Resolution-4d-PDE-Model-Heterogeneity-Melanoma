############################################################################
###### Module takes as inputs the melanoma, ECNE, chemical molecule, and
###### drug distributions for the simulation of the derivative of the ECNE
###### distribution w.r.t. time
###### #####################################################################

function partialv(c,v,m,p)

    phi_v,del_v,del_v0 = v_param()

    ## Source Terms

    ### Remodelling

    Rem = phi_v*v.*(1 .-v)

    ### Degradation

    Degrad = -del_v*m[:,:,2].*v - del_v0*v


    ## Total

    out = Rem + Degrad

    return out

end
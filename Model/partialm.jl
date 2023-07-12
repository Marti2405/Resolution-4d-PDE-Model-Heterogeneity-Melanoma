###########################################################################
##### Module takes as inputs the melanoma, ECNE, chemical molecule, and
##### drug distributions for the simulation of the derivative of the 
##### chemical molecule distribution w.r.t. time
##### #####################################################################

function partialm(c,v,m,p)

    Dm,phi_m,del_m = m_param()

    ## spatial dynamics

    ### Diffusion

    Coeff = zeros(size(m)...,2)
    Coeff[:,:,1,:] .= Dm[1]
    Coeff[:,:,2,:] .= Dm[2]
    Diff = divx2(Coeff.*gradx2(m))


    ## Source Terms

    ### Proliferation

    Prol = zeros(size(m)...)
    Prol[:,:,1] = phi_m[1]*(1 .-m[:,:,1]).*v
    Prol[:,:,2] = phi_m[2]*(1 .-m[:,:,2]).*inty(c)

    ### Degradation

    Degrad = zeros(size(m)...)
    Degrad[:,:,1] = -del_m[1]*m[:,:,1].*inty(c)
    Degrad[:,:,2] = -del_m[2]*m[:,:,2]


    ## Total

    out = Diff + Prol + Degrad

    return out

end
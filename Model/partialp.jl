###########################################################################
##### Module takes as inputs the melanoma, ECNE, chemical molecule, and
##### drug distributions for the simulation of the derivative of the drug
##### molecule distribution w.r.t. time
##### #####################################################################

function partialp(c,v,m,p)

    Dp,del_p,ndel_p = p_param()

    ## spatial dynamics

    ### Diffusion

    Coeff = zeros(size(p)...,2)
    Coeff[:,:,1,:] .= Dp[1]
    Coeff[:,:,2,:] .= Dp[2]
    Diff = divx2(Coeff.*gradx2(p))


    ## Source Terms

    ### Degradation

    Degrad = zeros(size(p)...)
    Degrad[:,:,1] = -(del_p[1]*inty(c) .+ ndel_p[1]).*p[:,:,1]
    Degrad[:,:,2] = -(del_p[2]*inty(c) .+ ndel_p[2]).*p[:,:,2]


    
    ## Total

    out = Diff + Degrad

    return out

end
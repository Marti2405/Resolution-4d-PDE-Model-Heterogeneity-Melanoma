###########################################################################
##### Module takes as inputs the melanoma, ECNE, chemical molecule, and
##### drug distributions for the simulation of the derivative of the 
##### melanoma cell distribution w.r.t. time
##### #####################################################################

function partialc(c,v,m,p)

    global Np, y, y1, y2

    
    Dx,chi_m,chi_v,Dy,r_mu,mLin,Phi,cPreg,del_c,pDreg = c_param()

    ## spatial dynamics

    ### Spatial Diffusion

    Coeff = zeros(size(c)...,2)
    Coeff[:,:,:,:,1] .= Dx[1] 
    Coeff[:,:,:,:,2] .= Dx[2]
    Diffx = divx(Coeff.*gradx(c))

    
    ### Spatial Advection (Chemo-/Haptotaxis)

    Cin = zeros(size(c)...,2)
    Cin[:,:,:,:,1] = c 
    Cin[:,:,:,:,2] = c
    Advecx = -divx(Cin.*(chi_m*gradx(projx(m[:,:,1])) +chi_v*gradx(projx(v))))

    ## metabolic dynamics

    ### Metabolic Diffusion

    Coeff = zeros(size(c)...,2)
    Coeff[:,:,:,:,1] = Dy[1]*projy(y2)
    Coeff[:,:,:,:,2] = 4*Dy[2]*projy((y1.-.5).^2)
    Diffy = divy(Coeff.*grady(c))

    ### Metabolic Advection

    Coeff[:,:,:,:,1] = projy(mLin[1])   
    Advecy = r_mu*divy(Coeff.*Cin)

    

    ## Source Terms

    ### Proliferation

    Prol = Phi*projy(cPreg).*projx(m[:,:,1].*(1 .-inty(c))).*c

    ### Degradation

    Degrad = zeros(size(c)...)
    for i in 1:Np
        Degrad = Degrad - projx(del_c[i]*p[:,:,i]).*projy(pDreg[i]) 
    end
    Degrad = Degrad.*c

    ## Total

    out = Diffx + Advecx
    out = out + Diffy + Advecy
    out = out + Prol + Degrad

    return out

end
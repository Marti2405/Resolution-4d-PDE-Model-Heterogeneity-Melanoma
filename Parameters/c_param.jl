###########################################################################
##### Code written by Arran Hodgkinson; University of Exeter; January 2022
##### Module generates a set of parameters for the simulation of the
##### derivative of the melanoma cell distribution w.r.t. time
##### #####################################################################

function  c_param()

    global Np,y1,y2
    
    xi_max = 0.15; xi_up = 0.5; xi_hyp = 0.5
    nu_max = 0.15; nu_up = 0.2; nu_hyp = 0.8
    
    ## Spatial Diffusion/Advection Parameters
    
    Dx = [1e-5 1e-5]'
    chi_m = 2e-4
    chi_v = 2e-4
    
    ## Metabolic Diffusion/Advection Parameters
    
    Dy = [1e-3 2e-3]'
    r_mu = 12e-3
    
    r_shft = .8; r_min = .17
    mLin = Array{Any}(undef, 2) # Create array that can hold multiple data types
    mLin[1] = (y1 - ((xi_max.-r_shft)/(nu_max.-r_shft)).*(y2 .-r_shft).-r_shft).*((1 .- r_min)*y2 .+ r_min)
    mLin[2] = zeros(size(y1)...)
    
    ## Proliferation Parameters
    
    Phi = 0.5
    
    phi_0 = 0.08; phi_up = 0.45; phi_max = 0.47
    alpha = 100; w_up_xi = 2; w_up_nu = 35; w_max = 20; w_hyp = 0.25
    cPreg = (phi_0
        .+ phi_up.*exp.(-(w_up_xi*(y1.-xi_up).^2+w_up_nu*(y2.-nu_up).^2))
        + phi_max.*exp.(-w_max*((y1.-xi_max).^2+(y2.-nu_max).^2)))
        *(1 ./(1 .+exp.(-alpha*((y1.-xi_hyp).^2+(y2.-nu_hyp).^2 .-w_hyp^2))))
    
    ## Degradation Parameters
    del_c = [1,1]'
    
    pDreg = Array{Any}(undef,Np)
    alp_dxi = 7; alp_dnu = 8; alp_dxn = 8; alp_y1d = .4; alp_y2d = .2
    pDreg[1] = exp.(-alp_dxi.*((y1.-alp_y1d).^2)
        - alp_dnu.*((y2.-alp_y2d).^2) - alp_dxn.*((y1.-alp_y1d).*(y2.-alp_y2d)))
    # alp_dxi = 10; alp_dnu = 2; alp_dxn = 0; alp_y1d = .9; alp_y2d = .5
    # pDreg[2] = exp(-alp_dxi*((y1-alp_y1d).^2) ... # 1. localised in y1-y2
    #     - alp_dnu*((y2-alp_y2d).^2) - alp_dxn*((y1-alp_y1d).*(y2-alp_y2d)))
    # alp_dxi = 200; alp_dnu = 200; alp_dxn = 0; alp_y1d = 0; alp_y2d = 0
    # pDreg[2] = pDreg[2] + exp(-alp_dxi*((y1-alp_y1d).^2) ...
    #     - alp_dnu*((y2-alp_y2d).^2) - alp_dxn*((y1-alp_y1d).*(y2-alp_y2d)))
    # pDreg[2] = 1-pDreg[1]; #2. complimentary & equal in scope
    # pDreg[2] = intx2[pDreg[1]]*pDreg[2]./intx2[pDreg[2]]
    # frc0 = .31; pDreg[2] = 1-pDreg[1]; # 4. .314 complmnt w bckgrnd
    # pDreg[2] = (1-frc0)*pDreg[2] + frc0*max(pDreg[2](:))
    # pDreg[2] = intx2[pDreg[1]]*ones(size(y1)); # 5. homogeneous
    # [40,4,40,4] almost; decline followed by rapid relapse
    alp_dxi = 8; alp_dnu = 1; alp_dxn = 0; alp_y1d = 0; alp_y2d = .5
    pDreg[2] = exp.(-alp_dxi.*((y1.-alp_y1d).^2) # 6. localised in south-west
        - alp_dnu.*((y2.-alp_y2d).^2) - alp_dxn.*((y1.-alp_y1d).*(y2.-alp_y2d)))
    alp_dxi = 8; alp_dnu = 1; alp_dxn = 0; alp_y1d = 1; alp_y2d = 0
    pDreg[2] = pDreg[2] + exp.(-alp_dxi.*((y1.-alp_y1d).^2)
        - alp_dnu.*((y2.-alp_y2d).^2) - alp_dxn.*((y1.-alp_y1d).*(y2.-alp_y2d)))
    
    #=
    figure(11); clf; hold on; set(gca,"FontSize",16)
      imagesc(y1[:,1],y2[1,:],pDreg[1]'); axis([0 1 0 1]); #colorbar()
      xlabel("y_1','FontSize',20); ylabel('y_2','FontSize",20)
    figure(12); clf; hold on; set(gca,"FontSize",16)
      imagesc(y1[:,1],y2[1,:],pDreg[2]'); axis([0 1 0 1]); colorbar()
      xlabel("y_1','FontSize',20); ylabel('y_2','FontSize",20)
    =#

    return [Dx,chi_m,chi_v,Dy,r_mu,mLin,Phi,cPreg,del_c,pDreg]
    
end
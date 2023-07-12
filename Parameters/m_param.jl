###########################################################################
##### Code written by Arran Hodgkinson; University of Exeter; January 2022
##### Module generates a set of parameters for the simulation of the
##### derivative of the chemical species distribution w.r.t. time
##### #####################################################################

function m_param()

    ## Spatial Diffusion Parameters
    
    Dm = [2e-4 2e-4]'
    
    ## Production Parameters
    
    phi_m = [.1 .1]'
    
    ## Degradation Parameters
    
    del_m = [.2 .02]'
    
    return [Dm,phi_m,del_m]
 end
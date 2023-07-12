###########################################################################
##### Code written by Arran Hodgkinson; University of Exeter; January 2022
##### Module generates a set of parameters for the simulation of the
##### derivative of the ECNE distribution w.r.t. time
##### #####################################################################

function v_param()

    ## Remodelling Parameters
    
    phi_v = .02
    
    ## Degradation Parameters
    
    del_v = .05; del_v0 = .01
    
    return [phi_v,del_v,del_v0]

end
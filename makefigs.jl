# Code written by Arran Hodgkinson, University of Exeter, January 2022
# Module takes input values from the melanoma treatment simulator and
# generates figures corresponding to the simulated results. Results are
# saved as images in the 'images' folder

import Plots

function makefigs(it, c, v, m, p, t, vC, svC, ptt, svfg, svpth)
    global Np, dt, x, y, y1, y2, dx, dy
    svfg = svfg==1
    sf = 0.3  # scaling factor

    

    # Figure Cy,t
        hdlx = round(2 * (length(x) - 1)) + 1
        if svfg
            F = Plots.plot(visible = false)
        else
            Plots.plot(5)
        end

        CC = zeros(hdlx, hdlx)
        CC[1:2:end, 1:2:end] = intx(c)'
        CC[2:2:end, :] = 0.5 .* (CC[1:2:end-2, :] + CC[3:2:end, :])
        CC[:, 2:2:end] = 0.5 .* (CC[:, 1:2:end-2] + CC[:, 3:2:end])

        

        # Adjust dimensions of y to match CC
        num_columns = div(hdlx, 2) + mod(hdlx, 2)
        y_adjusted = repeat(y, 1, num_columns)
  
        # Trim or pad y_adjusted to match the desired size
        yCC = sort([y_adjusted[1:num_columns, :]...][1:hdlx-1])
        Plots.heatmap!(yCC, yCC, CC[1:end-1,1:end-1])
        #plot!(size=(800,800)) # define the size of the plot default=(600,400)
        Plots.xlims!(0, 1)
        Plots.ylims!(0, 1)
        Plots.xlabel!("y_1")
        Plots.ylabel!("y_2")


        if svfg
            savefig(F, joinpath("$(svpth)/Cy,t", "Cy,t=$it.png"))
            #println("Image Cy,t=$it.png saved")
        end
    # end figure Cy,t
    

    
    # Figure Cx,t
        hdly = round(2 * (length(y) - 1)) + 1
        if svfg
            F = Plots.plot(visible = false)
        else
            Plots.plot(5)
        end

        CC = zeros(hdly, hdly)
        CC[1:2:end, 1:2:end] = inty(c)'
        CC[2:2:end, :] = 0.5 .* (CC[1:2:end-2, :] + CC[3:2:end, :])
        CC[:, 2:2:end] = 0.5 .* (CC[:, 1:2:end-2] + CC[:, 3:2:end])

        # Adjust dimensions of x to match CC
        num_columns = div(hdly, 2) + mod(hdly, 2)
        x_adjusted = repeat(x, 1, num_columns)
  
        # Trim or pad y_adjusted to match the desired size
        xCC = sort([x_adjusted[1:num_columns, :]...][1:hdly-1])
        Plots.heatmap!(xCC, xCC, CC[1:end-1,1:end-1])
        #plot!(size=(800,800)) # define the size of the plot default=(600,400)
        Plots.xlims!(0, 1)
        Plots.ylims!(0, 1)
        Plots.xlabel!("x_1")
        Plots.ylabel!("x_2")
        
        if svfg
            savefig(F, joinpath("$(svpth)/Cx,t", "Cx,t=$it.png"))
            #println("Image Cx,t=$it.png saved")
        end
    # end figure Cx,t


#    # Figure vC
#        ln = 0.1
#        fcol = [0.5 1 0.5; 1 1 0.5]
#        if svfg
#            F = Plots.plot(visible = false)
#        else
#            Plots.plot(5)
#        end
#
#        vcmax = findfirst((vC[2:end] .- vC[1:end-1]) .< 0)
#        
#        if isnothing(vcmax)
#            vcmax = vC[ppos]
#        else
#            vcmax = vC[vcmax]
#        end
#    
#        
#        pmx = 1800
#        ivcm = pmx / vcmax
#        println("size icvm : $(size(ivcm))")
#        if isnan.(ivcm)
#            ivcm[1]=1
#        end
#    
#        for i in 1:Np
#            
#            # lptt = size(ptt[i])
#            lptt = 0 # debug : size throws error for ptt array
#            
#            for j in 1:lptt
#                if isinf(ptt[i][j, 1])
#                    continue
#                else
#                    if isinf(ptt[i][j, 2])
#                        crr = [ptt[i][j, 1], maximum(t)]
#                    else
#                        crr = ptt[i][j, :]
#                    end
#                    fill(sf .* [crr, flip(crr)], ivcm .* [pmx, pmx, 0, 0], fcol[i, :], "EdgeAlpha", 0, "FaceAlpha", 0.5)
#                end
#                if j < lptt
#                    if isinf(ptt[i][j+1, 1])
#                        crr = [ptt[i][j, 2], maximum(t)]
#                    else
#                        crr = [ptt[i][j, 2], ptt[i][j+1, 1]]
#                    end
#                    fill(sf .* [crr, flip(crr)], ivcm .* [pmx, pmx, 0, 0], [1, 0.5, 0.5], "EdgeAlpha", 0, "FaceAlpha", 0.5)
#                end
#            end
#        end
#
#        x8 = sf .*t 
#        y8 = ivcm.*vC
#        plot(x8,y8)
#
#        if !isempty(svC)
#            plot(sf .* t, ivcm .* svC[1] .* ones(size(t)))
#            plot(sf .* t, ivcm .* svC[2] .* ones(size(t)))
#        end
#
#        
#        xlabel!("t (days)", FontSize = 20)
#        ylabel!("Volume (mm^3)", FontSize = 20)
#
#        if svfg
#            savefig(F, joinpath("$(svpth)", "$vC.png"))
#            println("Image vC.png saved")
#        end
#    # end figure vC


#    # Figure p,t 
#        Plots.heatmap(x, x, dropdims(sum(p;dims=3);dims=3))
#        #plot!(size=(800,800)) # define the size of the plot default=(600,400)
#        Plots.xlims!(0, 1)
#        Plots.ylims!(0, 1)
#        Plots.xlabel!("x_1")
#        Plots.ylabel!("x_2")
#        
#        if svfg
#            savefig(F, joinpath("$(svpth)/p,t", "p,t=$it.png"))
#            println("Image p,t=$it.png saved")
#        end
#
#    # end figure p,t


end

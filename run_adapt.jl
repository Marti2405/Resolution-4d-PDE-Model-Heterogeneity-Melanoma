"


Code description : 


    Number of treatments (`Np`), time step size (`dt`), maximum simulation time (`Tmx`), 
    arrays for time (`t`), spatial coordinates (`x`, `y`, `y1`, `y2`), and initial distributions (`c`, `v`, `m`, `p`)

    A figure is created to visualize the initial distribution `v0` and saved as an image

    The code enters a loop that iterates over each time step `i` in the time array `t`
        Within the loop, there is a treatment regimen condition that checks 
        if a treatment is ongoing or not based on specific criteria
        If a treatment is ongoing (`trt == 1`), the drug distribution `p` is 
        modified according to the treatment regimen. The drug concentration is 
        increased for a specified duration, and then decreased afterward. 
        Other drug distributions `p` not associated with the ongoing treatment are decreased
    
    Next, the code performs the predictor step, where the variables `cnew`, `vnew`, `mnew`, and `pnew` are updated 
    using the `partialc`, `partialv`, `partialm`, and `partialp` functions, respectively

    The corrector step follows, where the variables `c`, `v`, `m`, and `p` are updated by taking a 
    weighted average of the predictor step values and the values obtained from the partial differential equations

    Certain variable values are constrained to ensure they remain within valid ranges

    Figures are generated and saved periodically, and the tumor volume `vC` is updated based on the integral of `c`
    At the end of the simulation, the final results are saved and additional figures are created
    The total time taken for the simulation is printed
    

"


using Dates,Plots




# Display line
function line()
    println("----------------------------")
end




##-----------------------------------------------------------------------
##                  Import all the necessary functions


#folders path
folder_path = ["Operators","Model","Parameters"]

line()
println("Loading functions...\n")
println("Folder Main :")
fmain = ["init.jl", "makefigs.jl"]
for file in fmain
    include(file)
    println("   $(file) included")
end

# Loop through the files and include each file
for folder in folder_path
  # Get a list of files in the folder
  file_list = readdir(folder)
  println("Folder : ",folder)

  for file in file_list
    include(joinpath(folder, file))
    println("   $(file) included")
  end
end

println("\nAll Functions loaded successfully")
line()



# Create directory for saving images
svpth = "images/run_adapt/"
##-----------------------------------------------------------------------








##----------------------------------------------------------------------
##                             Run Adapt

# parameters
dt = 0.2
Tmx = 360 #360 Number of days
t = 0:dt:Tmx
Np = 2 # number of drugs
ptt = Vector{Any}(undef, 2)
stt = 40 # starting day of the treatment
svC = [0.055, 0.01] # thresholds for tumor volume
ttrt = -1
ptrt = 1 # type of drug (here we only apply drug 1)
ch = 0
trt = 0 # 0->no treatment,  1->treatement


dx = 0.02
dy = 0.02
global x = 0:dx:1
global y = 0:dy:1
y1 = repeat(y, 1, length(y))
y2 = repeat(y', length(y), 1)



# Initialize variables
(c, v, m, p) = init()
v0 = copy(v)
vC = zeros(length(t)) # volume of the tumor in function of time
treatment_log = arr = Array{Bool}(undef, length(t)) # to know when treatment was applied
println("Initialisation successfull")
line()






t_start = time() # save the time
tpos = 0 #init tpos
iter = 0 # track number of iterations

first_ptt = true # boolean to check the first time that we put elements in the Vector ptt
mean_time = 0 # compute mean time of iteration


for i in t 
    time_iteration = time() # timer to track performance
    global iter = i
    global tpos +=1

    #println("Computing PDE at t=$(iter)")
    
    
    # treatement regimen
    if (i>=stt) && rem(i, 1) < 1e-4 # if the treatment start has been reached
        
        if vC[tpos-1] >= svC[1] # if the volume of our tumor is over the top threshold
            println("Tumor over threshold, we apply treatment. Volume:$(vC[tpos-1])") #---------------------------------------------------------
            global treatment_log[tpos] = true
            global ch = 0
            global trt = 1 # we apply tratment
            if ttrt<0
                global ttrt = i
                
                if first_ptt
                    ptt[ptrt] = [i,Inf]'
                    global first_ptt = false
                else
                    ptt[ptrt] = vcat(ptt[ptrt], [i, Inf]')
                end
    
                
            end
        elseif vC[tpos-1] < svC[2] # if the volume of the tumor is under the bottom threshold
            println("Tumor under threshold, we don't apply treatment. Volume:$(vC[tpos-1])") #---------------------------------------------------------
            ttrt = -1
            trt = 0 # we don't apply treatment
            if ch == 0
                ptt[ptrt][end, 2] = i
                ch = 1
            end
        else # if we are between the thresholds
            if treatment_log[tpos-1] # if we where in treatment we stay in treatment
                global treatment_log[tpos] = true
            end
        end
    
    end
    

    if trt == 1 # if we apply treatment

        if i >= ttrt && i < ttrt + 1
            global p[:, :, ptrt] += 0.3 * dt * v0

        elseif i >= ttrt + 1
            global p[:, :, ptrt] += 0.3 * dt * (1 .- p[:, :, ptrt]) .* v

        end
        
        # Degradate all drugs that are not being used
        for j in 1:Np # for j in the number of drugds
            if j == ptrt # if we are applying this drug skip
                continue
            end
            global p[:, :, j] -= 0.5 * dt * p[:, :, j] # if we are not applying this drug then apply drug degradtion
        end

    else # if we don't apply treatment

        global p -= 0.5 * dt * p # drug degradtion

    end
    


    # predictor step
    cnew = c + dt .* partialc(c, v, m, p)
    vnew = v + dt .* partialv(c, v, m, p)
    mnew = m + dt .* partialm(c, v, m, p)
    pnew = p + dt .* partialp(c, v, m, p)
    # cnew, vnew, mnew, pnew = threaded_predictor_step(c, v, m, p, dt)
    
    
    # corrector step
    global c = 0.5 .* (c + cnew) .+ 0.5 .* dt .* partialc(cnew, vnew, mnew, pnew)
    global v = 0.5 .* (v + vnew) .+ 0.5 .* dt .* partialv(cnew, vnew, mnew, pnew)
    global m = 0.5 .* (m + mnew) .+ 0.5 .* dt .* partialm(cnew, vnew, mnew, pnew)
    global p = 0.5 .* (p + pnew) .+ 0.5 .* dt .* partialp(cnew, vnew, mnew, pnew)
    

    # remain in valid ranges
    global c[c .< 0] .= 0
    global v[v .< 0] .= 0
    global v[v .> 1] .= 1
    global m[m .< 0] .= 0
    global m[m .> 1] .= 1
    global p[p .< 0] .= 0
    global p[p .> 1] .= 1
    # c, v, m, p = threaded_update(c, v, m, p)

    
    

    
    vC[tpos] = intx2(inty(c)) # compute the number of cells (tumor volume)
    
    println("Tumor volume : $(vC[tpos])") # print the tumor volume
    
    

    # Figures
    if rem(i, 1) < 1e-4
        makefigs(i, c, v, m, p, t, vC, svC, ptt, 1, svpth) # create makefigs file

        time_current_iter = time()-time_iteration
        println("\nDone $(i/dt) in $(time_current_iter) seconds\n")
        global mean_time += time_current_iter
        sleep(0.01)
    end


end

mean_time = 1/(iter+5) * mean_time

line()
println("\nMean time for 1 iteration : $(mean_time) seconds\n")
line()
# save the evolution of the volume of the tumor
F = plot(t,vC)
savefig(F, joinpath("$(svpth)v", "vC.png"))
println("Image vC.png Saved")
line()

# save v0 
F = plot(size=(800, 800), aspect_ratio=:equal, xlim=(0, 1), ylim=(0, 1), colorbar=true)
heatmap!(x, x, v0', aspect_ratio=:equal, clim=(minimum(v0), maximum(v0)), color=:jet)
savefig(joinpath("$(svpth)v", "v0.png"))
println("Image v0.png Saved")
line()


# save v at final time
F = plot(size=(800, 800), aspect_ratio=:equal, xlim=(0, 1), ylim=(0, 1), colorbar=true)
heatmap!(x, x, v', aspect_ratio=:equal, clim=(minimum(v0), maximum(v0)), color=:jet)
savefig(joinpath("$(svpth)v", "vT.png"))
println("Image vT.png Saved")
line()



# show time of execution
line()
println("\nCompleted in $(time()-t_start) seconds\n")
line()
println("\n")
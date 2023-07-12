# Animation Creation

using Images, Plots



function create_visual(animation_name,image_folder,name_start,name_end,video_duration_sec=30)

    img_names = readdir(image_folder)
    number_images = size(img_names)[1]
    
    frames_per_second = Int(round(number_images/video_duration_sec))


    # Load the images into an array
    
    images = [load(joinpath(image_folder, "$(name_start)$(i)$(name_end)")) for i in 1:number_images-1]


    frames = length(images)

    # Creation of animation in function of the time
    println("Creating animation for the folder $image_folder")
    anim = @animate for i in 1:frames
        plot(images[i],xaxis=false,yaxis=false,title="$(animation_name[1:end-4])") # update plot
        
    end

    #save animation
    gif(anim,"$animation_name",fps=frames_per_second)
end


create_visual("Cy,t_Anim.mp4","images/run_adapt/Cy,t","Cy,t=",".0.png")
create_visual("Cx,t_Anim.mp4","images/run_adapt/Cx,t","Cx,t=",".0.png")
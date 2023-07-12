###########################################################################
##### Module generates initial conditions for the melanoma treatment
##### simulator; with a gaussian probability distribution for the cancer
##### cell population & an arbitrary form chosen for the extra-cellular
##### nutritional environment [ECNE], v. Chemical & drug species are 0
##### #####################################################################

function  init()

    global x,y1,y2

    lenx = length(x) #leny = length(y)
    x1 = repeat(x,1,lenx)
    x2 = repeat(x',lenx,1)
    X1 = projx(x1)
    X2 = projx(x2)
    Y1 = projy(y1)
    Y2 = projy(y2)
    
    # c = zeros(lenx,lenx,leny,leny)
    c = exp.(-250*((X1.-0.5).^2+(X2.-0.5).^2+(Y1.-0.2).^2+(Y2.-0.2).^2))
    
    
    # v = zeros(lenx,lenx); # v = .5*ones(lenx,lenx) - .5*inty[c]
    xx1 = (1/3)*(x1.+1.5); xx2 = (1/3)*(x2.+1.5); ze = 12*pi
    h = 0.5 .+ 0.5*sin.((ze.*xx1)./(xx2.+1)).*sin.(ze*xx1.*xx2).*sin.((ze*(1.0.-xx1))./(xx2.+1)).*sin.(ze*(xx1.-1).*(xx2.-1))

    v = min.(h,1.0.-inty(c))
    
    m = zeros(lenx,lenx,2)
    m[:,:,1] = v
    
    p = zeros(lenx,lenx,2)
    # p[:,:,1] = v; p[:,:,2] = v
    
    return [c,v,m,p]
end


module SimpleFlashLag

export flashlag_sim,
       stimulation!,
       propagate!

value(x) = x == 0 ? 0.4 : x == 1 ? 0.2 : x == 2 ? -0.2 : 0.0
J(network) = [value(abs(i - j)) for i in 1:size(network,1), j in 1:size(network, 1)]
overall_input(network, layer, λ) = layer == 1 ? 0.0 : J(network)*Θ(network[:,layer-1] - λ)

Θ(x) = int(x .>= 0)

function propagate!(network, Ω, λ)
    _, layers = size(network)
    for j in layers:-1:1
        network[:, j] = (1 - Ω)*network[:,j] + overall_input(network, j, λ)
    end
end

stimulation!(network, i, j = 1) = network[i,j] += 1.5

function flashlag_sim(;Ω = 0.6, λ = 0.65, lines = 18,layers = 3, maxtime = 13, plot = false)
    #Networks
    moving = zeros(Float64, lines, layers)
    flash  = zeros(Float64, lines, layers)

    #History
    hist_moving = zeros(Float64, lines, layers, maxtime)
    hist_flash  = zeros(Float64, lines, layers, maxtime)

    #Dynamic
    time = 2:maxtime
    for t in time
        propagate!(moving, Ω, λ)         #propagate moving
        propagate!(flash,   Ω, λ)        #propagate flash
        t >= 4 && stimulation!(flash, 7) #flash on position (7,1) from t=4 onwards
        stimulation!(moving, t+2)        #moving stimulus starting from 3rd position
        hist_moving[:,:,t], hist_flash[:,:,t] = moving[:,:], flash[:,:] #save on history
    end

    return hist_moving, hist_flash
end

end #module

plt.figute(figsize = (10,10))
# After read.jl
plot()
j=1
for i=6:48
    global spm
    if i==6
        spm=SPMinter(P2[1:i],F1k[:,j][1:i],)
        plot!(spm,xlims=(-10,1))
    else
        spm=SPMinter(P2[1:i],F1k[:,j][1:i],)
        plot!(spm,xlims=(-10,1))
    end #if
end #for
plot!(xlims=(-3,1.5),ylims=(-1,1))
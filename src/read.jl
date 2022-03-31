using JLD2
using Plots
using Gaussquad
using SPMinterpolation
Ppoints=50
# k=gausslegendremesh(10. ^-4,10. ^4,32,2);
P2=Array{Float64}(undef,Ppoints,1)
F1k=Array{Float64}(undef,Ppoints,32)
kname=Array{String}(undef,Ppoints,1)

@time for i=1:Ppoints
    local a, b
    global P2,F1k
    a,b=load("/Users/kjy/Desktop/program/julia/Gamma5/data/F1k_等间距_broken50/F1k$i.jld2","P2", "F1k")
    P2[i]=a
    # F1k[i,:]=1 ./b
    F1k[i,:]=b
end #for i

# for i=1:32
#     global kname
#     kn=k[i]
#     kname[i]="k=$kn"
# end

# scatter(P2,F1k,xlim=[0,1.3],lable=k)

# 总程序
print("++++++++++++++++++++++++Begin++++++++++++++++++++++++++\n")
include("./ini.jl")
include("./mkChebyshevD.jl")
# plist=[i/32. for i=1:32]
plist=[0. ]

lengthplist=length(plist)
@time for indexforp2=1:lengthplist
    global P2
    P2=plist[indexforp2]
    include("./mkvar.jl")
    include("./mkkernel.jl")
    include("./getsolu.jl")
    # jldsave("/Users/kjy/Desktop/program/julia/Gamma5/data/F1k/F1k$indexforp2.jld2";P2, F1k)
    print("$P2 for $indexforp2/$lengthplist done\n")
end

# include("./eva_kernel_new.jl")

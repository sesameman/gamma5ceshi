using Plots
kernel=[kernel11 kernel12 kernel13 kernel14;kernel21 kernel22 kernel23 kernel24;kernel31 kernel32 kernel33 kernel34;kernel41 kernel42 kernel43 kernel44;]
#δ=Diagonal(ones(4*dim))
# kernel_eva=Array{BigFloat}(undef, 4*dim, 4*dim);
right=[z4*ones(dim) ;zeros(3*dim)]
kernel_eva=I-kernel

solution=kernel_eva\right
F1=Array{Float64}(undef, kstep, zstep)
F2=Array{Float64}(undef, kstep, zstep)
F3=Array{Float64}(undef, kstep, zstep)
F4=Array{Float64}(undef, kstep, zstep)
FF=Array{Array}(undef,4,1)
FF[1]=F1
FF[2]=F2
FF[3]=F3
FF[4]=F4

Threads.@threads for i=1:4*dim
    f=Int((i-1)÷dim+1)
    k=getk((i-1)%dim+1)
    z=getz((i-1)%dim+1)
    FF[f][k,z]=solution[i]
end

F1k=zeros(kstep)
F2k=zeros(kstep)
F3k=zeros(kstep)
F4k=zeros(kstep)
for i=1:kstep
    for j=1:zstep
        F1k[i]+=F1[i,j]*weightz[j]*2/pi
        # F2k[i]+=F2[i,j]*weightz[j]*2/pi
        # F3k[i]+=F3[i,j]*weightz[j]*2/pi
        # F4k[i]+=F4[i,j]*weightz[j]*2/pi
    end
end
# plot(meshk,F1k,scale=:log10)
using Gaussquad
using FastGaussQuadrature
using Plots

const mt = 0.5;
const τ = ℯ^2-1;
const Λ = 0.234;
const ω = 0.5;
const dd = (0.82)^3/ω;
const Nf = 4;
const rm = 12/(33 - 2*Nf);
const m = 0.003;


cutdown=10^-4
cutup=10. ^4
kstep=128

k2=10000.
q2=10000.

sum1= Vector{Float64}(undef, kstep);
sum2= Vector{Float64}(undef, kstep);
F(x)=(1-exp(-x/(4*mt)^2))/x;
D(t)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*F(t)/log(τ+(1+t/Λ^2)^2));


meshq,weightq=gausslegendremesh(cutdown,cutup,kstep,2);

z2step=640
z1step=64
meshz2,weightz2=gausslegendre(z2step)
meshz1,weightz1=gausschebyshev(z1step, 2)
meshzk,weightzk=gausschebyshev(kstep, 2)

Threads.@threads for i=1:length(meshzk)
zk=meshzk[i]
sqr=sqrt(k2*q2)
for z1i=1:z1step
    for z2i=1:z2step
        weight1=weightz1[z1i]
        weight2=weightz2[z2i]
        z1=meshz1[z1i]
        z2=meshz2[z2i]
        w=weight1*weight2
        sum2[i]+=w*D(k2+q2-2*sqr*(zk*z1+sqrt((1-zk^2)*(1-z1^2))*z2))
    end 
end
end #th for


using Gaussquad
using FastGaussQuadrature

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
kstep=32

k2=1.

sum1= Vector{Float64}(undef, kstep);
sum2= Vector{Float64}(undef, kstep);
F(x)=(1-exp(-x/(4*mt)^2))/x;
D(t)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*F(t)/log(τ+(1+t/Λ^2)^2));


meshq,weightq=gausslegendremesh(cutdown,cutup,kstep,2);

meshz2,weightz2=gausslegendre(16)
meshz1,weightz1=gausschebyshev(16, 2)

Threads.@threads for i=1:32
q2=meshq[i]
sqr=sqrt(k2*q2)
sum1[i]=2*gausschebyshevint64(z->D(k2+q2-2*sqr*z))
sum2[i]=0.
for z1i=1:16
    for z2i=1:16
        weight1=weightz1[z1i]
        weight2=weightz2[z2i]
        z1=meshz1[z1i]
        z2=meshz2[z2i]
        w=weight1*weight2
        sum2[i]+=w*(2*z1^2-1)*gausschebyshevint64(zk->D(k2+q2-2*sqr*(zk*z1+sqrt((1-zk^2)*(1-z1^2))*z2)))*2/pi
    end 
end
end #th for

div=[sum1[i]/sum2[i] for i=1:32]
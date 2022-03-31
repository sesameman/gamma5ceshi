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


kstep=512
z2step=512
z1step=32

k2=10000.
q2=10000.

sum1= zeros(kstep);
sum2= zeros(kstep);
function F(x)
    if x<10^-7
        a=1/4
    else
        a=((-expm1(-x/(4*mt)^2))/x)::Float64;
    end #if
    a
end
D(t)=8*pi^2*(dd*exp(-t/(ω^2))/ω^4+rm*F(t)/log(τ+(1+t/Λ^2)^2));


meshq,weightq=gausslegendremesh(cutdown,cutup,kstep,2);

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
        kmiq=k2+q2-2*sqr*(zk*z1+sqrt((1-zk^2)*(1-z1^2))*z2)
        if kmiq<0
            print("$z2&$z1\n")
        end
        sum2[i]+=w*D(kmiq)
        sum1[i]+=w/kmiq
    end 
end
end #th for


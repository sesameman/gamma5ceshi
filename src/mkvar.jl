Threads.@threads for i=1:dim
    k[i]=kfunction(i)
    z[i]=zfunction(i)
    qPlus2[i]=qPlus2function(i)
    qSubt2[i]=qSubt2function(i)
    kdotp[i]=kdotpfunction(i)
    pdotq[i]=pdotqfunction(i)
    A1[i]=A1function(i)
    B1[i]=B1function(i)
    A2[i]=A2function(i)
    B2[i]=B2function(i)
    branchplus[i]=branchfunction(qPlus2[i])
    branchsubt[i]=branchfunction(qSubt2[i])
    branch[i]=branchplus[i]*branchsubt[i]
end
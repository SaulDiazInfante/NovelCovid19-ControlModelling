function SerologicalShields!(du,u,par,t)
    β_a,β_s,p,γ_e,γ_a,γ_s,γ_h,h,ξ,μ,α,S_ind,E_ind,Ia_ind,Is_ind,Ihsub_ind,Ihcrit_ind,R_ind,D_ind = par

    TOTI_a = sum(u[Ia_ind]) #Total asymptomatic infectious individuals
    TOTI_s = sum(u[Is_ind]) #Total symptomatic infectious individuals

    ShieldRTOT = sum(u[R_ind[3:6]]) #Total serological shields - those who are recovered and 20-60 yo.
    NTOT = sum(u[1:40])+sum(u[61:70]) #All individuals not in the hospitalised (Ihsub,Ihcrit) or dead (D) states
    SHIELD = ShieldRTOT

    for aa = 1:length(S_ind)#for each age class
        INDS = [S_ind[aa],E_ind[aa],Ia_ind[aa],Is_ind[aa],Ihsub_ind[aa],Ihcrit_ind[aa],R_ind[aa],D_ind[aa]]
        S,E,Ia,Is,Ihsub,Ihcrit,R,D = u[INDS]
        du[INDS[1]] = -β_a*S*TOTI_a/(NTOT+α*SHIELD) - β_s*S*TOTI_s/(NTOT+α*SHIELD) #dS/dt
        du[INDS[2]] = +β_a*S*TOTI_a/(NTOT+α*SHIELD) + β_s*S*TOTI_s/(NTOT+α*SHIELD) -  γ_e*E   #dE/dt
        du[INDS[3]] = p[aa]*γ_e*E - γ_a*Ia  #dI_a /dt
        du[INDS[4]] = (1-p[aa])*γ_e*E - γ_s*Is  #dI_s/dt
        du[INDS[5]] = h[aa]*(1-ξ[aa])*γ_s*Is - γ_h*Ihsub #dIhsub/dt
        du[INDS[6]] = h[aa]*ξ[aa]*γ_s*Is -  γ_h*Ihcrit #dIhcrit/dt
        du[INDS[7]] = γ_a*Ia + (1-h[aa])*γ_s*Is + γ_h*Ihsub + (1-μ[aa])*γ_h*Ihcrit  #dR/dt
        du[INDS[8]] = μ[aa]*γ_h*Ihcrit #dD/dt
    end

end

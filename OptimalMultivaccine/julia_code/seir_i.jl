function seir_i!(du,u,par,t)
    β_a,
    β_s,
    p,
    γ_e,
    γ_a,
    γ_s,
    γ_h,
    h,
    ξ,
    μ,
    α,
    S_ind,
    E_ind,
    Ia_ind,
    Is_ind,
    H_ind,
    R_ind,
    M_ind = par

    TOTI_a = sum(u[Ia_ind]) #Total asymptomatic infectious individuals
    TOTI_s = sum(u[Is_ind]) #Total symptomatic infectious individuals
    NTOT = sum(u[1:12])+sum(u[21:24]) #All individuals not in the hospitalised (Ihsub,Ihcrit) or dead (D) states

    for aa = 1:length(S_ind)#for each age class
        INDS = [
            S_ind[aa],
            E_ind[aa],
            Ia_ind[aa],
            Is_ind[aa],
            H_ind[aa],
            R_ind[aa],
            M_ind[aa]
        ]
        S, E, Ia, Is, H, R, M = u[INDS]
        du[INDS[1]] = -β_a * S
        du[INDS[2]] = β_a*S -  γ_e * E
        du[INDS[3]] = (1-p[aa]) * γ_e * E- γ_a * Ia
        du[INDS[4]] = p[aa] * γ_e * E - γ_s * Is
        du[INDS[5]] = h[aa]*ξ[aa]*γ_s*Is -  γ_h*Ihcrit
        du[INDS[6]] = γ_a*Ia + (1-h[aa])*γ_s*Is + γ_h*Ihsub + (1-μ[aa])*γ_h*Ihcrit  #dR/dt
        du[INDS[7]] = μ[aa]*γ_h*Ihcrit #dD/dt
    end
  end

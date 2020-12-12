function compute_a_0_i(prm)
    r_not_i = prm["r_not_i"]
    p_i = prm["p_i"]
    gamma_i = prm["gamma_i"]
    gamma_a_i = prm["gamma_a_i"]
    theta_i = prm["theta_i"]
    a_i = r_not_i ./ (p_i ./ gamma_i + theta_i .* (1 .- p_i) ./ gamma_a_i)
    return a_i
end

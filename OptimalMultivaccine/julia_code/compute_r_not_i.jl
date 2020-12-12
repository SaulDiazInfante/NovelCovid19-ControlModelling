function compute_r_not_i(prm)
    r_not_i = prm["r_not_i"]
    p_i = prm["p_i"]
    gamma_i = prm["gamma_i"]
    gamma_a_i = prm["gamma_a_i"]
    theta_i = prm["theta_i"]
    a_0_i = prm["a_0_i"]
    r_not_i = a_0_i .*  (p_i ./ gamma_i + theta_i .* (1 .- p_i) ./ gamma_a_i)
    return r_not_i
end

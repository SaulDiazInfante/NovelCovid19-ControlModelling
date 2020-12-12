function force_infection(u, par, t)
    T_0 = par["T_i"][1];
    T_1 = par["T_i"][2];
    T_2 = par["T_i"][3];
    T_end = par["T_i"][3];
    a_0 = par["a_0_i"]
    if t >= T_0 and t < T_1
        a_i_t = (1 -0.95 .* s_b) .* a_0
    elseif t >= T_1 and t <= T_2
        a_i_t = (1 -0.95 .* s_2) .* a_0
    elseif t >= T_2 and t <= T_end
        a_i_t = (1 -0.95 .* s_3) .* a_0
    else
        a_i_t = (1 -0.95 .* s_r) .* a_0
    return  a_i_t
end

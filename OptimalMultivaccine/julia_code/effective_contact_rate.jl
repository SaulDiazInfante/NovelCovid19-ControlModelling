import Pkg; Pkg.add("Pkg")
import JSON
function contact_rate!(t, parameters=prm_dict)
    T_0 = par["T_i"][1];
    T_1 = par["T_i"][2];
    T_2 = par["T_i"][3];
    T_end = par["T_i"][3];
    s_b = par["s_b"]
    s_r = par["s_r"]
    s_2_i = par["s_2_i"];
    s_3_i = par["s_3_i"];
    a_0 = par["a_0_i"]
    if t >= T_0 && t < T_1
        a_i_t = (1 .- 0.95 .* s_b) .* a_0
    elseif t >= T_1 && t <= T_2
        a_i_t = (1 .- 0.95 .* s_2_i) .* a_0
    elseif t >= T_2 && t <= T_end
        a_i_t = (1 .- 0.95 .* s_3_i) .* a_0
    else
        a_i_t = (1 .- 0.95 .* s_r) .* a_0
    end
    return  a_i_t
end
#=
t = collect(range(0.0, stop=500.0, length=500))
prm_dict = Dict()
# file_name = "parameters.json"
#file_name = "parameters_per_group.json"
file_name = "parameters_base_line_scene.json"
open(file_name, "r") do f
    global prm_dict
    dict_txt = read(f, String)
    prm_dict = JSON.parse(dict_txt)
end
a_i_t_ = contact_rate!.(t)
a_i_t_ =  permutedims(
    reshape(
        hcat(a_i_t_...),
        (length(a_i_t_[1]), length(a_i_t_))
    )
)
using Plots
plot(t, a_i_t_[:,1])
plot!(t, a_i_t_[:,2])
plot!(t, a_i_t_[:,3])
=#

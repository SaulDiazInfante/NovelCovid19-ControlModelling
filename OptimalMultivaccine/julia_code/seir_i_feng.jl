import Pkg; Pkg.add("Pkg")
import JSON
ENV["GRDIR"] = ""
using Pkg; Pkg.build("GR")
import Pkg; Pkg.add("DifferentialEquations")
using Pkg
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("StatsPlots")
Pkg.add("Gym")
using Plots, DifferentialEquations, LaTeXStrings, StatsPlots, Pkg, Gym
Pkg.status()
struct POPDEMOG
    AgeStructure #fraction of population in each age class
    PopulationSize #number of people in population
end
include("seir_i.jl")
include("compute_r_not_i.jl")
ExamplePOP = POPDEMOG(
    [0.12, 0.54, 0.34],
    10 * 10 ^ 6
)
MinAge = [40, 60, 80]
MiddleAge = [45, 65, 85]
plot(MiddleAge, ExamplePOP.AgeStructure * 100.0,
    xlabel="Age",
    ylims=(0,60),
    ylabel = "Pop.Demography (%)",
    legend = false,
    yguidefontsize=8,
    markershape = :auto,
    markersize = 2)
savefig("ExampleDemography.pdf")
current()
# Define state indexes
NumAgeClasses = 3;
NumEpiStates = 7;
S_ind = 1:NumAgeClasses
E_ind = (1 * NumAgeClasses + 1): 2 * NumAgeClasses
Ia_ind = (2 * NumAgeClasses + 1): 3 * NumAgeClasses
Is_ind = (3 * NumAgeClasses + 1): 4 * NumAgeClasses
H_ind = (4 * NumAgeClasses + 1): 5 * NumAgeClasses
R_ind = (5 * NumAgeClasses + 1): 6 * NumAgeClasses
M_ind = (6 * NumAgeClasses + 1): 7 * NumAgeClasses
# Load parameters
# TODO: Code T_j and implement this parameters
# in the signature of function a_i(t)
# TODO: implemetn the matrix s_ij
#
sol_init = solve(seir_i,
    Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
)
#
prm_dict = Dict()
# file_name = "parameters.json"
#file_name = "parameters_per_group.json"
file_name = "parameters_base_line_scene.json"
open(file_name, "r") do f
    global prm_dict
    dict_txt = read(f, String)
    prm_dict = JSON.parse(dict_txt)
end

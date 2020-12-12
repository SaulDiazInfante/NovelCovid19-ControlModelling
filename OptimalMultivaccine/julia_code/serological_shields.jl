import Pkg; Pkg.add("Pkg")
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
include("SerologicalShields.jl")
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
NumEpiStates = 8;
S_ind = 1:NumAgeClasses
E_ind = (1 * NumAgeClasses + 1): 2 * NumAgeClasses
Ia_ind = (2 * NumAgeClasses + 1): 3 * NumAgeClasses
Is_ind = (3 * NumAgeClasses + 1): 4 * NumAgeClasses
Ihsub_ind = (4 * NumAgeClasses + 1): 5 * NumAgeClasses
Ihcrit_ind = (5 * NumAgeClasses + 1):6 * NumAgeClasses
R_ind = (6 * NumAgeClasses + 1):7 * NumAgeClasses
D_ind = (7 * NumAgeClasses + 1):8 * NumAgeClasses
#Define model parameters

β_A = 0.4
β_S = 0.8
γ_e = (1.0/4.0)
γ_al = (1.0/4.0)
γ_sl = (1.0/4.0)
γ_ah = (1.0/6.0)
γ_sh = (1.0/6.0)
γ_h = (1.0/10.0)
p = [0.9, 0.8, 0.7]
# fraction of cases that are asymptomatic
h = [0.1, 0.3, 1.2] ./ 100;
# fraction of cases that require a hospital visit (by age class)
ξ = [5, 6.3, 12.2] ./ 100;
μ = 0.5 .* ones(NumAgeClasses)
#levels of serological shielding
α0 = 0.0
α1 = 2.0
α2 = 20.0
params0 = [
    β_A,
    β_S,
    p,
    γ_e,
    γ_al,
    γ_sl,
    γ_h,
    h,
    ξ,
    μ,
    α0,
    S_ind,
    E_ind,
    Ia_ind,
    Is_ind,
    Ihsub_ind,
    Ihcrit_ind,
    R_ind,
    D_ind
]
params2 = [
    β_A,
    β_S,
    p,
    γ_e,
    γ_al,
    γ_sl,
    γ_h,
    h,
    ξ,
    μ,
    α1,
    S_ind,
    E_ind,
    Ia_ind,
    Is_ind,
    Ihsub_ind,
    Ihcrit_ind,
    R_ind,
    D_ind
]
params20 = [β_A,
    β_S,
    p,
    γ_e,
    γ_al,
    γ_sl,
    γ_h,
    h,
    ξ,
    μ,
    α2,
    S_ind,
    E_ind,
    Ia_ind,
    Is_ind,
    Ihsub_ind,
    Ihcrit_ind,
    R_ind,
    D_ind
]
#
params0h = [
    β_A,
    β_S,
    p,
    γ_e,
    γ_ah,
    γ_sh,
    γ_h,
    h,
    ξ,
    μ,
    α0,
    S_ind,
    E_ind,
    Ia_ind,
    Is_ind,
    Ihsub_ind,
    Ihcrit_ind,
    R_ind,
    D_ind
]
#
params2h = [
    β_A,
    β_S,
    p,
    γ_e,
    γ_ah,
    γ_sh,
    γ_h,
    h,
    ξ,
    μ,
    α1,
    S_ind,
    E_ind,
    Ia_ind,
    Is_ind,
    Ihsub_ind,
    Ihcrit_ind,
    R_ind,
    D_ind
]
#
params20h = [
    β_A,
    β_S,
    p,
    γ_e,
    γ_ah,
    γ_sh,
    γ_h,
    h,
    ξ,
    μ,
    α2,
    S_ind,
    E_ind,
    Ia_ind,
    Is_ind,
    Ihsub_ind,
    Ihcrit_ind,
    R_ind,
    D_ind
]
POP = ExamplePOP #assign population statistics
# seed population with 500 symptomatic and 9500 asymptomatic cases
uss0 = zeros(NumAgeClasses * NumEpiStates) #vector of initial states
uss0[S_ind] = POP.AgeStructure .* POP.PopulationSize
#put everything in S
uss0[3] = uss0[3] - 1 #remove infected from S ages 20-29
uss0[13] = 1  #assign E cases ages 20-29
uss0 = uss0 ./ POP.PopulationSize;  #renormalise
#
#
# Low R_0 outbreak
#
tspan = (0.0, 400.0)
#run the simulations between 0 days and 400 days
#
#initial spin-up
prob_init = ODEProblem(
    SerologicalShields!,
    uss0, tspan, params0
)
#set up callback to stop integration when number of cases reaches a threshold
CaseThreshold = 10_000 / POP.PopulationSize
condition(u,t,integrator)  = sum(u[1:NumAgeClasses]) - (1 - CaseThreshold)
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)
#solve the initial spin up, checking for positivity and using callback definition
sol_init = solve(prob_init,
    Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
    saveat=1,
    isoutofdomain=(y,p,t)->any(x->x<0,y),
    callback = cb
)

#assign end spinup states as new initial condition.
uInit = sol_init.u[end];

#run scenarios
prob_alph0 = ODEProblem(SerologicalShields!,uInit,tspan,params0)
prob_alph2 = ODEProblem(SerologicalShields!,uInit,tspan,params2)
prob_alph20 = ODEProblem(SerologicalShields!,uInit,tspan,params20)

sol_alph0 = solve(prob_alph0,
    Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
    saveat=1,
    isoutofdomain=(y,p,t) -> any(x -> x < 0, y)
)
sol_alph2 = solve(prob_alph2,
    Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
    saveat=1,
    isoutofdomain=(y, p, t) -> any(x -> x < 0, y)
)
sol_alph20 = solve(
    prob_alph20,
    Tsit5(),
    abstol=1e-8,
    reltol=1e-8,
    saveat=1,
    isoutofdomain=(y,p,t) -> any(x -> x<0, y)
    );

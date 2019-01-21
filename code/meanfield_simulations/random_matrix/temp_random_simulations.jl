# This script runs the simulations to be approximated by the meanfeild approach.
# it should be run directly from the project directory (Bacteria_EcoFunc)

# #activate directory
using Pkg
Pkg.activate(".")

#Load simulation module
using EcoFuncSim; println("Loaded EcoFuncSim")
eco = EcoFuncSim

using Distributions; println("Loaded Distributions")
using LinearAlgebra; println("Loaded LinearAlgebra")
using DelimitedFiles; println("Loaded DelimitedFiles")

#include the helper functions
# include("functions.jl")
#Set up simulation parameters
nSp = 100 #number of species in each simulation
nSim = 10 #number of simulations per temperature
nTemp = 100 #number of temperatures to simuluate at
Tref = 285.0 #Reference temperature for all params
#pre-allocate the temperature vector and results array
T_vec = range(275.0,stop = 295.0,length = nTemp)
results = Array{Any,2}(undef,nTemp,nSim)

#Define distributions
#All temperature sensitvities are drawn from a normal distribution.
#U = uptake, R = Respiration , a = interactions
uEU = 0.32; uER = 0.65; uEa = 0.0 #means
sEU = 0.1 ; sER = 0.1 ; sEa = 0.1 # variences
U0  = 300.0 ; R0  = 100.0 ; a0  = 10.0/nSp; #B0

#export parameters
params = Dict(:uEU => uEU, :uER => uER, :uEa => uEa,
     :sEU => sEU, :sER => sER, :sEa => sEa,
     :U0 => U0, :R0 => R0, :a0 => a0)

writedlm("results/meanfield_simulations/random_matrix/temp_meanfield_sims/params.csv",params,',')

C0 = ones(nSp) #starting biomass
#Simulation Loop
for T = 1:nTemp
    #calculate Tdiff to simplify calculations
    Tdiff = (1/(eco.k * T_vec[T])) - (1/(eco.k * Tref))
    #define parameter distributions. they all folow the LogNormal
    dU = LogNormal(log(U0)- uEU * (Tdiff),abs(Tdiff) * sEU)
    dR = LogNormal(log(R0)- uER * (Tdiff),abs(Tdiff) * sER)
    da = LogNormal(log(a0)- uEa * (Tdiff),abs(Tdiff) * sEa)
    for sim = 1:nSim
        println("Temp: ",T_vec[T], " Sim: ",sim)
        #get actual rates
        U = rand(dU,nSp)
        R = rand(dR,nSp)
        a = rand(da,nSp,nSp)

        #create parameters
        P = make_params(U,R,a,intra = false)
        #do simulations
        sol = eco.simulate(C0,P,t_end = 5.0)
        # @assert all(mapslices(var,hcat(sol.u[end-10:end]...),dims=2) .< 1e-4) "Not reached equilibrium"
        dir = string("results/meanfield_simulations/random_matrix/temp_meanfield_sims/Temp_",string(T_vec[T]),"_sim_",string(sim),".csv")
        writedlm(dir,sol.u, ',')

        dir = string("results/meanfield_simulations/random_matrix/temp_meanfield_sims/Resp_Temp_",string(T_vec[T]),"_sim_",string(sim),".csv")
        writedlm(dir,R, ',')

    end
end

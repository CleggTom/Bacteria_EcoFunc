using DifferentialEquations, NamedTuples

#define boltzman equation
function boltz(B0,E,Tref,T)
    return(B0 * exp((-E/8.617e-5) * ((1/T)-(1/Tref)) ))
end

#define the parameters
p = @NT(p0 = 10.0, RA0 = 1.0, u0 = 100.0, RH0 = 10.0,
        Ep = 0.32, Era = 0.65, Eu = 0.65, Erh = 0.98,
        Tref = 293.15, T = 293.15,
        a_hh = 0.1, a_aa = 0.1, e_a = 0.3, e_h = 0.3,
        D_h = 0.1, D_a = 0.1)

#define the problem
function flux(du,u,p,t)
        du[1] = u[1] * ((p[:e_a] * boltz(p[:p0],p[:Ep],p[:Tref],p[:T])) -
                        (boltz(p[:RA0],p[:Era],p[:Tref],p[:T])) -
                        (p[:D_a]) -
                        (p[:a_aa] * u[1]))

        du[2] = u[2] * ((u[3] * p[:e_h] * boltz(p[:u0],p[:Eu],p[:Tref],p[:T])) -
                        (boltz(p[:RH0],p[:Erh],p[:Tref],p[:T])) -
                        (p[:D_h]) -
                        (p[:a_hh] * u[2]))

        du[3] = u[1] * ((p[:a_aa] * u[1]) + (p[:D_a])) +
                u[2] * ((p[:a_hh] * u[2]) + (p[:D_h]) -
               (u[3] * p[:e_h]*boltz(p[:u0],p[:Eu],p[:Tref],p[:T]) ) )
end

du = zeros(3)
u0 = [1.0,1.0,1.0]

#293 degrees
prob = ODEProblem(flux,u0,(0,1000.0),p)
sol = solve(prob)
result = hcat(sol.u...)'
writedlm("data/Sim_293.csv",result,',')

#275 degrees

#define the parameters
p = @NT(p0 = 10.0, RA0 = 1.0, u0 = 100.0, RH0 = 10.0,
        Ep = 0.32, Era = 0.65, Eu = 0.65, Erh = 0.98,
        Tref = 293.15, T = 275.15,
        a_hh = 0.1, a_aa = 0.1, e_a = 0.3, e_h = 0.3,
        D_h = 0.1, D_a = 0.1)

prob = ODEProblem(flux,u0,(0,1000.0),p)
sol = solve(prob)
result = hcat(sol.u...)'
writedlm("data/Sim_275.csv",result,',')

# 310
p = @NT(p0 = 10.0, RA0 = 1.0, u0 = 100.0, RH0 = 10.0,
        Ep = 0.32, Era = 0.65, Eu = 0.65, Erh = 0.98,
        Tref = 293.15, T = 310.15,
        a_hh = 0.1, a_aa = 0.1, e_a = 0.3, e_h = 0.3,
        D_h = 0.1, D_a = 0.1)

prob = ODEProblem(flux,u0,(0,5.0),p)
sol = solve(prob)
result = hcat(sol.u...)'
writedlm("data/Sim_310.csv",result,',')

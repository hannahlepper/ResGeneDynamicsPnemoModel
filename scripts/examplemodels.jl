using DifferentialEquations
using Plots

# A single-serotype model with rate of gain of resistance
function RGmod(du, u, p, t)
    
    #assign parameters to local names to make easier to read
    β, μ, c, τ, γ = p

    X = 1 - sum(u) #proportion of population uncolonised
    S, R = u

    du[1] = β*X*S - μ*S - τ*S - γ*S #changes to Sd
    du[2] = (1-c)*β*X*R - μ*R + γ*S #changes to Rd

end

β = 1.5
μ = 1/12
c = 0.3
τ = 0.01
γ = 0.01
p = [β, μ, c, τ, γ]
u0 = [0.1,0.1]
probRG = ODEProblem(RGmod, u0, (0.,10000.), p)
solRG = solve(probRG)
plot(solRG, xlab = "time in years", ylab = "prevalence", labels = ["S" "R"])

# A multi-serotype model which does not incorporate rate of gain of resistance
# To do: 
#keep same carriage duration for all serotypes, add different transmission rates for all serotypes
#Add rate of gain of resistance
function pdmod(du, u, p, t)
    
    #assign parameters to local names to make easier to read
    β, c, τ, a, D = p[1:5]
    D = Int64(D)
    mu = p[6:(5+D)]

    X = 1 - sum(u) #proportion of population uncolonised

    for d in 1:D #loop through serotypes
        
        Sd, Rd = [u[d], u[D+d]] #total individuals colonised by S and R versions of ST (easier to read)
        
        qd = (1 - (Sd+Rd)/(sum(u)) + 1/D)^a #NFDS term

        du[d] = qd*β*X*Sd - mu[d]*Sd - τ*Sd #changes to Sd
        du[D+d] = qd*(1-c)*β*X*Rd - mu[d]*Rd #changes to Rd

    end

end

#Carriage duration calculation
#To do: replace with a transmission rate calculation function
function carrdur(μ, δ, D)
    if (D == 1)
        return [μ]
    else
        mud = zeros(D)
        for i in 1:D
            mud[i] = μ * (1 + δ * (2 * ((i-1) / (D-1)) - 1) )
        end
        return mud
    end
end

#Parameter set up
β = 1.5 #transmission rate
c = 0.6 #resistance cost
τ = 0.1 #treatment rate
a = 2. #diversifying selection strength
D = 10 #number of serotypes

μaverage = 1/12
mus = carrdur(μaverage, 0.3, D) #calculate the carriage durations for each serotype
p = vcat([β, c, τ, a, D], mus) #group parameter values

#Initial conditions set up
u0 = zeros(D*2) .+ 0.01

#Model set up
prob = ODEProblem(pdmod, u0, (0.,100000.), p)

#Numerically solve the model and plot
sol = solve(prob)
plot(sol)
bar(sol(10000.))

#Now plot resistance frequency by serotype
Ss = sol(100000.)[1:D]
Rs = sol(100000.)[D+1:D*2]
resistancefrequency = Rs ./ (Ss .+ Rs)
bar(resistancefrequency, 
    xlab = "serotype", ylab = "resistance frequency",
    labels = nothing, xticks = [1:1:10;])
# run_gromet_demo.jl
# This file provides a demonstration of the `run_sim` interface

# Include necessary imports
using AlgebraicPetri
using JSON
include("src/ModelStratify.jl")
include("src/GrometInterop.jl")

# Run the simulation with the parameter json (can also pass Dict of parameters)
res = GrometInterop.run_sim("models/gromets/seird.json","models/parameters/seird.json")

# Save parameter results
open("gromet_res.json", "w") do f
  JSON.print(f, res, 2)
end


######################
# Functions as Rates #
######################

# Since we want to call a lower run_sim function, we have to hand-convert the
# gromet to a PetriNet

# Read in gromet json
pn = GrometInterop.gromet2petrinet("models/gromets/seird.json")

# Function that simulates social distancing
# if the infected population is greater than 10% of the uninfected population,
# enforce social distancing
function exposure_rate(u, t)
  if( u[:I] > 0.1*(u[:S]+u[:E]+u[:R]) )
    return 0.0001
  else
    return 0.001
  end
end

concs = Dict(:S=>10000, :E=>0, :I=>1, :R=>0, :D=>0)
rates = Dict{Symbol, Union{Real, Function}}(:inf=> 0.05, :exp => exposure_rate, :death=>0.001, :rec=>0.01)

res = GrometInterop.run_sim(pn, concs, rates, (0.0,20.0), collect(0.0:0.1:20.0))

open("gromet_sd_res.json", "w") do f
  JSON.print(f, res, 2)
end

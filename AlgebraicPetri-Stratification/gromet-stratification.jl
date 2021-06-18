using AlgebraicPetri
using JSON
using Catlab.Graphs.BasicGraphs
using Catlab.CategoricalAlgebra

include("src/ModelStratify.jl")
include("src/GrometInterop.jl")

# DEFINE EPI MODEL
# Read in gromet json
model = GrometInterop.gromet2petrinet("models/gromets/seird.json")
# Graph initial model
AlgebraicPetri.Graph(model)

# DEFINE CONNECTION GRAPH
# Read in json string
conn_json = open(f->read(f, String), "connection_graphs/scale_graph.json");
# Parse json to object
conn_graph = ModelStratify.deserialize(conn_json, ModelStratify.ScaleGraph{Number});
ModelStratify.show_graph(conn_graph)

# PERFORM DEMOGRAPHIC STRATIFICATION
demographic_model = ModelStratify.dem_strat(model, conn_graph, :S, :E, [:E,:I]);
# Display graph of model
AlgebraicPetri.Graph(demographic_model)
# Save svg and json files of model
open("demographic_gromet.json", "w") do f
   JSON.print(f, GrometInterop.petrinet2gromet(demographic_model, "SimpleDemoSIRD"), 2)
end

# PERFORM SPATIAL STRATIFICATION
spatial_model = ModelStratify.diff_strat(model, conn_graph);
# Display graph of model
AlgebraicPetri.Graph(spatial_model)
# Save svg and json files of model
open("spatial_gromet.json", "w") do f
  JSON.print(f, GrometInterop.petrinet2gromet(spatial_model, "SimpleSpatialSIRD"), 2)
end

using AlgebraicPetri
# using ArgParse
using Catlab.Graphs.BasicGraphs
using Catlab.CategoricalAlgebra
using JSON
using Genie
using Genie.Router
using Genie.Renderer.Json
 
include("ModelStratify.jl")

function load_model(file::String)
    model_json = open(f->read(f, String), file);
    model = ModelStratify.deserialize(model_json, LabelledPetriNet);
    model
end

function load_connection_graph(file::String)
    conn_json = open(f->read(f, String), file);
    connections = ModelStratify.deserialize(conn_json, BasicGraphs.Graph);
    connections
end

function load_states(file::String)
    states_json = open(f->read(f, String), file);
    states = JSON.parse(states_json);

    println(stderr, states)

    sus = Symbol(states["sus"]);
    exp = Symbol(states["exp"]);
    inf = [Symbol(i) for i = states["inf"]];

    sus, exp, inf
end

function demographic_stratification(topology::String, connections::String, states::String)
    println(stderr, "performing demographic stratification...")
    model = load_model(topology)
    connections = load_connection_graph(connections)
    sus, exp, inf = load_states(states)

    demographic_model = apex(ModelStratify.dem_strat(model, connections, sus, exp, inf));
    return ModelStratify.serializeToString(demographic_model);
end

function spatial_stratification(topology::String, connections::String)
    println(stderr, "performing spatial stratification...")
    model = load_model(topology)
    connections = load_connection_graph(connections)

    spatial_model = apex(ModelStratify.diff_strat(model, connections));
    return ModelStratify.serializeToString(spatial_model);
end

route("/") do 
    if haskey(@params, :type)
        if @params(:type) == "dem"
            Json.json(demographic_stratification(@params(:top), @params(:conn), @params(:states)))
        elseif @params(:type) == "spat"
            Json.json(spatial_stratification(@params(:top), @params(:conn)))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    up(8001, async=false)
end

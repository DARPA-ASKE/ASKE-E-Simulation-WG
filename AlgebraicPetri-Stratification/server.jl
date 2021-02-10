using AlgebraicPetri
using Catlab.Graphs.BasicGraphs
using Catlab.CategoricalAlgebra
using JSON
using Genie
using Genie.Router
using Genie.Requests
using Genie.Renderer.Json
 
include("ModelStratify.jl")

function load_model(model_json::Dict)
    model_str = JSON.json(model_json)
    model = ModelStratify.deserialize(model_str, LabelledPetriNet);
    model
end

function load_connection_graph(conn_json::Dict)
    conn_str = JSON.json(conn_json)
    connections = ModelStratify.deserialize(conn_str, BasicGraphs.Graph);
    connections
end

function load_states(states::Dict)
    sus = Symbol(states["sus"]);
    exp = Symbol(states["exp"]);
    inf = [Symbol(i) for i = states["inf"]];

    sus, exp, inf
end

function demographic_stratification(topology::Dict, connections::Dict, states::Dict)
    println(stderr, "performing demographic stratification...")
    model = load_model(topology)
    connections = load_connection_graph(connections)
    sus, exp, inf = load_states(states)

    demographic_model = apex(ModelStratify.dem_strat(model, connections, sus, exp, inf));
    return ModelStratify.serializeToString(demographic_model);
end

function spatial_stratification(topology::Dict, connections::Dict)
    println(stderr, "performing spatial stratification...")
    model = load_model(topology)
    connections = load_connection_graph(connections)

    spatial_model = apex(ModelStratify.diff_strat(model, connections));
    return ModelStratify.serializeToString(spatial_model);
end

route("/", method = POST) do 
    # @show jsonpayload()

    payload = jsonpayload()
    if haskey(payload, "type")
        if payload["type"] == "dem"
            Json.json(demographic_stratification(payload["top"], payload["conn"], payload["states"]))
        elseif payload["type"] == "spat"
            Json.json(spatial_stratification(payload["top"], payload["conn"]))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    up(8001, "0.0.0.0", async=false)
end

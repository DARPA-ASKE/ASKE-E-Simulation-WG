using AlgebraicPetri
using ArgParse
using Catlab.Graphs.BasicGraphs
using Catlab.CategoricalAlgebra
using JSON
 
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
    ModelStratify.serialize(demographic_model);
end

function spatial_stratification(topology::String, connections::String)
    println(stderr, "performing spatial stratification...")
    model = load_model(topology)
    connections = load_connection_graph(connections)

    spatial_model = apex(ModelStratify.diff_strat(model, connections));
    ModelStratify.serialize(spatial_model);
end

function parse_arguments()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--dem", "-d"
            help = "demographic stratification"
            action = :command
        "--spat", "-s"
            help = "spatial stratification"
            action = :command
        "--conn"
            help = "JSON file containing connection graph"
            required = true
        "--states"
            help =  """
                    JSON file containing extra state information:
                    {\"sus\": <str>, \"exp\": <str>, \"inf\": [<str>]}. 
                    Entries contain, respectively, symbol used for susceptible
                    states, symbol used for exposed states, and symbols
                    used for states that can expose susceptible states. REQUIRED
                    for demographic stratification
                    """
        "top"
            help = "JSON file containing original topology"
            required = true
    end

    return parse_args(s)
end

function main()
    args = parse_arguments()

    if haskey(args, "dem")
        demographic_stratification(args["top"], args["conn"], args["states"])
    elseif haskey(args, "spat")
        spatial_stratification(args["top"], args["conn"])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

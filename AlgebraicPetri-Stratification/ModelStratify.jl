module ModelStratify

# The available functions, all have necessary docstrings
export dem_strat, diff_strat, serialize, deserialize, save_petri, save_json, save_model

using AlgebraicPetri

using JSON

using Catlab
using Catlab.CategoricalAlgebra
using Catlab.Graphs.BasicGraphs
using Catlab.WiringDiagrams
using Catlab.Programs
using Catlab.Graphics
using Catlab.Graphics.Graphviz: run_graphviz

import Base.convert
Base.convert(::Type{Symbol}, str::String) = Symbol(str)

# Define helper functions for defining the two types of
# reactions in an epidemiology Model. Either a state
# spontaneously changes, or one state causes another to change

function cross_uwd(g::Catlab.Graphs.BasicGraphs.Graph)
    rel = RelationDiagram{Symbol}(0)

    # Add populations
    juncs = add_junctions!(rel, nv(g), variable=[Symbol("pop$i") for i in 1:nv(g)])
    srcs = subpart(g, :src)
    tgts = subpart(g, :tgt)
    # Add cross boxes
    add_parts!(rel, :Box, ne(g), name=[Symbol("cross_$(srcs[i])_$(tgts[i])") for i in 1:ne(g)])
    add_parts!(rel, :Port, ne(g), junction=srcs, box=1:ne(g))
    add_parts!(rel, :Port, ne(g), junction=tgts, box=1:ne(g))

    # Add epidemiology model boxes
    boxes = add_parts!(rel, :Box, nv(g), name=[Symbol("ep$i") for i in 1:nv(g)])
    add_parts!(rel, :Port, nv(g), junction=juncs, box=boxes)
    rel
end

function stratify(epi_petri::Function, connection_graph::Catlab.Graphs.BasicGraphs.Graph, diffusion_petri::Function)
    conn_uwd = cross_uwd(connection_graph)

    # Calls diffusion_petri for each edge as (src, tgt)
    ep_map = Dict{Symbol, OpenLabelledPetriNet}([Symbol("ep$i")=>epi_petri(i) for i in 1:nv(connection_graph)])
    srcs = subpart(connection_graph, :src)
    tgts = subpart(connection_graph, :tgt)
    for i in 1:ne(connection_graph)
        ep_map[Symbol("cross_$(srcs[i])_$(tgts[i])")] = diffusion_petri(srcs[i], tgts[i])
    end

    oapply(conn_uwd, ep_map)
end

dem_connection(epi_model::Function, sus_state::Symbol, exp_state::Symbol, inf_states::Array{Symbol}, x::Int, y::Int) = begin
    append_ind(x::Symbol, ind::Int) = Symbol("$(x)_$ind")

    ep1 = apex(epi_model(x))
    ep2 = apex(epi_model(y))

    sus1 = append_ind(sus_state, x)
    sus2 = append_ind(sus_state, y)
    exp1 = append_ind(exp_state, x)
    exp2 = append_ind(exp_state, y)
    inf1 = [append_ind(inf, x) for inf in inf_states]
    inf2 = [append_ind(inf, y) for inf in inf_states]

    LabelledPetriNet(vcat(subpart(ep1, :sname), subpart(ep2, :sname)),
                          [Symbol("crx_$(sus2)_$(inf)")=>((sus2, inf)=>(inf, exp2)) for inf in inf1]...)

end

diff_connection(epi_model::Function, x::Int, y::Int) = begin
    ep1 = apex(epi_model(x))
    ep2 = apex(epi_model(y))
    states1 = subpart(ep1, :sname)
    states2 = subpart(ep2, :sname)
    LabelledPetriNet(vcat(states1, states2),
                     [Symbol("diff_$(states1[i])_$(states2[i])")=>(states1[i]=>states2[i]) for i in 1:ns(ep1)]...)
end

add_index(epi_model::LabelledPetriNet, ind::Int) = begin
    new_petri = deepcopy(epi_model)
    snames = subpart(epi_model, :sname)
    tnames = subpart(epi_model, :tname)
    set_subpart!(new_petri, :sname, [Symbol("$(name)_$ind") for name in snames])
    set_subpart!(new_petri, :tname, [Symbol("$(name)_$ind") for name in tnames])
    new_petri
end



""" diff_strat(epi_model, connection_graph)

  This function takes in a LabelledPetriNet and a graph which describes
  geographical connections. It returns a LabelledPetriNet which models
  diffusion between geographic populations described by the given graph.
"""
function diff_strat(epi_model::LabelledPetriNet, connection_graph::Catlab.Graphs.BasicGraphs.Graph)
    epi_func(ind) = begin
        ind_epi = add_index(epi_model, ind)
        Open(ind_epi, subpart(ind_epi, :sname))
    end
    conn(x, y) = begin
        conn_epi = diff_connection(epi_func, x, y)
        Open(conn_epi, subpart(conn_epi, :sname)[1:ns(epi_model)],
                       subpart(conn_epi, :sname)[(ns(epi_model)+1):(2*ns(epi_model))])
    end
    stratify(epi_func, connection_graph, conn)
end

""" dem_strat(epi_model, connection_graph, sus_state, exp_state, inf_states)

  This function takes in a LabelledPetriNet and a graph which describes
  infection connections between populations. It also takes in the symbol used
  for susceptible states, the symbol used for the exposed state, and an array
  of symbols for states that can expose susceptible states. It returns a
  LabelledPetriNet which models diffusion between populations described by the
  given graph.
"""
function dem_strat(epi_model::LabelledPetriNet, connection_graph::Catlab.Graphs.BasicGraphs.Graph, sus_state::Symbol, exp_state::Symbol, inf_states::Array{Symbol})
    epi_func(ind) = begin
        ind_epi = add_index(epi_model, ind)
        Open(ind_epi, subpart(ind_epi, :sname))
    end
    conn(x, y) = begin
        conn_epi = dem_connection(epi_func, sus_state::Symbol,
                                  exp_state::Symbol, inf_states::Array{Symbol}, x, y)
        Open(conn_epi, subpart(conn_epi, :sname)[1:ns(epi_model)],
                       subpart(conn_epi, :sname)[(ns(epi_model)+1):(2*ns(epi_model))])
    end
    stratify(epi_func, connection_graph, conn)
end

""" Serialize an ACSet object to a JSON string
"""
serialize(x::ACSet; io=stdout) = JSON.print(io, x.tables)

""" Deserialize a dictionary from a parsed JSON string to an object of the given ACSet type
"""
function deserialize(input::Dict, type)
    out = type()
    for (k,v) ∈ input
        add_parts!(out, Symbol(k), length(v))
    end
    for l ∈ values(input)
        for (i, j) ∈ enumerate(l)
            for (k,v) ∈ j
                set_subpart!(out, i, Symbol(k), v)
            end
        end
    end
    out
end

""" Deserialize a JSON string to an object of the given ACSet type
"""
deserialize(input::String, type) = deserialize(JSON.parse(input), type)

""" Save Petri net as an svg image
"""
save_petri(g, fname::AbstractString, format::AbstractString) =
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, AlgebraicPetri.Graph(g), format=format)
    end

""" Save Graph as an svg image
"""
save_graph(g, fname::AbstractString, format::AbstractString) =
    open(string(fname, ".", format), "w") do io
        run_graphviz(io, g, format=format)
    end

""" Show Graph
"""
show_graph(g) = to_graphviz(g, node_labels=true)

""" Save serialization to json file
"""
save_json(C, fname::AbstractString) =
    open(string(fname, ".json"),"w") do f
        serialize(C; io=f)
    end

""" Save both Petri graph as svg and json file
"""
function save_model(model, fname::AbstractString)
    save_json(model, fname)
    save_petri(model, fname, "svg");
end


end

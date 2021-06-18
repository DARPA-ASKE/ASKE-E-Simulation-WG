module GrometInterop

using Catlab
using JSON: parsefile
using JSON
using Catlab.CategoricalAlgebra
using Catlab.Theories
using Catlab.Present
using LabelledArrays
using DifferentialEquations
using AlgebraicPetri

export gromet2petrinet, petrinet2gromet, run_sim

@present TheorySemagram(FreeSchema) begin
  (B,P,J,W,C,E,A)::Ob
  Value::Data

  (srcA,tgtA)::Hom(A, P)
  (srcC, tgtC)::Hom(C, B)
  (srcE, tgtE)::Hom(E, J)
  port::Hom(W,P)
  junction::Hom(W,J)
  box::Hom(P,B)

  bvalue::Attr(B, Value)
  pvalue::Attr(P, Value)
  jvalue::Attr(J, Value)
  wvalue::Attr(W, Value)
  cvalue::Attr(C, Value)
  evalue::Attr(E, Value)
  avalue::Attr(A, Value)
end

const AbstractSemagram = AbstractACSetType(TheorySemagram)
const Semagram = ACSetType(TheorySemagram, index=[:srcA, :tgtA,
                                            :srcC, :tgtC,
                                            :srcE, :tgtE,
                                            :port, :junction]){Any}

g2j_types = Dict{String, Type}("T:Integer"  => Int64,
                               "T:Number"   => Number,
                               "T:Real"     => Real,
                               "T:Float"    => Float64,
                               "T:Boolean"  => Bool,
                               "T:Nothing"  => Nothing,
                               "T:Character"=> Char,
                               "T:Symbol"   => Symbol)

j2g_types = Dict([g2j_types[k]=>k for k in keys(g2j_types)])

function add_uid!(uids, uid, obj::Tuple{Symbol, Int64})
  uid in keys(uids) && error("$uid is not unique in Gromet.")
  uids[uid] = obj
end

function gromet2semagram(g::Dict)
  s = Semagram()

  uids = Dict{String, Tuple{Symbol, Int64}}()
  # Import boxes and junctions
  for k in ["boxes", "junctions", "ports", "wires"]
    if isnothing(g[k])
      g[k] = []
    end
  end

  for b in g["boxes"]
    b_ind = add_part!(s, :B, bvalue=b)
    add_uid!(uids, b["uid"], (:B, b_ind))
  end

  for j in g["junctions"]
    j_ind = add_part!(s, :J, jvalue=j)
    add_uid!(uids, j["uid"], (:J, j_ind))
  end

  # Import ports
  for p in g["ports"]
    p_ind = add_part!(s, :P, pvalue=p, box=uids[p["box"]][2])
    add_uid!(uids, p["uid"], (:P, p_ind))
  end

  # Import wires
  for w in g["wires"]

    if "endpoints" in keys(w)
      s_uid, t_uid = w["endpoints"]
    else
      s_uid, t_uid = (w["src"], w["tgt"])
    end
    src, tgt = (uids[s_uid], uids[t_uid])

    if src[1] == :P
      if tgt[1] == :P
        w_ind = add_part!(s, :A, srcA=src[2], tgtA=tgt[2], avalue=w)
        add_uid!(uids, w["uid"], (:A, w_ind))
      elseif tgt[1] == :J
        w_ind = add_part!(s, :W, port=src[2], junction=tgt[2], wvalue=w)
        add_uid!(uids, w["uid"], (:W, w_ind))
      else
        error("Invalid gromet (wire from $(s_uid) to $(t_uid))")
      end
    elseif src[1] == :J
      if tgt[1] == :P
        w_ind = add_part!(s, :W, junction=src[2], port=tgt[2], wvalue=w)
        add_uid!(uids, w["uid"], (:W, w_ind))
      elseif tgt[1] == :J
        w_ind = add_part!(s, :E, srcE=src[2], tgtE=tgt[2], evalue=w)
        add_uid!(uids, w["uid"], (:E, w_ind))
      else
        error("Invalid gromet (wire from $(s_uid) to $(t_uid))")
      end
    else
      error("Invalid gromet (wire from $(s_uid) to $(t_uid))")
    end
  end

  s
end

function gromet2semagram(gfile::String)
  g = parsefile(gfile)
  gromet2semagram(g)
end

function semagram2petrinet(sg)
  if nparts(sg, :J) != 0
    md = sg[1, :jvalue]
    local pn
    labelled = !all(isnothing.([md["name"] for md in sg[:jvalue]]))
    valued = !all(isnothing.([md["value_type"] for md in sg[:jvalue]]))

    # Determine which type of Petrinet best fits
    if valued

      # We assume that all values of States are the same type, and all values of Rate are the same type
      r_type = g2j_types[sg[first(filter(i->sg[i, :jvalue]["type"] == "Rate", 1:nparts(sg, :J))), :jvalue]["value_type"]]
      c_type = g2j_types[sg[first(filter(i->sg[i, :jvalue]["type"] == "State", 1:nparts(sg, :J))), :jvalue]["value_type"]]
      if labelled
        pn = LabelledReactionNet{r_type, c_type}()
      else
        pn = ReactionNet{r_type, c_type}()
      end
    else
      if labelled
        pn = LabelledPetriNet()
      else
        pn = PetriNet()
      end
    end

    sg2pn = Array{Pair{Symbol, Int64}, 1}(undef, nparts(sg, :J))
    # Import junction data
    for i in 1:nparts(sg, :J)
      j = sg[i, :jvalue]
      if j["type"] == "State"
        j_ind = add_part!(pn, :S)
        sg2pn[i] = :S=>j_ind
        if labelled
          pn[j_ind, :sname] = Symbol(j["name"])
        end
        if valued
          pn[j_ind, :concentration] = isnothing(j["value"]) ? 0.0 : j["value"]
        end
      elseif j["type"] == "Rate"
        j_ind = add_part!(pn, :T)
        sg2pn[i] = :T=>j_ind
        if labelled
          pn[j_ind, :tname] = Symbol(j["name"])
        end
        if valued
          pn[j_ind, :rate] = isnothing(j["value"]) ? 0.0 : j["value"]
        end
      else
        error("$(j["type"]) is an invalid type for PetriNet junctions")
      end
    end
    # Import wire data
    for e in 1:nparts(sg, :E)
      src = sg2pn[sg[e, :srcE]]
      tgt = sg2pn[sg[e, :tgtE]]
      if src[1] == :S
        tgt[1] == :T || error("Invalid conection from state to state in wire $e")
        add_part!(pn, :I, is=src[2], it=tgt[2])
      else
        tgt[1] == :S || error("Invalid conection from transition to transition in wire $e")
        add_part!(pn, :O, ot=src[2], os=tgt[2])
      end
    end
    pn
  else
    PetriNet()
  end
end

# This function assumes that all of the relevant metadata is included in the
# `Value` attribute for each object
function semagram2gromet(sg, type::String, name=nothing)
  # Gromet with single box as root
  root_uid = "B:$name"
  gromet = Dict{String, Any}("syntax"=>"Gromet",
                "type" => type,
                "name" => name,
                "metadata" => nothing,
                "uid" => "$(name)$type",
                "root" => root_uid,
                "types"=> nothing,
                "literals" => nothing,
                "variables" => nothing)

  # Add junctions, wires, boxes, and ports
  gromet["junctions"] = isempty(sg[:jvalue]) ? nothing : sg[:jvalue]
  gromet["ports"] = isempty(sg[:pvalue]) ? nothing : sg[:pvalue]
  gromet["boxes"] = isempty(sg[:bvalue]) ? nothing : sg[:bvalue]
  wires = vcat(sg[:cvalue], sg[:evalue], sg[:avalue], sg[:wvalue])
  gromet["wires"] = isempty(wires) ? nothing : wires
  gromet
end

function make_uid!(uids::Set{String}, val::String)
  if val in uids
    pcs = split(val, "_")
    ind = pcs[end]
    uid = join(pcs[1:(end-1)], "_")

    cur_ind = tryparse(Int64,ind)
    if isnothing(cur_ind)
      cur_ind = 0
      uid = val
    end

    while "$(uid)_$cur_ind" in uids
      cur_ind += 1
      cur_ind < 1e4 || error("UID generation for $uid has reached index $cur_ind. Current indexing scheme is insufficient.")
    end
    uid_res = "$(uid)_$cur_ind"
    push!(uids, uid_res)
    uid_res
  else
    push!(uids, val)
    val
  end
end

function petrinet2semagram(pn::AbstractPetriNet, name::String)
  sg = Semagram()
  s2j = zeros(Int64, ns(pn))
  t2j = zeros(Int64, nt(pn))

  labelled = pn isa LabelledPetriNet || pn isa LabelledReactionNet
  valued = pn isa ReactionNet || pn isa LabelledReactionNet

  uids = Set{String}()
  for s in 1:ns(pn)
    s_uid = make_uid!(uids, "J:$(sname(pn, s))")
    s_val = Dict("syntax"=>"Junction",
                 "type"=>"State",
                 "name"=>labelled ? sname(pn, s) : nothing,
                 "metadata"=>nothing,
                 "value"=>nothing, # valued ? concentration(pn, s) : nothing,
                 "value_type"=>valued ? j2g_types[typeof(concentration(pn,s))] : nothing,
                 "uid"=>s_uid )
    s2j[s] = add_part!(sg, :J, jvalue=s_val)
  end

  for t in 1:nt(pn)
    t_uid = make_uid!(uids, "J:$(tname(pn, t))")
    t_val = Dict("syntax"=>"Junction",
                 "type"=>"Rate",
                 "name"=>labelled ? tname(pn, t) : nothing,
                 "metadata"=>nothing,
                 "value"=>nothing, #valued ? rate(pn, t) : nothing,
                 "value_type"=>valued ? j2g_types[typeof(rate(pn, t))] : nothing,
                 "uid"=>t_uid)
    t2j[t] = add_part!(sg, :J, jvalue=t_val)
  end

  for i in 1:ni(pn)
    w_uid = make_uid!(uids, "W:$(sname(pn, pn[i,:is])).$(tname(pn, pn[i, :it]))")
    w_src = sg[s2j[pn[i,:is]], :jvalue]["uid"]
    w_tgt = sg[t2j[pn[i,:it]], :jvalue]["uid"]
    w_val = Dict("syntax"=>"Wire",
                 "type"=>nothing,
                 "name"=>nothing,
                 "metadata"=>nothing,
                 "value"=>nothing,
                 "value_type"=>nothing,
                 "uid"=>w_uid,
                 "src"=>w_src,
                 "tgt"=>w_tgt)
    add_part!(sg, :E, srcE=s2j[pn[i,:is]], tgtE=t2j[pn[i,:it]], evalue=w_val)
  end

  for o in 1:no(pn)
    w_uid = make_uid!(uids, "W:$(tname(pn, pn[o,:ot])).$(sname(pn, pn[o, :os]))")
    w_src = sg[t2j[pn[o,:ot]], :jvalue]["uid"]
    w_tgt = sg[s2j[pn[o,:os]], :jvalue]["uid"]
    w_val = Dict("syntax"=>"Wire",
                 "type"=>nothing,
                 "name"=>nothing,
                 "metadata"=>nothing,
                 "value"=>nothing,
                 "value_type"=>nothing,
                 "uid"=>w_uid,
                 "src"=>w_src,
                 "tgt"=>w_tgt)
    add_part!(sg, :E, srcE=t2j[pn[o,:ot]], tgtE=s2j[pn[o,:os]], evalue=w_val)
  end

  wires = sg[:evalue]
  junctions = sg[:jvalue]
  b_val = Dict("wires"=> isempty(wires) ? nothing : [w["uid"] for w in wires],
               "boxes"=>nothing,
               "junctions"=> isempty(junctions) ? nothing : [j["uid"] for j in junctions],
               "syntax"=>"Relation",
               "type"=>nothing,
               "name"=>"$name",
               "metadata"=>nothing,
               "uid"=>"B:$name",
               "ports"=> nothing)
  add_part!(sg, :B, bvalue=b_val)
  sg
end

""" petrinet2gromet(pn::AbstractPetriNet, name::String)
  This function takes in any PetriNet and generates a gromet PNC object
"""
function petrinet2gromet(pn::AbstractPetriNet, name::String)
  semagram2gromet(petrinet2semagram(pn, name), "PetriNetClassic", name)
end

""" gromet2petrinet(gromet)
This function takes in a gromet PNC object (or filename) and generates the
appropriate flavor of PetriNet.
"""
function gromet2petrinet(gromet)
  semagram2petrinet(gromet2semagram(gromet))
end

function test()
  println("Loading Classic")
  gromet2semagram("SimpleSIR_gromet_PetriNetClassic.json")
  println("Loading Bilayer")
  gromet2semagram("SimpleSIR_gromet_Bilayer.json")
  println("Loading Diet")
  gromet2semagram("SimpleSIR_gromet_PrTNet.json")
  println("Loading FN")
  gromet2semagram("SimpleSIR_gromet_FunctionNetwork.json")

  for f in ["SimpleSIR_gromet_PetriNetClassic.json", "sir_no_labels_rates.json", "sir_no_labels_no_rates.json", "sir_with_rates.json", "sird.json", "seird.json", "seir.json"]
    println("$f to petrinet")
    sg = gromet2semagram(f)
    pn = semagram2petrinet(sg)
    println("Back to gromet")
    gr = petrinet2gromet(pn, "SimpleSIR")
  end
end

""" run_sim(pn, concs_d, rates_d, t_range, tsteps)
This function takes in a LabelledReactionNet, initial concentrations, reaction
rates, and the relevant time information, simulates, the system, and returns
a results json file.
"""
function run_sim(pn::LabelledReactionNet{R,C},
                 concs_d::Dict{Symbol, C}, rates_d::Dict{Symbol, <:Union{R, Function}},
                 t_range::Tuple{<:Real,<:Real}, tsteps::Array{<:Real,1}) where {R,C}
  sim = vectorfield(pn)
  concs = @LArray collect(values(concs_d)) Tuple(keys(concs_d))
  prob = ODEProblem(sim, concs, t_range, rates_d)

  # Add error handling around `sol` function to take advantage of Galois `status` field

  sol = solve(prob, Tsit5())
  res = Dict("times"=>tsteps, "values"=>Dict([s=>[sol(t)[s] for t in tsteps] for s in snames(pn)]))
  return Dict("result"=>res)
end

function run_sim(pn::LabelledReactionNet{R,C}, parameters,
                 start::Number, stop::Number, step::Number) where {R,C}
  rates = Dict{Symbol, R}()
  concs = Dict{Symbol, C}()

  for p in keys(parameters)
    if Symbol(p) ∈ snames(pn)
      concs[Symbol(p)] = parameters[p]
    elseif Symbol(p) ∈ tnames(pn)
      rates[Symbol(p)] = parameters[p]
    else
      error("$p not in States or Rates of Petri net")
    end
  end
  run_sim(pn, concs, rates, (start, stop), collect(range(start, stop, step=step)))
end

function run_sim(gromet::Dict, p::Dict)
  run_sim(gromet2petrinet(gromet), p["parameters"], p["start"], p["stop"], p["step"])
end

""" run_sim(gfile, param_file)
This file accepts two json files, the first of which contains a gromet PNC and
the second of which contains a parameter file. A dictionary with simulation
results is returned from this function.

TODO: Add options for more advanced time-stepping methods
"""
function run_sim(gfile::String, param_file::String)
  run_sim(parsefile(gfile), parsefile(param_file))
end

end

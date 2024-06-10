export TheoryCausalLoop, AbstractCausalLoop, CausalLoopUntyped, CausalLoop, CausalLoopF,
nn, ne, vname,
sedge, tedge, convertToCausalLoop, nnames, CausalLoopF, epol, epols,
Polarity, POL_ZERO, POL_REINFORCING, POL_BALANCING, POL_UNKNOWN, POL_NOT_WELL_DEFINED,
add_node!, add_nodes!, add_edge!, add_edges!, discard_zero_pol,
outgoing_edges, incoming_edges, extract_loops, is_walk, is_circuit, walk_polarity, cl_cycles


using MLStyle

# using DataMigrations

import Graphs: SimpleDiGraph, simplecycles, SimpleEdge

import Catlab.Graphs: nv


import Base: +, *, -, /

# Whenever we map from something with attributes to something without, can think
# of it like mapping the attribute to a type with a single instance, which should
# be equivalent to identity

# Need some more work on the theory to ensure it's functorial, though.

@present TheoryCausalLoopNameless(FreeSchema) begin
  
  V::Ob
  P::Ob
  M::Ob

  sp::Hom(P, V)
  tp::Hom(P, V)

  sm::Hom(M, V)
  tm::Hom(M, V)

end


@present TheoryCausalLoopPM <: TheoryCausalLoopNameless begin
  
  Name::AttrType
  vname::Attr(V, Name)

end


@present TheoryCausalLoopZ <: TheoryCausalLoopPM begin
  Z::Ob

  sz::Hom(Z, V)
  tz::Hom(Z, V)
end

@present TheoryCausalLoopFull <: TheoryCausalLoopZ begin
  NWD::Ob
  U::Ob

  snwd::Hom(NWD, V)
  tnwd::Hom(NWD, V)

  su::Hom(U, V)
  tu::Hom(U, V)
end


@abstract_acset_type AbstractCausalLoop
@acset_type CausalLoopNamelessUntyped(TheoryCausalLoopNameless, index=[:sp,:tp, :sn, :tn]) <: AbstractCausalLoop
@acset_type CausalLoopPMUntyped(TheoryCausalLoopPM, index=[:sp,:tp, :sn, :tn]) <: AbstractCausalLoop
@acset_type CausalLoopZUntyped(TheoryCausalLoopZ, index=[:sp,:tp, :sn, :tn, :sz, :tz]) <: AbstractCausalLoop
@acset_type CausalLoopFullUntyped(TheoryCausalLoopFull, index = [:sp,:tp, :sn, :tn, :sz, :tz]) <: AbstractCausalLoop

const CausalLoopNameless = CausalLoopNamelessUntyped
const CausalLoopPM = CausalLoopPMUntyped{Symbol}
const CausalLoopZ = CausalLoopZUntyped{Symbol}
const CausalLoopFull = CausalLoopFullUntyped{Symbol}








# const CausalLoop = CausalLoopUntyped{Symbol} 
# const CausalLoopF = CausalLoopFUntyped{Symbol, Polarity}


@enum Polarity begin
  POL_ZERO
  POL_REINFORCING
  POL_BALANCING
  POL_UNKNOWN
  POL_NOT_WELL_DEFINED
end


# @present TheoryCausalLoop(FreeSchema) begin
#   E::Ob
#   V::Ob

#   s::Hom(E,V)
#   t::Hom(E,V)

#   # Attributes:
#   Name::AttrType
  
#   vname::Attr(V, Name)
# end

@present TheoryCausalLoop <: SchGraph begin
  # E::Ob
  # V::Ob

  # s::Hom(E,V)
  # t::Hom(E,V)

  # Attributes:
  Name::AttrType
  
  vname::Attr(V, Name)
end


@present TheoryCausalLoopPol <: TheoryCausalLoop begin
  Polarity::AttrType
  epol::Attr(E, Polarity)
end

# TODO: Make subtyping a bit more sensible
@abstract_acset_type AbstractSimpleCausalLoop
@acset_type CausalLoopUntyped(TheoryCausalLoop, index=[:src,:tgt]) <: AbstractSimpleCausalLoop
@acset_type CausalLoopPolUntyped(TheoryCausalLoopPol, index=[:src,:tgt]) <: AbstractSimpleCausalLoop
const CausalLoop = CausalLoopUntyped{Symbol}
const CausalLoopPol = CausalLoopPolUntyped{Symbol, Polarity}

function from_clp(cl::CausalLoopPol)
  pols = subpart(cl, :epol)
  names = Dict([i => x for (i, x) in enumerate(subpart(cl, :vname))])
  src = map(x -> names[x], subpart(cl, :src))
  tgt = map(x -> names[x], subpart(cl, :tgt))
  st = Vector{Pair{Symbol, Symbol}}(map(((x,y),) -> x => y, zip(src, tgt)))
  CausalLoopF(subpart(cl, :vname), st, pols)
end


function to_clp(nodes::Vector{Symbol}, reinf::Vector{Pair{Int, Int}}, 
  bal::Vector{Pair{Int, Int}}, zero::Vector{Pair{Int, Int}}, unknown::Vector{Pair{Int, Int}},
  nwd::Vector{Pair{Int, Int}})

  ne = length(reinf) + length(bal) + length(zero) + length(unknown) + length(nwd)
  pols = vcat(
    repeat([POL_REINFORCING], length(reinf)),
    repeat([POL_BALANCING], length(bal)),
    repeat([POL_ZERO], length(zero)),
    repeat([POL_UNKNOWN], length(unknown)),
    repeat([POL_NOT_WELL_DEFINED], length(nwd)),
  )

  edges = vcat(reinf, bal, zero, unknown, nwd)
  src = map(first, edges)
  tgt = map(last, edges)

  clp = CausalLoopPol()
  add_parts!(clp, :V, length(nodes); vname=nodes)
  add_parts!(clp, :E, ne; src = src, tgt = tgt, epol = pols)

  clp

end

function to_clp(cl::CausalLoopPM)
  to_clp(
    Vector{Symbol}(subpart(cl, :vname)), 
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sp), subpart(cl, :tp)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sm), subpart(cl, :tm)))),
    Vector{Pair{Int,Int}}(),
    Vector{Pair{Int,Int}}(),
    Vector{Pair{Int,Int}}()
    )
end

function to_clp(cl::CausalLoopZ)
  to_clp(
    Vector{Symbol}(subpart(cl, :vname)), 
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sp), subpart(cl, :tp)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sm), subpart(cl, :tm)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sz), subpart(cl, :tz)))),
    Vector{Pair{Int,Int}}(),
    Vector{Pair{Int,Int}}()
    )
end

function to_clp(cl::CausalLoopFull)
  to_clp(
    Vector{Symbol}(subpart(cl, :vname)), 
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sp), subpart(cl, :tp)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sm), subpart(cl, :tm)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :sz), subpart(cl, :tz)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :su), subpart(cl, :tu)))),
    Vector{Pair{Int,Int}}(map(((x,y),) -> x => y, zip(subpart(cl, :snwd), subpart(cl, :tnwd)))),
    )
end

# function to_clp(cl::CausalLoopZ)
#   clp = CausalLoopPol()

#   p = np(cl)
#   m = nm(cl)
#   z = nz(cl)

#   add_parts!(clp, :V, nparts(cl, :V))
#   add_parts!(clp, :E, p; s = subpart(cl, :sp), t = subpart(cl, :tp), epol=repeat([POL_REINFORCING], p))
#   add_parts!(clp, :E, m; s = subpart(cl, :sm), t = subpart(cl, :tm), epol=repeat([POL_BALANCING], m))
#   add_parts!(clp, :Z, m; s = subpart(cl, :sz), t = subpart(cl, :tz), epol=repeat([POL_ZERO], z))

# end

# function to_clp(cl::CausalLoopFull)
#   clp = CausalLoopPol()

#   p = np(cl)
#   m = nm(cl)
#   z = nz(cl)
#   nwd = nnwd(cl)
#   u = nu(cl)


#   add_parts!(clp, :V, nparts(cl, :V))
#   add_parts!(clp, :E, p; s = subpart(cl, :sp), t = subpart(cl, :tp), epol=repeat([POL_REINFORCING], p))
#   add_parts!(clp, :E, m; s = subpart(cl, :sm), t = subpart(cl, :tm), epol=repeat([POL_BALANCING], m))
#   add_parts!(clp, :E, m; s = subpart(cl, :sz), t = subpart(cl, :tz), epol=repeat([POL_ZERO], z))
#   add_parts!(clp, :E, m; s = subpart(cl, :snwd), t = subpart(cl, :tnwd), epol=repeat([POL_NOT_WELL_DEFINED], nwd))
#   add_parts!(clp, :E, m; s = subpart(cl, :su), t = subpart(cl, :tu), epol=repeat([POL_UNKNOWN], u))
# end


  
function *(p1::Polarity, p2::Polarity)
  if (p1 == POL_ZERO || p2 == POL_ZERO) return POL_ZERO end
  if (p1 == POL_UNKNOWN || p2 == POL_UNKNOWN) return POL_UNKNOWN end
  if (p1 == POL_NOT_WELL_DEFINED || p2 == POL_NOT_WELL_DEFINED) return POL_NOT_WELL_DEFINED end
  if ((p1 == POL_BALANCING && p2 == POL_BALANCING) || (p1 == POL_REINFORCING && p2 == POL_REINFORCING)) return POL_REINFORCING end
  return POL_BALANCING
end

function +(p1::Polarity, p2::Polarity)
  @match (p1, p2) begin
    (POL_ZERO, _) => p2
    (_, POL_ZERO) => p1

    (POL_UNKNOWN, _) || (_, POL_UNKNOWN) => POL_UNKNOWN
    
    (POL_NOT_WELL_DEFINED, _) || (_, POL_NOT_WELL_DEFINED) => POL_NOT_WELL_DEFINED
    (POL_REINFORCING, POL_BALANCING) || (POL_BALANCING, POL_REINFORCING) => POL_NOT_WELL_DEFINED

    (POL_REINFORCING, POL_REINFORCING) => POL_REINFORCING
    (POL_BALANCING, POL_BALANCING) => POL_BALANCING
  end
end

# -(p::Polarity) = p == POL_ZERO ? POL_ZERO : DomainError(p)

# function -(p1::Polarity, p2::Polarity)
#   @match (p1, p2) begin
#     (_, POL_ZERO) => p1
#     (POL_ZERO, POL_ZERO) => POL_ZERO
#     (POL_ZERO, _) => DomainError((p1, p2))

#     (POL_UNKNOWN, _) => POL_UNKNOWN

#     (POL_NOT_WELL_DEFINED, POL_UNKNOWN) => DomainError((p1, p2))
#     (POL_NOT_WELL_DEFINED, _) => POL_NOT_WELL_DEFINED

#     (POL_BALANCING, POL_BALANCING) => POL_BALANCING
#     (POL_BALANCING, _) => DomainError((p1, p2))

#     (POL_REINFORCING, POL_REINFORCING) => POL_REINFORCING
#     (POL_REINFORCING, _) => DomainError((p1, p2))

#   end
# end

# """p1 is a path of concatenated pols, p2 is what's being removed from the end """
# function /(p1::Polarity, p2::Polarity)
#   @match (p1, p2) begin
#     (POL_BALANCING, POL_BALANCING) || (POL_REINFORCING, POL_REINFORCING) => POL_REINFORCING
#     (POL_REINFORCING, POL_BALANCING) || (POL_BALANCING, POL_REINFORCING) => POL_BALANCING
#     (POL_ZERO, POL_ZERO) => POL_UNKNOWN # it's possible the previous value was the only zero
    
#     (_, POL_ZERO) => DivideError()

#     (POL_UNKNOWN, POL_UNKNOWN) => POL_UNKNOWN
#     (_, POL_UNKNOWN) => DivideError()

    
#     (POL_NOT_WELL_DEFINED, _) => POL_NOT_WELL_DEFINED
#     (_, POL_NOT_WELL_DEFINED) => DivideError()
#   end
# end




# @present TheoryCausalLoopF <: TheoryCausalLoop begin
#   Polarity::AttrType
#   epolarity::Attr(E, Polarity)
# end

# @abstract_acset_type AbstractCausalLoop
# @acset_type CausalLoopUntyped(TheoryCausalLoop, index=[:s,:t]) <: AbstractCausalLoop
# @acset_type CausalLoopFUntyped(TheoryCausalLoopF, index=[:s,:t]) <: AbstractCausalLoop
# const CausalLoop = CausalLoopUntyped{Symbol} 
# const CausalLoopF = CausalLoopFUntyped{Symbol, Polarity}

add_vertex!(c::AbstractCausalLoop;kw...) = add_part!(c,:V;kw...) 
add_vertices!(c::AbstractCausalLoop,n;kw...) = add_parts!(c,:V,n;kw...)

add_plus!(c::AbstractCausalLoop;kw) = add_part!(c, :P; kw...)
add_pluses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :P, n; kw...)

add_minus!(c::AbstractCausalLoop;kw) = add_part!(c, :M; kw...)
add_minuses!(c::AbstractCausalLoop,n;kw) = add_part!(c, :M, n; kw...)

# add_edge!(c::AbstractCausalLoop,s,t;kw...) = add_part!(c,:E,s=s,t=t;kw...) 
# add_edges!(c::AbstractCausalLoop,n,s,t;kw...) = add_parts!(c,:E,n,s=s,t=t;kw...)

"""
    CausalLoop(ns,es)
Create causal loop diagram from collection of nodes and collection of edges.
"""
# CausalLoop(ns,es) = begin
#     c = CausalLoop()
#     ns = vectorify(ns)
#     es = vectorify(es)
    
#     ns_idx=state_dict(ns)
#     add_nodes!(c, length(ns), nname=ns)

#     s=map(first,es)
#     t=map(last,es)
#     add_edges!(c, length(es), map(x->ns_idx[x], s), map(x->ns_idx[x], t))

#     c
# end



CausalLoopF() = CausalLoopPM()
CausalLoopF(ns::Vector{Symbol}, es::Vector{Pair{Symbol, Symbol}}, pols::Vector{Polarity}) = begin
  @assert length(pols) == length(es)

  if (POL_NOT_WELL_DEFINED ∈ pols || POL_UNKNOWN ∈ pols)
    c = CausalLoopFull()
  elseif (POL_ZERO ∈ pols)
    c = CausalLoopZ()
  else
    c = CausalLoopPM()
  end


  # c = CausalLoopF()
  ns = vectorify(ns)
  es = vectorify(es)
  
  ns_idx=state_dict(ns)
  add_vertices!(c, length(ns), vname=ns)

  s=map(first,es)
  t=map(last,es)


  for i in eachindex(pols)
    src = s[i]
    tgt = t[i]

    if pols[i] == POL_REINFORCING
      add_part!(c, :P; sp = ns_idx[src], tp = ns_idx[tgt])
    elseif pols[i] == POL_BALANCING
      add_part!(c, :M; sm =  ns_idx[src], tm = ns_idx[tgt])
    elseif pols[i] == POL_ZERO
      add_part!(c, :Z; sz =  ns_idx[src], tz =ns_idx[tgt])
    elseif pols[i] == POL_NOT_WELL_DEFINED
      add_part!(c, :NWD; snwd =  ns_idx[src], tnwd = ns_idx[tgt])
    elseif pols[i] == POL_UNKNOWN
      add_part!(c, :U; su =  ns_idx[src], tu = ns_idx[tgt])
    end
  end

  c
end



# # type piracy, hooray
# # return the count of each components
# """ return count of nodes of CLD """
# nv(c::AbstractCausalLoop) = 

# """ return count of edges of CLD """
# TODO: This should be unnecessary now, since it subtypes graph.  Just use ne and nv like normal.
nedges(c::AbstractSimpleCausalLoop) = nparts(c,:E) #edges
# nv(c::AbstractSimpleCausalLoop) = nparts(c,:V) #vertices




np(c::AbstractCausalLoop) = nparts(c, :P)
nm(c::AbstractCausalLoop) = nparts(c, :M)
nz(c::CausalLoopZ) = nparts(c, :Z)
nnwd(c::CausalLoopFull) = nparts(c, :NWD)
nu(c::CausalLoopFull) = nparts(c, :U)

nedges(c::AbstractCausalLoop) = np(c) + nm(c)
nedges(c::CausalLoopZ) = np(c) + nm(c) + nz(c)
nedges(c::CausalLoopFull) = np(c) + nm(c) + nz(c) + nnwd(c) + nu(c)





sedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :sp) : subpart(c, e - np(c), :sm)
end



tedge(c::CausalLoopPM, e) = begin
  @assert e <= nedges(c)
  e <= np(c) ? subpart(c, e, :tp) : subpart(c, e - np(c), :tm)
end



""" return node's name with index n """
# nname(c::AbstractCausalLoop,n) = subpart(c,n,:nname) # return the node's name with index of s
""" return edge's name with target number t """
sedge(c::CausalLoopPol,e) = subpart(c,e,:src)
""" return edge's name with edge number e """
tedge(c::CausalLoopPol,e) = subpart(c,e,:tgt)

""" return node names of CLD """
# nnames(c::AbstractCausalLoop) = [nname(c, n) for n in 1:nn(c)]

epol(c::CausalLoopPol,e) = subpart(c,e,:epol)

# epols(c::CausalLoopPol) = [epol(c, n) for n in 1:ne(c)]


outgoing_edges(c::CausalLoopPol, n) = collect(filter(i -> sedge(c,i) == n, 1:ne(c)))
incoming_edges(c::CausalLoopPol, n) = collect(filter(i -> tedge(c,i) == n, 1:ne(c)))



# function convertToCausalLoop(p::AbstractStockAndFlowStructure)
    
#     sns=snames(p)
#     fns=fnames(p)
#     svns=svnames(p)
#     flowVariableIndexs=[flowVariableIndex(p,f) for f in 1:nf(p)]
#     vNotf=setdiff(1:nvb(p),flowVariableIndexs)
#     vNotfns=[vname(p,v) for v in vNotf]
    
#     ns=vcat(sns,fns,svns,vNotfns)

#     lses=[sname(p,subpart(p,ls,:lss))=>svname(p,subpart(p,ls,:lssv)) for ls in 1:nls(p)]
#     lsvfes=[svname(p,subpart(p,lsv,:lsvsv))=>subpart(p,lsv,:lsvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lsv,:lsvv),:fv))) : vname(p,subpart(p,lsv,:lsvv)) for lsv in 1:nlsv(p)]
#     lfves=[sname(p,subpart(p,lv,:lvs))=>subpart(p,lv,:lvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvv),:fv))) : vname(p,subpart(p,lv,:lvv)) for lv in 1:nlv(p)]
#     fies=[fname(p,subpart(p,i,:ifn))=>sname(p,subpart(p,i,:is)) for i in 1:ni(p)]
#     foes=[fname(p,subpart(p,o,:ofn))=>sname(p,subpart(p,o,:os)) for o in 1:no(p)]

#     es=vcat(lses,lsvfes,lfves,fies,foes)

#     return CausalLoop(ns,es)
# end

"""
Convert StockFlow to CLD.
Nodes: stocks, flows, sum variables, parameters, nonflow dynamic variables
Edges: morphisms in stock flow
"""
# function convertToCausalLoop(p::AbstractStockAndFlowStructureF)
    
#     sns=snames(p)
#     fns=fnames(p)
#     svns=svnames(p)
#     pns=pnames(p)
#     flowVariableIndexs=[flowVariableIndex(p,f) for f in 1:nf(p)]
#     vNotf=setdiff(1:nvb(p),flowVariableIndexs)
#     vNotfns=[vname(p,v) for v in vNotf]
    
#     ns=vcat(sns,fns,svns,vNotfns,pns)

#     lses=[sname(p,subpart(p,ls,:lss))=>svname(p,subpart(p,ls,:lssv)) for ls in 1:nls(p)]
#     lsvfes=[svname(p,subpart(p,lsv,:lsvsv))=>subpart(p,lsv,:lsvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lsv,:lsvv),:fv))) : vname(p,subpart(p,lsv,:lsvv)) for lsv in 1:nlsv(p)]
#     lfves=[sname(p,subpart(p,lv,:lvs))=>subpart(p,lv,:lvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvv),:fv))) : vname(p,subpart(p,lv,:lvv)) for lv in 1:nlv(p)]
#     fies=[fname(p,subpart(p,i,:ifn))=>sname(p,subpart(p,i,:is)) for i in 1:ni(p)]
#     foes=[fname(p,subpart(p,o,:ofn))=>sname(p,subpart(p,o,:os)) for o in 1:no(p)]
#     lpvs=[pname(p,subpart(p,lp,:lpvp))=>subpart(p,lp,:lpvv) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lp,:lpvv),:fv))) : vname(p,subpart(p,lp,:lpvv)) for lp in 1:nlpv(p)]
#     lvvs=[subpart(p,lv,:lvsrc) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvsrc),:fv))) : vname(p,subpart(p,lv,:lvsrc))=>subpart(p,lv,:lvtgt) in flowVariableIndexs ? fname(p,only(incident(p,subpart(p,lv,:lvtgt),:fv))) : vname(p,subpart(p,lv,:lvtgt)) for lv in 1:nlvv(p)]


#     es=vcat(lses,lsvfes,lfves,fies,foes,lpvs,lvvs)

#     return CausalLoop(ns,es)
# end

# function from_catlab_graph(g::Catlab.Graph, p::Vector{Polarity})
#   cl = CausalLoopF()
#   add_parts!(cl, :N, Catlab.nn(g))
#   add_parts!(cl, :E, Catlab.ne(g); s = subpart(g, :src), t = subpart(g, :tgt), epolarity = p)
#   cl
# end

# function from_graphs_graph(g::Graphs.Graph, p::Vector{Polarity})
#   cl = CausalLoopF()
#   add_parts!(cl, :N, Graphs.nn(g))
#   add_parts!(cl, :E, Graphs.ne(g); s = subpart(g, :src), t = subpart(g, :tgt), epolarity = p)
# end


# function simple_cycles(c::)
# end

# function to_graphs_graph(cl::CausalLoopPM)
#   edges = collect(zip(vcat(subpart(cl, :sp), subpart(cl, :sm)), vcat(subpart(cl, :tp), subpart(cl, :tm))))
#   g = SimpleDiGraph(SimpleEdge.(edges))
#   g
#   # return (g, np(cl))
# end

function to_graphs_graph(cl::CausalLoopPol)
  g = SimpleDiGraph(SimpleEdge.(zip(subpart(cl, :src), subpart(cl, :tgt))))
  g
end

# function cl_cycles(cl::CausalLoopPM)
#   g, last_pos = to_graphs_graph(cl)
#   for cycle ∈ simplecycles(g)
#     cycle_length = length(cycle)
#     node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)

#   end
# end

function cl_cycles(cl::K) where K <: AbstractCausalLoop
  cl_cycles(to_clp(cl))
end

function cl_cycles(cl::CausalLoopPol)
  edges = collect(zip(subpart(cl, :src), subpart(cl, :tgt)))
  # Unique are sufficient for making simple graph.
  g = SimpleDiGraph(SimpleEdge.(edges))

  all_cycles = Vector{Vector{Int}}()
  # Edges => Polarity
  for cycle ∈ simplecycles(g)
    cycle_length = length(cycle)
    # Last pair is cycle[end], cycle[1]
    node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)
    nonsimple_cycles = Vector{Vector{Int}}()
    for (p1, p2) in node_pairs
      # grab all edges with p1 as start and p2 as end
      push!(nonsimple_cycles, Vector{Int}(intersect(incident(cl, p1, :src), incident(cl, p2, :tgt))))
    end
    # generated_cycles = Vector{Vector{Int}}(Base.Iterators.product(nonsimple_cycles...))
    # For loop instead of comprehension to get around product being multidimensional
    for c in Base.Iterators.product(nonsimple_cycles...)
      push!(all_cycles, collect(c))
    end
  end

  all_cycles

end


# function cl_cycles(cl::CausalLoopPol)

#   edges = collect(zip(subpart(cl, :s), subpart(cl, :t)))
#   pair_to_edge = state_dict(edges) # Note, this is unique pairs, not all.
#   # Unique are sufficient for making simple graph.
#   g = SimpleDiGraph(SimpleEdge.(edges))

#   # Edges => Polarity
#   cycle_pol = Vector{Pair{Vector{Int}, Polarity}}()
#   for cycle ∈ simplecycles(g)
#     cycle_length = length(cycle)
#     # Last pair is cycle[end], cycle[1]
#     node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)
#     nonsimple_cycles = Vector{Vector{Int}}()
#     for (p1, p2) in node_pairs
#       matching_edges = Vector{Int}()
#       # grab all edges with p1 as start and p2 as end
#       # This is the bit that could be made more efficient, it loops over all edges every time
#       for cl_node in 1:ne(cl)
#         if sedge(cl, cl_node) == p1 && tedge(cl, cl_node) == p2
#           push!(matching_edges, cl_node)
#         end
#       end
#       push!(nonsimple_cycles, matching_edges)
#     end

#     for cycle_instance in Base.Iterators.product(nonsimple_cycles...)
#       balancing_count = 0
#       is_unknown = false
#       is_zero = false
#       is_not_well_defined = false
#       for node_number in cycle_instance
#         current_pol = epol(cl, node_number)
#         if current_pol == POL_BALANCING
#           balancing_count += 1
#         elseif current_pol == POL_UNKNOWN
#           is_unknown = true
#         elseif current_pol == POL_ZERO
#           is_zero = true
#           break
#         elseif current_pol == POL_NOT_WELL_DEFINED
#           is_not_well_defined = true
#         end
#       end

#       collected_cycle = collect(cycle_instance)

#       if is_zero
#         push!(cycle_pol, collected_cycle => POL_ZERO)
#       elseif is_unknown
#         push!(cycle_pol, collected_cycle => POL_UNKNOWN)
#       elseif is_not_well_defined
#         push!(cycle_pol, collected_cycle => POL_NOT_WELL_DEFINED)
#       elseif iseven(balancing_count)
#         push!(cycle_pol, collected_cycle => POL_REINFORCING)
#       else
#         push!(cycle_pol, collected_cycle => POL_BALANCING)
#       end
#     end
#   end

#   cycle_pol
          
# end



# function cl_cycles(cl::CausalLoopPM) 
#   edges = collect(zip(vcat(subpart(cl, :sp), subpart(cl, :sm)), vcat(subpart(cl, :tp), subpart(cl, :tm))))
#   # Unique are sufficient for making simple graph.
#   g = SimpleDiGraph(SimpleEdge.(edges))

#   all_cycles = Vector{Vector{Int}}()

#   # Edges => Polarity
#   for cycle ∈ simplecycles(g)
#     cycle_length = length(cycle)
#     # Last pair is cycle[end], cycle[1]
#     node_pairs = ((cycle[node_index], cycle[(node_index % cycle_length) + 1]) for node_index in 1:cycle_length)
#     nonsimple_cycles = Vector{Vector{Int}}()
#     for (p1, p2) in node_pairs
#       # grab all edges with p1 as start and p2 as end
#       push!(nonsimple_cycles, Vector{Int}(
#         union(
#           intersect(incident(cl, p1, :sp), incident(cl, p2, :tp)),
#           intersect(incident(cl, p1, :sm), incident(cl, p2, :tm)) .+ np(cl)
#         )
#       )
#       )
#     end
#     # generated_cycles = Vector{Vector{Int}}(Base.Iterators.product(nonsimple_cycles...))
#     # For loop instead of comprehension to get around product being multidimensional
#     for c in Base.Iterators.product(nonsimple_cycles...)
#       push!(all_cycles, collect(c))
#     end
#   end

#   all_cycles

# end



function extract_loops(cl::K) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  cycles = cl_cycles(cl)
  map(x -> x => walk_polarity(cl, x), cycles)
end

# TODO: terrible
epol(cl::CausalLoopPM, e) = begin
  @assert e <= nedges(cl)
  e <= np(cl) ? POL_REINFORCING : POL_BALANCING
end

epol(cl::CausalLoopZ, e) = begin
  @assert e <= nedges(cl)
  if e <= np(cl)
    POL_REINFORCING
  elseif e <= (np(cl) + nm(cl))
    POL_BALANCING
  else
    POL_ZERO
  end
end

epol(cl::CausalLoopFull, e) = begin
  @assert e <= nedges(cl)
  if e <= np(cl)
    POL_REINFORCING
  elseif e <= (np(cl) + nm(cl))
    POL_BALANCING
  elseif e <= (np(cl) + nm(cl) + nz(cl))
    POL_ZERO
  elseif e <= (np(cl) + nm(cl) + nz(cl) + nu(cl))
    POL_UNKNOWN
  else
    POL_NOT_WELL_DEFINED
  end
end


function walk_polarity(cl::K, edges::Vector{Int}) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  foldl(*, map(x -> epol(cl, x), edges); init = POL_REINFORCING)
end


function extract_all_nonduplicate_paths(clp::CausalLoopPol)


  function rec_search!(path, nodes, paths)
    target = tedge(clp, path[end])
    outgoing = Vector{Int}(incident(clp, target, :src)) # all connected edges from target vertex

    for o in outgoing
      outgoing_tgt = tedge(clp, o) # note, clp is defined in outer func
      if outgoing_tgt ∈ nodes
        continue
      end

      new_nodes = Set{Int}([nodes..., outgoing_tgt])

      new_path = [path..., o]
      push!(paths, new_path => paths[path] * epol(clp, o))
      rec_search!(new_path, new_nodes, paths)
    end
  end


  paths = Dict{Vector{Int}, Polarity}(Vector{Int}() => POL_REINFORCING)
  for e in 1:nedges(clp)
    nodes = Set{Int}([sedge(clp, e)])
    push!(paths, [e] => epol(clp, e))
    rec_search!([e], nodes, paths)
  end

  paths

end

function extract_all_nonduplicate_paths(cl::K) where K <: AbstractCausalLoop
  clp = to_clp(cl)
  extract_all_nonduplicate_paths(clp)
end

# function walk_polarity(cl::CausalLoopF, edges::Vector{Int})
#   foldl(*, map(x -> epol(cl, x), edges); init = POL_REINFORCING)
# end

""" 
Cycles are uniquely characterized by sets of edges, not sets of nodes
We construct a simple graph, then for each edge, we check if there exist
multiple edges with the same start and end node
We then product all of them, to get every cycle.

This could be made more efficient, but it should be fine for now.
"""
# function extract_loops(cl::CausalLoopF)

#   cycle_pol = Vector{Pair{Vector{Int}, Polarity}}()


#   for cycle_instance in cl_cycles(cl)
#     collected_cycle = collect(cycle_instance)
#     push!(cycle_pol, collected_cycle => walk_polarity(cl, collected_cycle))
#   end

#   cycle_pol
          
# end

function is_walk(cl::CausalLoopPM, edges::Vector{Int})
  all(x -> tedge(cl, edges[x]) == sedge(cl, edges[x+1]), eachindex(edges[1:end-1]))
end

function is_circuit(cl::CausalLoopPM, edges::Vector{Int})
  is_path(cl, edges) && sedge(cl, edges[1]) == tedge(cl, edges[end])
end

# TODO: How, pray tell, is this a functor
function betweenness(cl::K) where K <: AbstractCausalLoop
  g = to_graphs_graph(cl)
  Graphs.betweenness_centrality(g)
end


# function num_loops(cl::CausalLoopF, name::Symbol)
#   el = cl_cycles(cl)
#   node_num = only(incident(cl, :nname, name))
#   return count(x -> node_num ∈ x, el)
# end






# # Graphs.betweenness_centrality
# # ^ Use this if we don't care about polarities
# function betweenness(cl::CausalLoopF; max_edges = typemax(Int))
#   g = to_graphs_graph(cl)
#   num_graph_vert = vertices(g)[end] # ensures we don't go over if the conversion deleted the end 
#   # deleting the end can only happen if there are no edges which come from or leave the node, so the betweenness centrality value for it is 0


#   # past a certain threshold, we start getting infinity confused with very large.

#   # Though we're using typemax(Int) for this so why does this matter
#   @assert length(edges(g)) < max_edges

#   betweenness_values = zeroes(nn(cl))

#   for node in nn(cl)
#     if node > num_graph_vert
#       break
#     end

#     dij = dijkstra_shortest_paths(g, node ; maxdist=max_edges)
#     preds = dij.predecessors

#     # cache the 




#   end



  
# end
# NOTE: simplecycles returns NODES!!!

#TODO: Deal with nameless AbstractCausalLoop, or create a subtype without nameless

function to_graphs_graph(cl::K) where K <: AbstractCausalLoop
  to_graphs_graph(to_clp(cl))
  # edges = collect(zip(subpart(cl, :s), subpart(cl, :t)))
  # g = SimpleDiGraph(SimpleEdge.(edges)) # Note, this will discard the final nodes if they have no edges
  # g 
end

# Acts as if there can be more than one edge between nodes
function num_loops_var_on(c::K, name::Symbol) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  name_index = only(incident(c, name, :vname))
  node_cycles = map(x -> map(y -> tedge(c, y), x), cl_cycles(c)) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), node_cycles)
end

# Acts as if there is at most one edge between nodes
function num_indep_loops_var_on(c::K, name::Symbol) where K <: Union{AbstractCausalLoop, CausalLoopPol}
  g = to_graphs_graph(c)
  sc = simplecycles(g)
  # @show sc
  name_index = only(incident(c, name, :vname))
  # retuer

  # node_cycles = map(x -> map(y -> tedge(c, y), x),sc) # sedge is equivalent here.
  # if a node is in a cycle, of course, it will be in both the src set and the tgt set
  return count(∋(name_index), sc)
end


# function to_catlab_graph(cl::K) where K <: AbstractCausalLoop
#   g = Catlab.Graph(nn(cl))
#   add_parts!(g, :E, ne(cl) ; src = subpart(cl, :s), tgt = subpart(cl, :t))
#   g
# end



# ! This works, but we have version conflicts right now
# ! I added DataMigrations to TOML; presumably, just need to wait a few days for updated requirements in used packages

# function to_catlab_graph(cl::CausalLoopNameless)
#   mig = @migration SchGraph TheoryCausalLoopNameless begin
#     E => @cases (p::P; m::M)
#     V => V
#     src => begin
#       p => sp
#       m => sm
#     end
#     tgt => begin
#       p => tp
#       m => tm
#     end 
#   end
#   return migrate(Graph, cl, mig)
# end

# function to_graphs_graph(cl::)
# end


# function discard_zero_pol(c)
#   cl = CausalLoopF()
#   add_vertices!(cl, nn(c) ; nname = nnames(c))
#   for edge in 1:ne(c)
#     pol = epol(c, edge)
#     if pol != POL_ZERO
#       add_edge!(cl, sedge(c, edge), tedge(c, edge) ; epolarity = pol)
#     end
#   end
#   cl
# end




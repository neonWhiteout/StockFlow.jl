module UsefulFunctions

using StockFlow
using StockFlow.Syntax
using MLStyle


export get_link_index, generate_variable_names, infer_links



function get_link_index(sf, part1::String, part2::String)


    stocks = [String(k) for k in values((getfield(sf, :subparts)[:sname]).m)]
    
    flows = [String(k) for k in values((getfield(sf, :subparts)[:fname]).m)]
    # sort!(flows, by = key -> ((getfield(sf, :subparts)[:fv]).m)[:key]) # sorts flows by fv

    sums = [String(k) for k in values((getfield(sf, :subparts)[:svname]).m)]
    dynamic_variables = [String(k) for k in values((getfield(sf, :subparts)[:vname]).m)]
    params = [String(k) for k in values((getfield(sf, :subparts)[:pname]).m)]

    @assert allunique(vcat(stocks, flows, sums, dynamic_variables, params)) # make sure all the variable names are unique!

    # TODO: Include an assert that all the elements of the links are unique?

    @match part1 begin
        if part1 ∈ stocks && part2 ∈ sums end => begin
            v1 = findfirst(item -> item == part1, stocks)
            v2 = findfirst(item -> item == part2, sums)
            it = collect(zip([k[2] for k in (getfield(sf, :subparts)[:lss]).m], [k[2] for k in (getfield(sf, :subparts)[:lssv]).m]))

            return Dict("LS" => findfirst(item -> item == (v1, v2), it))
        end
        if part1 ∈ stocks && part2 ∈ dynamic_variables end => begin
            v1 = findfirst(item -> item == part1, stocks)
            v2 = findfirst(item -> item == part2, dynamic_variables)
            it = collect(zip([k[2] for k in (getfield(sf, :subparts)[:lvs]).m], [k[2] for k in (getfield(sf, :subparts)[:lvv]).m]))

            return Dict("LV" => findfirst(item -> item == (v1, v2), it))
        end
        if part1 ∈ flows && part2 ∈ stocks end => begin
            v1 = findfirst(item -> item == part1, flows)
            v2 = findfirst(item -> item == part2, stocks)
            it1 = collect(zip([k[2] for k in (getfield(sf, :subparts)[:ifn]).m], [k[2] for k in (getfield(sf, :subparts)[:is]).m]))
            it2 = collect(zip([k[2] for k in (getfield(sf, :subparts)[:ofn]).m], [k[2] for k in (getfield(sf, :subparts)[:os]).m]))

            # println(it1)
            # println(values(it1))
            # println(it2)
            # println(v1 => v2)

            
            

            flowdict = Dict("I" => findfirst(item -> item == (v1, v2), it1), "O" => findfirst(item -> item == (v1, v2), it2))
            return filter(x -> ~(isnothing(x[2])), flowdict)
        end
        if part1 ∈ sums && part2 ∈ dynamic_variables end => begin
            v1 = findfirst(item -> item == part1, sums)
            v2 = findfirst(item -> item == part2, dynamic_variables)
            it = collect(zip([k[2] for k in (getfield(sf, :subparts)[:lsvsv]).m], [k[2] for k in (getfield(sf, :subparts)[:lsvv]).m]))

            return Dict("LSV" => findfirst(item -> item == (v1, v2), it))
        end
        if part1 ∈ dynamic_variables && part2 ∈ dynamic_variables end => begin
            v1 = findfirst(item -> item == part1, dynamic_variables)
            v2 = findfirst(item -> item == part2, dynamic_variables)
            it = collect(zip([k[2] for k in (getfield(sf, :subparts)[:lvsrc]).m], [k[2] for k in (getfield(sf, :subparts)[:lvtgt]).m]))

            return Dict("LVV" => findfirst(item -> item == "I" => findfirst(i(v1, v2), it)))
        end
        if part1 ∈ params && part2 ∈ dynamic_variables end => begin
            v1 = findfirst(item -> item == part1, params)
            v2 = findfirst(item -> item == part2, dynamic_variables)
            it = collect(zip([k[2] for k in (getfield(sf, :subparts)[:lpvp]).m], [k[2] for k in (getfield(sf, :subparts)[:lpvv]).m]))

            return Dict("LPV" => findfirst(item -> item == (v1, v2), it))
        end
        _ => begin 
            println("Failed to find a link between provided pair!  Perhaps you provided the arguments in the wrong order?")
            return Dict()
        end
    end




end


function generate_variable_names(sf)
    stocks = [String(k) for k in values((getfield(sf, :subparts)[:sname]).m)]
    
    flows = [String(k) for k in values((getfield(sf, :subparts)[:fname]).m)]
    # sort!(flows, by = key -> ((getfield(sf, :subparts)[:fv]).m)[:key]) # sorts flows by fv

    sums = [String(k) for k in values((getfield(sf, :subparts)[:svname]).m)]
    dynamic_variables = [String(k) for k in values((getfield(sf, :subparts)[:vname]).m)]
    params = [String(k) for k in values((getfield(sf, :subparts)[:pname]).m)]

    ops = [String(k) for k in values((getfield(sf, :subparts)[:vop]).m)]

    dv_op = ["dvop_" * String(dynamic_variables[i]) * "_" * String(ops[i]) for i in 1:length(dynamic_variables)]
    
    lvsposition = [k for k in values((getfield(sf, :subparts)[:lvsposition]).m)]
    lsvsvposition = [k for k in values((getfield(sf, :subparts)[:lsvsvposition]).m)]
    lvsrcposition = [k for k in values((getfield(sf, :subparts)[:lvsrcposition]).m)]
    lpvpposition = [k for k in values((getfield(sf, :subparts)[:lpvpposition]).m)]


    fv = ["fv_" * String(flows[i]) * "_" * String(dynamic_variables[j]) for (i, j) in zip(1:length(flows), values((getfield(sf, :subparts)[:fv]).m))]

    ls = ["ls_" * String(stocks[i]) * "_" * String(sums[j]) for (i,j) in zip(values((getfield(sf, :subparts)[:lss]).m), values((getfield(sf, :subparts)[:lssv]).m))]
    inflows = ["i_" * String(flows[i]) * "_" * String(stocks[j]) for (i,j) in zip(values((getfield(sf, :subparts)[:ifn]).m), values((getfield(sf, :subparts)[:is]).m))]
    outflows = ["o_" * String(flows[i]) * "_" * String(stocks[j]) for (i,j) in zip(values((getfield(sf, :subparts)[:ofn]).m), values((getfield(sf, :subparts)[:os]).m))]
    lv = ["lv_" * String(stocks[i]) * "_" * String(dynamic_variables[j]) * "_" * string(lvsposition[k]) for (i,j,k) in zip(values((getfield(sf, :subparts)[:lvs]).m), values((getfield(sf, :subparts)[:lvv]).m), values((getfield(sf, :subparts)[:lvsposition]).m))]
    lsv = ["lsv_" * String(sums[i]) * "_" * String(dynamic_variables[j]) * "_" * string(lsvsvposition[k]) for (i,j,k) in zip(values((getfield(sf, :subparts)[:lsvsv]).m), values((getfield(sf, :subparts)[:lsvv]).m),  values((getfield(sf, :subparts)[:lsvsvposition]).m))]
    lvv = ["lvv_" * String(dynamic_variables[i]) * "_" * String(dynamic_variables[j]) * "_" * string(lvsrcposition[k]) for (i,j,k) in zip(values((getfield(sf, :subparts)[:lvsrc]).m), values((getfield(sf, :subparts)[:lvtgt]).m), values((getfield(sf, :subparts)[:lvsrcposition]).m))]

    lpv = ["o_" * String(params[i]) * "_" * String(dynamic_variables[j]) * "_" * string(lpvpposition[k]) for (i,j,k) in zip(values((getfield(sf, :subparts)[:lpvp]).m), values((getfield(sf, :subparts)[:lpvv]).m), values((getfield(sf, :subparts)[:lpvpposition]).m))]


    elements = ["stocks" => stocks, "flows" => flows, "sums" => sums, "dynamic_variables" => dynamic_variables, "params" => params, "ops" => ops, "dv_op" => dv_op, "fv" => fv, "ls" => ls, "i" => inflows, 
    "o" => outflows, "lv" => lv, "lsv" => lsv, "lvv" => lvv, "lpv" => lpv]

    
end



function extract_subpart(sf, subpart::String)
    return extract_subpart(sf, Symbol(subpart))
end

function extract_subpart(sf, subpart::Symbol)
    return (getfield(sf, :subparts)[subpart]).m
end

function extract_subpart_ordered(sf, subpart::String)
    return extract_subpart_ordered(sf, Symbol(subpart))
end

function extract_subpart_ordered(sf, subpart::Symbol)
    return [k[2] for k in extract_subpart(sf, subpart)] # Here's hoping it's ordered, heh.  It should be.
    # As in, the pairs should be 1 => _, 2 => _, ..., rather than 5 => _, 14 => _, ...
end

function extract_subpart_ordered_tuple(sf, subpart::String...)
    return extract_subpart_ordered_tuple(sf, [Symbol(k) for k in subpart]) # TODO: fix lmao
end

function extract_subpart_ordered_tuple(sf, subpart::Symbol...)
    if isempty(subpart)
        return []
    end

    subparts = [extract_subpart_ordered(sf, sp) for sp in subpart]
    # println(subparts)
    # println("HERE")
    collected_tuples = [tuple(vec...) for vec in zip(subparts...)]
    # println(collected_tuples)
    return collected_tuples
    # k = collect([extract_subpart_ordered(sf, sp) for sp in subpart])
    # println(k)
    # return k
end



"""
In a lot of cases, we don't need to explicitly define where links map to; 
if A <- C -> B, and we have A -> A' and B -> B' and a unique C' such that A' <- C' -> B', we can assume C -> C'
This can extend to variables and parameters, though they're a bit more complicated; if v = op(a, b), a -> a', b -> b', op -> op',
and there exists a v' such that v' = op'(a', b'), we can assume v -> v'.  May require tracing a tree to disambiguate. 


I'm *thinking* that in a lot of cases, we only need to define the mappings of stocks, flows, and potentially sum variables.  The rest can be inferred.
If we have the mappings of stocks, flows, sv, vars and params, then everything else should be able to be inferred.

:S => [2,4,1,3], ...
necMaps must contain keys S, F, SV, P, V
"""
function infer_links(sfsrc, sftgt, necMaps::Dict{Symbol, Vector{Int64}})
    LS = :LS => infer_particular_link(sfsrc, sftgt, :lss => necMaps[:S], :lssv => necMaps[:SV])
    I = :I => infer_particular_link(sfsrc, sftgt, :ifn => necMaps[:F], :is => necMaps[:S])
    O = :O => infer_particular_link(sfsrc, sftgt, :ofn => necMaps[:F], :os => necMaps[:S])
    LV = :LV => infer_particular_link(sfsrc, sftgt, :lvs => necMaps[:S], :lvv => necMaps[:V])
    LSV = :LSV => infer_particular_link(sfsrc, sftgt, :lsvsv => necMaps[:SV], :lsvv => necMaps[:V])
    LVV = :LVV => infer_particular_link(sfsrc, sftgt, :lvsrc => necMaps[:V], :lvtgt => necMaps[:V])
    LPV = :LPV => infer_particular_link(sfsrc, sftgt, :lpvp => necMaps[:P], :lpvv => necMaps[:V])

    return Dict(LS, I, O, LV, LSV, LVV, LPV)
end


# function infer_particular_link(sfsrc, sftgt, linkname1::Symbol, linkname2::Symbol, map1::Vector{Int64}, map2::Vector{Int64}) # linkname isn't necessary but means we don't need to check all subparts


function infer_particular_link(sfsrc, sftgt, link1::Pair{String, Vector{Int64}}, link2::Pair{String, Vector{Int64}})
    linkname1, map1 = link1
    linkname2, map2 = link2
    return infer_particular_link(sfsrc, sftgt, Symbol(linkname1) => map1, Symbol(linkname2) => map2)
end

"""
This isn't foolproof; may need to ensure position is preserved.
"""
function infer_particular_link(sfsrc, sftgt, link1::Pair{Symbol, Vector{Int64}}, link2::Pair{Symbol, Vector{Int64}}) # linkname isn't necessary but means we don't need to check all subparts
    # maps go index -> index.  EG, [1,3,2] means 1 -> 1, 2 -> 3, 3 -> 2
    # linkmap = collect(zip(map1, map2))::Vector{Tuple{Int64, Int64}}

    linkname1, map1 = link1
    linkname2, map2 = link2

    link_domain = extract_subpart_ordered_tuple(sfsrc, linkname1, linkname2)::Vector{Tuple{Int64, Int64}}
    # Corresponds to column 2 and 3 when printing the sf
    # Note, does NOT currently include position, which could lead to ambiguity
    link_codomain = extract_subpart_ordered_tuple(sftgt, linkname1, linkname2)::Vector{Tuple{Int64, Int64}}
    # link_codomain = collect(zip(extract_subpart_ordered(sftgt, linkname1), extract_subpart_ordered(sftgt, linkname2)))

    # link_maps = [map1, map2]

    inferred_links = Vector{Int64}() # points to index should map to

    for (i, (x1, x2)) in enumerate(link_domain) # eg, v1 is a stock, v2 is a sum variable, and they have a link.  Let's say stock A -> A′, sv N -> N′.
        # We're checking that there exists a unique link between A′ and N′.
        # We do this by grabbing all the links L in the domain and L′ in the codomain, converting L based on what its constituents map to, checking that
        # there's a unique match, and appending the unique match's index to a vector (of int64, since everything is just an int)
        
        m1, m2 = map1[x1], map2[x2] 
        
        potential_values = Vector{Int64}()
        for (j, (y1, y2)) in enumerate(link_codomain)
            if (m1 == y1) && (m2 == y2)
                push!(potential_values, j)
            end
        end
        # potential_values = filter(((index, (i, j)),) -> i == m1 && j == m2, [(index, (i, j)) for (index, (i,j)) in enumerate(link_codomain)])
        @assert length(potential_values) == 1 "Didn't find one value to map ($v1, $v2) to for $linkname1, $linkname2 !  Found: $(length(potential_values))"
        push!(inferred_links, potential_values[1])
    end

    return inferred_links
end


# function infer_ls(sfsrc, sftgt, smap::Dict{Int64, Int64}, svmap::Dict{Int64, Int64})

# end




# function infer_links(sfsrc, sftgt, necMaps::Dict{String, Vector{String}})
# end


end
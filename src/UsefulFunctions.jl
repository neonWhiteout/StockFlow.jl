module UsefulFunctions

using StockFlow
using StockFlow.Syntax
using MLStyle
using DataStructures # new requirement maybe

export get_link_index, generate_variable_names, infer_links, infer_from_stocks_and_flows



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

function extract_subpart_ordered_tuple(sf, subpart::Symbol...)::Vector{Vector{Any}} # usually Int64 but works with names too.  Might just be Int64 and Symbol
    if isempty(subpart)
        return []
    end

    subparts = [extract_subpart_ordered(sf, sp) for sp in subpart]
    # println(subparts)
    # println("HERE")
    # collected_tuples = [tuple(vec...) for vec in zip(subparts...)]
    collected_tuples = [collect(vec) for vec in zip(subparts...)]

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
    # could probably include vop by just giving vop as link1 and link2 arguments?
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

    link_domain = extract_subpart_ordered_tuple(sfsrc, linkname1, linkname2) # ::Vector{Vector{Int64}} (at present, it's just Vector{Vector{Any}})
    # Corresponds to column 2 and 3 when printing the sf
    # Note, does NOT currently include position, which could lead to ambiguity
    link_codomain = extract_subpart_ordered_tuple(sftgt, linkname1, linkname2) # ::Vector{Vector{Int64}}
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
        @assert length(potential_values) == 1 "Didn't find one value to map ($x1, $x2) to for $linkname1, $linkname2 !  Found: $(potential_values)"
        push!(inferred_links, potential_values[1])
    end

    return inferred_links
end


# function infer_ls(sfsrc, sftgt, smap::Dict{Int64, Int64}, svmap::Dict{Int64, Int64})

# end




# function infer_links(sfsrc, sftgt, necMaps::Dict{String, Vector{String}})
# end


# function infer_sums(sfsrc, sftgt, )

# function extract_nth_mappings(sf, i, fields...)
#     return extract_subpart_ordered_tuple(sf, fields)[i]
# end

# We're dealing with finset which means we can totally brute force everything
# We shouldn't, but we can
"""
I'm not sure if this will be possible, but worth a try, I suppose.

We just need to infer sums, dynamic variables and parameters, then run infer_links

necMaps needs only :S => [3, 5, ...] and :F => [4, 6, 1, ...]
"""
function infer_from_stocks_and_flows(sfsrc, sftgt, necMaps) # TODO: Make the variable names not absolutely garbage
    # Like bloody hell I have no idea what's going on
    # part of that is a consequence of defining mappings between mappings but part of it is just me being dumb
    # maybe replace instances of 'var' with 'dyvar' or 'dv'?


    stockmaps = necMaps[:S]
    flowmaps = necMaps[:F]










    var_domain = 1:nvb(sfsrc) # literally just 1:n
    var_codomain = 1:nvb(sftgt) # 1:m

    flowvar_domain = extract_subpart_ordered(sfsrc, :fv)
    flowvar_codomain = extract_subpart_ordered(sftgt, :fv)

    flowvar_maps = zeros(Int64, length(var_domain)) # TODO: just use nvb instead of length
    flowvar_reverse_maps = [Set{Int64}() for _ in var_codomain]

    unmapped_vars = Set{Int64}(var_domain)


    # println(flowmaps)

    for (srcflow, tgtflow) in enumerate(flowmaps)
        srcflow_var = flowvar_domain[srcflow]
        tgtflow_var = flowvar_codomain[tgtflow]

        flowvar_maps[srcflow_var] = tgtflow_var

        push!(flowvar_reverse_maps[tgtflow_var], srcflow_var)

        # push!(flowvar_maps, (srcflow_var, tgtflow_var))
        delete!(unmapped_vars, srcflow_var)
    end # ok now we actually have the maps flowvar -> flowvar





    # var_mappings = zeros(Int64, length(var_domain)) # initialize array for var mappings.
    # var_mappings_successful_set = Set(var_domain)
    # # Working off the assumption that the flowvar won't change when mapping.  
    # println(flowvar_domain, var_mappings, var_codomain, flowvar_codomain)

    # for (i, x) in enumerate(flowvar_domain) # i is index in domain fv, x is actual entry.  x points to a variable index in domain.
    #     var_mappings[x] = var_codomain[flowvar_codomain[i]]
    #     pop!(var_mappings_successful_set, x)
    #     # (fv -> v) = (corresponding fv -> v)
    # end

    # So by this point we should have all the flow variables correctly mapped.  We now need to look at non-flow variables.
    # We can check if there exist unique variable/stock or 

    lv_domain = extract_subpart_ordered_tuple(sfsrc, :lvs, :lvv)
    lv_codomain = extract_subpart_ordered_tuple(sftgt, :lvs, :lvv)

    lsv_domain = extract_subpart_ordered_tuple(sfsrc, :lsvsv, :lsvv)
    lsv_codomain = extract_subpart_ordered_tuple(sftgt, :lsvsv, :lsvv)

    lvv_domain = extract_subpart_ordered_tuple(sfsrc, :lvsrc, :lvtgt)
    lvv_codomain = extract_subpart_ordered_tuple(sftgt, :lvsrc, :lvtgt)

    lpv_domain = extract_subpart_ordered_tuple(sfsrc, :lpvp, :lpvv)
    lpv_codomain = extract_subpart_ordered_tuple(sftgt, :lpvp, :lpvv)


    var_maps = copy(flowvar_maps)::Vector{Int64} # TODO: make this a nonsucky data structure.  Vector{Int64} initialized with 0 is probably fine.
    # k fixed should be ok now just
    reverse_maps = copy(flowvar_reverse_maps)

    while !(isempty(unmapped_vars)) # if there exist variables which aren't just flow variables
        # things get weird here.  Now, we're going to check for variable links with known variables and see if they're unique.
        # possibly requires recursion, eg if v_x = v_a + v_b, v_x′ = v_a + v_c
        # We will also potentially need to check for stock or sum variable links

        println("UNMAPPED: $(unmapped_vars)")


        vars_to_check_vector = collect(unmapped_vars)


        # for vars, src is on right side, tgt is on left.  So, we're looking for src; we want to go backwards, trying to find flowvars made of unknown vars.
        # Unless we have a weird case where a variable has no relation to a flow, we should find all of them

        # known: lvtgt
        # unknown: lvsrc

        for var::Int64 in vars_to_check_vector
            varlinks = filter((((lvsrc, lvtgt),) -> lvsrc == var), lvv_domain) # all variable -> variable links with var

            println("varlinks:$varlinks")

            known_varlinks = filter(((lvsrc, lvtgt),) -> !(lvtgt in unmapped_vars), varlinks) # all variable links such that we know what tgt maps to
            # specifically, all links which contain var, and some other variable


            println("known_variables:$known_varlinks")
            # If our known variable links to multiple unknown, then we're screwed; can't figure out which link this one should go with.  So, need to filter those out.
            # This will only occur in cases where a variable is result of two other variables.  Should be able to fix this by refactoring.
            # So there are some cases where this just won't work at all.
            known_counter = DefaultDict{Int64,Int64}(0)

            # known_varlinkset = Set{Int64}([k[2] for k in known_varlinks])

            # for (lvsrc, lvtgt) in lvv_codomain # uh I think this can be done outside the loop?  
            #     # TODO: figure out if I can move it without any trouble
            #     if length(flowvar_reverse_maps[lvtgt]) 
            #     known_counter[lvtgt] += 1
            # end
            println("known_counter: $known_counter")
            println("reverse_maps: $reverse_maps")

            known_varlinks_with_unique_known = filter(((lvsrc, lvtgt),) -> length(filter(x -> !(x in unmapped_vars), 
            reverse_maps[lvsrc])) < 2, known_varlinks)
                

            println("known_varlinks_with_unique_known: $known_varlinks_with_unique_known")
            
            println("YEYEY")
            if length(known_varlinks_with_unique_known) == 1 # Hooray!  Unique link! 
                # that is, at this point, we have a dynamic variable which we know what it maps to, and we know it has one variable link
                # so, what we can do is find that variable link in the codomain and use that to figure out what this variable should map to
                known_varlink = known_varlinks_with_unique_known[1]
                known_varlink_tgt = known_varlink[2]
                println(known_varlink_tgt)
                println(var_maps)

                known_var_image = var_maps[known_varlink_tgt]

                # known_var_image = var_maps[src]

                # known_var_image = filter(((src, tgt),) -> src == known_varlink_tgt, var_maps)[1][2] # if there are more than one of these something has gone horribly wrong.
                # TODO: Change datatype of flowvar_maps/var_maps.  Probably just make it a straight Vector{Int64} of length var_domain and default vals 0.  Standard index -> value.
                println("known_var_image: $known_var_image")
                println("lvv_codomain: $lvv_codomain")

                known_var_image_varlink = filter(((lvsrc, lvtgt),) -> lvtgt == known_var_image && lvsrc in unmapped_vars, lvv_codomain) # once again, if there are more than one, something has gone horribly wrong
                # known_varlinks_with_unique_known should be forcing this to be length 1
                if !(length(known_var_image_varlink) == 1)
                    println("didn't get one varlink!  $(known_var_image_varlink), continuing!")
                    continue
                end
                # @assert length(known_var_image_varlink) == 1 "didn't get one varlink!  $(known_var_image_varlink)"

                known_var_image_varlink_src = known_var_image_varlink[1][1]

                var_maps[var] = known_var_image_varlink_src


                push!(reverse_maps[known_var_image_varlink_src], var)
                # push!(var_maps, (var, known_var_image_varlink_src))
                delete!(unmapped_vars, var)
                println("unmapped_vars: $unmapped_vars")
                println("DSAD")
            end


            # println(varlinks)
        end

        if length(unmapped_vars) == length(vars_to_check_vector) # that is, we iterated over all the elements in the set and discovered no new links.
            # if this happens, we do not have enough information to determine the remaining dynamic_variable links from flowlinks alone.
            # TODO: Add checks with stocks and sum variables.
            println("ah crap")
            println("Set contains: $(unmapped_vars)")
            @assert false

        end


    println("var_maps: $var_maps")
    println("Unmapped: $unmapped_vars")



    end

end


    # I suppose we can start with looking at flow variables?






end
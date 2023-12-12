using ..StockFlow.Syntax
using MLStyle
using MLStyle.Modules.AST

@present TheoryStockAndFlowG <: TheoryStockAndFlowF begin
  
  Polarity::AttrType
  
  lspol::Attr(LS, Polarity)
  opol::Attr(O, Polarity)
  lsvpol::Attr(LSV, Polarity)
  lvvpol::Attr(LVV, Polarity)
  lpvpol::Attr(LPV, Polarity)
  
  Exponent::AttrType

  U::Ob
  CU::Ob
  LUCU::Ob


  uname::Attr(U, Name)
  exp::Attr(LUCU, Exponent)
  

  lucuu::Hom(LUCU, U)
  lucucu::Hom(LUCU, CU)

  scu::Hom(S, CU)
  svcu::Hom(SV, CU)
  fcu::Hom(F, CU)
  vcu::Hom(V, CU)
  pcu::Hom(P, CU)

end


@abstract_acset_type AbstractStockAndFlowG <: AbstractStockAndFlowF
@acset_type StockAndFlowGUntyped(TheoryStockAndFlowG, index=[:is,:os,:ifn,:ofn,:fv,:lvs,:lvv,:lsvsv,:lsvv,:lss,:lssv,:lvsrc,:lvtgt,:lpvp,:lpvv, :lucuu, :lucucu, :scu, :svcu, :fcu, :vcu, :pcu]) <: AbstractStockAndFlowG
const StockAndFlowG = StockAndFlowGUntyped{Symbol,Symbol,Int8,Int8,Int8}

add_unit!(p::AbstractStockAndFlowG;kw...) = add_part!(p,:U;kw...)
add_units!(p::AbstractStockAndFlowG,n;kw...) = add_parts!(p,:U,n;kw...)

add_cunit!(p::AbstractStockAndFlowG;kw...) = add_part!(p,:CU;kw...)
add_cunits!(p::AbstractStockAndFlowG,n;kw...) = add_parts!(p,:CU,n;kw...)

add_UCUlink!(p::AbstractStockAndFlowG,u,cu;kw...) = add_part!(p,:LUCU;lucuu=u,lucucu=cu,kw...)
add_UCUlinks!(p::AbstractStockAndFlowG,n,u,cu;kw...) = add_parts!(p,:LUCU,n;lucuu=u,lucucu=cu, kw...)

uname(p::AbstractStockAndFlowG,u) = subpart(p,u,:uname) 
unames(p::AbstractStockAndFlowG,u) = subpart(p,:uname) 

set_unames!(p::AbstractStockAndFlowG, names) = set_subpart!(p, :uname, names)
set_exps!(p::AbstractStockAndFlowG, exps) = set_subpart!(p, :exp, exps)

set_lspols!(p::AbstractStockAndFlowG, pols) = set_subpart!(p, :lspol, pols)
set_opols!(p::AbstractStockAndFlowG, pols) = set_subpart!(p, :opol, pols)
set_lsvpols!(p::AbstractStockAndFlowG, pols) = set_subpart!(p, :lsvpol, pols)
set_lvvpols!(p::AbstractStockAndFlowG, pols) = set_subpart!(p, :lvvpol, pols)
set_lpvpols!(p::AbstractStockAndFlowG, pols) = set_subpart!(p, :lpvpol, pols)




function convert(::Type{StockAndFlowG}, sff::AbstractStockAndFlowF)
    sfg = StockAndFlowG()
    
    add_stocks!(sfg, ns(sff), sname=snames(sff))
    add_svariables!(sfg, nsv(sff), svname=svnames(sff))
    add_variables!(sfg, nvb(sff), vname=vnames(sff))
    add_flows!(sfg, (fv(sff, i) for i in 1:nf(sff)), nf(sff), fname=fnames(sff))
    add_parameters!(sfg, np(sff), pname=pnames(sff))
    
    add_inflows!(sfg, ni(sff), subpart(sff, :is), subpart(sff, :ifn))
    add_outflows!(sfg, no(sff), subpart(sff, :os), subpart(sff, :ofn))
    add_Vlinks!(sfg, nlv(sff), subpart(sff, :lvs), subpart(sff, :lvv))
    add_Slinks!(sfg, nls(sff), subpart(sff, :lss), subpart(sff, :lssv))
    add_SVlinks!(sfg, nsv(sff), subpart(sff, :lsvsv), subpart(sff, :lsvv))
    add_VVlinks!(sfg, nlvv(sff), subpart(sff, :lvsrc), subpart(sff, :lvtgt))
    add_Plinks!(sfg, nlpv(sff), subpart(sff, :lpvp), subpart(sff, :lpvv))
    
    
    if (nvb(sff) > 0)
        vop = subpart(sff, :vop)
        set_subpart!(sfg, :vop, vop)
    end


    if (nlv(sff) > 0)
        lvv = subpart(sff, :lvsposition)
        set_subpart!(sfg, :lvsposition, lvv)
    end
    
    
    if (nlvv(sff) > 0)
        lvtgt = subpart(sff, :lvsrcposition)
        set_subpart!(sfg, :lvsrcposition, lvtgt)
    end
    

    if (nlsv(sff) > 0)
        lsvv = subpart(sff, :lsvsvposition)
        set_subpart!(sfg, :lsvsvposition, lsvv)
    end
    

    if (nlpv(sff) > 0)
        lpvv = subpart(sff, :lpvpposition)
        set_subpart!(sfg, :lpvpposition, lpvv)
    end
    
    return sfg

end


StockAndFlowG(sff::AbstractStockAndFlowF) = convert(StockAndFlowG, sff)





function parse_stock_units!(stocks, s::Union{Expr, Symbol})
    s_dict = @capture ($Stock:$units) s
    push!(stocks, (s_dict[:Stock], s_dict[:units]))
end

function parse_param_units!(params, p::Union{Expr, Symbol})
    p_dict = @capture ($Param:$units) p
    push!(params, (p_dict[:Param], p_dict[:units]))
end



function extract_exponents(exp)
    @match exp begin
        ::Symbol => 
            Dict(exp => 1)
        :(1/$B) => Dict(k => -v for (k, v) in extract_exponents(B))
        Expr(:call, :*, A, B...) =>
            merge!(+, extract_exponents(A), [extract_exponents(b) for b ∈ B]...)
        Expr(:call, :/, A, B...) =>
            merge!(+, extract_exponents(A), [Dict(k => -v for (k, v) in extract_exponents(b)) for b ∈ B]...)
        :($A^$n) =>
            Dict(A => n)
    end
end




function parse_stock_and_flow_units_syntax(statements::Vector{Any})
    stocks::Vector{Tuple{Symbol, Union{Symbol, Expr}}} = []
    params::Vector{Tuple{Symbol, Union{Symbol, Expr}}} = []
    dyvars::Vector{Tuple{Symbol,Expr}} = []
    flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
    sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []
    current_phase = (_, _) -> ()
    for statement in statements
        @match statement begin
            QuoteNode(:stocks) => begin
                current_phase = s -> parse_stock_units!(stocks, s)
            end
            QuoteNode(:parameters) => begin
                current_phase = p -> parse_param_units!(params, p)
            end
            QuoteNode(:dynamic_variables) => begin
                current_phase = d -> parse_dyvar!(dyvars, d)
            end
            QuoteNode(:flows) => begin
                current_phase = f -> parse_flow!(flows, f)
            end
            QuoteNode(:sums) => begin
                current_phase = s -> parse_sum!(sums, s)
            end
            QuoteNode(kw) =>
                error("Unknown block type for Stock and Flow syntax: " * String(kw))
            _ => current_phase(statement)
        end
    end

    s = StockAndFlowBlock(map(((s, u),) -> s,  stocks), map(((p, u),) -> p,  params), dyvars, flows, sums)
    sunits = map(((s, u),) -> u,  stocks)
    punits = map(((p, u),) -> u,  params)
    return s, sunits, punits
end















macro stock_and_flow_units(block)
    Base.remove_linenums!(block)
    block_args = block.args
    return quote
      local syntax_lines, sunits, punits = parse_stock_and_flow_units_syntax($block_args)
      local saff_args = stock_and_flow_syntax_to_arguments(syntax_lines)
      new_dyvars = map(kv -> kv.first => get(kv.second), saff_args.dyvars)
      sfu = StockAndFlowFUnits(
          saff_args.stocks,
          saff_args.params,
          new_dyvars,
          saff_args.flows,
          saff_args.sums,
      )
    #   println(new_dyvars)
        units = Dict{Symbol, Int}()
        for expr in sunits ∪ punits
            if typeof(expr) == Symbol
                if expr ∉ keys(units)
                    push!(units, expr => length(units) +1)
                end
            else
                for arg in expr.args[2:end]
                    if typeof(arg) == Expr
                        for subarg in arg.args[2:end]
                            if subarg ∉ keys(units)
                                push!(units, subarg => length(units)+1)
                            end
                        end
                    elseif arg ∉ keys(units)
                        push!(units, arg => length(units)+1)
                    end
                end
            end
        end


    #   units = unique([x.args[2:end]... for x in sunits] ∪ [x.args[2:end]... for x in punits])
        # units = unique([x.args[]])
        non_unitless = filter(x -> x != :UNITLESS, keys(units))
        add_parts!(sfu, :U, length(non_unitless); uname=non_unitless)
        # add_parts!(sfu, :LCUS, ns(sfu));

        add_parts!(sfu, :CU, ns(sfu) + nsv(sfu) + nf(sfu) + nvb(sfu) + np(sfu))
        
        stock_end = ns(sfu)
        sum_end = stock_end + nsv(sfu)
        flow_end = sum_end + nf(sfu)
        dyvar_end = flow_end + nvb(sfu)
        param_end = dyvar_end + np(sfu)

        

        set_subpart!(sfu, :scu, 1:stock_end)
        set_subpart!(sfu, :svcu, stock_end+1:sum_end)
        set_subpart!(sfu, :fcu, sum_end+1:flow_end)
        set_subpart!(sfu, :vcu, flow_end+1:dyvar_end)
        set_subpart!(sfu, :pcu, dyvar_end+1:param_end)

        # lucuu, lucucu, exponent
        # LUCUs = Vector{Tuple{Int, Int, Int}}()
        lucuu = Vector{Int}()
        lucucu = Vector{Int}()
        exp = Vector{Int}()
        LUCU_count = 0


        exp_dict = Dict{Symbol, Dict{Symbol, Int}}()

        for (stock, expr) in enumerate(sunits)
            exponent_count = filter!(((k,v),) -> (v != 0), extract_exponents(expr))
            for (unit, val) in exponent_count
                if unit == :UNITLESS
                    exponent_count = Dict{Symbol, Int}()
                    break
                end

                push!(lucuu, units[unit])
                push!(lucucu, stock)
                push!(exp, val)
                LUCU_count += 1
                # push!(LUCUs, (units[unit], stock, exponent_count))
            end
            push!(exp_dict, sname(sfu, stock) => exponent_count)
        end

        for (param, expr) in enumerate(punits)
            exponent_count = filter!(((k,v),) -> v != 0, extract_exponents(expr))
            for (unit, val) in exponent_count
                if unit == :UNITLESS
                    exponent_count = Dict{Symbol, Int}()
                    break
                end
                push!(lucuu, units[unit])
                push!(lucucu, param + dyvar_end)
                push!(exp, val)

                LUCU_count += 1

                # push!(LUCUs, (units[unit], dyvar_end + 1+ param, exponent_count))
            end
            push!(exp_dict, pname(sfu, param) => exponent_count)

        end


        for sum ∈ 1:nsv(sfu)
            incoming_stocks = stockssv(sfu, sum)
            stocks = map(x -> sname(sfu, x), incoming_stocks)
            stock_units = map(x -> exp_dict[x], stocks)

    
            # stock_units = map( x -> subpart(sfu,x,:scu), incoming_stocks)
            @assert allequal(stock_units)
            if length(stock_units) > 0
                for (s, v) in stock_units[1]
                    push!(lucuu, units[s])
                    push!(lucucu, sum + stock_end)
                    push!(exp, v)
                    LUCU_count += 1
                end
                push!(exp_dict, svname(sfu, sum) => stock_units[1])
            else
                @assert false "Don't know how to deal with empty sums for units!"
            #     push!(sum_units, :())
            end
        end


        completed_dyvars = Set{Symbol}()
        while length(completed_dyvars) != nvb(sfu)
            for (i, dv) in enumerate(new_dyvars)
                dyvar, pop = dv
                p, op = pop
                if dyvar in completed_dyvars
                    continue
                end
                if typeof(p) == Symbol && p in keys(exp_dict)
                    for (unit, exponent) in exp_dict[p]
                        push!(lucuu, units[unit])
                        push!(lucucu, i + flow_end)
                        push!(exp, exponent)
                        LUCU_count += 1
                    end
                    push!(exp_dict, dyvar => exp_dict[p])
                    push!(completed_dyvars, dyvar)
                end
                if length(p) == 2 && p[1] in keys(exp_dict) && p[2] in keys(exp_dict)
                    p1, p2 = p
                    if (op == :+) || (op == :-)
                        @assert exp_dict[p1] == exp_dict[p2]
                        for (unit, exponent) in exp_dict[p1]
                            push!(lucuu, units[unit])
                            push!(lucucu, i + flow_end)
                            push!(exp, exponent)
                            LUCU_count += 1
                        end
                        push!(exp_dict, dyvar => exp_dict[p])
                        push!(completed_dyvars, dyvar)

                    elseif (op == :*)
                        pair_units = merge(+, exp_dict[p1], exp_dict[p2])
                        for (unit, exponent) in pair_units
                            if exponent == 0
                                continue
                            end
                            push!(lucuu, units[unit])
                            push!(lucucu, i + flow_end)
                            push!(exp, exponent)
                            LUCU_count += 1
                        end
                        push!(exp_dict, dyvar => pair_units)
                        push!(completed_dyvars, dyvar)

                    elseif (op == :/)
                        pair_units = merge(+, exp_dict[p1], Dict(k => -v for (k, v) in exp_dict[p2]))
                        for (unit, exponent) in pair_units
                            if exponent == 0
                                continue
                            end
                            push!(lucuu, units[unit])
                            push!(lucucu, i + flow_end)
                            push!(exp, exponent)
                            LUCU_count += 1
                        end
                        push!(exp_dict, dyvar => pair_units)
                        push!(completed_dyvars, dyvar)

                    end
                    #ignoring other operations for now
                end
            end
        end



        for flow in 1:nf(sfu)
            for (unit, exponent) in exp_dict[vname(sfu, fv(sfu, flow))]
                if exponent == 0
                    continue
                end
                push!(lucuu, units[unit])

                push!(lucucu, flow + sum_end)

                push!(exp, exponent)
                LUCU_count += 1
            end
            push!(exp_dict, fname(sfu, flow) => exp_dict[vname(sfu, fv(sfu, flow))])
        end



        # println(units)
        # println(exp)
        # println(lucucu)
        # println(lucuu)

        add_parts!(sfu, :LUCU, LUCU_count; exp=exp, lucucu=lucucu, lucuu=lucuu)






        # add_parts!(sfu, :LCUS, ns(sfu))
        # add_parts!(sfu, :LCUV, nv(sfu))
        # add_parts!(sfu, :LCUSV, nsv(sfu))
        # add_parts!(sfu, :LCUF, n)
      return sfu
    end
end








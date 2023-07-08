module Syntax2


export @stock_and_flow_string

# export @stock_and_flow_manual


using MLStyle
using StockFlow
using StockFlow.Syntax
using Catlab
using Catlab.CategoricalAlgebra
# using Syntax

state_dict(n) = Dict(s=>i for (i, s) in enumerate(n))




# macro stock_and_flow_manual(block)

#     sf = StockAndFlowF()

#     Base.remove_linenums!(block)

#     LS::Vector{Symbol} = []
#     SV::Vector{Symbol} = []
#     LSV::Vector{Symbol} = []
#     S::Vector{Symbol} = []
#     LV::Vector{Symbol} = []
#     I::Vector{Symbol} = []
#     O::Vector{Symbol} = []
#     F::Vector{Symbol} = []
#     V::Vector{Symbol} = []

#     LPV::Vector{Symbol} = []
#     P::Vector{Symbol} = []
#     LVV::Vector{Symbol} = []

#     vop::Vector{Symbol} = []

#     current_phase = (_, _) -> ()

#     for statement in block.args
#         @match statement begin
#             QuoteNode(:LS) => begin
#                 current_phase = ls -> push!(LS, ls)
#             end
#             QuoteNode(:SV) => begin
#                 current_phase = sv -> push!(SV, sv)
#             end
#             QuoteNode(:LSV) => begin
#                 current_phase = lsv -> push!(LSV, lsv)
#             end
#             QuoteNode(:S) => begin
#                 current_phase = s -> push!(S, s)
#             end
#             QuoteNode(:LV) => begin
#                 current_phase = lv -> push!(LV, lv)
#             end
#             QuoteNode(:I) => begin
#                 current_phase = i -> push!(I, i)
#             end
#             QuoteNode(:O) => begin
#                 current_phase = o -> push!(O, o)
#             end
#             QuoteNode(:F) => begin
#                 current_phase = f -> push!(F, f)
#             end
#             QuoteNode(:V) => begin
#                 current_phase = v -> push!(V, v)
#             end
#             QuoteNode(:LPV) => begin
#                 current_phase = lpv -> push!(LPV, lpv)
#             end
#             QuoteNode(:P) => begin
#                 current_phase = p -> push!(P, p)
#             end
#             QuoteNode(:LVV) => begin
#                 current_phase = lvv -> push!(LVV, lvv)
#             end
#             QuoteNode(:vop) => begin
#                 current_phase = vop_var -> push!(vop, vop_var)
#             end
            
            
#             QuoteNode(kw) =>
#                 throw("Unknown block type for Stock and Flow syntax: " * String(kw))
#             _ => current_phase(statement)
#         end
#     end

#     str_to_vector = Dict("LS" => LS, "SV" => SV, "LSV" => LSV, "S" => S, "LV" => LV, "I" => I, "O" => O, "F" => F, "V" => V, "LPV" => LPV, "P" => P, "LVV" => LVV)

#     # for part in ["LS", "SV", "LSV", "S", "LV", "I", "O", "F", "V", "LPV", "P", "LVV"]
#     for part in keys(str_to_vector)
#         # println("TYO")
#         n = length(str_to_vector[part])
#         names = str_to_vector[part]
#         # names = map(string, str_to_vector[part])
#         if (part == "S") 
#             # println("DaANC")
#             add_stocks!(sf, n, sname=names)
#         elseif (part == "SV")
#             add_parts!(sf, Symbol(part), n, svname=names)
#         elseif (part == "F")
#             add_parts!(sf, Symbol(part), n, fname=names)
#         elseif (part == "V")
#             add_parts!(sf, Symbol(part), n, vname=names, vop=vop) #vop needs to be in right order
#         else
#             add_parts!(sf, Symbol(part), n)
#         end
#     end



#     return sf

# end


# function foo() println("bar") end
# foo()





function open_feet(stockflow, block)
    Base.remove_linenums!(block)
    feet::Vector{StockAndFlow0} = parse_open_feet_syntax(block.args)
    return Open(stockflow, feet...)
end



macro open_feet(stockflow, block)
    Base.remove_linenums!(block)
    feet::Vector{StockAndFlow0} = parse_open_feet_syntax(block.args)
    return Open(eval(stockflow), feet...)
 end
 
     
 
 
 function parse_open_feet_syntax(statements::Vector{Any})
     # stockflows::Vector{Union{StockAndFlowStructure,StockAndFlow,StockAndFlowF}} = [] # should only be 1
     feet::Vector{StockAndFlow0} = []
     current_phase = (_,_) -> ()
     for statement in statements
         @match statement begin
             QuoteNode(:feet) => begin
                 current_phase = fe -> parse_foot!(feet, fe)
             end
             QuoteNode(kw) =>
                  throw("Unknown block type for open feet syntax: " * String(kw))
             _ => current_phase(statement)
         end
     end
     # # TODO: We only want there to be one stock and flow, right?
     # # s = OpenFeetBlock(first(stockflows), feet)
     return feet
     # return s 
 end
 
 # struct OpenFeetBlock
 #     stockflow::Union{StockAndFlowStructure,StockAndFlow,StockAndFlowF}
 #     feet::Vector{StockAndFlow0}
 # end
 
 
 # function parse_stockflow!(stockflows::Vector{Union{StockAndFlowStructure,StockAndFlow,StockAndFlowF}}, stockflow::Symbol) #TODO: Confirm stockflow should be a symbol rather than an Expr
 #     push!(stockflows, eval(stockflow)) #TODO: Confirm if this is the correct way to do this
 #     # In Python, it was drilled into my head that eval is evil; not sure if this extends to Julia
 # end
 
 function parse_foot!(feet, fo::Expr) # ::Vector{StockAndFlow0}
     println(fo)
     @match fo begin
         :($f1 => $f2) => begin
             println(f1)
             println(f2)
             println(:f1)
             println(:f2)
             push!(feet, foot(Symbol(f1), Symbol(f2), Symbol(f1) => Symbol(f2)))
         end
         _ => println("Fail!")
     end
     
 end
 
         # Expr(c, _, _) || Expr(c, _, _, _) =>
         # throw("Unhandled expression in foot definition " * String(c))
     
 
 # macro stock_and_flow_feet(block)
 #     Base.remove_linenums!(block)
 #     syntax_lines = parse_stock_and_flow_syntax(block.args)
 #     saff_args = stock_and_flow_syntax_to_arguments(syntax_lines)
 #     return Open(
 #         saff_args.stocks,
 #         saff_args.params,
 #         saff_args.dyvars,
 #         saff_args.flows,
 #         saff_args.sums,
 #     )
 # end
 
 
 # function parse_stock_and_flow_feet_syntax(statements::Vector{Any})
 #     stocks::Vector{Symbol} = []
 #     params::Vector{Symbol} = []
 #     dyvars::Vector{Tuple{Symbol,Expr}} = []
 #     flows::Vector{Tuple{Symbol,Expr,Symbol}} = []
 #     sums::Vector{Tuple{Symbol,Vector{Symbol}}} = []
 #     feet::Vector{Tuple{Symbol, Symbol}} = []
 #     current_phase = (_, _) -> ()
 #     for statement in statements
 #         @match statement begin
 #             QuoteNode(:stocks) => begin
 #                 current_phase = s -> parse_stock!(stocks, s)
 #             end
 #             QuoteNode(:parameters) => begin
 #                 current_phase = p -> parse_param!(params, p)
 #             end
 #             QuoteNode(:dynamic_variables) => begin
 #                 current_phase = d -> parse_dyvar!(dyvars, d)
 #             end
 #             QuoteNode(:flows) => begin
 #                 current_phase = f -> parse_flow!(flows, f)
 #             end
 #             QuoteNode(:sums) => begin
 #                 current_phase = s -> parse_sum!(sums, s)
 #             end
 #             QuoteNode(:feet) => begin
 #                 # current_phase = fe -> parse_feet!(feet, fe)
 #             end
 #             QuoteNode(kw) =>
 #                 throw("Unknown block type for Stock and Flow syntax: " * String(kw))
 #             _ => current_phase(statement)
 #         end
 #     end
 
 #     s = StockAndFlowBlock(stocks, params, dyvars, flows, sums)
 #     # leg = 
 #     return s
 # end
 
 
 
 # function baz() end
 
 
 
 """
Macro similar to stock_and_flow, can be additionally given a vector of strings to replace \$v1, ..., \$vn in the block argument
"""
macro stock_and_flow_string(replacestrings, block)
    Base.remove_linenums!(block)


    replacestrings = replacestrings.args
    # println(block.args)



    # println(string(block.args))


    # for i in block.args
    #     println(typeof(i))
    # end

    # expression_string::String = string(block.args)
    
    

    n = length(replacestrings)
    pairs = ["v$i" => replacestrings[n-i+1] for i in 1:n] # replacing vn with final string, ..., v1 with first.
    # Doing this in reverse order so v10 matches before v1.


    new_expression_vector = []
    for line in block.args
        type_of_current_line = typeof(line)
        current_line = string(line)
        replaced_expr = replace(current_line, pairs...)
        println(replaced_expr)
        new_expression = type_of_current_line(replaced_expr)
        push!(new_expression_vector, new_expression)
        # println(line)
        # println(typeof(line))
        # push!(new_expression_vector, Expr(replace(string(line), pairs...)))
    end



    syntax_lines = parse_stock_and_flow_syntax(new_expression_vector)
    saff_args = stock_and_flow_syntax_to_arguments(syntax_lines)
    return StockAndFlowF(
        saff_args.stocks,
        saff_args.params,
        saff_args.dyvars,
        saff_args.flows,
        saff_args.sums,
    )
end
    

 



end
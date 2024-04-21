module Rewrite

export sfrewrite, @rewrite

# regex to search for uncommented print:  [^#]\s*@show

using ...StockFlow

using ..Syntax
import ..Syntax: parse_dyvar, parse_flow, parse_sum, parse_stock, 
parse_param, sf_to_block, StockAndFlowBlock, stock_and_flow_syntax_to_arguments,
StockArgT, ParamArgT, create_flow_definitions, Binop, Ref, get


using Catlab.ACSets.ACSetInterface
using Catlab.CategoricalAlgebra
using AlgebraicRewriting
using AlgebraicRewriting: rewrite


using MLStyle
using MLStyle.Modules.AST

"""
Convert a stockflow block to a stockflow
"""
function block_to_sf(block)
  args = stock_and_flow_syntax_to_arguments(block)
  sf = StockAndFlowF(args.stocks, args.params,  
    map(kv -> kv.first => StockFlow.Syntax.get(kv.second), args.dyvars), 
    args.flows, args.sums)
  return sf
end


struct SFPointer
  type::Symbol
  index::Int
end





function add_operand_link!(sf, operand, operand_type, dyvar_index, operand_index)
  # @show sf, operand, operand_type, dyvar_index, operand_index
  if operand_type == :S
    # @show sf, operand
    stock_index = findfirst(==(operand), snames(sf))
    add_part!(sf, :LV ; lvs = stock_index, lvv = dyvar_index, lvsposition = operand_index)
  elseif operand_type == :P
    param_index = findfirst(==(operand), pnames(sf))
    add_part!(sf, :LPV ; lpvp = param_index, lpvv = dyvar_index, lpvpposition = operand_index)
  elseif operand_type == :V
    dyvar2_index = findfirst(==(operand), vnames(sf))
    add_part!(sf, :LVV ; lvsrc = dyvar2_index, lvtgt = dyvar_index, lvsrcposition = operand_index)
  elseif operand_type == :SV
    sum_index = findfirst(==(operand), svnames(sf))
    add_part!(sf, :LSV ; lsvsv = sum_index, lsvv = dyvar_index, lsvsvposition = operand_index)
  end
end

function add_object_of_type!(sf, object, object_type, object_index, sf_block)
  if object_type == :S
    add_stock!(sf ; sname = object)
  elseif object_type == :P
    add_parameter!(sf ; pname = object)
  elseif object_type == :V
    dyvar_op = sf_block.dyvars[object_index][2].args[1]
    add_variable!(sf ; vname = object, vop = dyvar_op)
  elseif object_type == :SV
    add_svariable!(sf ; svname = object)
  elseif object_type == :F
    flow_dyvar = sf_block.flows[object_index][2].args[2]
    flow_dyvar_index = findfirst(==(flow_dyvar), vnames(sf))
    @assert !(isnothing(flow_dyvar_index)) "Tried adding a flow before its dyvar!"
    add_flow!(sf, flow_dyvar_index ; fname = object)
  end
end


function remove_from_sums2!(stock_name, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  for sum ∈ sf_block.sums
    sum_name = sum[1]
    sum_stocks = sum[2]
    # for loop because it's technically possible for a stock
    # to be linked to a sum twice.
    for sum_stock in sum_stocks
      if sum_stock == stock_name
        stock_index = name_dict[stock_name].index
        if sum_name ∉ L_set
          add_part!(L, :SV ; svname = sum_name)
          push!(L_set, sum_name)
          # add_part!(L, :LS ; lvs = stock_index, lvv = nsv(L))
          sum_index = nsv(L)
        else
          sum_index = findfirst(==(sum_name), svnames(L))
        end
        # add_part!(L, :LS ; lvs = stock_index, lvv = sum_index)
        # ! WE'RE ASSUMING THERE IS AT MAX ONE LINK BETWEEN A SUM AND A STOCK
        link = (stock_name, sum_name)
        # link = (stock_index, sum_index)
        if !(link in L_connect_dict[:LS])
          push!(L_connect_dict[:LS], link)
        end
        if !(link in remove_connect_dict[:LS])
          push!(remove_connect_dict[:LS], link)
        end
      end
    end
  end
end


function dyvar_link_name_from_object_type(object_type)
  if object_type == :S
    return :LV
  elseif object_type == :SV
    return :LSV
  elseif object_type == :P
    return :LPV
  elseif object_type == :V
    return :LVV
  end
end

function add_links_from_dict!(sf, connect_dict)
  # @show connect_dict
  for ((lvs, lvv, lvsposition)) in connect_dict[:LV]
    # @show sf
    lvs = findfirst(==(lvs), snames(sf))
    lvv = findfirst(==(lvv), vnames(sf))

    add_part!(sf, :LV; lvs = lvs, lvv = lvv, lvsposition = lvsposition)
  end
  for ((lpvp, lpvv, lpvpposition)) in connect_dict[:LPV]
    lpvp = findfirst(==(lpvp), pnames(sf))
    lpvv = findfirst(==(lpvv), vnames(sf))

    # @show sf

    add_part!(sf, :LPV; lpvp = lpvp, lpvv = lpvv, lpvpposition = lpvpposition)
  end
  for ((lsvsv, lsvv, lsvsvposition)) in connect_dict[:LSV]
    lsvsv = findfirst(==(lsvsv), svnames(sf))
    lsvv = findfirst(==(lsvv), vnames(sf))

    add_part!(sf, :LSV; lsvsv = lsvsv, lsvv = lsvv, lsvsvposition = lsvsvposition)
  end
  for ((lvsrc, lvtgt, lvsrcposition)) in connect_dict[:LVV]

    lvsrc = findfirst(==(lvsrc), vnames(sf))
    lvtgt = findfirst(==(lvtgt), vnames(sf))

    add_part!(sf, :LVV; lvsrc = lvsrc, lvtgt = lvtgt, lvsrcposition = lvsrcposition)
  end
  for ((lss, lssv)) in connect_dict[:LS]
    # @show sf
    # @show lss, lssv
    lss = findfirst(==(lss), snames(sf))
    lssv = findfirst(==(lssv), svnames(sf))

    add_part!(sf, :LS; lss = lss, lssv = lssv)
  end
  for ((is, ifn)) in connect_dict[:I]
    is = findfirst(==(is), snames(sf))
    ifn = findfirst(==(ifn), fnames(sf))

    add_part!(sf, :I; is = is, ifn = ifn)
  end
  for ((os, ofn)) in connect_dict[:O]
    os = findfirst(==(os), snames(sf))
    ofn = findfirst(==(ofn), fnames(sf))
    add_part!(sf, :O; os = os, ofn = ofn)
  end
end

function remove_from_dyvars2!(object_name, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  for dyvar ∈ sf_block.dyvars
    dyvar_name = dyvar[1]
    dyvar_expression = dyvar[2].args
    dyvar_op = dyvar_expression[1]
    for (operand_index, operand) in enumerate(dyvar_expression[2:end])
      if operand == object_name

        object_type = name_dict[object_name].type
      

        if dyvar_name ∉ L_set
          add_part!(L, :V ; vname = dyvar_name, vop = dyvar_op)
          push!(L_set, dyvar_name)
          dyvar_index = nvb(L)
        else
          dyvar_index = findfirst(==(dyvar_name), vnames(L))
        end

        # add_operand_link!(L, object_name, object_type, dyvar_index, operand_index)
        link_name = dyvar_link_name_from_object_type(object_type)
        link = (object_name, dyvar_name, operand_index)
        if !(link in L_connect_dict[link_name])
          push!(L_connect_dict[link_name], link)
        end

        if !(link in remove_connect_dict[link_name])
          push!(remove_connect_dict[link_name], link)
        end


      end # if
    end # for oper
  end # for dyvar
end #func



function add_redefintions!(L, L_redef_queue, R_dyvar_queue, R_sum_queue, R_flow_queue, name_dict, L_set, sf_block, L_connect_dict, remove_connect_dict, removed_set)
  for object in L_redef_queue
    if length(object.args) == 2
      object_name = object.args[1]
    elseif length(object.args) == 3 # flow
      object_name = object.args[3].args[2].args[1]
    end
    object_pointer = name_dict[object_name]
    object_type = object_pointer.type
    object_index = object_pointer.index

    if object_type == :V
      # @show L_set

      dyvar_operands = object.args[2].args
      dyvar_op = dyvar_operands[1]
      push!(R_dyvar_queue, object)
      if !(object_name in L_set) 
        push!(L_set, object_name)
        original_dyvar_op = sf_block.dyvars[name_dict[object_name].index][2].args[1]
        add_variable!(L ; vname = object_name, vop = original_dyvar_op)
        # dyvar_index = nvb(L)
      end
      dyvar_index = findfirst(==(object_name), vnames(L))

      # @show "SDSA"
      for (operand_index, operand) in enumerate(sf_block.dyvars[object_index][2].args[2:end])

        if operand in keys(name_dict)
          # need to add all links to L
          operand_pointer = name_dict[operand]
          operand_type = operand_pointer.type
          operand_original_index = operand_pointer.index

          # @show L_set
          # @show L
          if !(operand in L_set)
            # note, IT IS NOT REMOVED!
            # It could be removed separately, but as is, this should
            # place it in I

            # TODO: User's responsibility to add the values in the redef under
            # the correct header.  If we have a redef N = [S, I],
            # user needs to add S and I under a :stocks
            # Is this reasonable?

            push!(L_set, operand)

            add_object_of_type!(L, operand, operand_type, operand_original_index, sf_block)


          end # if !operand

          # only want to add link if there was one in that position originally,
          # hence the continue.

          original_def = sf_block.dyvars[name_dict[object_name].index]
          # if original_def[2].args[operand_index + 1] != operand
          #   continue
          # end

          # NEW
          # we dont want to add a duplicate if the link is already going to be added from the connect dict
          dyvar_link_type = dyvar_link_name_from_object_type(operand_type)
          # if operand_type == :S
          #   op_real_index = findfirst(==(operand), snames(L))
          # elseif operand_type == :P
          #   op_real_index = findfirst(==(operand), pnames(L))
          # elseif operand_type == :V
          #   op_real_index = findfirst(==(operand), vnames(L))
          # elseif operand_type == :SV
          #   op_real_index = findfirst(==(operand), svnames(L))
          # end
          # link = (op_real_index => dyvar_index, operand_index)
          link = (operand, object_name, operand_index)
          # @show "here", link
          # @show L_connect_dict
          if !(link in L_connect_dict[dyvar_link_type])
            # @show "HDASDASDSA"
            push!(L_connect_dict[dyvar_link_type], link)
          end
          # we want to remove all links such that they aren't in the redefinition
          # operand is in old, we want to check if it's in the new

          # TODO: wrong
          if (!(link in remove_connect_dict[dyvar_link_type])) && (length(object.args[2].args[2:end]) < operand_index || object.args[2].args[operand_index + 1] != operand )
            push!(remove_connect_dict[dyvar_link_type], link)
          end
          # @show L_connect_dict, remove_connect_dict
          # END NEW
          # add_operand_link!(L, operand, operand_type, dyvar_index, operand_index)

        end # if operand in namedict
      end # for each operand (except operator)

    elseif object_type == :SV
      if !(object_name in L_set)
        push!(L_set, object_name)
        add_svariable!(L ; svname = object_name)
      end
      push!(R_sum_queue, object)

      sum_index = findfirst(==(object_name), svnames(L))
      for operand in object.args[2].args
        if operand in keys(name_dict)
          if !(operand in L_set)
            push!(L_set, operand)
            add_stock!(L ; sname = operand)
          end
          
          stock_index = findfirst(==(operand), snames(L))
          link = (operand, object_name)
          if !(link in L_connect_dict[:LS])
            push!(L_connect_dict[:LS], link)
          end
          # Don't remove link if it's retained
          # all_old_links = zip(get_lss())
          original_links = sf_block.sums[name_dict[object_name].index][2]

          if !(link in remove_connect_dict[:LS]) && !(operand in original_links)
            push!(remove_connect_dict[:LS], link)
          end          
          # add_part!(L, :LS ; lss = stock_index, lssv = sum_index)
        end
      end
    elseif object_type == :F
      # If the dyvar has changed, need to have it in L, nothing in I, re-add new to R
      original_flow_var = sf_block.flows[name_dict[object_name].index][2].args[2]    # flow_var = object.args[3].args[2].args[2]
      # @show 
      if !(original_flow_var in L_set)
        push!(L_set, original_flow_var)
        dyvar_op = sf_block.dyvars[name_dict[original_flow_var].index][2].args[1]
        add_variable!(L ; vname = original_flow_var, vop = dyvar_op)
      end
      if !(object_name in L_set)
        # @show vnames(L), original_flow_var
        flow_var_index = findfirst(==(original_flow_var), vnames(L))
        push!(L_set, object_name)
        add_flow!(L, flow_var_index ; fname = object_name)
      end

      flow_index = findfirst(==(object_name), fnames(L))

      original_inflow = sf_block.flows[flow_index][3]
      # @show original_inflow
      original_outflow = sf_block.flows[flow_index][1]


      # original_inflow = object.args[2]
      # original_outflow = object.args[3].args[3]

      if (original_inflow != :F_NONE)
        # @show original_outflow
        if !(original_inflow in L_set)
          push!(L_set, original_inflow)
          add_stock!(L ; sname = original_inflow)
        end
        link = (original_inflow, object_name)
        if !(link in L_connect_dict[:I])
          # @show "YES"
          push!(L_connect_dict[:I], link)
        end
        if !(link in remove_connect_dict[:I])
          push!(remove_connect_dict[:I], link)
        end         
        # add_part!(L, :I ; ifn = flow_index, is = findfirst(==(original_inflow), snames(L)))
      end 

      if (original_outflow != :F_NONE)
        if !(original_outflow in L_set)
          push!(L_set, original_outflow)
          add_stock!(L ; sname = original_outflow)
        end
        # @show original_outflow, snames(L)
        link = (original_outflow, object_name)
        if !(link in L_connect_dict[:O])
          push!(L_connect_dict[:O], link)
        end
        if !(link in remove_connect_dict[:O])
          push!(remove_connect_dict[:O], link)
        end         
        # add_part!(L, :O ; ofn = flow_index, os = findfirst(==(original_outflow), snames(L)))
      end 

      # different dynamic variables
      if !(original_flow_var == object.args[3].args[2].args[2])
        # @show L_connect_dict
        # @show L
        push!(removed_set, object_name)

        if (original_inflow != :F_NONE)
          link = (original_inflow, object_name)
          # @show link, L_connect_dict[:I]
          if !(link in L_connect_dict[:I])
            push!(L_connect_dict[:I], link)
          end
          if !(link in remove_connect_dict[:I])
            push!(remove_connect_dict[:I], link)
          end
        end
        if (original_outflow != :F_NONE)
          link = (original_outflow, object_name)
          if !(link in L_connect_dict[:O])
            push!(L_connect_dict[:O], link)
          end
          if !(link in remove_connect_dict[:O])
            push!(remove_connect_dict[:O], link)
          end
        end
      end

      # needs to have dyvar added...
      push!(R_flow_queue, object)
      # end


    end 


  end
end




# TODO
function remove_from_flows!(object_name, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  for flow in sf_block.flows
    inflow_stock = flow[3]
    outflow_stock = flow[1]
    flow_name = flow[2].args[1]
    flow_dyvar = flow[2].args[2]

    link = (object_name, flow_name)

    if (inflow_stock == object_name) || (outflow_stock == object_name)
      if !(flow_name in L_set)
        push!(L_set, flow_name)
        if !(flow_dyvar in L_set)
          dyvar_op = sf_block.dyvars[name_dict[flow_dyvar].index].args[1]
          push!(L_set, flow_dyvar)
          add_variable!(L ; vname = flow_dyvar, vop = dyvar_op)
        end
        add_flow!(L, findfirst(==(flow_dyvar), vnames(L)) ; fname = flow_name, )
      end
    end

    if (inflow_stock == object_name)
      if !(link in L_connect_dict[:I])
        push!(L_connect_dict[:I], link)
      end
      if !(link in remove_connect_dict[:I])
        push!(remove_connect_dict[:I], link)
      end
    end

    if (outflow_stock == object_name)
      if !(link in L_connect_dict[:O])
        push!(L_connect_dict[:O], link)
      end
      if !(link in remove_connect_dict[:O])
        push!(remove_connect_dict[:O], link)
      end
    end


  end

end


function remove_block!(removed, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
  
  object_pointer = name_dict[removed]
  object_type = object_pointer.type
  push!(removed_set, removed)

  if !(removed in L_set)
    push!(L_set, removed)
    add_object_of_type!(L, removed, object_type, object_pointer.index, sf_block)
  end





  #TODO: factor out

  # @show "u r here", object_type


  if object_type == :S
    object_definition = sf_block.stocks[object_pointer.index]
    # add_part!(L, :S ; sname = object_definition.val)
    remove_from_sums2!(object_definition.val, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
    remove_from_dyvars2!(object_definition.val, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
    remove_from_flows!(removed, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
    

  elseif object_type == :P
    object_definition = sf_block.params[object_pointer.index]
    # add_part!(L, :P ; pname = object_definition.val)
    remove_from_dyvars2!(object_definition.val, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  elseif object_type == :V
    # @show "I can do the calculus."
    object_definition = sf_block.dyvars[object_pointer.index]
    dyvar_op = object_definition[2].args[1]
    vname = object_definition[1]

    for (operand_index, operand) in enumerate(object_definition[2].args[2:end])
      if !(operand in L_set)
        push!(L_set, operand)
        add_object_of_type!(L, operand, name_dict[operand].type, name_dict[operand].index, sf_block)
      end
      operand_link_type = dyvar_link_name_from_object_type(name_dict[operand].type)
      link = (operand, removed, operand_index)
      if !(link in L_connect_dict[operand_link_type])
        push!(L_connect_dict[operand_link_type], link)
      end
      if !(link in remove_connect_dict[operand_link_type])
        push!(remove_connect_dict[operand_link_type], link)
      end
    end

    # add_part!(L, :V ; vname = vname, vop = dyvar_op)
    # @show L
    remove_from_dyvars2!(vname, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
    # @show L
  elseif object_type == :SV
    object_definition = sf_block.sums[object_pointer.index]
    svname = object_definition[1]

    for sum_stock in object_definition[2]
      if !(sum_stock in L_set)
        push!(L_set, sum_stock)
        add_stock!(L ; sname = sum_stock)
      end
      link = (sum_stock, removed)
      if !(link in L_connect_dict[:LS])
        push!(L_connect_dict[:LS], link)
      end
      if !(link in remove_connect_dict[:LS])
        push!(remove_connect_dict[:LS], link)
      end
    end

    # add_part!(L, :SV ; svname = svname)
    remove_from_dyvars2!(svname, sf_block, L, L_set, name_dict, L_connect_dict, remove_connect_dict)
  elseif object_type == :F
    flow = sf_block.flows[object_pointer.index]
    inflow_stock = flow[3]
    outflow_stock = flow[1]
    flow_name = flow[2].args[1]
    flow_dyvar = flow[2].args[2]

    inflow_link = (inflow_stock, flow_name)
    outflow_link = (outflow_stock, flow_name)

    if !(flow_dyvar in L_set)
      fv = findfirst(==(flow_dyvar), vnames(L)).args[1]
      push!(L_set, flow_dyvar)
      add_variable(L ; vname = flow_dyvar, fv = fv)
    end

    if !(inflow_stock in L_set)
      push!(L_set, inflow_stock)
      add_stock(L ; sname = inflow_stock)
    end


    if !(outflow_stock in L_set)
      push!(L_set, inflow_stock)
      add_stock(L ; sname = inflow_stock)
    end


    if !(inflow_link in L_connect_dict[:I])
      push!(L_connect_dict[:I], inflow_link)
    end


    if !(inflow_link in remove_connect_dict[:I])
      push!(remove_connect_dict[:I], inflow_link)
    end

    if !(outflow_link in L_connect_dict[:O])
      push!(L_connect_dict[:O], outflow_link)
    end

    if !(outflow_link in remove_connect_dict[:O])
      push!(remove_connect_dict[:O], outflow_link)
    end


    # push!(L_connect_dict)


    # link = (object_name, flow_name)
    

  end
end



function remove_links_block!(l, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
  # src = nothing
  # tgt = nothing
  # pos = 0
  @match l begin
    
    :($(src::Symbol) => $(tgt::Symbol)) => begin
      tgt_type = name_dict[tgt].type

      if !(tgt in L_set)
        push!(L_set, tgt)
        add_object_of_type!(L, tgt, tgt_type, name_dict[tgt].index, sf_block)
      end
      if !(src in L_set)
        push!(L_set, src)
        add_object_of_type!(L, src, name_dict[src].type, name_dict[src].index, sf_block)
      end

      if tgt_type == :SV
        # link_number = count(x -> x == src, sf_block[name_dict[tgt].index][2])
        if !((src, tgt) in L_connect_dict[:LS])
          push(L_connect_dict[:LS], (src, tgt))
        end
        push!(remove_connect_dict[:LS], (src, tgt))
      elseif tgt_type == :V # tgt must be dyvar
        positions = findall(==(src), sf_block.dyvars[name_dict[tgt].index][2].args[2:end])
        link_type = dyvar_link_name_from_object_type(name_dict[src].type)
        for pos in positions
          if !((src, tgt, pos) in L_connect_dict[link_type])
            push!(L_connect_dict[link_type], (src, tgt, pos))
          end
          push!(remove_connect_dict[link_type], (src, tgt, pos))
        end
      elseif tgt_type == :F
        tgt_pointer = name_dict[tgt]
        flow_definintion = sf_block.flows[tgt_pointer.index]
        flow_dyvar = flow_definintion[2].args[2]
        if !(flow_dyvar in L_set)
          flow_dyvar_op = sf_block.dyvars[name_dict[flow_dyvar].index][2].args[1]
          push!(L_set, flow_dyvar)
          add_variable!(L ; vname = flow_dyvar, vop = flow_dyvar_op)
        end
        if flow_definintion[3] == src
          if !((src, tgt) in L_connect_dict[:I])
            push!(L_connect_dict[:I], (src, tgt))
          end
          push!(remove_connect_dict[:I], (src, tgt))
        end
        if flow_definintion[1] == src
          if !((src, tgt) in L_connect_dict[:O])
            push!(L_connect_dict[:O], (src, tgt))
          end
          push!(remove_connect_dict[:O], (src, tgt))
        end
      end # if tgt_type == 
    end # :($(src:: ... => begin

    :($(src::Symbol) => $(tgt::Symbol), $(position::Int)) => begin
      tgt_type = name_dict[tgt].type

      if !(tgt in L_set)
        push!(L_set, tgt)
        add_object_of_type!(L, tgt, tgt_type, name_dict[tgt].index, sf_block)
      end
      if !(src in L_set)
        push!(L_set, src)
        add_object_of_type!(L, src, name_dict[src].type, name_dict[src].index, sf_block)
      end

      if tgt_type == :SV
        if !((src, tgt) in L_connect_dict[:LS])
          push(L_connect_dict[:LS], (src, tgt))
        end
        push(remove_connect_dict[:LS], (src, tgt))
      
      elseif tgt_type == :V
        if !((src, tgt, position) in L_connect_dict[link_type])
          push!(L_connect_dict[link_type], (src, tgt, position))
        end
        push!(remove_connect_dict[link_type], (src, tgt, position))
      elseif tgt_type == :F
        if position == 2
          if !((src, tgt) in L_connect_dict[:I])
            push!(L_connect_dict[:I], (src, tgt))
          end
          push!(remove_connect_dict[:I], (src, tgt))

        elseif position == 1
          if !((src, tgt) in L_connect_dict[:O])
            push!(L_connect_dict[:O], (src, tgt))
          end
          push!(remove_connect_dict[:O], (src, tgt))
        end
      end #if tgt_type              

      # TODO: error case

    end # ... $(position::Int)) => begin
  end # match
end # current_phase



function add_links_block!(l, L_set, name_dict, L, sf_block, R_link_vector)
  src = nothing
  tgt = nothing
  position = nothing
  @match l begin 
    :($(srca::Symbol) => $(tgta::Symbol), $(positiona::Int)) => begin
      src = srca
      tgt = tgta
      position = positiona
    end
    :($(srca::Symbol) => $(tgta::Symbol)) => begin
      src = srca
      tgt = tgta
      position = 0
    end
  end
  if (src in keys(name_dict)) && !(src in L_set)
    add_object_of_type!(L, src, name_dict[src].type, name_dict[src].index, sf_block)
    push!(L_set, src)
  end
  if (tgt in keys(name_dict)) && !(tgt in L_set)
    add_object_of_type!(L, tgt, name_dict[tgt].type, name_dict[tgt].index, sf_block)
    push!(L_set, tgt)
  end
  push!(R_link_vector, (src => tgt, position))      
end 


function dyvar_swaps_block!(dw, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict, R_link_vector)
  capture_dict = @capture $old > $new dw
  new = capture_dict[:new]
  old = capture_dict[:old]

  if !(old in L_set)
    push!(L_set, old)
    add_object_of_type!(L, old, name_dict[old].type, name_dict[old].index, sf_block)
  end

  for ((dyvar_name, dyvar_expr)) in sf_block.dyvars
    dyvar_op = dyvar_expr.args[1]
    dyvar_operands = dyvar_expr.args[2:end]
    matching_indices = findall(==(old), dyvar_operands)
    if !(isempty(matching_indices))
      if !(dyvar_name in L_set)
        push!(L_set, dyvar_name)
        add_variable!(L ; vname = dyvar_name, vop = dyvar_op)
      end

      for index in matching_indices
        old_link = (old, dyvar_name, index)
        # new_link = (new, dyvar_name, index)

        old_link_type = dyvar_link_name_from_object_type(name_dict[old].type)
        if !(old_link in L_connect_dict[old_link_type])
          push!(L_connect_dict[old_link_type], old_link)
        end
        if !(old_link in remove_connect_dict[old_link_type])
          push!(remove_connect_dict[old_link_type], old_link)
        end


        
        if (new in keys(name_dict)) && !(new in L_set)
          push!(L_set, new)
          add_object_of_type!(L, new, name_dict[new].type, name_dict[new].index, sf_block)
        end

        push!(R_link_vector, (new => dyvar_name, index))

      end

      


    end
    # if src in dyvar_expr.args[2:end]

  end
end





function sfrewrite2(sf::K, block::Expr) where {K <: AbstractStockAndFlowF}

  Base.remove_linenums!(block)
  name_vector = [snames(sf) ; svnames(sf) ; vnames(sf) ; fnames(sf) ; pnames(sf)]
  @assert allunique(name_vector) "Not all names are unique!  $(filter(x -> count(y -> y == x, name_vector) >= 2, name_vector))"

  # Unfortunately, this makes it un-extendable to general schemas...it'll only work for StockAndFlowF...
  sf_block::StockAndFlowBlock = sf_to_block(sf)

  L = StockAndFlowF()

  name_dict = Dict{Symbol, SFPointer}([
    [s => SFPointer(:S, i) for (i, s) ∈ enumerate(snames(sf))] ;
    [sv => SFPointer(:SV, i) for (i, sv) ∈ enumerate(svnames(sf))] ;
    [v => SFPointer(:V, i) for (i, v) ∈ enumerate(vnames(sf))] ;
    [f => SFPointer(:F, i) for (i, f) ∈ enumerate(fnames(sf))] ;
    [p => SFPointer(:P, i) for (i, p) ∈ enumerate(pnames(sf))]
  ])

  # symbols added to L
  L_set = Set{Symbol}()

  L_connect_dict = Dict{Symbol, Vector{Tuple}}(
    :LV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LPV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LVV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LSV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LS => Vector{Tuple{Symbol, Symbol}}(),
    :I => Vector{Tuple{Symbol, Symbol}}(),
    :O => Vector{Tuple{Symbol, Symbol}}(),
  )

  remove_connect_dict = Dict{Symbol, Vector{Tuple}}(
    :LV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LPV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LVV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LSV => Vector{Tuple{Symbol, Symbol, Int}}(),
    :LS => Vector{Tuple{Symbol, Symbol}}(),
    :I => Vector{Tuple{Symbol, Symbol}}(),
    :O => Vector{Tuple{Symbol, Symbol}}(),
  )

  L_redef_queue = Vector{Expr}()

  # used to construct I.  Take L_set - removed_set
  removed_set = Set{Symbol}()

  R_stock_queue = Vector{Symbol}()
  R_param_queue = Vector{Symbol}()
  R_sum_queue = Vector{Expr}()
  R_dyvar_queue = Vector{Expr}()
  R_flow_queue = Vector{Expr}()

  # R_name_dict = Dict{Symbol, SFPointer}()

  R_link_vector = Vector{Tuple}()


  current_phase = (_, _) -> ()
  for statement in block.args
    @match statement begin
      QuoteNode(:removes) => begin
        current_phase = removed -> remove_block!(removed, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
      end

      # TODO: figure out if I remove links at this stage.  If I do, then remove_links needs to come after
      # all the parts which are defined in it
      QuoteNode(:remove_links) => begin
       current_phase = l -> remove_links_block!(l, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict)
      end # quotenode

      QuoteNode(:add_links) => begin
        current_phase = l -> add_links_block!(l, L_set, name_dict, L, sf_block, R_link_vector)
      end
      QuoteNode(:dyvar_swaps) => begin
        current_phase = dw -> dyvar_swaps_block!(dw, removed_set, L_set, name_dict, L, sf_block, L_connect_dict, remove_connect_dict, R_link_vector)
      end



      QuoteNode(:redefs) => begin
        current_phase = r -> push!(L_redef_queue, r)
      end
      
      QuoteNode(:stocks) => begin
        current_phase = s -> push!(R_stock_queue, s)
      end
      QuoteNode(:parameters) => begin
        current_phase = p -> push!(R_param_queue, p)
      end
      QuoteNode(:dynamic_variables) => begin
        current_phase = d -> push!(R_dyvar_queue, d)
      end
      QuoteNode(:flows) => begin
        current_phase = f -> push!(R_flow_queue, f)
      end
      QuoteNode(:sums) => begin
        current_phase = s -> push!(R_sum_queue, s)
      end
      QuoteNode(kw) =>
        error("Unknown block type for rewrite syntax: " * String(kw))
      _ => current_phase(statement)

    end # match
  end # for


  # extract links from redefinitions

  #TODO: both redefinitions and other objects can cause objects to be added; need to make sure things are only added once.
  #See current issue with sir -> seir
  add_redefintions!(L, L_redef_queue, R_dyvar_queue, R_sum_queue, R_flow_queue, name_dict, L_set, sf_block, L_connect_dict, remove_connect_dict, removed_set)

  # @show "OK"
  # @show L

  add_links_from_dict!(L, L_connect_dict)
  # @show L


  # @show allunique(L_connect_dict[:LPV]), L_connect_dict[:LPV]
  # @show remove_connect_dict[:LPV]
  # @show removed_set
  # @show L

  # @show "OK"
  # @show L

  for dyvar in R_dyvar_queue
    # @show dyvar
    for operand in dyvar.args[2].args[2:end]
      if !(operand in L_set) && operand in keys(name_dict)
        push!(L_set, operand)
        add_object_of_type!(L, operand, name_dict[operand].type, name_dict[operand].index, sf_block)
      end
    end
  end


  for sum in R_sum_queue
    for stock_sum in sum.args[2].args
      if !(stock_sum in L_set) && stock_sum in keys(name_dict)
        push!(L_set, stock_sum)
        add_object_of_type!(L, stock_sum, name_dict[stock_sum].type, name_dict[stock_sum].index, sf_block)
      end
    end
  end
  for flow in R_flow_queue
    outflow = flow.args[2]
    inflow = flow.args[3].args[3]
    dyvar = flow.args[3].args[2].args[2]
    # @show dyvar
    if dyvar isa Expr
      stack = [dyvar]
      completed = Vector{Symbol}()
      while !(isempty(stack))
        dyvar_it = stack[1]
        for val in dyvar_it.args[2:end]
          if val isa Symbol
            push!(completed, val)
          elseif val isa Expr
            
            push!(stack, val)
          end
          # @show stack
        end
        popfirst!(stack)

      end
      # @show completed

      for e_value in completed
        # @show e_value
        if !(e_value in L_set) && e_value in keys(name_dict)
          push!(L_set, e_value)
          add_object_of_type!(L, e_value, name_dict[e_value].type, name_dict[e_value].index, sf_block)
        end
      end
    end
    for value in [outflow, inflow, dyvar]
      if !(value in L_set) && value in keys(name_dict)
        push!(L_set, value)
        add_object_of_type!(L, value, name_dict[value].type, name_dict[value].index, sf_block)
      end
    end
  end

  # @show L
  # L should now be complete
  
  # @show L_set, removed_set

  I_set = setdiff(L_set, removed_set)
  I_connect_dict = Dict(key => [value for value in filter(v -> !(v in remove_connect_dict[key]), L_connect_dict[key])] for key in keys(L_connect_dict))
  # @show I_connect_dict L_connect_dict remove_connect_dict
  # I_connect_dict = Dict(key => collect(filter(!isnothing, indexin(L_connect_dict[key], remove_connect_dict[key]))) for key in keys(L_connect_dict))
  I = StockAndFlowF()
  
  I_set_flows = filter(x -> name_dict[x].type == :F, I_set)


  # TODO: This isn't preserving any links.  Verify if I need to...
  # probably need to preserve inflows and outflows, at least...
  for object in I_set
    if object in I_set_flows
      continue
    end
    object_pointer = name_dict[object]
    object_type = object_pointer.type
    object_index = object_pointer.index

    add_object_of_type!(I, object, object_type, object_index, sf_block)

  end

  # @show I
  for object in I_set_flows
    add_object_of_type!(I, object, :F, name_dict[object].index, sf_block)
  end

  # @show I_connect_dict
  # @show I
  add_links_from_dict!(I, I_connect_dict)

  # @show I

  R = deepcopy(I)
  R_name_dict = Dict{Symbol, SFPointer}([
    [s => SFPointer(:S, i) for (i, s) ∈ enumerate(snames(R))] ;
    [sv => SFPointer(:SV, i) for (i, sv) ∈ enumerate(svnames(R))] ;
    [v => SFPointer(:V, i) for (i, v) ∈ enumerate(vnames(R))] ;
    [f => SFPointer(:F, i) for (i, f) ∈ enumerate(fnames(R))] ;
    [p => SFPointer(:P, i) for (i, p) ∈ enumerate(pnames(R))]
  ])

  not_yet_added_stocks = collect(filter(x -> !(x in snames(R)), R_stock_queue))
  not_yet_added_params = collect(filter(x -> !(x in pnames(R)), R_param_queue))

  
  # @show R_name_dict

  type_index_counter = ns(R) + 1
  (x -> (push!(R_name_dict, x => SFPointer(:S, type_index_counter)); type_index_counter+=1)).(not_yet_added_stocks)

  type_index_counter = np(R) + 1
  (x -> (push!(R_name_dict, x => SFPointer(:P, type_index_counter)); type_index_counter+=1)).(not_yet_added_params)
  # @show R_name_dict



  add_stocks!(R, length(not_yet_added_stocks) ; sname = not_yet_added_stocks)
  add_parameters!(R, length(not_yet_added_params) ; pname = not_yet_added_params)



  if (length(R_sum_queue) > 0)
    sums = parse_sum.( R_sum_queue)
    sum_names = map(first, sums)
    sum_lists = map(last, sums)
  else
    sum_names = Vector{Symbol}()
    sum_lists = Vector{Expr}()
  end

  if (length(R_dyvar_queue) > 0)
    dyvars = parse_dyvar.(R_dyvar_queue)
    # filtered_dyvars = collect(filter(x -> !(x.args[1] in vnames(R)), R_dyvar_queue))
    # @show filtered_dyvars
    dyvar_names = map(x -> getindex(x, 1), dyvars)
    dyvar_definitions = map(x -> getindex(x, 2), dyvars)
    dyvar_ops = (x -> x.args[1]).(dyvar_definitions)
  else
    dyvar_names = Vector{Symbol}()
    dyvar_definitions = Vector{Expr}()
    dyvar_ops = Vector{Symbol}()
  end

  # flow_definitions = map(x -> (x.args[2], x.args[3].args[2], x.args[3].args[3]), R_flow_queue)
  # flow
  

  # x = collect(filter(x -> (x.args[3].args[1] in fnames(R)), R_flow_queue ) )
  # @show x

  if (length(R_flow_queue) > 0)
    parsed_flows = parse_flow.(R_flow_queue)
    # @show parsed_flows
    updated_flows_and_dyvars = create_flow_definitions(parsed_flows, vcat(dyvar_names, vnames(I)))
    updated_flows = updated_flows_and_dyvars[1]
    new_dyvars = updated_flows_and_dyvars[2]
    append!(dyvar_names, map(x -> getindex(x, 1), new_dyvars))
    new_dyvar_names = map(x -> getindex(x, 1), new_dyvars)
    new_dyvar_definitions = map(x -> getindex(x, 2), new_dyvars)
    new_dyvar_unwrapped = collect(zip(new_dyvar_names, map(Syntax.get, new_dyvar_definitions)))
    # @show new_dyvar_unwrapped
    # @show new_dyvar_unwrapped[1][1]
    # @show new_dyvar_unwrapped[1][2]
    new_dyvar_definitions = map(kv -> kv[2][1] isa Tuple ? Expr(:(=), kv[1], Expr(:call, kv[2][2], kv[2][1][1], (kv[2][1][2]))) : :Expr(:(=), kv[1], Expr(:call, kv[2][2], (kv[2][1]))),  new_dyvar_unwrapped)

    # new_dyvar_expressions = map(kv -> kv[1] isa Tuple ? Expr(:call, last(kv), kv[1][1], kv[1][2]) : Expr(:call, last(kv), kv[1]), new_dyvar_unwrapped)

    append!(dyvar_definitions, new_dyvar_definitions)
    append!(dyvar_ops, (x -> x.args[2].args[1]).(new_dyvar_definitions))
    append!(R_dyvar_queue, new_dyvar_definitions)

    # @show updated_flows
    flow_stocks_out = map(x -> x[1], parsed_flows)
    flows = map(x -> x[2], parsed_flows)
    flow_stocks_in = map(x -> x[3], parsed_flows)
    flow_names = (x -> (x.args[1])).(flows)
    flow_variables = last.(updated_flows)
    # flow_variables = (x -> (x.args[2])).(flows)
  else
    flow_stocks_out = Vector{Symbol}()
    flow_stocks_in = Vector{Symbol}()
    flows = Vector{Expr}()
    flow_names = Vector{Symbol}()
    flow_variables = Vector{Symbol}()
  end


  # add dyvars and sums to R_name_dict, so we know their index when we need to
  # link them to a dyvar
  type_index_counter = nvb(R) + 1
  (x -> (push!(R_name_dict, x => (SFPointer(:V, type_index_counter))); type_index_counter+=1)).(dyvar_names)

  type_index_counter = nsv(R) + 1
  (x -> (push!(R_name_dict, x => (SFPointer(:SV, type_index_counter))); type_index_counter+=1)).(sum_names)
  # (x -> push!(R_name_dict, x -> SFPointer(:F, nf(R)))).flow_names


  # filtered_dyvars = collect(filter(x -> !(x.args[1] in vnames(R)), R_dyvar_queue))
  filtered_dyvars = collect(filter(x -> !(x.args[1] in vnames(R)), R_dyvar_queue))
  
  filtered_dyvar_names = map(x -> first(x.args), filtered_dyvars)

  # @show filtered_dyvar_names
  filtered_dyvar_ops = map(x -> x.args[2].args[1], filtered_dyvars)

  #TODO filter svariables

  filtered_sum_names = map(x -> x.args[1], collect(filter(x -> !(x.args[1] in svnames(R)), R_sum_queue  )    ))

  add_svariables!(R, length(filtered_sum_names), svname = filtered_sum_names)
  add_variables!(R, length(filtered_dyvar_names), vname = filtered_dyvar_names, vop = filtered_dyvar_ops)


  filtered_flows = filter(x -> !(x.args[1] in fnames(R)), flows)
  # @show updated_flows
  flow_variables_indices = map(x -> findfirst(==(x), vnames(R)), map(x -> x[2], updated_flows))
  filtered_flow_names = map(x -> x.args[1], filtered_flows)
  # filtered_flows = map(x -> x.args[3].args[2].args[1], collect(filter(x -> !(x.args[3].args[2].args[1] in fnames(R)) ,R_flow_queue)))
  # flow_variable_indices = map(x -> findfirst(==(x), vnames(R)), flow_variables)
  # @show flow_variables, filtered_flows
  # @show vnames(R)
  # @show flow_variable_indices
  add_flows!(R, flow_variables_indices, length(filtered_flow_names), fname = filtered_flow_names)

  # @show R

  # @show I_connect_dict[:LPV]

  for sum_index in eachindex(R_sum_queue)
    real_sum_index = findfirst(==(R_sum_queue[sum_index].args[1]), svnames(R))
    for stock in sum_lists[sum_index]
      stock_index = findfirst(==(stock), snames(R))
      link = (stock, R_sum_queue[sum_index].args[1])
      if !(link in I_connect_dict[:LS])
        add_part!(R, :LS ; lss = stock_index, lssv = real_sum_index)
      end
    end
  end

  # @show R_dyvar_queue

  for dyvar_index in eachindex(R_dyvar_queue)
    real_dyvar_index = findfirst(==(R_dyvar_queue[dyvar_index].args[1]), vnames(R))
    dyvar_definition = R_dyvar_queue[dyvar_index].args[2]
    # @show dyvar_definition
    for (operand_index, operand) in enumerate(dyvar_definition.args[2:end])
      link = (operand, R_dyvar_queue[dyvar_index].args[1], operand_index)
      if !(link in I_connect_dict[dyvar_link_name_from_object_type(R_name_dict[operand].type)])
        # @show "triggered!"
        add_operand_link!(R, operand, R_name_dict[operand].type, real_dyvar_index, operand_index)
      end
    end
  end

  for flow_index in eachindex(R_flow_queue)
    stock_in = flow_stocks_in[flow_index]
    if (stock_in != :F_NONE)
      stock_in_index = R_name_dict[stock_in].index
      link = (stock_in, R_flow_queue[flow_index].args[3].args[2].args[1])
      if !(link in I_connect_dict[:I])
        add_part!(R, :I ; is = stock_in_index, ifn = flow_index)
      end
    end

    stock_out = flow_stocks_out[flow_index]
    if (stock_out != :F_NONE)
      stock_out_index = R_name_dict[stock_out].index
      link = (stock_out, R_flow_queue[flow_index].args[3].args[2].args[1])
      if !(link in I_connect_dict[:O])
        add_part!(R, :O ; os = stock_out_index, ofn = flow_index)
      end
    end 
  end

  # R_name_dict = Dict{Symbol, SFPointer}([
  #   [s => SFPointer(:S, i) for (i, s) ∈ enumerate(snames(sf))] ;
  #   [sv => SFPointer(:SV, i) for (i, sv) ∈ enumerate(svnames(sf))] ;
  #   [v => SFPointer(:V, i) for (i, v) ∈ enumerate(vnames(sf))] ;
  #   [f => SFPointer(:F, i) for (i, f) ∈ enumerate(fnames(sf))] ;
  #   [p => SFPointer(:P, i) for (i, p) ∈ enumerate(pnames(sf))]
  # ])
  # @show R_name_dict
  # @show R_link_vector
  for link in R_link_vector
    src = link[1][1]
    tgt = link[1][2]
    src_pointer = R_name_dict[src]
    tgt_pointer = R_name_dict[tgt]

    src_index = src_pointer.index
    tgt_index = tgt_pointer.index

    tgt_type = tgt_pointer.type

    if tgt_type == :SV
      add_part!(R, :LS ; lss = src_index, lssv = tgt_index)
    elseif tgt_type == :F
      position = link[2]
      if position == 1
        add_part!(R, :O ; is = src_index, ifn = tgt_index)
      elseif position == 2
        add_part!(R, :I ; os = src_index, ofn = tgt_index)
      end
    elseif tgt_type == :V
      add_operand_link!(R, src, src_pointer.type, tgt_index, link[2])
    end


  end

  # @show R
  # @show R.subparts.lpvp.m
  # @show R


  # @show L
  # @show I
  # @show R

  hom1 = homomorphism(I,L)
  hom2 = homomorphism(I,R)


  rule = Rule(hom1, hom2)


  sf_rewritten = rewrite(rule, sf)


  if isnothing(sf_rewritten)
    @show sf
    @show homomorphism(L, sf)
    @show hom1
    @show hom2

    @show L
    @show I
    @show R
    # println(L)
    # println(I)
    # println(R)
    return error("Failed to apply rule!")
  end
  return sf_rewritten






end





"""
Given a stockflow and a syntax block, create a new stockflow by using a rewrite
rule to apply modified homomorphisms to it.

Define three new blocks L, I and R, such that sf ⊇ L ⊇ I and R ⊇ I.  Define
homomorphisms I -> L and I -> R.  Apply this change to the original stockflow.

Use :stocks, :flows, :sums, :dynamic_variables and :parameters to add those
particular objects.  Same format as @stock_and_flow.

Use :dyvar_swaps to replace all instances of an object in a dynamic variable
with another object.  Does not delete the original object.

Use :redefs to provide a new definition for an existing sum or dyvar.

Use :removes to delete objects.

```julia
@rewrite aged_sir begin

  :redefs
  v_meanInfectiousContactsPerSv_cINC = cc_C * v_prevalencev_INC_post
  v_meanInfectiousContactsPerSv_cINA = cc_A * v_prevalencev_INA_post


  :parameters
  fcc
  fca
  fac
  faa

  :dynamic_variables

  v_CCContacts = fcc * v_prevalencev_INC
  v_CAContacts = fca * v_prevalencev_INA
  
  v_ACContacts = fac * v_prevalencev_INC
  v_AAContacts = faa * v_prevalencev_INA
  
  v_prevalencev_INC_post = v_CCContacts + v_CAContacts
  v_prevalencev_INA_post = v_ACContacts + v_AAContacts

end
```
"""
macro rewrite(sf, block)
  escaped_block = Expr(:quote, block)
  sf = esc(sf)
  quote
    sfrewrite2($sf, $escaped_block)
  end
end
  


end 
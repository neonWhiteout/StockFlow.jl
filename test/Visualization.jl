@testset "Causal Loop F" begin
  @testset "Graphing Empty Causal Loop F" begin
    empty = CausalLoopPM()
    empty_graph = GraphSF(empty)
    @test isempty(empty_graph.stmts)
  end

  node_info(statement) = (statement.name, statement.attrs[:label])
  edge_info(statement) = (statement.path[1].name, statement.path[2].name, statement.attrs[:label])

  @testset "Graphing standard Causal Loop F" begin
    ABA = CausalLoopPM([:A, :B], [:A => :B, :B => :A], [POL_NEGATIVE, POL_POSITIVE])
    ABA_graph = GraphSF(ABA)
    
    n1 = node_info(ABA_graph.stmts[1])
    n2 = node_info(ABA_graph.stmts[2])

    e1 = edge_info(ABA_graph.stmts[3])
    e2 = edge_info(ABA_graph.stmts[4])

    # Orders + before -
    @test n1 == ("n1", "A") && n2 == ("n2", "B") && e2 == ("n1", "n2", "-") && e1 == ("n2", "n1", "+")


    C = CausalLoopPM([:C], [:C => :C for _ in 1:2], [POL_NEGATIVE, POL_POSITIVE])
    C_graph = GraphSF(C)

    c_n1 = node_info(C_graph.stmts[1])

    c_e1 = edge_info(C_graph.stmts[2])
    c_e2 = edge_info(C_graph.stmts[3])

    
    @test (c_n1 == ("n1", "C") && c_e2 == ("n1", "n1", "-") && c_e1 == ("n1", "n1", "+"))

    DE = CausalLoopPM([:D, :E], [:D => :E], [POL_POSITIVE])
    DE_graph = GraphSF(DE ; schema = "C0")

    de_n1 = node_info(DE_graph.stmts[1])
    de_n2 = node_info(DE_graph.stmts[2])

    de_e1 = DE_graph.stmts[3].path[1].name, DE_graph.stmts[3].path[2].name

    @test (de_n1 == ("n1", "D") && de_n2 == ("n2", "E")
      && ! hasproperty(DE_graph.stmts[3].attrs, :label) && de_e1 == ("n1", "n2"))


  end


  @testset "Graph_RB for Causal Loop F" begin
    
    @test isempty(Graph_RB(CausalLoopPM()).stmts)

    C2 = CausalLoopPM([:C2], [:C2 => :C2 for _ in 1:5], [POL_NEGATIVE, POL_POSITIVE, POL_NEGATIVE, POL_POSITIVE, POL_NEGATIVE])
  #   C2_RB = Graph_RB(C2)


  #   # number nodes + 2 * number edges
  #   C_nodes = [node_info(C2_RB.stmts[i]) for i in 1:11]
  #   # number edges * 3 
  #   C_edges = [(C2_RB.stmts[i].path[1].name, C2_RB.stmts[i].path[2].name) for i in 12:26]

  #   @show C_nodes, C_edges
  #   @test (C_nodes == [("n1", "C2"), 
  #   ("n2", "+"), ("n3", "+"), ("n4", "-"), ("n5", "+"), ("n6", "-"),
  #   ("n7", "B"), ("n8", "R"), ("n9", "±"), ("n10", "?"), ("n11", "0")]

  #   && C_edges == [

  #   ("n1", "n2"), ("n2", "n1"),
  #   ("n1", "n3"), ("n3", "n1"),
  #   ("n1", "n4"), ("n4", "n1"),
  #   ("n1", "n5"), ("n5", "n1"), 
  #   ("n1", "n6"), ("n6", "n1"),

  #   ("n2", "n7"), ("n3", "n8"), ("n4", "n9"), ("n5", "n10"), ("n6", "n11")
  #   ])


  end


 
end
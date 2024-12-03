using NLopt, PhyloNetworks;

function edge2taxonAgeGrad(net::PhyloNetworks.HybridNetwork; checkPreorder::Bool=true):: Matrix{Float64}
  checkPreorder && preorder!(net)
  node_pos = Dict((n.number => i for (i, n) in enumerate(net.nodes_changed)))
  out = zeros(Float64, (length(net.nodes_changed), length(net.edge)))
  for (j, e) in enumerate(net.edge)
    weight = (e.hybrid ? 2.0 * e.gamma : 2.0)
    out[node_pos[e.node[2 - Int(e.isChild1)].number], j] -= weight
    out[node_pos[e.node[Int(e.isChild1) + 1].number], j] += weight
  end
  out
end
"""
    calibrateFromPairwiseDistances!(net, distances::Matrix{Float64},
        taxon_names::Vector{<:AbstractString})

function from PhyloNetworks package: https://github.com/JuliaPhylo/PhyloNetworks.jl
Originally published: https://doi.org/10.1093/molbev/msx235
Customized/extended to add also cost for changing the internal edge length from SNaQ, 
not only minimizing difference between patristic distances and distance matrix.

Calibrate the network to match (as best as possible) input
pairwise distances between taxa, such as observed from sequence data.
`taxon_names` should provide the list of taxa, in the same order
in which they they are considered in the `distances` matrix.
The optimization criterion is the sum of squares between the
observed distances, and the distances from the network
(weighted average of tree distances, weighted by γ's).
The network's edge lengths are modified.

Warning: for many networks, mutiple calibrations can fit the pairwise
distance data equally well (lack of identifiability).
This function will output *one* of these equally good calibrations.

optional arguments (default):
- checkPreorder (true)
- forceMinorLength0 (false) to force minor hybrid edges to have a length of 0
- NLoptMethod (`:LD_MMA`) for the optimization algorithm.
  Other options include `:LN_COBYLA` (derivative-free); see NLopt package.
- tolerance values to control when the optimization is stopped:
  ftolRel (1e-12), ftolAbs (1e-10) on the criterion, and
  xtolRel (1e-10), xtolAbs (1e-10) on branch lengths / divergence times.
- verbose (false)
- max_iter (1000) maximum number of interations NLopt  
- lower_bound (0.0) - lower bound for edge_length
- upper_bound (0.0) - upper bound for edge_length/node_age 
                      if upper_bound <= lower_boud no upper bound is set
- node_weight (0.0) - multiplier for total cost of internal edge length differences: 
                      0.0 = automatically estimated according to the node_leaf_ratio_type
- edge_weight (1.0) - multiplier for each internal edge length difference: 
                      1.0 = fixed number for each edge;
                      []=automatically estimated from how many times is the edge used for connecting pairs of leafs;
                      otherwise vector with weights for each edge in net.edge order
- node_leaf_ratio_type ("mean") - node_leaf_ratio to be used if node_weight=0.0:
                                  "mean" - mean distance between leafs / mean length of internal edges; 
                                  "sum" - sum of distances between leafs / sum of internal edge length)
"""
function calibrateFromPairwiseDistances_ext!(net::PhyloNetworks.HybridNetwork,
      D::Array{Float64,2}, taxNames::Vector{<:AbstractString};
      checkPreorder::Bool=true, forceMinorLength0::Bool=false, verbose::Bool=false,
      ultrametric::Bool=true, NLoptMethod::Symbol=:LD_MMA,
      ftolRel::Float64=PhyloNetworks.fRelBL, ftolAbs::Float64=PhyloNetworks.fAbsBL,
      xtolRel::Float64=PhyloNetworks.xRelBL, xtolAbs::Float64=PhyloNetworks.xAbsBL, 
      node_weight::Float64=0.0, lower_bound::Float64=0.0, upper_bound::Float64=0.0, max_iter::Int=1000, 
      edge_weight::Union{Float64, AbstractVector{<:Float64}} = 1.0, node_leaf_ratio_type::AbstractString="mean")

    checkPreorder && preorder!(net)
    # fixit: remove root node if of degree 2, and if ultrametric=false
    defaultedgelength = median(D)/(length(net.edge)/2)
    for e in net.edge
        if e.length == -1.0 e.length=defaultedgelength; end
        # get smarter starting values: NJ? fast dating?
    end
    if ultrametric # get all node ages in pre-order
        na = getNodeAges(net)
    else 
      na = Float64[]; 
    end
    internal_edge_filter = BitVector(map(e->(e.isMajor || !forceMinorLength0) && all(map(n->!n.leaf, e.node)), net.edge))
    edge_gamma = [(e.hybrid ? e.gamma : 1.0) for e in net.edge[internal_edge_filter]]
    ori_edge_length = [e.length for e in net.edge][internal_edge_filter]
    # get number and indices of edge/nodes to be optimized
    if forceMinorLength0 && !ultrametric
        parind = filter(i -> net.edge[i].isMajor, 1:net.numEdges)
        nparams = length(parind) # edges to be optimized: indices
        par = [e.length for e in net.edge][parind] # and lengths
        for i in 1:net.numEdges
            if !net.edge[i].isMajor net.edge[i].length=0.0; end
        end
    elseif !forceMinorLength0 && !ultrametric
        nparams = length(net.edge) # number of parameters to optimize
        par = [e.length for e in net.edge]
        parind = 1:nparams
    elseif !forceMinorLength0 && ultrametric
        nparams = net.numNodes - net.numTaxa # internal nodes to be optimized
        parind = filter(i -> !net.nodes_changed[i].leaf, 1:net.numNodes)
        par = na[parind]
    else # forceMinorLength0 && ultrametric
        nparams = net.numNodes - net.numTaxa - net.numHybrids
        parind = filter(i -> !(net.nodes_changed[i].leaf || net.nodes_changed[i].hybrid),
                        1:net.numNodes)
        par = na[parind]
        hybInd = filter(i -> net.nodes_changed[i].hybrid, 1:net.numNodes)
        hybParentInd = Int[] # index in nodes_changed of minor parent
        hybGParentI = Int[] # index in 1:nparams of minor (grand-)parent in param list
        for i in hybInd
            n = net.nodes_changed[i]
            p = getparentminor(n)
            pi = findfirst(n -> n===p, net.nodes_changed)
            push!(hybParentInd, pi)
            pii = findfirst(isequal(pi), parind)
            while pii===nothing # in case minor parent of n is also hybrid node
                p = getparentminor(p)
                pi = findfirst(n -> n===p, net.nodes_changed)
                pii = findfirst(isequal(pi), parind)
            end
            push!(hybGParentI, pii)
        end
    end
    # match order of leaves in input matrix, versus pre-order
    nodenames = [n.name for n in net.nodes_changed] # pre-ordered
    ntax = length(taxNames)
    tipind = Int[] # pre-order index for leaf #i in dna distances
    for l in taxNames
        i = findfirst(isequal(l), nodenames)
        i !== nothing || error("taxon $l not found in network")
        push!(tipind, i)
    end
    # initialize M=dist b/w all nodes, G=gradient (constant)
    M = PhyloNetworks.pairwiseTaxonDistanceMatrix(net, keepInternal=true,
            checkPreorder=false, nodeAges=na)
    # Mori = deepcopy(M)
    if !ultrametric && sort([e.number for e in net.edge]) != collect(1:net.numEdges)
        for i in 1:net.numEdges # renumber edges, needed for G
            net.edge[i].number = i
        end
    end # G assumes edges numbered 1:#edges, if optim edge lengths
    # edges will be weighted by how many leaf-taxon pairs they connect
    if edge_weight isa AbstractVector
      if isempty(edge_weight)
        edge_weight = sum(PhyloNetworks.pairwiseTaxonDistanceGrad(net, checkEdgeNumber=false)[tipind, tipind, internal_edge_filter], dims=(1,2))[1, 1, :] ./ 2.0
      else
        edge_weight = edge_weight[internal_edge_filter]
      end
    end
    G = PhyloNetworks.pairwiseTaxonDistanceGrad(net, checkEdgeNumber=false, nodeAges=na) .* 2
    if ultrametric
      node_grad = edge2taxonAgeGrad(net; checkPreorder=false)[parind, internal_edge_filter]
    else
      parind_internal = par_ind[internal_edge_filter[par_ind]]
    end
    node_leaf_ratio = (node_leaf_ratio_type == "sum" ? (sum(D) / 2.0 / sum(edge_gamma .* edge_weight .* ori_edge_length)) : mean(D) * sum(edge_gamma) / sum(edge_gamma .* edge_weight .* ori_edge_length))
    if node_weight == 0.0
      node_weight = node_leaf_ratio
      # println("node_weight=", node_weight)
    end
    # contraints: to force a parent to be older than its child
    if ultrametric
      numConstraints = length(parind) -1 + net.numHybrids
      # calculate indices in param list of child & (grand-)parent once
      chii = Int[] # non-root internal node, can be repeated: once per constraint
      anii = Int[] # closest ancestor in param list
      ci = 1 # index of constraint
      for i in 2:length(net.nodes_changed) # 1=root, can skip
        n = net.nodes_changed[i]
        if n.leaf continue; end # node ages already bounded by 0
        if n.hybrid && forceMinorLength0          # get index in param list of
          nii = hybGParentI[findfirst(isequal(i), hybInd)] # minor grand-parent (same age)
        else
          nii = findfirst(isequal(i), parind)
        end
        for e in n.edge
          if getchild(e) == n # n child of e
            p = getparent(e)  # parent of n
            if forceMinorLength0 && n.hybrid && !e.isMajor
                continue; end # p and n at same age already
            pi = findfirst(isequal(p.number), [no.number for no in net.nodes_changed])
            if forceMinorLength0 && p.hybrid
              pii = hybGParentI[findfirst(isequal(pi), hybInd)]
            else
              pii = findfirst(isequal(pi), parind)
            end
            push!(chii, nii)
            push!(anii, pii)
            verbose && println("node $(net.nodes_changed[parind[nii]].number) constrained by age of parent $(net.nodes_changed[parind[pii]].number)")
            ci += 1
          end
        end
      end
      length(chii) == numConstraints ||
        error("incorrect number of node age constraints: $numConstraints")
      function ageConstraints(result, nodeage, grad)
        if length(grad) > 0 # grad: nparams x nConstraints: ∂cj/∂xi = grad[i,j]
            fill!(grad, 0.0)
        end
        for j in 1:numConstraints
          nii = chii[j]; pii = anii[j]
          result[j] = nodeage[nii] - nodeage[pii] - lower_bound # jth constraint: cj ≤ 0
          if length(grad) > 0 # nparams x nConstraints
            grad[nii, j] =  1.
            grad[pii, j] = -1.
          end
        end
      end
    end
    opt = NLopt.Opt(NLoptMethod,nparams) # :LD_MMA to use gradient
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    NLopt.maxeval!(opt, max_iter) # max iterations
    NLopt.ftol_rel!(opt,ftolRel)
    NLopt.ftol_abs!(opt,ftolAbs)
    NLopt.xtol_rel!(opt,xtolRel)
    NLopt.xtol_abs!(opt,xtolAbs)
    # NLopt.maxtime!(opt, t::Real)
    NLopt.lower_bounds!(opt, fill(lower_bound, nparams))
    if upper_bound > lower_bound
      NLopt.upper_bounds!(opt, fill(upper_bound, nparams))
    end
    if ultrametric
      NLopt.inequality_constraint!(opt,ageConstraints,fill(0.0,numConstraints))
    end
    counter = [0]
    function obj(x::Vector{Float64}, grad::Vector{Float64})
        verbose && println("mismatch objective, BL = $(x)")
        counter[1] += 1
        # update edge lengths or node ages
        if ultrametric # update na using x, in place
            for i in 1:nparams # na=0 at leaves already
                na[parind[i]] = x[i]
            end
            if forceMinorLength0 # set hybrid age to minor parent age
                for i in 1:net.numHybrids
                    na[hybInd[i]] = na[hybParentInd[i]] # pre-order important
                end
            end
        else # not ultrametric: optimize branch lengths
            for i in 1:nparams # update network
                net.edge[parind[i]].length = x[i]
            end
        end
        # update distances in M, in place
        PhyloNetworks.pairwiseTaxonDistanceMatrix!(M,net,na)
        # ss = 0.0 # sum of squares between M and observed distances
        # ss = sum((M[node_filter, node_filter] .- Mori[node_filter, node_filter]) .^ 2) ./ 2.0 # adding distance difference of all internal nodes
        edge_diff = [e.length for e in net.edge[internal_edge_filter]] .- ori_edge_length
        ss = node_weight * sum(edge_gamma .* edge_weight .* edge_diff .* edge_diff)
        for i in 2:ntax; for j in 1:(i-1)
            ss += (M[tipind[i],tipind[j]]-D[i,j])^2 
        end; end
        if length(grad) > 0 # sum_ij 2 * dM_ij/dx_t * (M_ij-D_ij)
          if ultrametric
            grad .= (2.0 * node_weight) .* (node_grad * (edge_weight .* edge_diff))
          else
            grad .= 0.0 # for t in 1:nparams grad[t] = 0.0; end;
            grad[parind_internal] .= (2.0 * node_weight) .* (edge_gamma .* edge_weight .* edge_diff)
          end
          for i in 2:ntax; for j in 1:(i-1);
            for t in 1:nparams
              G[tipind[i],tipind[j],parind[t]] != 0.0 || continue
              grad[t] += G[tipind[i],tipind[j],parind[t]] *
                        (M[tipind[i],tipind[j]]-D[i,j])
            end
            if ultrametric && forceMinorLength0
              for i in 1:net.numHybrids # na[hybrid] set to na[minor parent]
                grad[hybGParentI[i]] += G[tipind[i],tipind[j],hybInd[i]] *
                        (M[tipind[i],tipind[j]]-D[i,j])
              end
            end
          end; end
        end
        return ss
    end
    NLopt.min_objective!(opt,obj)
    fmin, xmin, ret = NLopt.optimize(opt,par) # optimization here!
    verbose && println("got $(round(fmin, digits=5)) at $(round.(xmin, digits=5)) after $(counter[1]) iterations (return code $(ret))")
    return fmin, xmin, ret, counter[1], node_leaf_ratio
end

using NLopt, PhyloNetworks;

const fAbsBL = 1e-15
const fRelBL = 1e-12
const xAbsBL = 1e-15
const xRelBL = 1e-8

function edge2taxonAgeGrad(net::PhyloNetworks.HybridNetwork; checkPreorder::Bool=true, lower_bound::Float64=0.0):: Matrix{Float64}
  checkPreorder && preorder!(net)
  node_pos = Dict((n.number => i for (i, n) in enumerate(net.nodes_changed)))
  out = zeros(Float64, (length(net.nodes_changed), length(net.edge)))
  for (j, e) in enumerate(net.edge)
    # weight = (e.hybrid ? e.gamma : 1.0)
    out[node_pos[e.node[2 - Int(e.isChild1)].number], j] = -1.0  # -= weight  child
    out[node_pos[e.node[Int(e.isChild1) + 1].number], j] = 1.0  # += weight  parent
  end
  out
end

#####
# returns vector of node ages in the order of net.nodes_changed
# metric      - :major (only age from major tree is retrieved; minor edge is only considered if it is the only child edge i.e. in case of chained hybrid nodes)
#               :majormean node ages is calculated as average of its children with major edges (unless there is only minor edge e.g. in case of chained hybrid nodes), but all children ages (including from minor edges) are considered as minimal age for parent to prevent negative edges
#               :max (maximum node age is returned)
# lower_bound - minimal edge length for non-hybrid edges
# lower_bound_hyb - minimal edge length for hybrid edges
function getNodeAges(net, metric::Symbol; checkPreorder::Bool=true, lower_bound::Float64=0.0, lower_bound_hyb::Float64=0.0)
  checkPreorder && preorder!(net)
  # println("getNodeAges edges=$(map(x->x.length, net.edge))")
  out = Vector{Float64}(undef, size(net.nodes_changed, 1))
  node2i = Dict{PhyloNetworks.Node, Int}((Pair(n, i) for (i, n) in enumerate(net.nodes_changed)))
  for i in reverse(axes(net.nodes_changed, 1))
    node = net.nodes_changed[i]
    if node.leaf
      out[i] = 0.0
      continue
    end
    age = 0.0
    age_min = 0.0
    age_cnt = 0
    major = false
    for e in node.edge
      (e.node[1 + Int(e.isChild1)] === node) || continue  # only edge where the node is parent
      child_i = node2i[e.node[2 - Int(e.isChild1)]]
      e_nolen = (e.length < 0.0) && !isapprox(e.length, 0.0)
      e_age = out[child_i] + max(e.length, (e.hybrid ? lower_bound_hyb : lower_bound))
      e_age_min = out[child_i] + (e.hybrid ? lower_bound_hyb : lower_bound)
      (age_min < e_age_min) && (age_min = e_age_min)
      e_major = e.isMajor
      # println("#$(i) ($(node.number)) -> $(child_i) ($(e.node[2 - Int(e.isChild1)].number)) age=$(out[child_i]) + $(e.length) gamma=$(e.gamma)")
      if metric === :major
        if major < e_major
          major = e_major
          age = e_age
        elseif (major == e_major) && (age < e_age)
          age = e_age
        end
      elseif metric === :max
        (age < e_age) && (age = e_age)
      elseif !e_nolen  # :majormean
        if major < e_major
          major = e_major
          age = e_age
          age_cnt = 1
        elseif major == e_major
          age += e_age
          age_cnt += 1
        end
      end
    end
    out[i] = (metric === :majormean ? (age_cnt == 0 ? age_min : max(age / age_cnt, age_min)) : age)
  end
  # println("getNodeAges ages=$(out)")
  out
end

function getParentChildIndex(net; checkPreorder::Bool=true)::Tuple{Vector{Int}, Vector{Int}}
  checkPreorder && preorder!(net)
  node2i = Dict((node=>i for (i, node) in enumerate(net.nodes_changed)))
  [node2i[e.node[1 + Int(e.isChild1)]] for e in net.edge], [node2i[e.node[2 - Int(e.isChild1)]] for e in net.edge]
end

# returns vector of node heigh (distance from the root) in the order of net.nodes_changed
# metric      - :major (only height from major tree is retrieved; minor edge is only considered if it is the only child edge i.e. in case of chained hybrid nodes)
#               :weighted (weights the major and minor edge with gamma, edges with no length are ignored), moreover the node height must be greater or equal to heights of all its parent nodes
#               :max (maximum node height is returned)
# lower_bound - minimal edge length for non-hybrid edges
# lower_bound_hyb - minimal edge length for hybrid edges
function getNodeHeights(net, metric::Symbol; checkPreorder::Bool=true, lower_bound::Float64=0.0, lower_bound_hyb::Float64=0.0)
  checkPreorder && preorder!(net)
  out = Vector{Float64}(undef, size(net.nodes_changed, 1))
  node2i = Dict{PhyloNetworks.Node, Int}((Pair(n, i) for (i, n) in enumerate(net.nodes_changed)))
  for i in axes(net.nodes_changed, 1)
    node = net.nodes_changed[i]
    height = 0.0
    height_min = 0.0
    major = false
    gamma = 0.0
    for e in node.edge
      (e.node[2 - Int(e.isChild1)] === node) || continue  # only edge where the node is child
      parent_i = node2i[e.node[1 + Int(e.isChild1)]]
      e_nolen = (e.length < 0.0) && !isapprox(e.length, 0.0)
      e_height = out[parent_i] + max(e.length, (e.hybrid ? lower_bound_hyb : lower_bound))
      e_height_min = out[parent_i] + (e.hybrid ? lower_bound_hyb : lower_bound)
      (height_min < e_height_min) && (height_min = e_height_min)
      e_major = e.isMajor
      e_gamma = (e.hybrid ? e.gamma : 1.0)
      # println("#$(i) ($(node.number)) -> $(child_i) ($(e.node[2 - Int(e.isChild1)].number)) age=$(out[child_i]) + $(e.length) gamma=$(e.gamma)")
      if (metric === :major)  # there can be only one major parent
        if major < e_major
          major = e_major
          height = e_height
        elseif (major == e_major)  # there can be at most only one minor parent
          height = e_height
        end
      elseif metric === :max
        (height < e_height) && (height = e_height)
      elseif !e_nolen  
        if metric === :majormean
          if major < e_major  # there can be only one major parent
            major = e_major
            height = e_height
          elseif (major == e_major)  # there can be at most only one minor parent
            height = e_height
          end
        else  # :weighted
          height += e_height * e_gamma
          gamma += e_gamma
        end
      end
    end
    out[i] = (metric === :weighted ? (gamma == 0.0 ? height_min : max(height / gamma, height_min)) : 
              (metric === :majormean ? max(height, e_height_min) : height))
  end
  out
end

function getNetWithCorrectedEdgeLength(net, fun_length::Function=identity; ignore_nolength::Bool=true)
  out = deepcopy(net)
  for e in out.edge
    ignore_nolength && (e.length < 0.0) && !isapprox(e.length, 0.0) && continue
    e.length = fun_length(e.length)
  end
  out
end

"""
    calibrateFromPairwiseDistances!(net, distances::Matrix{Float64},
        taxon_names::Vector{<:AbstractString})

function from PhyloNetworks package: https://github.com/JuliaPhylo/PhyloNetworks.jl
Originally published: https://doi.org/10.1093/molbev/msx235
Copyright (c) 2014-2018: Claudia Solis-Lemus and Cecile Ane.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

############################
Customized/extended to add cost for changing the internal branch length estimated by SNaQ after applying power function to it (as SNaQ branches does not correlate well with calendar time), 
which was enabler to search for global instead of local optima.
Square difference between patristic distances and distance matrix are summed with weigthed square differences of internal branches after applying the power function
The calibrateFromPairwiseDistances-ext.jl script is licensed under the MIT "Expat" License
Copyright (c) 2024: Ivan Rurik
############################

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
- forceMinorLength0 (false) - to force minor hybrid edges to have a length of 0
- NLoptMethod ([:GN_ISRES, :LD_MMA]) - either just a local optimizer or a vector of optimizers that should be run sequentially 
  optimizers must support inequality constraint if ultrametric == true
  for global optimizers, upper bound must be also specified
  local optimizers: see NLopt package
    :LD_MMA - gradient based local optimizer
    :LN_COBYLA - derivative-free local optimizer; 
  global optimizers: see NLopt package
    :GN_ISRES
    :GN_ORIG_DIRECT
    :GN_ORIG_DIRECT_L
    :GN_AGS
- tolerance values to control when the optimization is stopped: either a value or vector of values (different value for each optimization in NLoptMethod) 
  ftolRel ([1e-11, 1e-15]), ftolAbs ([1e-12, 1e-12]) on the criterion, and
  xtolRel ([1e-13, 1e-15]), xtolAbs ([1e-8, 1e-8]) on branch lengths / divergence times.
- verbose (false)
- max_iter ([5000000, 10000]) maximum number of interations for each NLoptMethod optimized: either a value or vector of values (different value for each optimization in NLoptMethod) 
- lower_bound (0.0) - lower bound for edge_length
- upper_bound (0.0) - upper bound for edge_length/node_age (needed for global optimizer)
                      if upper_bound <= lower_bound, it will be automatically set based on D matrix and net edge lengths
- node_power (1.0) - initial exponent for the power function applied to internal branch length
- node_power_calib (1.0) - determines the interval for calibrating the initial exponent [node_power / node_power_calib, node_power * node_power_calib] (e.g. 1.0= means initial exponent is fixed, 10= means the exponent will be fitted within range [node_power / 10.0, node_power * 10.0])
- node_multi (1.0) - the pairwise distances will be multiplied by the reciprocal value of this coeficient (1/node_multi)
                     (0.0= means it will be fitted automatically; allowed when node_power_calib > 1.0 and node_weight == 0.0)
- node_weight (0.0) - multiplier for total cost of internal edge length differences: 
                      0.0 = automatically estimated according to the node_leaf_ratio_type
- edge_weight (1.0) - multiplier for each internal edge length difference: 
                      1.0 = fixed number for each edge;
                      []=automatically estimated from how many times is the edge used for connecting pairs of leafs;
                      otherwise vector with weights for each edge in net.edge order
- node_leaf_ratio_type (:mean) - node_leaf_ratio to be used if node_weight=0.0:
                                  :mean - mean distance between leafs / mean length of internal edges; 
                                  :sum - sum of distances between leafs / sum of internal edge length
"""
function calibrateFromPairwiseDistances_ext!(net::PhyloNetworks.HybridNetwork,
      D::Array{Float64,2}, taxNames::Vector{<:AbstractString};
      checkPreorder::Bool=true, forceMinorLength0::Bool=false, verbose::Bool=false,
      ultrametric::Bool=true, lower_bound::Float64=0.0, upper_bound::Float64=0.0, node_power::Float64=1.0, node_power_calib::Float64=1.0, node_multi::Float64=1.0, node_age_metric::Symbol=:weighted,
      node_weight::Float64=0.0, edge_weight::Union{Float64, AbstractVector{<:Float64}} = 1.0, node_leaf_ratio_type::Symbol=:mean,
      NLoptMethod::Union{Symbol, AbstractVector{Symbol}}=[:GN_ISRES, :LD_MMA], max_iter::Union{Integer, AbstractVector{<:Integer}}=[5000000, 10000], neg_edge_tol::Float64=0.0,
      ftolRel::Union{Real, AbstractVector{<:Real}}=[fRelBL * 10000.0, fRelBL], ftolAbs::Union{Real, AbstractVector{<:Real}}=[fAbsBL, fAbsBL],
      xtolRel::Union{Real, AbstractVector{<:Real}}=[xRelBL * 100.0, xRelBL], xtolAbs::Union{Real, AbstractVector{<:Real}}=[xAbsBL, xAbsBL]
)::Tuple{Symbol, Float64, Float64, Float64, Vector{Tuple{Symbol, Symbol, Float64, Vector{Float64}, Int, Vector{Float64}}}}
  (NLoptMethod isa Symbol) && (NLoptMethod = [NLoptMethod])
  (max_iter isa Integer) && (max_iter = [max_iter])
  (ftolRel isa Real) && (ftolRel = [ftolRel])
  (ftolAbs isa Real) && (ftolAbs = [ftolAbs])
  (xtolRel isa Real) && (xtolRel = [xtolRel])
  (xtolAbs isa Real) && (xtolAbs = [xtolAbs])
  power_range = (node_power / node_power_calib, node_power * node_power_calib)
  checkPreorder && preorder!(net)
  # edge_parent_ind, edge_child_ind = getParentChildIndex(net, checkPreorder=false)
  # fixit: remove root node if of degree 2, and if ultrametric=false
  # match order of leaves in input matrix, versus pre-order
  Dfltr = tril(trues(size(D)), -1)
  nodenames = [n.name for n in net.nodes_changed] # pre-ordered
  ntax = length(taxNames)
  tipind = Int[] # pre-order index for leaf #i in dna distances
  for l in taxNames
      i = findfirst(isequal(l), nodenames)
      i !== nothing || error("taxon $l not found in network")
      push!(tipind, i)
  end
  # calc node_leaf_ratio
  internal_edge_filter = BitVector(map(e->(e.isMajor || !forceMinorLength0) && all(map(n->!n.leaf, e.node)) && (e.length != -1.0), net.edge))  # edges without edge.length are ignored
  edge_gamma = [(e.hybrid ? e.gamma : 1.0) for e in net.edge[internal_edge_filter]]
  # edges will be weighted by how many leaf-taxon pairs they connect
  if edge_weight isa AbstractVector
    if isempty(edge_weight)
      edge_weight = sum(PhyloNetworks.pairwiseTaxonDistanceGrad(net, checkEdgeNumber=false)[tipind, tipind, internal_edge_filter], dims=(1,2))[1, 1, :] ./ 2.0
    else
      edge_weight = edge_weight[internal_edge_filter]
    end
  end
  # get number and indices of edge/nodes to be optimized
  if forceMinorLength0 && !ultrametric
      parind = filter(i -> net.edge[i].isMajor, 1:net.numEdges)
      npar = length(parind) # edges to be optimized: indices
  elseif !forceMinorLength0 && !ultrametric
      npar = length(net.edge) # number of parameters to optimize
      parind = 1:npar
  elseif !forceMinorLength0 && ultrametric
      npar = net.numNodes - net.numTaxa # internal nodes to be optimized
      parind = filter(i -> !net.nodes_changed[i].leaf, 1:net.numNodes)
  else # forceMinorLength0 && ultrametric
      # hybrid nodes are the nodes where the hybridization "flows in" and have 2 edges to parents (major and minor); if minor edge length==0, age of hybridization node is already determined by the are of the minor parent node (from which the hybridization "flows out")
      npar = net.numNodes - net.numTaxa - net.numHybrids
      parind = filter(i -> !(net.nodes_changed[i].leaf || net.nodes_changed[i].hybrid),  
                      1:net.numNodes)
      hybInd = filter(i -> net.nodes_changed[i].hybrid, 1:net.numNodes)
      hybParentInd = Int[] # index in nodes_changed of minor parent
      hybGParentI = Int[] # index in 1:npar of minor (grand-)parent in param list
      for i in hybInd
          n = net.nodes_changed[i]
          p = getparentminor(n)
          p_i = findfirst(n -> n===p, net.nodes_changed)
          push!(hybParentInd, p_i)
          pii = findfirst(isequal(p_i), parind)
          while pii===nothing # in case minor parent of n is also hybrid node
              p = getparentminor(p)
              p_i = findfirst(n -> n===p, net.nodes_changed)
              pii = findfirst(isequal(p_i), parind)
          end
          push!(hybGParentI, pii)
      end
  end
  nparall = npar
  lower_bound_ori = 0.0
  # additional params
  if node_power_calib > 1.0
    nparall += 1
    power_ind = nparall
    lower_bound_ori = 1e-15
    # node_power = sqrt(power_range[1] * power_range[2])
  else
    power_ind = 0
    # node_power = power_range[1]
  end
  (node_multi > 0.0) || ((power_ind != 0) && (node_weight == 0.0)) || error("node_multi==0 is allowed only with node_weight==0.0 and node_power_calib > 1.0")
  total_edge_weight = edge_gamma .* edge_weight
  ori_edge_length = [max(e.length, lower_bound_ori) for e in net.edge][internal_edge_filter]
  log_ori_edge_length = log.(ori_edge_length)
  # contraints: to force a parent to be older than its child
  # number of constraints == number of edges between internal nodes
  if ultrametric
    numConstraints = length(parind) -1 + net.numHybrids
    # calculate indices in param list of child & (grand-)parent once
    chii = Int[] # non-root internal node, can be repeated: once per constraint
    anii = Int[] # closest ancestor in param list
    hybridii = Bool[]
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
          forceMinorLength0 && n.hybrid && !e.isMajor && continue  # p and n at same age already
          push!(hybridii, n.hybrid)  # hybrid node i.e. node that has minor and major hybrid parent edges
          p_i = findfirst(isequal(p.number), [no.number for no in net.nodes_changed])
          if forceMinorLength0 && p.hybrid
            pii = hybGParentI[findfirst(isequal(p_i), hybInd)]
          else
            pii = findfirst(isequal(p_i), parind)
          end
          push!(chii, nii)
          push!(anii, pii)
          verbose && println("node $(net.nodes_changed[parind[nii]].number) constrained by age of parent $(net.nodes_changed[parind[pii]].number)")
          ci += 1
        end
      end
    end
    length(chii) == numConstraints || error("incorrect number of node age constraints: $numConstraints")
    function ageConstraints(result, nodeage, grad)
      if length(grad) > 0 # grad: npar x nConstraints: ∂cj/∂xi = grad[i,j] (i.e. grad[i,j] = ∂result(j)/∂nodeage(i))
          fill!(grad, 0.0)
      end
      for j in 1:numConstraints
        nii = chii[j]; pii = anii[j]
        result[j] = nodeage[nii] - nodeage[pii] + (hybridii[j] ? 0.0 : lower_bound)  # jth constraint: cj ≤ 0  (i.e. negative edge length: edge_length(-lower_bound)=parent_age-child_age(-lower_bound))
        if length(grad) > 0 # npar x nConstraints
          grad[nii, j] =  1.
          grad[pii, j] = -1.
        end
      end
      verbose && (counter % counter_interval == 1) && println("- node age constraint: $(count(result .> 0.0)) negative edge lengths")
    end
  end
  # verbose && println(map(x->x.length, net.edge))
  # now empty edge lengths will be prefilled with default value
  defaultedgelength = max(median(D[Dfltr])/(length(net.edge)/2), lower_bound) ^ node_power
  for (i, e) in enumerate(net.edge)
      e.length = ((e.length < 0.0) && !isapprox(e.length, 0.0) ? e.length=defaultedgelength : e.length ^ node_power)  # get smarter starting values: NJ? fast dating?
      e.number = i  # pairwiseTaxonDistanceGrad assumes edges are numbered 1:#edges as it creates a cube; if optim edge lengths (i.e. !ultrametric)
  end
  # verbose && println(map(x->x.length, net.edge[internal_edge_filter]))
  function get_par(prev_par=nothing)
    if !isnothing(prev_par)
      update_par(prev_par)
    end
    if ultrametric # get all node ages in pre-order
      lna = getNodeAges(net, node_age_metric, checkPreorder=false, lower_bound=lower_bound, lower_bound_hyb=0.0)
      lpar = lna[parind]
    else 
      lna = Float64[]; 
      lpar = [(forceMinorLength0 && !e.isMajor ? 0.0 : max(e.length, (e.hybrid ? 0.0 : lower_bound))) for e in net.edge][parind]
    end
    (power_ind == 0) || push!(lpar, (isnothing(prev_par) ? node_power : prev_par[power_ind]))  # node calibration parameter
    return (lna, lpar)
  end
  # initialize M=dist b/w all nodes, G=gradient (constant)
  (na, par) = get_par()
  M = PhyloNetworks.pairwiseTaxonDistanceMatrix(net, keepInternal=true, checkPreorder=false, nodeAges=na)
  indM = view(M, tipind, tipind)
  # Mori = deepcopy(M)
  G = PhyloNetworks.pairwiseTaxonDistanceGrad(net, checkEdgeNumber=false, nodeAges=na) .* 2  # not depending on actual nodeAges (it only checks if is empty or not)
  indG = view(G, tipind, tipind, :)
  if ultrametric
    node_grad = edge2taxonAgeGrad(net; checkPreorder=false)[parind, internal_edge_filter]
  else
    parind_internal = par_ind[internal_edge_filter[par_ind]]
  end
  # initialize all fixed params
  # verbose && println(map(x->x.length, net.edge[internal_edge_filter]))
  edge_ori = (ori_edge_length .^ node_power)
  weighted_edge_ori = total_edge_weight .* (ori_edge_length .^ node_power)
  sum_weighted_edge_ori = sum(weighted_edge_ori)
  sum_D2 = sum(D[Dfltr] .^ 2)
  # node_weight
  if node_weight == 0.0
    edge_diff = [e.length for e in net.edge[internal_edge_filter]] .- edge_ori
    weighted_edge_diff = total_edge_weight .* edge_diff
    node_weight_numerator = (node_leaf_ratio_type === :sum ? sum(D[Dfltr]) : mean(D[Dfltr]) * sum(edge_gamma))
    node_multi_var = (node_multi == 0 ? 2.0 * sum_D2 / (2.0 * sum(indM[Dfltr] .* D[Dfltr]) - node_weight_numerator * sum(weighted_edge_diff .* edge_diff) / sum_weighted_edge_ori) : node_multi)
    node_weight_fix = (power_ind == 0 ? node_weight_numerator / sum_weighted_edge_ori / node_multi_var : 0.0)
  else 
    node_weight_fix = node_weight
    node_multi_var = node_multi
  end
  verbose && println("node_multi=$(round(node_multi_var, digits=5)) ($(round(node_multi, digits=5))) @ node_power=$(round(node_power, digits=5)) node_weight=$(round(node_weight, digits=5))")
  # node_leaf_ratio = (node_leaf_ratio_type === :sum ? (sum(D) / 2.0 / sum(edge_gamma .* edge_weight .* ori_edge_length)) : mean(D) * sum(edge_gamma) / sum(edge_gamma .* edge_weight .* ori_edge_length))
  # if node_weight == 0.0
  #   node_weight = node_leaf_ratio
  #   # println("node_weight=", node_weight)
  # end
  #
  # bounds
  max_par = max(maximum(D[Dfltr]) / node_multi_var, 2.0 * maximum(par))
  if upper_bound < max_par
    upper_bound = max_par
  end
  if node_power != 1.0
    max_par = maximum(getNodeHeights(getNetWithCorrectedEdgeLength(net, x->x^node_power), :max; checkPreorder=false, lower_bound=lower_bound))
    if upper_bound < max_par
      upper_bound = max_par
    end
  end
  if node_power_calib > 1.0
    max_par = maximum((maximum(getNodeHeights(getNetWithCorrectedEdgeLength(net, x->x^(node_power*cal)), :max; checkPreorder=false, lower_bound=lower_bound)) for cal in (1.0/node_power_calib, node_power_calib)))
    if upper_bound < max_par
      upper_bound = max_par
    end
  end
  if ultrametric
    nparx = net.numNodes - net.numTaxa # internal nodes to be optimized
    parindx = filter(i -> !net.nodes_changed[i].leaf, 1:net.numNodes)
    tmp = Vector{Float64}(undef, nparx)
    node2ix = Dict{PhyloNetworks.Node, Int}((Pair(net.nodes_changed[parindx[i]], i) for i in axes(parindx,1)))
    for i in nparx:-1:1
      parix = parindx[i]
      node = net.nodes_changed[parix]
      tmp[i] = maximum(map(x->(x[1] == 0 ? 0.0 : (tmp[x[1]] - (x[2] ? lower_bound : 0.0))), ((get(node2ix, e.node[2 - Int(e.isChild1)], 0), e.hybrid) for e in node.edge if e.node[1 + Int(e.isChild1)] === node)), init=0.0) + lower_bound  # taking maximum of children + lower_bound (all hybrid edges can be 0.0)
    end
    if forceMinorLength0
      lower_bounds = tmp[map(x->!x.hybrid, net.nodes_changed[parindx])]
    else
      lower_bounds = copy(tmp)
    end
    for i in 1:nparx
      parix = parindx[i]
      node = net.nodes_changed[parix]
      tmp[i] = minimum(map(x->(x[1] == 0 ? (upper_bound + lower_bound) : (tmp[x[1]] + (x[2] ? lower_bound : 0.0))), ((get(node2ix, e.node[1 + Int(e.isChild1)], 0), e.hybrid) for e in node.edge if e.node[2 - Int(e.isChild1)] === node)), init=upper_bound + lower_bound) - lower_bound  # taking minimum of parents - lower_bound (all hybrid edges can be 0.0)
    end
    if forceMinorLength0
      upper_bounds = tmp[map(x->!x.hybrid, net.nodes_changed[parindx])]
    else
      upper_bounds = copy(tmp)
    end
  else
    lower_bounds = [(e.hybrid ? 0.0 : lower_bound) for e in net.edge[parind]]
    upper_bounds = fill(upper_bound, npar)
  end
  # power_range parameter
  if power_ind != 0
    push!(lower_bounds, node_power / node_power_calib)  # power_range[1])
    push!(upper_bounds, node_power * node_power_calib)  # power_range[2])
  end
  counter = 0
  node_weight_var = 0.0
  function update_par(x::Vector{Float64})
    # update edge lengths or node ages
    if ultrametric # update na using x, in place
      # for i in 1:npar # na=0 at leaves already
      #     na[parind[i]] = x[i]
      # end
      na[parind] .= x[begin:npar]
      if forceMinorLength0 # set hybrid age to minor parent age
          for i in 1:net.numHybrids
              na[hybInd[i]] = na[hybParentInd[i]] # pre-order important
          end
      end
      # for (i, e) in enumerate(net.edge)
      #   e.length = na[edge_parent_ind[i]] - na[edge_child_ind[i]]
      # end
    else # not ultrametric: optimize branch lengths
        for i in 1:npar # update network
            net.edge[parind[i]].length = x[i]
        end
    end
    # update distances in M, in place
    PhyloNetworks.pairwiseTaxonDistanceMatrix!(M,net,na)
  end
  function obj(x::Vector{Float64}, grad::Vector{Float64})
    # verbose && println("mismatch objective, BL = $(x)")
    counter += 1
    update_par(x)
    # ss = 0.0 # sum of squares between M and observed distances
    # ss = sum((M[node_filter, node_filter] .- Mori[node_filter, node_filter]) .^ 2) ./ 2.0 # adding distance difference of all internal nodes
    if power_ind != 0
      node_power = x[power_ind]
      edge_ori = ori_edge_length .^ node_power
      weighted_edge_ori = total_edge_weight .* edge_ori
      sum_weighted_edge_ori = sum(weighted_edge_ori)
    end
    edge_diff = [e.length for e in net.edge[internal_edge_filter]] .- edge_ori
    weighted_edge_diff = total_edge_weight .* edge_diff
    node_weight_var = (node_weight_fix == 0.0 ? node_weight_numerator / sum_weighted_edge_ori / node_multi_var : node_weight_fix)
    ss_node = node_weight_var * sum(weighted_edge_diff .* edge_diff)
    verbose && (counter % counter_interval == 1) && println("counter=", counter, ", ss_node=", round(ss_node, digits=5), ", node_multi=", round(node_multi_var, digits=5), ", node_power=", round(node_power, digits=5))
    # ss = ss_node = (node_weight / x[end]) * sum(edge_gamma .* edge_weight .* edge_diff .* edge_diff)
    # edge_diff = [e.length for e in net.edge[internal_edge_filter]] .- (ori_edge_length * x[end])
    # ss = ss_node = (node_weight / x[end]) * sum(edge_gamma .* edge_weight .* edge_diff .* edge_diff)
    # koef2 = x[end] * x[end]
    # ss = ss_node = (node_weight / koef2) * sum(edge_gamma .* edge_weight .* edge_diff .* edge_diff)
    # for i in 2:ntax; for j in 1:(i-1)
    #     ss += (M[tipind[i],tipind[j]]-D[i,j]/node_multi)^2 
    # end; end
    tipsD = D[Dfltr] ./ node_multi_var
    tipsDiff = indM[Dfltr] .- tipsD
    # println("size(tipsDiff)=$(size(tipsDiff))")
    # println("size(indG[Dfltr, parind])=$(size(indG[Dfltr, parind]))")
    ss = ss_node + sum(tipsDiff .^ 2)
    if length(grad) > 0 # sum_ij 2 * dM_ij/dx_t * (M_ij-D_ij)
      if power_ind != 0
        grad[power_ind] = -2.0 * node_weight_var * sum(weighted_edge_diff .* edge_ori .* log_ori_edge_length) - ss_node * sum(weighted_edge_ori .* log_ori_edge_length) / sum_weighted_edge_ori
      end
      verbose && (counter % counter_interval == 1) && (power_ind != 0) && println("- grad power_ind=$(round(grad[power_ind], digits=5))")
      if ultrametric
        grad[1:npar] .= (2.0 * node_weight_fix) .* (node_grad * weighted_edge_diff)
        # grad .= [(2.0 * node_weight / x[end]) .* (node_grad * (edge_weight .* edge_diff)); - (2.0 * node_weight / x[end]) * sum(ori_edge_length .* edge_gamma .* edge_weight .* edge_diff) - (ss_node / x[end])] 
        # grad .= [(2.0 * node_weight / koef2) .* (node_grad * (edge_weight .* edge_diff)); - (2.0 * node_weight / koef2) * sum(ori_edge_length .* edge_gamma .* edge_weight .* edge_diff) - 2 * (ss_node / x[end])] 
      else
        grad .= 0.0 # for t in 1:npar grad[t] = 0.0; end;
        grad[parind_internal] .= (2.0 * node_weight_fix) .* weighted_edge_diff
        # grad[parind_internal] .= [(2.0 * node_weight / x[end]) .* (edge_gamma .* edge_weight .* edge_diff); - (2.0 * node_weight / x[end]) * sum(ori_edge_length .* edge_gamma .* edge_weight .* edge_diff) - (ss_node / x[end])]
        # grad[parind_internal] .= [(2.0 * node_weight / koef2) .* (edge_gamma .* edge_weight .* edge_diff); - (2.0 * node_weight / koef2) * sum(ori_edge_length .* edge_gamma .* edge_weight .* edge_diff) - 2 * (ss_node / x[end])] 
      end
      # for i in 2:ntax; for j in 1:(i-1);
      #   for t in 1:npar
      #     G[tipind[i],tipind[j],parind[t]] != 0.0 || continue
      #     grad[t] += G[tipind[i],tipind[j],parind[t]] *
      #               (M[tipind[i],tipind[j]]-D[i,j])
      #   end
      #   if ultrametric && forceMinorLength0
      #     for h in 1:net.numHybrids # na[hybrid] set to na[minor parent]
      #       grad[hybGParentI[h]] += G[tipind[i],tipind[j],hybInd[h]] *
      #               (M[tipind[i],tipind[j]]-D[i,j])
      #     end
      #   end
      # end; end
      grad[1:npar] .+= vec(sum(indG[Dfltr, parind] .* tipsDiff, dims=1))
      ultrametric && forceMinorLength0 && (grad[hybGParentI] .+= vec(sum(indG[Dfltr, hybind] .* tipsDiff, dims=1)))
    end
    return ss
  end
  # optimization
  out = Vector{Tuple{Symbol, Symbol, Float64, Vector{Float64}, Int, Vector{Float64}}}()
  last_rec = nothing
  prev_par = par
  counter_interval = 10000
  for (opt_i, opt_method) in enumerate(NLoptMethod)
    # println("*** Running $(opt_method)...")
    counter = 0
    opt = NLopt.Opt(opt_method, nparall) # :LD_MMA to use gradient
    # :LN_COBYLA for (non)linear constraits, :LN_BOBYQA for bound constraints
    NLopt.maxeval!(opt, max_iter[min(opt_i, end)]) # max iterations
    NLopt.ftol_rel!(opt, ftolRel[min(opt_i, end)])
    NLopt.ftol_abs!(opt, ftolAbs[min(opt_i, end)])
    NLopt.xtol_rel!(opt, xtolRel[min(opt_i, end)])
    NLopt.xtol_abs!(opt, xtolAbs[min(opt_i, end)])
    # NLopt.maxtime!(opt, t::Real)
    if ultrametric
      NLopt.inequality_constraint!(opt, ageConstraints, fill(neg_edge_tol, numConstraints))
    end
    NLopt.lower_bounds!(opt, lower_bounds)
    NLopt.upper_bounds!(opt, upper_bounds)
    NLopt.min_objective!(opt, obj)
    (na, par) = get_par(prev_par)
    fmin, xmin, ret = NLopt.optimize(opt, par) # optimization here!    
    verbose && println("return code=$(ret) fmin=$(round(fmin, digits=5)) xmin=$(round.(xmin, digits=5)) node_multi=$(round(node_multi_var, digits=5)) iter=$(counter)")
    # FORCED_STOP = -5, 
    # ROUNDOFF_LIMITED = -4
    # OUT_OF_MEMORY = -3
    # INVALID_ARGS = -2
    # FAILURE = -1
    # SUCCESS = 1
    # STOPVAL_REACHED = 2
    # FTOL_REACHED = 3
    # XTOL_REACHED = 4
    # MAXEVAL_REACHED = 5
    # MAXTIME_REACHED = 6
    if ret in (:FAILURE, :INVALID_ARGS, :OUT_OF_MEMORY, :ROUNDOFF_LIMITED, :FORCED_STOP)
      error("$ret exception in NLopt.$(opt_method): fmin=$(round(fmin, digits=5)) xmin=$(round.(xmin, digits=5)) node_multi=$(round(node_multi_var, digits=5)) iter=$(counter) xinit=$(round.(prev_par, digits=5))")
    end
    prev_par = xmin
    push!(out, (opt_method, ret, fmin, xmin, counter, prev_par))
    last_rec = (ret, node_weight_var, node_multi_var, node_power, out)
    counter_interval = 100
  end
  return last_rec
end

## -----------------------------------------------------------------------------
## numerical
## -----------------------------------------------------------------------------
function jac_x(d::JuMP.NLPEvaluator, x::Array{Float64,1})
    nconstr = MathProgBase.numconstr(d.m)
    nvar = MathProgBase.numvar(d.m)
    J = spzeros(nconstr, nvar)
    j!(J, x, model=d)
    return J
end
function dFdy_dFdx(d::JuMP.NLPEvaluator, x::Array{Float64,1}, opfdata::OPFData)
    ## evaluate
    J = jac_x(d, x)

    ## indices
    lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
    busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;
    nbus = length(buses); nline = length(lines); ngen = length(generators); nload = length(findall(buses.bustype .== 1))
    Pg_idx_offset = 0; Qg_idx_offset = ngen
    Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
    Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
    # ref_id = opfdata.bus_ref
    # gen_id = filter(x -> x ∉ ref_id, opfdata.generators.bus)
    # load_id = filter(x -> x ∉ Set([gen_id; ref_id]), opfdata.buses.bus_i)  ## purely load buses
    ref_idx = busIdx
    ref_idx = opfdata.bus_ref
    gen_idx = filter(x -> x ∉ ref_idx, opfdata.generators.bus)
    load_idx = filter(x -> x ∉ Set([gen_idx; ref_idx]), opfdata.buses.bus_i)  ## purely load buses
    ref_idx = opfdata.bus_ref

    #       Qg^{G ∪ R}                           Vm^L                       Va^{L ∪ G}
    yidx = [Qg_idx_offset .+ [gen_idx; ref_idx]; Vm_idx_offset .+ load_idx; Va_idx_offset .+ [gen_idx; load_idx]]
    #       Pd^{R ∪ G ∪ L}                    Qd^{R ∪ G ∪ L}
    xidx = [Pd_idx_offset .+ collect(1:nbus); Qd_idx_offset .+ collect(1:nbus)]
    #       P^{G ∪ L}          Q^{R ∪ G ∪ L}
    fidx = [gen_idx; load_idx; nbus .+ [ref_idx; gen_idx; load_idx]]

    ## return
    dFdy = J[fidx, yidx]
    dFdx = J[fidx, xidx]
    return dFdy, dFdx
end


# function get_idx_id(M::OPFModel, opfdata::OPFData)
#   @assert(M.kind == :S)
#   nbus = length(opfdata.buses)
#   ngen = length(opfdata.generators)
#   Pg_idx_offset = 0; Qg_idx_offset = ngen
#   Vm_idx_offset = 2*ngen; Va_idx_offset = 2*ngen+nbus
#   Pd_idx_offset = 2*ngen+2*nbus; Qd_idx_offset = 2*ngen+3*nbus
#   refid = opfdata.bus_ref
#   genid = filter(x -> x ∉ refid, opfdata.generators.bus)
#   loadid = filter(x -> x ∉ Set([genid; refid]), opfdata.buses.bus_i)  ## purely load buses
#
#   BUS_idx_to_id = mapIdxToBusId(opfdata)
#   GEN_idx_to_id = mapGenersToBuses(opfdata)
#
#   ## variable / parameters
#   ## ---------------------------------------------------------------------------
#   full = Dict{Symbol, NamedTuple}()
#   full[:Pg] = get_idx_and_id(M.m, :Pg, GEN_idx_to_id, Pg_idx_offset)
#   full[:Qg] = get_idx_and_id(M.m, :Qg, GEN_idx_to_id, Qg_idx_offset)
#   full[:Vm] = get_idx_and_id(M.m, :Vm, BUS_idx_to_id, Vm_idx_offset)
#   full[:Va] = get_idx_and_id(M.m, :Va, BUS_idx_to_id, Va_idx_offset)
#   full[:Pd] = get_idx_and_id(M.m, :Pd, BUS_idx_to_id, Pd_idx_offset)
#   full[:Qd] = get_idx_and_id(M.m, :Qd, BUS_idx_to_id, Qd_idx_offset)
#
#   part = Dict{Symbol, NamedTuple}()
#   part[:Pg_G] = filter_namedtuple(full[:Pg], genid)
#   part[:Pg_R] = filter_namedtuple(full[:Pg], refid)
#   part[:Pg_L] = (id=:NA, idx=:NA)
#   part[:Qg_G] = filter_namedtuple(full[:Qg], genid)
#   part[:Qg_R] = filter_namedtuple(full[:Qg], refid)
#   part[:Qg_L] = (id=:NA, idx=:NA)
#   part[:Vm_G] = filter_namedtuple(full[:Vm], genid)
#   part[:Vm_R] = filter_namedtuple(full[:Vm], refid)
#   part[:Vm_L] = filter_namedtuple(full[:Vm], loadid)
#   part[:Va_G] = filter_namedtuple(full[:Va], genid)
#   part[:Va_R] = filter_namedtuple(full[:Va], refid)
#   part[:Va_L] = filter_namedtuple(full[:Va], loadid)
#   part[:Pd_G] = filter_namedtuple(full[:Pd], genid)
#   part[:Pd_R] = filter_namedtuple(full[:Pd], refid)
#   part[:Pd_L] = filter_namedtuple(full[:Pd], loadid)
#   part[:Qd_G] = filter_namedtuple(full[:Qd], genid)
#   part[:Qd_R] = filter_namedtuple(full[:Qd], refid)
#   part[:Qd_L] = filter_namedtuple(full[:Qd], loadid)
#
#   ## pf equations
#   ## ---------------------------------------------------------------------------
#   pfe = Dict{Symbol, Array{Int,1}}()
#   pfe[:Pidx] = [linearindex(M.m[:P][i]) for i in eachindex(M.m[:P])]
#   pfe[:Qidx] = [linearindex(M.m[:Q][i]) for i in eachindex(M.m[:Q])]
#   pfe[:Prefidx] = [refid]
#   pfe[:Qrefidx] = [nbus + refid]
#   return full, part, pfe
# end
#
# function get_y(; part::Dict{Symbol, NamedTuple}, rem_ref::Bool=true)
#   out = Dict()
#   #                  3                1                1                1                4                1
#   if rem_ref == false
#     out[:idx] = [part[:Va_G].idx; part[:Va_R].idx; part[:Va_L].idx; part[:Vm_L].idx; part[:Qg_G].idx; part[:Qg_R].idx]
#     out[:id]  = [part[:Va_G].id; part[:Va_R].id; part[:Va_L].id; part[:Vm_L].id; part[:Qg_G].id; part[:Qg_R].id]
#   else
#     out[:idx] = [part[:Va_G].idx; part[:Va_L].idx; part[:Vm_L].idx; part[:Qg_G].idx; part[:Qg_R].idx]
#     out[:id]  = [part[:Va_G].id; part[:Va_L].id; part[:Vm_L].id; part[:Qg_G].id; part[:Qg_R].id]
#   end
#   out[:idx_to_sum] = part[:Qg_G].idx[findall(part[:Qg_G].id .∈ Ref(nonunique(part[:Qg_G].id)))]
#   out[:Jidx_to_sum] = findall(out[:idx] .∈ Ref(out[:idx_to_sum]))
#   return out
# end
#
# function get_F(; pfe::Dict{Symbol, Array{Int,1}}, rem_ref::Bool=true)
#   out = Dict()
#   if rem_ref == false
#     out[:idx] = [pfe[:Pidx]; pfe[:Qidx]]
#   else
#     out[:idx] = [filter(x -> x ∉ pfe[:Prefidx], pfe[:Pidx]); pfe[:Qidx]]
#   end
#   return out
# end



## -----------------------------------------------------------------------------
## analytic
## -----------------------------------------------------------------------------
function get_values(M::OPFModel)
  @assert(M.kind == :S)
  v = Dict{Symbol, Array{Float64,1}}()
  v[:Pg] = getvalue(M.m[:Pg])
  v[:Qg] = getvalue(M.m[:Qg])
  v[:Vm] = getvalue(M.m[:Vm])
  v[:Va] = getvalue(M.m[:Va])
  v[:Pd] = getvalue(M.m[:Pd])
  v[:Qd] = getvalue(M.m[:Qd])
  v[:x] = [v[:Pg]; v[:Qg]; v[:Vm]; v[:Va]; v[:Pd]; v[:Qd]]
  return v
end
# v = get_values(M)
# function get_data(opfdata::OPFData)
#   YffR, YffI, YttR, YttI, YftR, YftI, YtfR, YtfI, YshR, YshI = computeAdmitances(opfdata.lines, opfdata.buses, opfdata.baseMVA)
#
# function Fp(vals::Dict, data::Dict, b::Int64)
#
#   ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
#   + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
#   + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
#   - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - buses[b].Pd ) / baseMVA      # Sbus part
#
#   ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] )
#   + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
#   + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
#   - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - sum(baseMVA*Pd[l] for l in busIdx[b]) ) / baseMVA      # Sbus part

# J = J[F[:idx], y[:idx]]
# if ~isempty(y[:Jidx_to_sum])
#     v = sum(J[:, y[:Jidx_to_sum]], dims=2)
#     J = [J[:, (1):(y[:Jidx_to_sum][1])]  J[:, (y[:Jidx_to_sum][end]+1):(end)]]
# end

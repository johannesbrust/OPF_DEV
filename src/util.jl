function acopf_solve(opfmodel::JuMP.Model, opfdata::OPFData)

  #
  # Initial point - needed especially for pegase cases
  #
  Pg0,Qg0,Vm0,Va0 = acopf_initialPt_IPOPT(opfdata)

  # 04/01/19, J.B., modification for JuMP versions
  # if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
  # 	set_value(getindex(opfmodel, :Pg), Pg0)
  # 	set_value(getindex(opfmodel, :Qg), Qg0)
  # 	set_value(getindex(opfmodel, :Vm), Vm0)
  # 	set_value(getindex(opfmodel, :Va), Va0)
  # else
	setvalue(getindex(opfmodel, :Pg), Pg0)
    setvalue(getindex(opfmodel, :Qg), Qg0)
    setvalue(getindex(opfmodel, :Vm), Vm0)
    setvalue(getindex(opfmodel, :Va), Va0)
  # end
  status = :IpoptInit
  status = solve(opfmodel)

  if status != :Optimal
    println("Could not solve the model to optimality.")
  end
  return opfmodel, status
end
function acopf_solve(M::OPFModel, opfdata::OPFData); return OPFModel(acopf_solve(M.m, opfdata)..., M.kind); end

function acopf_outputAll(opfmodel::JuMP.Model, kind::Symbol, opfdata::MPCCases.OPFData)
  #shortcuts for compactness
  lines = opfdata.lines; buses = opfdata.buses; generators = opfdata.generators; baseMVA = opfdata.baseMVA
  busIdx = opfdata.BusIdx; FromLines = opfdata.FromLines; ToLines = opfdata.ToLines; BusGeners = opfdata.BusGenerators;

  nbus  = length(buses); nline = length(lines); ngen  = length(generators)

  # OUTPUTING
  println("Objective value: ", getobjectivevalue(opfmodel), "USD/hr")
  VM=getvalue(getindex(opfmodel,:Vm)); VA=getvalue(getindex(opfmodel,:Va))
  PG=getvalue(getindex(opfmodel,:Pg)); QG=getvalue(getindex(opfmodel,:Qg))
  if kind == :S
    PD=getvalue(getindex(opfmodel,:Pd)); QD=getvalue(getindex(opfmodel,:Qd))
  end

  println("============================= BUSES ==================================")
  println("  BUS    Vm     Va    |   Pg (MW)    Qg(MVAr) ")   # |    P (MW)     Q (MVAr)")  #|         (load)   ")

  println("                      |     (generation)      ")
  println("----------------------------------------------------------------------")
  for i in 1:nbus
    @printf("%4d | %6.2f  %6.2f | %s  | \n",
	    buses[i].bus_i, VM[i], VA[i]*180/pi,
	    (length(BusGeners[i])==0) ? "   --          --  " : @sprintf("%7.2f     %7.2f", baseMVA*PG[BusGeners[i][1]], baseMVA*QG[BusGeners[i][1]]))
  end
  println("\n")

  within=20 # percentage close to the limits


  nflowlim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nflowlim += 1
    end
  end

  if nflowlim>0
    println("Number of lines with flow limits: ", nflowlim)
    if kind == :D
      optvec=zeros(2*nbus+2*ngen)
    elseif kind == :S
      optvec=zeros(4*nbus+2*ngen)
    end
    optvec[1:ngen]=PG
    optvec[ngen+1:2*ngen]=QG
    optvec[2*ngen+1:2*ngen+nbus]=VM
    optvec[2*ngen+nbus+1:2*ngen+2*nbus]=VA
    if kind == :S
      optvec[2*ngen+2*nbus+1:2*ngen+3*nbus]=PD
      optvec[2*ngen+3*nbus+1:2*ngen+4*nbus]=QD
    end

    d = JuMP.NLPEvaluator(opfmodel)
    MathProgBase.initialize(d, [:Jac])
    if kind == :D
      consRhs = zeros(2*nbus+2*nflowlim)
    elseif kind ==:S
      consRhs = zeros(4*nbus+2*nflowlim)
    end
    MathProgBase.eval_g(d, consRhs, optvec)


    #println(consRhs)

    @printf("================ Lines within %d pct of flow capacity ===================\n", within)
    println("Line   From Bus    To Bus    At capacity")

    nlim=1
    for l in 1:nline
      if lines[l].rateA!=0 && lines[l].rateA<1.0e10
        flowmax=(lines[l].rateA/baseMVA)^2
        idx = 2*nbus+nlim

        if( (consRhs[idx]+flowmax)  >= (1-within/100)^2*flowmax )
          @printf("%3d      %3d      %3d        %5.3f pct\n", l, lines[l].from, lines[l].to, 100*sqrt((consRhs[idx]+flowmax)/flowmax))
          #@printf("%7.4f   %7.4f    %7.4f \n", consRhs[idx], consRhs[idx]+flowmax,  flowmax)
        end
        nlim += 1
      end
    end
  end

  return
end
function acopf_outputAll(M::OPFModel, opfdata::OPFData); return acopf_outputAll(M.m, M.kind, opfdata); end

# Compute initial point for IPOPT based on the values provided in the case data
function acopf_initialPt_IPOPT(opfdata::MPCCases.OPFData)
  Pg=zeros(length(opfdata.generators)); Qg=zeros(length(opfdata.generators)); i=1
  for g in opfdata.generators
    # set the power levels in in between the bounds as suggested by matpower
    # (case data also contains initial values in .Pg and .Qg - not used with IPOPT)
    Pg[i]=0.5*(g.Pmax+g.Pmin)
    Qg[i]=0.5*(g.Qmax+g.Qmin)
    i=i+1
  end
  @assert i-1==length(opfdata.generators)

  Vm=zeros(length(opfdata.buses)); i=1;
  for b in opfdata.buses
    # set the ini val for voltage magnitude in between the bounds
    # (case data contains initials values in Vm - not used with IPOPT)
    Vm[i]=0.5*(b.Vmax+b.Vmin);
    i=i+1
  end
  @assert i-1==length(opfdata.buses)

  # set all angles to the angle of the reference bus
  Va = opfdata.buses[opfdata.bus_ref].Va * ones(length(opfdata.buses))

  return Pg,Qg,Vm,Va
end

function get_idx_and_id(m::JuMP.Model, k::Symbol, idx_to_id::Dict, offset::Int64)
  IDX = [linearindex(e) for e in getindex(m, k)]
  ID = IDX[[idx_to_id[i] for i in eachindex(IDX)]] .- offset
  return (idx=IDX, id=ID)
end
function filter_namedtuple(NT::NamedTuple, ids)
  idmask = findall(x -> x ∈ ids, NT.id)
  return (idx = NT.idx[idmask], id = NT.id[idmask])
end

function nonunique(x::AbstractArray{T}) where T
    xs = sort(x)
    duplicatedvector = T[]
    for i=2:length(xs)
        if (isequal(xs[i],xs[i-1]) && (length(duplicatedvector)==0 || !isequal(duplicatedvector[end], xs[i])))
            push!(duplicatedvector,xs[i])
        end
    end
    duplicatedvector
end

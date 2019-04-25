################## PCD (Probability Constrained Dispatch) ################
# Formulation of probability constrained optimal power flow (OPF) problems.
#
# cc_acopf_modelB_v1.jl is an initial implementation to compute a JuMP model that uses
# equality constraints to represent sensitivities. This method is related to
# cc_acopf_model.jl (model A), and uses the probability constraints from that function. The
# sensitivity conditions are
#
# [partial F / partial y] [dy / dx] = -[partial F / partial x].
#
# Here F represents the power-flow equations and
# y = [Qg^L Vm^L Va], x = [Pd Qd]. Superscript '^L' refers to load busses. New
# variables for this model are contained in Gamma := [dy / dx].
#
# This function uses nonlinear expressions to test the algebraic representation
# of columns.
###########################################################################
# 04/24/19, J.B.

function cc_acopf_modelB_v1(opf_data, options::Dict=Dict(), data::Dict=Dict()) # , options::Dict=Dict(), data::Dict=Dict()

  ## parse options
  #print_level 		= 0;
 #lossless       = haskey(options, :lossless)       ? options[:lossless]       : false
 #current_rating = haskey(options, :current_rating) ? options[:current_rating] : false
 #epsilon_Vm     = haskey(options, :epsilon_Vm)     ? options[:epsilon_Vm]     : 0.05
 epsilon_Va     = haskey(options, :epsilon_Va)     ? options[:epsilon_Va]     : 0.05
 #epsilon_Qg     = haskey(options, :epsilon_Qg)     ? options[:epsilon_Qg]     : 0.05
 #γ              = haskey(options, :gamma)          ? options[:gamma]          : 1.0
 #relax_Gamma    = haskey(options, :relax_Gamma)    ? options[:relax_Gamma]    : false
 #print_level    = haskey(options, :print_level)    ? options[:print_level]    : 0
 #xtilde         = haskey(options, :xtilde)         ? options[:xtilde]         : true
 #Gamma_type     = haskey(options, :Gamma_type)     ? options[:Gamma_type]     : :d
 # if lossless && !current_rating
 #   println("warning: lossless assumption requires `current_rating` instead of `power_rating`\n")
 #   current_rating = true
 # end

  # shortcuts for compactness
  lines = opf_data.lines; buses = opf_data.buses; generators = opf_data.generators; baseMVA = opf_data.baseMVA
  busIdx = opf_data.BusIdx; FromLines = opf_data.FromLines; ToLines = opf_data.ToLines; BusGeners = opf_data.BusGenerators;

  nbus = length(buses); nline = length(lines); ngen = length(generators); nload = length(findall(buses.bustype .== 1))
  # @assert(nload + ngen == nbus); NOT assert bc some buses can have more than one generator...

  # branch admitances
  YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA)

  ## parse data
  # Sig_d    	= Matrix(Diagonal(ones(2nbus)))
  # Va_min 	= -pi * ones(nbus)
  # Va_max 	= pi * ones(nbus)

  Sig_d    	= haskey(data, :Sigma_d) ? data[:Sigma_d] : Matrix(Diagonal(ones(2nbus)))
  Va_min 	= haskey(data, :Va_min)  ? data[:Va_min]  : -pi * ones(nbus)
  Va_max 	= haskey(data, :Va_max)  ? data[:Va_max]  :  pi * ones(nbus)
  Z 		= Normal(0,1)

  #
  # model
  #
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    opfmodel = Model(with_optimizer(Ipopt.Optimizer))
  else
    opfmodel = Model(solver=IpoptSolver(print_level=1))
  end

  @variable(opfmodel, generators[i].Pmin <= Pg[i=1:ngen] <= generators[i].Pmax)
  @variable(opfmodel, generators[i].Qmin <= Qg[i=1:ngen] <= generators[i].Qmax)

  @variable(opfmodel, buses[i].Vmin <= Vm[i=1:nbus] <= buses[i].Vmax)
  @variable(opfmodel, Va_min[i]		<= Va[i=1:nbus] <= Va_max[i])
  #@variable(opfmodel, Va[1:nbus])

  ## assumes no buses have generator and load
  @variable(opfmodel, buses[i].Pd/baseMVA <= Pd[i=1:nbus] <= buses[i].Pd/baseMVA)
  @variable(opfmodel, buses[i].Qd/baseMVA <= Qd[i=1:nbus] <= buses[i].Qd/baseMVA)

  # Additional (array) variable that represents sensitivities.
  # GM := [dy / dx]
  @variable(opfmodel, GM[1:2*(nbus-1),1:(2*nload)]);

  #fix the voltage angle at the reference bus
  if "19" ∈ split(string(Pkg.installed()["JuMP"]), ".")
    set_lower_bound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)
    set_upper_bound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)
  else
    setlowerbound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)
    setupperbound(Va[opf_data.bus_ref], buses[opf_data.bus_ref].Va)
  end

  @NLobjective(opfmodel, Min, sum( generators[i].coeff[generators[i].n-2]*(baseMVA*Pg[i])^2
			             +generators[i].coeff[generators[i].n-1]*(baseMVA*Pg[i])
				     +generators[i].coeff[generators[i].n  ] for i=1:ngen))

  #
  # power flow balance
  #
  #real part
  @NLconstraint(opfmodel, P[b=1:nbus],
    ( sum( YffR[l] for l in FromLines[b]) + sum( YttR[l] for l in ToLines[b]) + YshR[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[busIdx[lines[l].from]]*( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Pg[g] for g in BusGeners[b]) - sum(baseMVA*Pd[l] for l in busIdx[b]) ) / baseMVA      # Sbus part
    ==0)
  #imaginary part
  @NLconstraint(opfmodel, Q[b=1:nbus],
    ( sum(-YffI[l] for l in FromLines[b]) + sum(-YttI[l] for l in ToLines[b]) - YshI[b] ) * Vm[b]^2
    + sum( Vm[b]*Vm[busIdx[lines[l].to]]  *(-YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] )
    + sum( Vm[b]*Vm[busIdx[lines[l].from]]*(-YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b]   )
    - ( sum(baseMVA*Qg[g] for g in BusGeners[b]) - sum(baseMVA*Qd[l] for l in busIdx[b]) ) / baseMVA      #Sbus part
    ==0)
  #
  # branch/lines flow limits
  #
  # 04/01/19, J.B., Comment function calls
  @constraintref F_fr[1:nline]  ## from bus, TODO: won't work in JuMP v0.19
  @constraintref F_to[1:nline]  ## to bus, TODO: won't work in JuMP v0.19
  nlinelim=0
  for l in 1:nline
    if lines[l].rateA!=0 && lines[l].rateA<1.0e10
      nlinelim += 1
      flowmax=(lines[l].rateA/baseMVA)^2

      #branch apparent power limits (from bus)
      Yff_abs2=YffR[l]^2+YffI[l]^2; Yft_abs2=YftR[l]^2+YftI[l]^2
      Yre=YffR[l]*YftR[l]+YffI[l]*YftI[l]; Yim=-YffR[l]*YftI[l]+YffI[l]*YftR[l]
      F_fr[l] = @NLconstraint(opfmodel,
	              Vm[busIdx[lines[l].from]]^2 *
              	( Yff_abs2*Vm[busIdx[lines[l].from]]^2 + Yft_abs2*Vm[busIdx[lines[l].to]]^2
              	  + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
              	)
                - flowmax <=0)

      #branch apparent power limits (to bus)
      Ytf_abs2=YtfR[l]^2+YtfI[l]^2; Ytt_abs2=YttR[l]^2+YttI[l]^2
      Yre=YtfR[l]*YttR[l]+YtfI[l]*YttI[l]; Yim=-YtfR[l]*YttI[l]+YtfI[l]*YttR[l]
      F_to[l] = @NLconstraint(opfmodel,
        	      Vm[busIdx[lines[l].to]]^2 *
                ( Ytf_abs2*Vm[busIdx[lines[l].from]]^2 + Ytt_abs2*Vm[busIdx[lines[l].to]]^2
                  + 2*Vm[busIdx[lines[l].from]]*Vm[busIdx[lines[l].to]]*(Yre*cos(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]])-Yim*sin(Va[busIdx[lines[l].from]]-Va[busIdx[lines[l].to]]))
                )
                - flowmax <=0)
    end
  end
  JuMP.registercon(opfmodel, :F_fr, F_fr)
  JuMP.registercon(opfmodel, :F_to, F_to)

  ## Jacobian columns w.r.t. Vm, Va
  # jac = 	| dFP/dVm dFP/dVa | = 	| jacVmP jacVaP |
  #			| dFQ/dVm dFQ/dVa | 	| jacVmQ jacVaQ |

	# Real parts
  # VmP (Voltage magnitudes)
  @NLexpression(opfmodel,jacVmP[b=1:nbus,j=1:nbus],
	sum( 2*(YshR[b] * Vm[b]) for l in 1:1 if b == j) +
	sum( 2*(YffR[l] * Vm[b]) + Vm[busIdx[lines[l].to]] * ( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]])) for l in FromLines[b] if b == j) +
	sum( Vm[b] * ( YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) + YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] if busIdx[lines[l].to] == j ) +
	sum( 2*(YttR[l] * Vm[b]) + Vm[busIdx[lines[l].from]] * ( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]])) for l in ToLines[b] if b == j ) +
	sum( Vm[b] * ( YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) + YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]]  )) for l in ToLines[b] if busIdx[lines[l].from] == j ))

	# VaP, Voltage angles for real parts
	@NLexpression(opfmodel,jacVaP[b=1:nbus,j=1:nbus],
	sum( Vm[b]*Vm[busIdx[lines[l].to]] * ( -YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) + YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] if b == j) +
	sum( Vm[b] * Vm[busIdx[lines[l].to]] * ( YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) -YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] if busIdx[lines[l].to] == j) +
	sum( Vm[b]*Vm[busIdx[lines[l].from]] * ( -YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) + YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )) for l in ToLines[b] if b == j ) +
	sum( Vm[b] * Vm[busIdx[lines[l].from]] * ( YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) -YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )) for l in ToLines[b] if busIdx[lines[l].from] == j) )

	# Imaginary parts
	# VmQ (Voltage magnitudes)
	@NLexpression(opfmodel,jacVmQ[b=1:nbus,j=1:nbus],
	sum( 2*(-YshI[b])*Vm[b] for l in 1:1 if b == j ) +
	sum( 2*(-YffI[l]*Vm[b]) + Vm[busIdx[lines[l].to]] * ( -YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] if b == j) +
	sum( Vm[b] * ( -YftI[l]*cos(Va[b]-Va[busIdx[lines[l].to]] ) + YftR[l]*sin(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] if busIdx[lines[l].to] == j ) +
	sum( 2*(-YttI[l])*Vm[b] + Vm[busIdx[lines[l].from]] * ( -YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]]  )) for l in ToLines[b] if b == j) +
	sum( Vm[b] * ( -YtfI[l]*cos(Va[b]-Va[busIdx[lines[l].from]] ) + YtfR[l]*sin(Va[b]-Va[busIdx[lines[l].from]]  )) for l in ToLines[b] if busIdx[lines[l].from] == j ) )

	# VaQ
	@NLexpression(opfmodel,jacVaQ[b=1:nbus,j=1:nbus],
	sum( Vm[b]*Vm[busIdx[lines[l].to]] * ( YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) + YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] if b == j ) +
	sum( Vm[b] * Vm[busIdx[lines[l].to]]* ( -YftI[l]*sin(Va[b]-Va[busIdx[lines[l].to]] ) -YftR[l]*cos(Va[b]-Va[busIdx[lines[l].to]]  )) for l in FromLines[b] if busIdx[lines[l].to] == j ) +
	sum( Vm[b]*Vm[busIdx[lines[l].from]] * ( YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) + YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )) for l in ToLines[b] if b == j) +
	sum( Vm[b] * Vm[busIdx[lines[l].from]] * ( -YtfI[l]*sin(Va[b]-Va[busIdx[lines[l].from]] ) -YtfR[l]*cos(Va[b]-Va[busIdx[lines[l].from]]  )) for l in ToLines[b] if busIdx[lines[l].from] == j))

	# Assume model is correct until this point

	# PF: ... - (sum(Qg) - Qd)
	# Equality constraints on sensitivities
	busidx 		= 1:nbus; #ldidx = Array{Int}(undef,nload);
	ldidx 		= busidx[buses.bustype.==1];
	gidx 		= busidx[buses.bustype.==2];
	refidx 		= opf_data.bus_ref;
	nrefidx 	= busidx[buses.bustype.!=3];

	o2nb 	= ones(2*nbus,1);
	I2nb 	= diagm(0 => o2nb[:]);

	#rhs1_ 	= I2nb[nrefidx,ldidx];
	rhs1_ 	= I2nb[1:nbus,ldidx];
	#rhs2_ 	= I2nb[ngen.+nrefidx,ngen.+ldidx];
	rhs2_ 	= I2nb[(nbus+1):(2*nbus),nbus.+ldidx];

	#zl 		= zeros(length(nrefidx),nload);
	zl 		= zeros(nbus,nload);

	rhs1 	= cat(rhs1_,zl,dims=2);
	rhs2 	= cat(zl,rhs2_,dims=2);

	dFdQ 	= I2nb[(nbus+1):2*nbus,nbus.+gidx];


	# Dictionaries to map Jacobian columns to Gamma rows
	QidxGam 	= Dict{Int}{Int}();
	VmidxGam 	= Dict{Int}{Int}();
	VaidxGam 	= Dict{Int}{Int}();

	for i in 1:length(gidx)
		QidxGam[gidx[i]] = i;
	end
	for i in 1:nload
		VmidxGam[ldidx[i]] 	= length(gidx) + i;
	end
	for i in 1:length(nrefidx)
		VaidxGam[nrefidx[i]] = nload + length(gidx) + i;
	end

	# Real sensitivity constraints
	for i in 1:nbus
		if i != refidx
			for j in 2*nload

				@NLconstraint(opfmodel,
				sum(jacVmP[i,l] * GM[VmidxGam[l],j] for l in ldidx ) +
				sum(jacVaP[i,l] * GM[VaidxGam[l],j] for l in nrefidx ) + rhs1[i,j] == 0)

			end
		end
	end

	# Imaginary sensitivity constraints
	for i in 1:nbus
		if i != refidx
			for j in 2*nload

				@NLconstraint(opfmodel,
				sum(-dFdQ[i,QidxGam[l]] * GM[QidxGam[l],j] for l in gidx ) +
				sum(jacVmQ[i,l] * GM[VmidxGam[l],j] for l in ldidx ) +
				sum(jacVaQ[i,l] * GM[VaidxGam[l],j] for l in nrefidx ) + rhs2[i,j] == 0)

			end
		end
	end

	sigidx 						= Array{Int}(undef,2*nload);
	sigidx[1:nload] 			= ldidx[1:nload];
	sigidx[(nload+1):2*nload] 	= nbus .+ ldidx[1:nload];

	Sig_r 						= diag(Sig_d[sigidx,sigidx], 0);

	# Chance constraints
	# Va max (cf. cc_acopf_model.jl)
	#@constraintref cc_Va_max[1:length(nrefix)]
	gamma = 1.0;
	for l in nrefidx

		Val 	= Va[l];
		vareps 	= 1.0 - (gamma * epsilon_Va) / (nbus-1);
		q 		= quantile(Z,vareps);

		@NLconstraint(opfmodel, sum( GM[VaidxGam[l],k]*GM[VaidxGam[l],k]*Sig_r[k] for k in 1:(2*nload) ) <= (Va_max[l] - Val) / q )
		# cc_Va_max[VaidxGam[l]] = @NLconstraint(opfmodel, sum(GM(VaidxGam[l],k) for k in 1:2*nload ) <= (Va_max[l] - Val) / q )

	end
	#JuMP.registercon(opfmodel, :cc_Va_max, cc_Va_max)

	# Va min (cf. cc_acopf_model.jl)
	gamma = 1.0;
	for l in nrefidx

		Val 	= Va[l];
		vareps 	= (gamma * epsilon_Va) / (nbus-1);
		q 		= quantile(Z,vareps);

		@NLconstraint(opfmodel, sum( GM[VaidxGam[l],k]*GM[VaidxGam[l],k]*Sig_r[k] for k in 1:(2*nload) ) <= (Va_min[l] - Val) / q )

	end

  @printf("Buses: %d  Lines: %d  Generators: %d\n", nbus, nline, ngen)
  println("Lines with limits  ", nlinelim)

  return OPFModel(opfmodel, :InitData, :S)
end

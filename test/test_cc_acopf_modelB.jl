################## PCD (Probability Constrained Dispatch) ################
# Formulation of probability constrained optimal power flow (OPF) problems.
#
# test_cc_acopf_modelB.jl is intended as an alternative approach for a
# chance constrained model.
###########################################################################
# 04/24/19, J.B.


using MPCCases
using JuMP, JuMPUtil, Ipopt, MathProgBase, Distributions
using SparseArrays, LinearAlgebra

include("../src/OPF.jl")
#include("../src/sopfmodel_jac_test.jl")
include("../src/sopfmodel_jac_testg.jl")
include("../src/cc_acopf_modelB_v1.jl")
include("jac_analytic.jl")
include("jac_analytic_sumif.jl")

const path = "/Users/johannesbrust/Dropbox/ANL/projects/LOAD_DYNAMIC/code/OPF_DEV/test/cases"
const case = "case9"

opfdata     = load_case(case, path, other=false);

lines       = opfdata.lines;
buses       = opfdata.buses;
generators  = opfdata.generators;
baseMVA     = opfdata.baseMVA;
busIdx      = opfdata.BusIdx;
FromLines   = opfdata.FromLines;
ToLines     = opfdata.ToLines;
BusGeners   = opfdata.BusGenerators;

nbus        = length(buses);
nline       = length(lines);
ngen        = length(generators);

jac         = zeros(2*nbus,2*ngen+4*nbus);

YffR,YffI,YttR,YttI,YftR,YftI,YtfR,YtfI,YshR,YshI = computeAdmitances(lines, buses, baseMVA)

# Jacobian using MPB
sm      = sacopf_model(opfdata)
sm      = acopf_solve(sm, opfdata)

sj      = sacopf_model_jac_testg(opfdata)
sj      = acopf_solve(sj,opfdata)

# CC Model
cm 		= cc_acopf_modelB_v1(opfdata)
cm 		= acopf_solve(cm,opfdata)


# nvars   = 2*ngen+4*nbus
# ncons   = 2*nbus
#
# # Jump model and computation of Jacobian
# # Uses MathProgBase
# jmdl    = sm.m
# d       = JuMP.NLPEvaluator(jmdl)
#
# MathProgBase.initialize(d,[:Jac])
#
# # Variables
# # Pg, Qg, Vm, Va, Pd, Qd
# xopt              = zeros(nvars)
# PG                = getvalue(getindex(jmdl,:Pg))
# QG                = getvalue(getindex(jmdl,:Qg))
# VM                = getvalue(getindex(jmdl,:Vm))
# VA                = getvalue(getindex(jmdl,:Va))
# PD                = getvalue(getindex(jmdl,:Pd))
# QD                = getvalue(getindex(jmdl,:Qd))
#
# idxs              = 1
# idxe              = ngen
# xopt[idxs:idxe]   = PG; idxs = idxe + 1; idxe = idxs + ngen - 1
# xopt[idxs:idxe]   = QG; idxs = idxe + 1; idxe = idxs + nbus - 1
# xopt[idxs:idxe]   = VM; idxs = idxe + 1; idxe = idxs + nbus - 1
# xopt[idxs:idxe]   = VA; idxs = idxe + 1; idxe = idxs + nbus - 1
# xopt[idxs:idxe]   = PD; idxs = idxe + 1; idxe = idxs + nbus - 1
# xopt[idxs:idxe]   = QD
#
# # Indexing
# busidx 	= 1:nbus;
# ldidx 	= busidx[buses.bustype.==1];
# gidx 	= busidx[buses.bustype.==2];
# refidx 	= opfdata.bus_ref;
# nrefidx = busidx[buses.bustype.!=3];
#
# o2nb 	= ones(2*nbus,1);
# I2nb 	= diagm(0 => o2nb[:]);
#
# 	#rhs1_ 	= I2nb[nrefidx,ldidx];
# rhs1_ 	= I2nb[1:nbus,ldidx];
# 	#rhs2_ 	= I2nb[ngen.+nrefidx,ngen.+ldidx];
# rhs2_ 	= I2nb[(nbus+1):(2*nbus),nbus.+ldidx];
#
# 	#zl 		= zeros(length(nrefidx),nload);
# nload 	= length(ldidx);
# zl 		= zeros(nbus,nload);
#
# rhs1 	= cat(rhs1_,zl,dims=2);
# rhs2 	= cat(zl,rhs2_,dims=2);
#
# dFdQ 	= I2nb[(nbus+1):2*nbus,nbus.+gidx];



# Jacobian computations
#full_jac          = zeros(2*ncons,nvars)
# (Idxi,Idxj)       = MathProgBase.jac_structure(d)
# jac_v             = zeros(length(Idxi))
#
# MathProgBase.eval_jac_g(d,jac_v,xopt)
#
# sjac              = sparse(Idxi,Idxj,jac_v,2*ncons,nvars) # sparse Jacobian
# jac_mpb           = Array(sjac) # dense Jacobian
#
# # Analytic Jacobian
# jac_a             = jac_analytic(opfdata,PG,QG,VM,VA,PD,QD);
#
# # Analytic Jacobian (sum-if implementation)
# jac_a_sif         = jac_analytic_sumif(opfdata,PG,QG,VM,VA,PD,QD);
#
# # Column extractions Vm, Va
#
# colscmp           = (2*ngen+1):(2*(ngen+nbus));
# cols_mpb          = jac_mpb[1:ncons,colscmp];
#
# cols_a            = jac_a[:,colscmp];
#
# err               = cols_mpb-cols_a;
#
# cols_asif         = jac_a_sif[:,colscmp];
#
# # Jacobian values
# #sjm               = sj.m;
# #jac_v_p           = getvalue(getindex(sjm,:jacp));
#
# jmdl    = sj.m
# d_      = JuMP.NLPEvaluator(jmdl)
#
# MathProgBase.initialize(d_,[:Jac])
#
# g               = zeros(2*ncons+4*nbus*nbus,1);
#
# MathProgBase.eval_g(d_,g,xopt);
#
# jacVmP              = zeros(nbus,nbus);
# idxs                =0;
# idxe                =0;
#
# jacVaP              = zeros(nbus,nbus);
# idxsa               =0;
# idxea               =0;
#
# jacVmQ              = zeros(nbus,nbus);
# idxsmq              =0;
# idxemq              =0;
#
# jacVaQ              = zeros(nbus,nbus);
# idxsaq              =0;
# idxeaq              =0;
#
# for k=1:nbus
#
#     idxs = 2*ncons +(k-1)*nbus + 1;
#     idxe = 2*ncons +k*nbus;
#
#     #jacJP[:,k] = g[idxs:idxe,1];
#     jacVmP[k,:] = g[idxs:idxe,1];
#
#     idxsa = 2*ncons + nbus*nbus +(k-1)*nbus + 1;
#     idxea = 2*ncons + nbus*nbus + k*nbus;
#
#     #jacJP[:,k] = g[idxs:idxe,1];
#     jacVaP[k,:] = g[idxsa:idxea,1];
#
#     idxsmq = 2*(ncons + nbus*nbus) + (k-1)*nbus + 1;
#     idxemq = 2*(ncons + nbus*nbus) + k*nbus;
#
#     #jacJP[:,k] = g[idxs:idxe,1];
#     jacVmQ[k,:] = g[idxsmq:idxemq,1];
#
#     idxsaq = 2*(ncons + nbus*nbus) + nbus*nbus + (k-1)*nbus + 1;
#     idxeaq = 2*(ncons + nbus*nbus) + nbus*nbus + k*nbus;
#
#     #jacJP[:,k] = g[idxs:idxe,1];
#     jacVaQ[k,:] = g[idxsaq:idxeaq,1];
#
# end
#
# colscmp1            = (2*ngen+1):(2*ngen+nbus);
# jac_a_cmp           = jac_a[1:nbus,colscmp1];
# jac_mpb_cmp         = jac_mpb[1:nbus,colscmp1];
#
# colscmp2            = (2*ngen+nbus+1):(2*(ngen+nbus));
# jac_a_cmp2          = jac_a[1:nbus,colscmp2];
# jac_mpb_cmp2        = jac_mpb[1:nbus,colscmp2];
#
# #colscmp1             = (2*ngen+1):(2*ngen+nbus);
# jac_a_cmpq           = jac_a[(nbus+1):2*nbus,colscmp1];
# jac_mpb_cmpq         = jac_mpb[(nbus+1):2*nbus,colscmp1];
#
# #colscmp2             = (2*ngen+nbus+1):(2*(ngen+nbus));
# jac_a_cmp2q          = jac_a[(nbus+1):2*nbus,colscmp2];
# jac_mpb_cmp2q        = jac_mpb[(nbus+1):2*nbus,colscmp2];

#obj_    = MathProgBase.eval_f(d_,xopt)

#MathProgBase.initialize(d,[:Jac])


# Calling test

#Pg = ones(ngen,1);
#Qg = ones(ngen,1);
#Vm = ones(nbus,1);
#Va = ones(nbus,1);
#Pd = ones(nbus,1);
#Qd = ones(nbus,1);

# Debugging Juno.@enter
#jac_analytic(opfdata,Pg,Qg,Vm,Va,Pd,Qd)

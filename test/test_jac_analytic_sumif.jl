################## PCD (Probability Constrained Dispatch) ################
# Formulation of probability constrained optimal power flow (OPF) problems.
#
# test_jac_analytic_sumif.jl is an initial script to compute constraint Jacobians using
# derivatives of the power flow equations based on conditions that resemble
# 'sum if' statements.
###########################################################################
# 04/19/19, J.B.


using MPCCases
using JuMP, JuMPUtil, Ipopt, MathProgBase
using SparseArrays, LinearAlgebra

include("../src/OPF.jl")
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

nvars   = 2*ngen+4*nbus
ncons   = 2*nbus

# Jump model and computation of Jacobian
# Uses MathProgBase
jmdl    = sm.m
d       = JuMP.NLPEvaluator(jmdl)

MathProgBase.initialize(d,[:Jac])

# Variables
# Pg, Qg, Vm, Va, Pd, Qd
xopt              = zeros(nvars)
PG                = getvalue(getindex(jmdl,:Pg))
QG                = getvalue(getindex(jmdl,:Qg))
VM                = getvalue(getindex(jmdl,:Vm))
VA                = getvalue(getindex(jmdl,:Va))
PD                = getvalue(getindex(jmdl,:Pd))
QD                = getvalue(getindex(jmdl,:Qd))

idxs              = 1
idxe              = ngen
xopt[idxs:idxe]   = PG; idxs = idxe + 1; idxe = idxs + ngen - 1
xopt[idxs:idxe]   = QG; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = VM; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = VA; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = PD; idxs = idxe + 1; idxe = idxs + nbus - 1
xopt[idxs:idxe]   = QD

# Jacobian computations
#full_jac          = zeros(2*ncons,nvars)
(Idxi,Idxj)       = MathProgBase.jac_structure(d)
jac_v             = zeros(length(Idxi))

MathProgBase.eval_jac_g(d,jac_v,xopt)

sjac              = sparse(Idxi,Idxj,jac_v,2*ncons,nvars) # sparse Jacobian
jac_mpb           = Array(sjac) # dense Jacobian

# Analytic Jacobian
jac_a             = jac_analytic(opfdata,PG,QG,VM,VA,PD,QD);

# Analytic Jacobian (sum-if implementation)
jac_a_sif         = jac_analytic_sumif(opfdata,PG,QG,VM,VA,PD,QD);

# Column extractions Vm, Va

colscmp           = (2*ngen+1):(2*(ngen+nbus));
cols_mpb          = jac_mpb[1:ncons,colscmp];

cols_a            = jac_a[:,colscmp];

err               = cols_mpb-cols_a;

cols_asif         = jac_a_sif[:,colscmp];


# Calling test

#Pg = ones(ngen,1);
#Qg = ones(ngen,1);
#Vm = ones(nbus,1);
#Va = ones(nbus,1);
#Pd = ones(nbus,1);
#Qd = ones(nbus,1);

# Debugging Juno.@enter
#jac_analytic(opfdata,Pg,Qg,Vm,Va,Pd,Qd)

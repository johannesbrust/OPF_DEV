################## PCD (Probability Constrained Dispatch) ################
# Formulation of probability constrained optimal power flow (OPF) problems.
#
# test_jul_syntax.jl
#
# Script to test Julia syntax on 'opf' data structures.
##########################################################################
# 12/04/19, J.B.

include("../src/OPF.jl")

const path = "/Users/johannesbrust/Dropbox/ANL/projects/LOAD_DYNAMIC/code/OPF/test/cases"
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

# Syntax tests

Vm = ones(nbus,1);
Va = ones(nbus,1);

b  = 1;

tl = ToLines[b];
fl = FromLines[b];

#e1 = 2*(sum(YffR[fl[:]]) + sum(YttR[tl[:]]) + YshR[b])*Vm[b] +
#    sum(Vm[busIdx[lines[fl[:]].to]] .*
#     ( YftR[fl[:]].*cos(Va[b].-Va[busIdx[lines[fl[:]].to]] ) .+
#        YftI[fl[:]].*sin(Va[b].-Va[busIdx[lines[fl[:]].to]]  )));

jac_    = 0;
jac_    = 2*(sum(YffR[fl[:]]) + sum(YttR[tl[:]]) + YshR[b])*Vm[b];

nfl     = length(fl);

for l in 1:nfl

        jac[b,b] +=  Vm[busIdx[lines[fl[l]].to]] *
        ( YftR[fl[l]]*cos(Va[b]-Va[busIdx[lines[fl[l]].to]] ) +
         YftI[fl[l]].*sin(Va[b]-Va[busIdx[lines[fl[l]].to]]  ));

end

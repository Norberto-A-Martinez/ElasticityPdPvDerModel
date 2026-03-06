/*
Author: Norberto Abrante Martinez

mixed integer non-linear programming model for the PV planning

*/
model;
/*DECLARAÇÃO DE CONJUNTOS*/
set N within {0..4000}; # EDS buses/nodes
set T within {1..24} ; # set of time periods
set S; # set of irradiation scenarios
set D; # set of power demand scenarios
set C := 1..3; # set of customer classes

/*DECLARAÇÃO DE PARÂMETROS*/
param delt_t := 1; # duration of period t in hours
param prob_S {S}; # probability of occurrence of scenario
param prob_D {D}; # probability of occurrence of scenario

param de_curve_type {T,C}; # load curve
param de_curve_scenarios {T,D}; # load curve
param de_curve {N,T,D}; # load curve

param Pd {N}; # active power demand
param Qd {N}; # reactive power demand
param PV_0 {N}; # initial PV generation at each node
param region {N}; # region identifier for each node
param irrad{T,S}; # distributed generation generation profile
param ET := 0.30988; # energy tariff
param SUT := 0.34389; # system usage tariff

param i := 0.03; # interest rate
param theta_pf := 0.97; # PV inverter power factor
param pv_panel_capacity := 0.7; # PV panel capacity kwp
param pv_inverter_capacity := 5; # PV inverter capacity kw
param pv_panel_cost := 100; # PV panel cost
param pv_inverter_cost := 1000; # PV inverter cost
param pv_installation_cost := 2500; # PV installation cost
param qqcz {T};

var ex_co {N} >= 0;
var ex_in {N} >= 0;
var ex_net {N} >= 0;
var an_savings {N} >= 0; # annual savings for the prosumer
var investment_pv {N} >=0;
var PV_panels {N} integer >= 0;
var PV_inverter {N} integer >= 0;
var PV_allocation {N} binary; # variable that indicates if the node n allocation

set M := 1..12 ordered; # set of months for energy credit accumulation
var E_cr {N,M} >= 0; # energy credits at each node for month
var cost_prosumer {N,M};
var cost_consumer {N,M};
var PV_gen {N,S,T} >= 0; # PV generation at each node
var PV_cut {N,S,T} >= 0; # PV curtailed power at each node
var qPV_gen {N,S,T}; # potência reativa fornecida pela geração distribuída no nó i
var E_in {N,S,D,T} >= 0;
var E_co {N,S,D,T} >= 0;
var pv_panel_cost_total = sum {n in N} (PV_panels[n]*pv_panel_cost);
var pv_inverter_cost_total = sum {n in N} (PV_inverter[n]*pv_inverter_cost);
var pv_installation_cost_total = sum {n in N} (PV_allocation[n]*pv_installation_cost);

# subject to investment_pv_limit {n in N}: 
#     investment_pv[n] <= 50e3;

subject to investment_pv_def {n in N}:
    investment_pv[n] = PV_panels[n]*pv_panel_cost + PV_inverter[n]*pv_inverter_cost + PV_allocation[n]*pv_installation_cost;

# tarifarry model
subject to liquid_power {n in N, s in S, d in D, t in T}:
    E_in[n,s,d,t] - E_co[n,s,d,t] = (PV_gen[n,s,t] - Pd[n]*de_curve[n,t,d])*delt_t;

subject to expected_net_consumption {n in N}:
    ex_net[n] = sum {t in T, s in S, d in D} (E_co[n,s,d,t] - E_in[n,s,d,t])*prob_S[s]*prob_D[d];

subject to expected_injection {n in N}:
    ex_in[n] = sum {t in T, s in S, d in D} E_in[n,s,d,t]*prob_S[s]*prob_D[d];

subject to expected_consumption {n in N}:
    ex_co[n] = sum {t in T, s in S, d in D} E_co[n,s,d,t]*prob_S[s]*prob_D[d];

subject to def_cost_pro {n in N, m in M}:
    cost_prosumer[n,m] = ET*max(0,ex_net[n] - E_cr[n,m]) + SUT*(0.9*ex_in[n] + ex_co[n]);

subject to def_cost_con {n in N, m in M}:
    cost_consumer[n,m] = (ET + SUT)*ex_co[n];

subject to def_energy_credit_ini {n in N, m in M: m == 1}:
    E_cr[n,m] = 0;

subject to def_energy_credit {n in N, m in M: m > 1}:
    E_cr[n,m] = E_cr[n,m-1] + max(0, ex_net[n]) - min(E_cr[n,m-1], max(0, ex_net[n]));

subject to savings_def {n in N}:
    an_savings[n] = sum {m in M} (cost_consumer[n,m] - cost_prosumer[n,m]);

subject to lpv_a {n in N}:
    investment_pv[n] <= sum{y in 1..5}an_savings[n]/(1 + i)^y;

subject to PV_panels_allocation {n in N}:
    PV_panels[n] <= 80*PV_allocation[n];

subject to pv_inverter_allocation {n in N}:
    PV_inverter[n] <= 20*PV_allocation[n];

subject to panels_generation {n in N, s in S, t in T}:
    PV_gen[n,s,t] = min(PV_panels[n]*pv_panel_capacity*irrad[t,s], PV_inverter[n]*pv_inverter_capacity);

subject to GD_reactive_1 {n in N, s in S, t in T}:
    - PV_gen[n,s,t]*tan(acos(theta_pf)) <= qPV_gen[n,s,t];

subject to GD_reactive_2 {n in N,  s in S, t in T}:
    qPV_gen[n,s,t] <= PV_gen[n,s,t]*tan(acos(theta_pf));

minimize FO: sum {n in N, m in M} cost_prosumer[n,m];

data;
param:	
T:	qqcz   :=
1	1
2	1
3	1
4	1
5	1
6	1
7	1
8	1
9	1
10	1
11	1
12	1
13	1
14	1
15	1
16	1
17	1
18	1
19	1
20	1
21	1
22	1
23	1
24	1
;

param: S: prob_S:=
	1	0.4
	2	0.3
	3	0.3
	# 4	0.35
;

param irrad :
    1	    2	    3:=
1   0.00	0.00	0.00
2   0.00	0.00	0.00
3   0.00	0.00	0.00
4   0.00	0.00	0.00
5   0.00	0.00	0.00
6   0.04	0.00	0.00
7   0.18	0.05	0.02
8   0.37	0.17	0.11
9   0.55	0.23	0.21
10  0.71	0.27	0.07
11  0.84	0.39	0.09
12  0.94	0.39	0.09
13  0.99	0.54	0.06
14  0.99	0.58	0.12
15  0.95	0.66	0.21
16  0.83	0.45	0.17
17  0.54	0.39	0.09
18  0.38	0.20	0.13
19  0.09	0.09	0.05
20  0.00	0.00	0.00
21  0.00	0.00	0.00
22  0.00	0.00	0.00
23  0.00	0.00	0.00
24  0.00	0.00	0.00
;

param: D: prob_D:=
	1	0.4
	2	0.3
	3	0.3
;

param de_curve_scenarios:	
    1           2         3 :=
1	0.5924      0.4668    0.88758
2	0.6536      0.4383    0.82992
3	0.7245      0.3459    0.65400
4	0.8101      0.3195    0.58590
5	0.8618      0.2848    0.48837
6	0.8579      0.2848    0.47964
7	0.6207      0.2848    0.56738
8	0.5542      0.2848    0.58853
9	0.4444      0.2848    0.56856
10	0.3365      0.2848    0.57760
11	0.2848      0.3048    0.57920
12	0.2848      0.3382    0.51521
13	0.3034      0.3992    0.53652
14	0.3506      0.4788    0.57943
15	0.4481      0.4686    0.60439
16	0.4384      0.3861    0.59670
17	0.5283      0.4117    0.68099
18	0.5455      0.4373    0.68766
19	0.5771      0.4263    0.70855
20	0.5711      0.4867    0.84150
21	0.6722      0.5527    0.89919
22	0.6561      0.5968    0.89718
23	0.6278      0.5538    0.78955
24	0.5944      0.4990    0.74354
;

param:
N: 	Pd    Qd    region :=
0	0     0     1
1	100   60    1
2	90    40    2
3	120   80    1
4	60    30    3
5	60    20    3
6	200   100   1
7	200   100   1
8	60    20    3
9	60    20    3
10	45    30    3
11	60    35    3
12	60    35    3
13	120   80    1
14	60    10    3
15	60    20    3
16	60    20    3
17	90    40    3
18	90    40    2
19	90    40    2
20	90    40    2
21	90    40    2
22	90    50    2
23	420   200   1
24	420   200   1
25	60    25    3
26	60    25    3
27	60    20    3
28	120   70    1
29	200   600   1
30	150   70    1
31	210   100   1
32	60    40    3
;

load amplcsv.dll;
table tab_d1 IN "amplcsv" "input-data/profiles_long.csv": [T,C], de_curve_type;
read table tab_d1;

# recalculation and normalization of the load curve
for {t in T, n in N, d in D}{
    let de_curve[n,t,d] := de_curve_type[t,region[n]]*de_curve_scenarios[t,d];
}

param curve_max := max {t in T, n in N, d in D} de_curve[n,t,d];
param curve_min := min {t in T, n in N, d in D} de_curve[n,t,d];

for {t in T, n in N, d in D}{
    let de_curve[n,t,d] := (de_curve[n,t,d] - curve_min)/(curve_max - curve_min);
}

display de_curve;

option solver knitro;
option cplex_options "mipgap = 1 mip_multistart=1 gapabs=0.01 outlev=2 threads=12";

suffix iis OUT;
option cplex_options 'names=1 iis=1 conflictalg=4';

solve FO;

# display P_inj, P_cons;
# infisibility diagnosis
if solve_result = "infeasible" then {
    print "";
    print "==================================================";
    print "   IRREDUCIBLE INCONSISTENT SUBSYSTEM (IIS) FOUND ";
    print "==================================================";
    print "";
    
    print "--- CONFLICTING CONSTRAINTS ---";
    # Iterate through all constraints
    for {a in 1.._ncons} {
        if _con[a].iis != 0 then {
            # If IIS status is different from 0 (non-zero), show the name
            printf "Conflict in Constraint: %s (Status: %d)\n", _conname[a], _con[a].iis;
        }
    }
    
    print "";
    print "--- VARIABLES WITH CONFLICTING BOUNDS ---";
    # Iterate through all variables
    for {j in 1.._nvars} {
        if _var[j].iis != 0 then {
            # If IIS status is different from 0, show the name
            printf "Problem with Variable: %s (Status: %d)\n", _varname[j], _var[j].iis;
        }
    }
    print "==================================================";

print "--- BOUND VIOLATION DIAGNOSIS ---";

# Check LOWER bound violations (Var < Lower Bound - 0.001)
for {j in 1.._nvars} {
    if _var[j] < _var[j].lb - 0.001 then {
        printf "MIN VIOLATION: %s = %f (Minimum was %f)\n", _varname[j], _var[j], _var[j].lb;
    }
}

# Check UPPER bound violations (Var > Upper Bound + 0.001)
for {j in 1.._nvars} {
    if _var[j] > _var[j].ub + 0.001 then {
        printf "MAX VIOLATION: %s = %f (Maximum was %f)\n", _varname[j], _var[j], _var[j].ub;
    }
}
}


display PV_allocation, PV_panels, PV_inverter;
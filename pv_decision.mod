/*
Author: Norberto Abrante Martinez
mixed integer non-linear programming model for the PV planning
*/

model;
set N within {0..4000}; # EDS buses/nodes
set T within {1..24} ; # set of time periods
set S := 1..3; # set of irradiation scenarios
set D := 1..3; # set of power demand scenarios
set M := 1..12 ordered; # set of months for energy credit accumulation

param delt_t := 1; # duration of period t in hours
param prob_S {M,S}; # probability of occurrence of scenario
param prob_D {M,D}; # probability of occurrence of scenario
param curve_load {M,D,T}; # load curve for each month, demand scenario and time period
param curve_pv{M,S,T}; # distributed generation generation profile
param Pd {N}; # active power demand
param ET; # energy tariff
param SUT; # system usage tariff
param i := (1 + 0.05)^(1/12) - 1; # interest rate
param theta_pf := 0.97; # PV inverter power factor
param pv_panel_capacity := 0.7; # PV panel capacity kwp
param pv_inverter_capacity := 5; # PV inverter capacity kw
param pv_panel_cost := 500; # PV panel cost
param pv_inverter_cost := 7000; # PV inverter cost
param pv_installation_cost := 4000; # PV installation cost
param max_inverters_factor := 0.2; # panels per Pd unit allowed
param qqcz {T};
param taxa_fioB := 0.9;

var ex_co {N,M} >= 0;
var ex_in {N,M} >= 0;
var credit_used {N,M} >= 0;    # Créditos sacados do banco neste mês
var E_net {N,M,S,D,T}; # net energy balance
var E_cr {N,M} >= 0; # energy credits at each node for month
var E_paid {N,M} >= 0; # energy efectively paid
var e_inst_injet {N,M,S,D,T} >= 0;
var e_inst_consu {N,M,S,D,T} >= 0;
var an_savings {N}; # annual savings for the prosumer
var investment_pv {N};
var PV_panels {N} integer >= 0;
var PV_inverter {N} integer >= 0;
var PV_allocation {N} binary; # variable that indicates if the node n allocation
var cost_prosumer {N,M};
var cost_consumer {N,M};
var PV_gen {N,M,S,T} >= 0; # PV generation at each node

subject to investment_pv_def {n in N}:
    investment_pv[n] = PV_panels[n]*pv_panel_cost 
                     + PV_inverter[n]*pv_inverter_cost 
                     + PV_allocation[n]*pv_installation_cost;

subject to investmentminit {n in N}:
    investment_pv[n] <= 30e3;

# subject to PV_panels_allocation {n in N}:
#     PV_inverter[n] <= (max_inverters_factor * Pd[n]) * PV_allocation[n];

subject to panels_gen_limit_panel {n in N, m in M, s in S, t in T}:
    PV_gen[n,m,s,t] <= PV_panels[n]*pv_panel_capacity*curve_pv[m,s,t];

subject to panels_gen_limit_inverter {n in N, m in M, s in S, t in T}:
    PV_gen[n,m,s,t] <= PV_inverter[n]*pv_inverter_capacity;

subject to liquid_power {m in M, n in N, s in S, d in D, t in T}:
     e_inst_injet[n,m,s,d,t] - e_inst_consu[n,m,s,d,t] 
     = (PV_gen[n,m,s,t] - Pd[n]*curve_load[m,d,t])*delt_t;

var bin {N,M} binary; # for the monthly
var bin2 {N,M,S,D,T} binary; # for the instantaneous
var ex_liqu {N,M};
var inj_liqu {N,M} >= 0;
var con_liqu {N,M} >= 0;

subject to guarantee_split_a1 {m in M, n in N, s in S, d in D, t in T}:
    e_inst_consu[n,m,s,d,t] <= 1e7*bin2[n,m,s,d,t];

subject to guarantee_split_b1 {m in M, n in N, s in S, d in D, t in T}:
    e_inst_injet[n,m,s,d,t] <= 1e7*(1-bin2[n,m,s,d,t]);

subject to expected_injection {n in N, m in M}:
    ex_in[n,m] = 30.4 * sum {t in T, s in S, d in D} e_inst_injet[n,m,s,d,t]*prob_S[m,s]*prob_D[m,d];

subject to expected_consumption {n in N, m in M}:
    ex_co[n,m] = 30.4 * sum {t in T, s in S, d in D} e_inst_consu[n,m,s,d,t]*prob_S[m,s]*prob_D[m,d];

subject to expected_liquid_calc {n in N, m in M}:
    ex_liqu[n,m] = ex_in[n,m] - ex_co[n,m];

subject to injected_consumed_liqueid {n in N, m in M}:
    ex_liqu[n,m] = inj_liqu[n,m] - con_liqu[n,m];

subject to guarantee_split_a {n in N, m in M}:
    inj_liqu[n,m] <= 1e7*bin[n,m];

subject to guarantee_split_b {n in N, m in M}:
    con_liqu[n,m] <= 1e7*(1-bin[n,m]);

subject to first_month_credt_usage {n in N}:
    E_cr[n,first(M)] = E_cr[n,last(M)];

subject to def_energy_credit {n in N, m in M: m > 1}:
    E_cr[n,m] = E_cr[n,m-1] + inj_liqu[n,m] - credit_used[n,m];

subject to credid_usage_a {n in N, m in M}:
   credit_used[n,m] <= con_liqu[n,m];

subject to credid_usage_b {n in N, m in M: m > 1}:
   credit_used[n,m] <= E_cr[n,m-1];

subject to def_cost_prosumer {n in N, m in M}:
    cost_prosumer[n,m] = (ET + SUT) * (con_liqu[n,m] - credit_used[n,m]) 
            + SUT * taxa_fioB * ex_in[n,m] + con_liqu[n,m]*; # wire b on what was compensated by the prosumer

subject to def_cost_consumer {n in N, m in M}:
    cost_consumer[n,m] = (ET + SUT) * 30.4 * sum {t in T, s in S, d in D} Pd[n]*curve_load[m,d,t]*prob_S[m,s]*prob_D[m,d];

subject to def_savings {n in N}:
    an_savings[n] = sum {m in M} (cost_consumer[n,m] - cost_prosumer[n,m]);

subject to investiment_limt_by_PresentValue {n in N}:
    investment_pv[n] <= sum{m in 1..48} (cost_consumer[n,(m-1) mod 12+1] - cost_prosumer[n,(m-1) mod 12+1])/(1 + i)^m;

maximize FO: sum{n in N} (sum{m in 1..300} (cost_consumer[n,(m-1) mod 12+1] - cost_prosumer[n,(m-1) mod 12+1])/(1 + i)^m - investment_pv[n]);

data;

load amplcsv.dll;
table _T IN "amplcsv" "input-data/profiles.csv": T <- [T];
read table _T;

table _Pd IN "amplcsv" "input-data/nodes.csv": N <- [N], Pd;
read table _Pd;

table load_table IN "amplcsv" "input-data/parametros_load_curvas.csv": [M,D,T], curve_load;
read table load_table;

table pv_table IN "amplcsv" "input-data/parametros_pv_curvas.csv": [M,S,T], curve_pv;
read table pv_table;

table pv_prob IN "amplcsv" "input-data/parametros_pv_probs.csv": [M,S], prob_S;
read table pv_prob;

table load_prob IN "amplcsv" "input-data/parametros_load_probs.csv": [M,D], prob_D;
read table load_prob;

option solver cplex;
option cplex_options "mipgap = 0.1 threads=12";

let ET:= 0.00988;
let SUT:= 0.14389;

# for {fator in 1..20}{
#     display ET, SUT;
#     solve FO;
#     printf "\n\n";
#     # put in a csv file the results
#     printf "N,PV_allocation,PV_panels,PV_inverter\n" > ("results\pv-solution_"& fator &".csv");
#     for {n in N}{
#         printf "%d,%d,%d,%d\n", n, PV_allocation[n], PV_panels[n], PV_inverter[n] > ("results\pv-solution_"& fator &".csv");
#     }
#     let ET := ET + 0.01;
#     let SUT := SUT + 0.01;
# display cost_consumer - cost_prosumer;
# }

solve FO;

printf "N,PV_allocation,PV_panels,PV_inverter\n";
for {n in N}{
    printf "%d,%d,%d,%d\n", n, PV_allocation[n], PV_panels[n], PV_inverter[n];
}

display an_savings;


printf "month, \tprosumer, \tconsumer\n";
for {m in M}{
    printf "%d\t%10.4f\n",m ,cost_consumer[1,m] - cost_prosumer[1,m] ;
}

printf'\n';
printf "month,\t E_cr,\t credit_used,\t ex_co,\t ex_in,\t ex_liqu,\t inj_liqu,\t con_liqu\n";
for {m in M}{
    printf "%d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n",
    m , E_cr[1,m], credit_used[1,m], ex_co[1,m], ex_in[1,m], ex_liqu[1,m], inj_liqu[1,m], con_liqu[1,m];
}
printf "\n";
display ex_co[1,1], ex_in[1,1];

# printf '\tPV, \tPD, \tNet energy, \tnet injected, \tnet consumed\n';
# for {s in S, d in D}{
#     printf "s-%d d-%d\n", s,d;
#     for {t in T}{
#         printf "%d:\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\n", 
#                 t,PV_gen[1,1,s,t], Pd[1]*curve_load[1,d,t], E_net[1,1,s,d,t], e_inst_injet[1,1,s,d,t], e_inst_consu[1,1,s,d,t];
#     }
# }

# printf "power\n";
# for {t in T}{
#     printf "%d, %8.3f %8.3f\n", t, PV_gen[1,1,1,t], Pd[1]*curve_load[1,1,t];
# }
# printf "\n";

printf "\nn: invest, VPL 10\n";
for {n in N}{
    printf "%d:\t%8.4f\t%8.4f\n", n, investment_pv[n], sum{y in 1..10}an_savings[n]/(1 + i)^y;
}

display max{m in M, n in N, s in S, d in D, t in T} e_inst_consu[n,m,s,d,t];
display max{m in M, n in N, s in S, d in D, t in T} e_inst_injet[n,m,s,d,t];
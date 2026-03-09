/*
Author: Norberto Abrante Martinez

mixed integer non-linear programming model for the PV planning

*/
model;
/*DECLARAÇÃO DE CONJUNTOS*/
set N within {0..4000}; # EDS buses/nodes
set T within {1..24} ; # set of time periods
set S := 1..3; # set of irradiation scenarios
set D := 1..3; # set of power demand scenarios
set M := 1..12 ordered; # set of months for energy credit accumulation

/*DECLARAÇÃO DE PARÂMETROS*/
param delt_t := 1; # duration of period t in hours
param prob_S {M,S}; # probability of occurrence of scenario
param prob_D {M,D}; # probability of occurrence of scenario

param curve_load {M,D,T}; # load curve for each month, demand scenario and time period
param curve_pv{M,S,T}; # distributed generation generation profile

param Pd {N}; # active power demand
param Qd {N}; # reactive power demand
param PV_0 {N}; # initial PV generation at each node
param region {N}; # region identifier for each node
param ET := 0.30988; # energy tariff
param SUT := 0.44389; # system usage tariff

param i := 0.03; # interest rate
param theta_pf := 0.97; # PV inverter power factor
param pv_panel_capacity := 0.7; # PV panel capacity kwp
param pv_inverter_capacity := 5; # PV inverter capacity kw
param pv_panel_cost := 0; # PV panel cost
param pv_inverter_cost := 0; # PV inverter cost
param pv_installation_cost := 0; # PV installation cost
param qqcz {T};

var ex_co {N,M} >= 0;
var ex_in {N,M} >= 0;
# var ex_net {N,M};
var deficit {N,M} >= 0;        # Energia que precisou da rede no mês
var surplus {N,M} >= 0;        # Energia que sobrou no mês
var credit_used {N,M} >= 0;    # Créditos sacados do banco neste mês
var E_paid {N,M} >= 0;         # Energia que será faturada na tarifa ET
var an_savings {N} >= 0; # annual savings for the prosumer
var investment_pv {N} >=0;
var PV_panels {N} integer >= 0;
var PV_inverter {N} integer >= 0;
var PV_allocation {N} binary; # variable that indicates if the node n allocation

var E_cr {N,M} >= 0; # energy credits at each node for month
var cost_prosumer {N,M};
var cost_consumer {N,M};
var PV_gen {N,M,S,T} >= 0; # PV generation at each node
var PV_cut {N,M,S,T} >= 0; # PV curtailed power at each node
var qPV_gen {N,M,S,T}; # potência reativa fornecida pela geração distribuída no nó i
var E_in {N,M,S,D,T} >= 0;
var E_co {N,M,S,D,T} >= 0;
var pv_panel_cost_total = sum {n in N} (PV_panels[n]*pv_panel_cost);
var pv_inverter_cost_total = sum {n in N} (PV_inverter[n]*pv_inverter_cost);
var pv_installation_cost_total = sum {n in N} (PV_allocation[n]*pv_installation_cost);

# subject to investment_pv_limit {n in N}: 
#     investment_pv[n] <= 50e3;

# subject to expected_net_consumption {n in N, m in M}:
#     ex_net[n,m] = sum {t in T, s in S, d in D} (E_co[n,m,s,d,t] - E_in[n,m,s,d,t])*prob_S[m,s]*prob_D[m,d];

# subject to def_cost_pro {n in N, m in M}:
#     cost_prosumer[n,m] = ET*max(0,ex_net[n,m] - E_cr[n,m]) + SUT*(0.9*ex_in[n,m] + ex_co[n,m]);

# subject to def_energy_credit_ini {n in N, m in M}:
#     E_cr[n,m] = 0;

# subject to def_energy_credit {n in N, m in M: m > 1}:
#     E_cr[n,m] = E_cr[n,m-1] + max(0, ex_net[n,m]) - min(E_cr[n,m-1], max(0, ex_net[n,m]));

subject to PV_panels_allocation {n in N}:
    PV_panels[n] <= 80*PV_allocation[n];

subject to pv_inverter_allocation {n in N}:
    PV_inverter[n] <= 40*PV_allocation[n];

# subject to panels_generation {n in N, m in M, s in S, t in T}:
#     PV_gen[n,m,s,t] = min(PV_panels[n]*pv_panel_capacity*curve_pv[m,s,t], PV_inverter[n]*pv_inverter_capacity);

subject to panels_gen_limit_panel {n in N, m in M, s in S, t in T}:
    PV_gen[n,m,s,t] <= PV_panels[n]*pv_panel_capacity*curve_pv[m,s,t];

subject to panels_gen_limit_inverter {n in N, m in M, s in S, t in T}:
    PV_gen[n,m,s,t] <= PV_inverter[n]*pv_inverter_capacity;

subject to GD_reactive_1 {n in N, m in M, s in S, t in T}:
    - PV_gen[n,m,s,t]*tan(acos(theta_pf)) <= qPV_gen[n,m,s,t];

subject to GD_reactive_2 {n in N, m in M,  s in S, t in T}:
    qPV_gen[n,m,s,t] <= PV_gen[n,m,s,t]*tan(acos(theta_pf));

subject to investment_pv_def {n in N}:
    investment_pv[n] = PV_panels[n]*pv_panel_cost + PV_inverter[n]*pv_inverter_cost + PV_allocation[n]*pv_installation_cost;

# tarifarry model
subject to liquid_power {m in M, n in N, s in S, d in D, t in T}:
    E_in[n,m,s,d,t] - E_co[n,m,s,d,t] = (PV_gen[n,m,s,t] - Pd[n]*curve_load[m,d,t])*delt_t;

subject to expected_injection {n in N, m in M}:
    ex_in[n,m] = sum {t in T, s in S, d in D} E_in[n,m,s,d,t]*prob_S[m,s]*prob_D[m,d];

subject to expected_consumption {n in N, m in M}:
    ex_co[n,m] = sum {t in T, s in S, d in D} E_co[n,m,s,d,t]*prob_S[m,s]*prob_D[m,d];

subject to def_cost_con {n in N, m in M}:
    cost_consumer[n,m] = (ET + SUT)*ex_co[n,m];

# Balanço Mensal Esperado
subject to expected_net_balancea {n in N, m in M}:
    deficit[n,m] <= ex_co[n,m];

subject to expected_net_balanceb {n in N, m in M}:
    surplus[n,m] <= ex_in[n,m];

subject to expected_net_balance {n in N, m in M}:
    deficit[n,m] - surplus[n,m] = ex_co[n,m] - ex_in[n,m];

# Restrições de Uso de Créditos 
subject to limit_credit_bank_m1 {n in N, m in M: m == 1}:
    credit_used[n,m] <= 0; # Mês 1 não tem saldo anterior

subject to limit_credit_bank {n in N, m in M: m > 1}:
    credit_used[n,m] <= E_cr[n,m-1];

subject to limit_credit_need {n in N, m in M}:
    credit_used[n,m] <= deficit[n,m];

# Atualização do Banco de Créditos
subject to def_energy_credit_m1 {n in N, m in M: m == 1}:
    E_cr[n,m] = surplus[n,m] - credit_used[n,m];

subject to def_energy_credit {n in N, m in M: m > 1}:
    E_cr[n,m] = E_cr[n,m-1] + surplus[n,m] - credit_used[n,m];

# Cálculo da Energia Paga
subject to def_E_paid {n in N, m in M}:
    E_paid[n,m] = deficit[n,m] - credit_used[n,m];

# Custo do Prosumidor
subject to def_cost_pro {n in N, m in M}:
    cost_prosumer[n,m] = ET * E_paid[n,m] + SUT * (0.9 * ex_in[n,m] + ex_co[n,m]);

subject to savings_def {n in N}:
    an_savings[n] = sum {m in M} (cost_consumer[n,m] - cost_prosumer[n,m]);

# subject to lpv_a {n in N}:
#     investment_pv[n] <= sum{y in 1..5}an_savings[n]/(1 + i)^y;

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
table load_table IN "amplcsv" "input-data/parametros_load_curvas.csv": [M,D,T], curve_load;
read table load_table;

table pv_table IN "amplcsv" "input-data/parametros_pv_curvas.csv": [M,S,T], curve_pv;
read table pv_table;

table pv_prob IN "amplcsv" "input-data/parametros_pv_probs.csv": [M,S], prob_S;
read table pv_prob;

table load_prob IN "amplcsv" "input-data/parametros_load_probs.csv": [M,D], prob_D;
read table load_prob;

# option solver knitro;
# option knitro_options "mip_multistart=1 ms_enable=1 outlev=2 threads=12";
option solver cplex;
option cplex_options 'names=1 iis=1 conflictalg=4';
option cplex_options "mipgap = 0.1 outlev=2 threads=12";

solve FO;
# printf "node deficit surplus ex_co ex_in\n";
# for {m in M}{
#     printf "%d %6.2f %6.2f %6.2f %6.2f\n",m ,deficit[1,m], surplus[1,m], ex_co[1,m], ex_in[1,m];

# }

display cost_prosumer, cost_consumer, E_cr, credit_used, deficit, surplus;

display PV_allocation, PV_panels, PV_inverter;
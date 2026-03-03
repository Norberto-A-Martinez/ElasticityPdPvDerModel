/*
Author: Norberto Abrante Martinez

mixed integer linear programming model for the distribution system planning 
with energy storage systems and capacitor banks allocation

*/

/*DECLARAÇÃO DE CONJUNTOS*/
set N within {0..4000}; # EDS buses/nodes
set L within {0..4000} cross N cross N; # EDS lines/branches
set B; # set of battery technologies
set T ordered; # set of time periods
set S; # set of irradiation scenarios
set Y ordered; # set of years
set D; # set of power demand scenarios

set GD; # set of distributed generation units
set C := 1..3; # set of customer classes

param Wmax := 100; # maximum number of pieces for the current linearization
set W := 1..Wmax ; # set of linearization pieces for the current calculation

/*DECLARAÇÃO DE PARÂMETROS*/
param tag; # usado para controlar qual modelo para o calculo da corrente utilizara
param d_Y {Y} ; # length of each year in days
param delt_t := 1; # duration of period t in hours
param prob_S {S}; # probability of occurrence of scenario
param prob_D {D}; # probability of occurrence of scenario
# current linearization calculation
param V2est {N,S,D,T}; # estimated voltage squared at each node
param mS {W}; # coefficient for the current linearization
param WS; # maximum value for the current linearization
param vb; # nominal voltage magnitude
param vmin ; # minimum voltage
param vmax ; # maximum voltage
param SE; # bus type that identifies the substation 0: load, 1: SE
param Pd_0 {N}; # Initial active power demand
param Qd_0 {N}; # Initial reactive power demand
param S_max := 8e3; # maximum apparent power of the substation
param Sd_0 {N}; # initial active power demand
param pf {N}; # power factor at each node
param de_curve_type {T,C}; # load curve
param de_curve_scenarios {T,D}; # load curve
param de_curve {N,T,D}; # load curve
param Smax {L}; # maximum apparent power in the circuit ij
param R {L}; # resistance in the circuit ij
param X {L}; # reactance in the circuit ij
param Z2 {L}; # impedance in the circuit ij squared
param Imax {L}; # maximum current magnitude limit in the circuit ij
param r3 := 1.7320508076; # square root of 3

# economic parameters
param i; # interest rate
param u; # equipment lifespan in years
param AP ; # AP factor to obtain the equivalent annualized value of a capital under interest i and time
param region {N}; # region identifier for each node
param allowed_busses {N}; # buses allowed for allocation of resources
param bat_inverter_capacity ; # battery inverter capacity
param bat_module_capacity ; # battery module capacity
param Rtef {B}; # round-trip efficiency
param DoD {B}; # depth of discharge
param beta {B}; # battery degradation rate
param life_Cy {B}; # battery life in cycles
param life_yr {B}; # battery life in years
param cosB {B}; # battery storage cost
param cosI {B}; # inverter cost
param man; # maintenance cost percentage
param M := 100; # maximum number of modules in the node
param irrad{T,S}; # distributed generation generation profile
param load {N}; # load at each node
param cp_size; # capacitor size in kVAr
param cp_module_cost; # module cost
param battery_instal_cost; # instalation cost
param capacitor_instal_cost; # instalation cost
param cp_Max; # maximum number of modules in the node
param ET {T}; # energy tariff
param SUT {T}; # system usage tariff
param theta_pf; # PV inverter power factor
param pv_panel_capacity ; # PV panel capacity kwp
param pv_inverter_capacity ; # PV inverter capacity kw
param pv_panel_cost ; # PV panel cost
param pv_inverter_cost ; # PV inverter cost
param pv_installation_cost ; # PV installation cost
param Pd {N}; # active power demand
param Qd {N}; # reactive power demand
param PV_0 {N}; # initial PV generation at each node

# EDS variables
var Pg {N,S,D,T}; # potência ativa fornecida pela subestação no nó i
var Qg {N,S,D,T}; # potência reativa fornecida pela subestação no nó i
var Qcp {N,S} >= 0; # potência reativa do capacitor no nó i
var V2 {N,S,D,T} >= 0; # variável que representa o quadrado de V[i]
var I2 {L,S,D,T} >= 0; # variável que representa o quadrado de I[i,j]
var P {L,S,D,T}; # fluxo de potência ativa no circuito ij
var Q {L,S,D,T}; # fluxo de potência ativa no circuito ij
var capacitor_allocation {N} binary; # variable that indicates if the node n allocation
var battery_allocation {N} binary; # variable that indicates if the node n allocation
var cp_modules {N} integer >= 0; # number of capacitor modules
var bat_inverters {N} integer; # number of battery inverters
var cp_operation {N,S} integer >= 0; # number of active capacitor modules
var bat_opt {N,S,D,T} binary; # battery operation - avoid simultaneous charge and discharge
var EoC_ini {N,S} >= 0; # battery initial state of charge
var bat_disch {N,S,D,T} >= 0; # battery charging power
var bat_charg {N,S,D,T} >= 0; # battery discharging power
var EoC {N,S,D,T} >= 0; # battery state of charge
var PV_allocation {N} binary; # variable that indicates if the node n allocation
var PV_inverter {N} integer >= 0;
var PV_panels {N} integer >= 0;
var PV_gen {N,S,T} >= 0; # PV generation at each node
var PV_cut {N,S,T} >= 0; # PV curtailed power at each node
var qPV_gen {N,S,T}; # potência reativa fornecida pela geração distribuída no nó i
var Pp {L,S,D,T} >= 0; # positive part of the active power flow
var Pn {L,S,D,T} >= 0; # negative part of the active power flow
var Qp {L,S,D,T} >= 0; # positive part of the reactive power flow
var Qn {L,S,D,T} >= 0; # negative part of the reactive power flow
var Dp {L,S,D,T,W} >= 0; # linearization variable for the active power flow
var Dq {L,S,D,T,W} >= 0; # linearization variable for the reactive power flow

# OBJECTIVE FUNCTION
# social vulnerability index pre-calculated in pros_data.dat

var power_loss_I2 = sum {(l,m,n) in L,  t in T, s in S, d in D} prob_D[d] * prob_S[s] * I2[l,m,n,s,d,t] * R[l,m,n];
var cost_p_loss_I2 = sum {(l,m,n) in L,  t in T, s in S, d in D} prob_D[d] * prob_S[s] * I2[l,m,n,s,d,t] * R[l,m,n]*(ET[t] + SUT[t]);


# capacitor costs
var investment_capacitors_modules = sum {n in N} (cp_modules[n] * cp_module_cost * cp_size);
var investment_capacitors_installation = sum {n in N} (capacitor_allocation[n] * capacitor_instal_cost);
var investment_capacitors = investment_capacitors_modules + investment_capacitors_installation;
#batery costs
var investment_battery_modules = sum {b in B, n in N} (bat_inverters[n] * cosB[b] * bat_module_capacity);
var investment_battery_inverter = sum {b in B, n in N} (bat_inverters[n] * cosI[b] * bat_inverter_capacity);
var investment_battery_installation = sum {n in N} (battery_allocation[n] * battery_instal_cost);
var investment_battery = investment_battery_modules + investment_battery_inverter + investment_battery_installation;
# PV system costs
var pv_panel_cost_total = sum {n in N} (PV_panels[n] * pv_panel_cost);
var pv_inverter_cost_total = sum {n in N} (PV_inverter[n] * pv_inverter_cost);
var pv_installation_cost_total = sum {n in N} (PV_allocation[n] * pv_installation_cost);

# total investment
var investment = investment_capacitors + investment_battery;

subject to total_investment: investment <= 500e3;

var investment_pv {N} >=0;

subject to investment_pv_limit {n in N}: investment_pv[n] <= 50e3;

subject to investment_pv_def {n in N}:
    investment_pv[n] = PV_panels[n] * pv_panel_cost + PV_inverter[n] * pv_inverter_cost + PV_allocation[n] * pv_installation_cost;

# # liquid present value of the PV investment.
# subject to lpv_a {n in N}:
#     0 <= investment_pv[n]/sum {t in T, s in S, d in D}(Pd[n]*de_curve[n,t,d] - PV_gen[n,s,t])/(1 + 0.4)^30;

# subject to lpv_b {n in N}:
#     investment_pv[n]/sum {t in T, s in S, d in D}(Pd[n]*de_curve[n,t,d] - PV_gen[n,s,t])/(1 + 0.4)^30 <= 5;

var P_inj {N,S,D,T} >= 0;
var P_cons {N,S,D,T} >= 0;

subject to liquid_power {n in N, s in S, d in D, t in T}:
    P_inj[n,s,d,t] - P_cons[n,s,d,t] = PV_gen[n,s,t] - Pd[n]*de_curve[n,t,d];

var liquid_tariff = sum {n in N, s in S, d in D, t in T} prob_S[s] * prob_D[d] * (ET[t]*(P_inj[n,s,d,t] - P_cons[n,s,d,t]) - SUT[t]*(0.9*P_inj[n,s,d,t] + P_cons[n,s,d,t]));

minimize FO1 : investment + cost_p_loss_I2;
minimize FO2 : liquid_tariff;
minimize FO3 : investment + cost_p_loss_I2 - liquid_tariff;

/* DISTRIBUTION SYSTEM RESTRICTIONS */
subject to substation_1 {d in D, s in S, t in T}:
    Pg[SE,s,d,t] <= S_max * 0.8;

subject to substation_2a {d in D, s in S, t in T}:
    -S_max <= Qg[SE,s,d,t];

subject to substation_2b {d in D, s in S, t in T}:
    Qg[SE,s,d,t] <= S_max;

subject to substation_3a {d in D, s in S, t in T}:
    Qg[SE,s,d,t] <= sqrt(2)*S_max - Pg[SE,s,d,t];

subject to substation_3b {d in D, s in S, t in T}:
    -sqrt(2)*S_max + Pg[SE,s,d,t] <= Qg[SE,s,d,t];

subject to balanco_potencia_ativa {n in N, t in T, s in S, d in D}:
    sum {(l,m,n) in L} P[l,m,n,s,d,t] 
    - sum {(l,n,m) in L} (P[l,n,m,s,d,t] + I2[l,n,m,s,d,t]*R[l,n,m]) 
    - bat_charg[n,s,d,t] + bat_disch[n,s,d,t] 
    + PV_gen[n,s,t] 
    + Pg[n,s,d,t] = Pd[n]*de_curve[n,t,d];

subject to balanco_potencia_reativa {n in N,  t in T,s in S, d in D}:
    sum {(l,m,n) in L} Q[l,m,n,s,d,t] 
    - sum {(l,n,m) in L} (Q[l,n,m,s,d,t] + I2[l,n,m,s,d,t]*X[l,n,m]) 
    + qPV_gen[n,s,t]
    + Qcp[n,s] + Qg[n,s,d,t] = Qd[n]*de_curve[n,t,d];

subject to queda_magnitude_tensao {(l,m,n) in L,  t in T,s in S, d in D}:
    V2[m,s,d,t] - V2[n,s,d,t] = 2*(P[l,m,n,s,d,t] * R[l,m,n] + Q[l,m,n,s,d,t] * X[l,m,n]) + I2[l,m,n,s,d,t] * Z2[l,m,n];

subject to calculo_magnitude_corrente_conico {(l,m,n) in L,  t in T, s in S, d in D: tag == 2}:
    V2[n,s,d,t] * I2[l,m,n,s,d,t] >= P[l,m,n,s,d,t]^2 + Q[l,m,n,s,d,t]^2;

subject to calculo_magnitude_corrente_nl {(l,m,n) in L,  t in T, s in S, d in D: tag == 1}:
    V2[n,s,d,t] * I2[l,m,n,s,d,t] = P[l,m,n,s,d,t]^2 + Q[l,m,n,s,d,t]^2;

subject to CURRENT_FLOW_SQUARE{(l,m,n) in L,  t in T, s in S, d in D: tag == 0} :
    V2est[n,s,d,t] * I2[l,m,n,s,d,t] = sum{w in W} mS[w] * (Dp[l,m,n,s,d,t,w] + Dq[l,m,n,s,d,t,w]);
    
subject to WELTA_CURRENT_FLOW_re{(l,m,n) in L,  t in T, s in S, d in D: tag == 0} :
    Pp[l,m,n,s,d,t] + Pn[l,m,n,s,d,t] = sum{w in W} Dp[l,m,n,s,d,t,w];

subject to WELTA_CURRENT_FLOW_im{(l,m,n) in L,  t in T, s in S, d in D: tag == 0} :
    Qp[l,m,n,s,d,t] + Qn[l,m,n,s,d,t] = sum{w in W} Dq[l,m,n,s,d,t,w];

subject to WELTA_CURRENT_FLOW_re_soma{(l,m,n) in L,  t in T, s in S, d in D: tag == 0} :
    P[l,m,n,s,d,t] = Pp[l,m,n,s,d,t] - Pn[l,m,n,s,d,t];

subject to WELTA_CURRENT_FLOW_im_soma{(l,m,n) in L,  t in T, s in S, d in D: tag == 0} :
    Q[l,m,n,s,d,t] = Qp[l,m,n,s,d,t] - Qn[l,m,n,s,d,t];

subject to MAXIMUN_WELTA_CURRENT_FLOW_re {(l,m,n) in L,  w in W, t in T, s in S, d in D: tag == 0}:
    Dp[l,m,n,s,d,t,w] <= WS;

subject to MAXIMUN_WELTA_CURRENT_FLOW_im {(l,m,n) in L,  w in W, t in T, s in S, d in D: tag == 0}:
    Dq[l,m,n,s,d,t,w] <= WS;

# current limit
subject to corrente {(l,m,n) in L,  t in T,s in S, d in D}:
    -Imax[l,m,n]^2 <= I2[l,m,n,s,d,t] <= Imax[l,m,n]^2;

# upper voltage limit
subject to tens_max {n in N,  t in T,s in S, d in D}:
    V2[n,s,d,t] <= (vmax*vb)^2;

# lower voltage limit
subject to tens_min {n in N,  t in T,s in S, d in D}:
    V2[n,s,d,t] >= (vmin*vb)^2;

# capacitors
subject to CapacitorAllocation {n in N}:
    capacitor_allocation[n] <= allowed_busses[n];

subject to capacitor_modules_decision {n in N}:
    cp_modules[n] <= M*capacitor_allocation[n];

subject to capacitor_operation_decision {n in N,  s in S}:
    cp_operation[n,s] <= cp_modules[n];

subject to reactive_power_delivered {n in N,  s in S}:
    Qcp[n,s] = cp_size*cp_operation[n,s];

# batteries
subject to BatteryAllocation {n in N}:
    battery_allocation[n] <= allowed_busses[n];

subject to eq_15a {n in N}:
    bat_inverters[n] <= M*battery_allocation[n];

subject to BatteryOpt_charg {n in N,  s in S, t in T, d in D}:
    bat_charg[n,s,d,t] <= M*bat_inverter_capacity*(1 - bat_opt[n,s,d,t]);

subject to BatteryOpt_disch {n in N,  s in S, t in T, d in D}:
    bat_disch[n,s,d,t] <= M*bat_inverter_capacity*bat_opt[n,s,d,t];

subject to eq_23 {n in N,  s in S, t in T, d in D}:
    bat_charg[n,s,d,t] <= bat_inverter_capacity*bat_inverters[n];

subject to eq_24 {n in N,  s in S, t in T, d in D}:
    bat_disch[n,s,d,t] <= bat_inverter_capacity*bat_inverters[n];

subject to eq_25 {n in N,  s in S, t in T, d in D, b in B: t == 1}:
    EoC[n,s,d,t] - EoC_ini[n,s] = delt_t*(Rtef[b]*bat_charg[n,s,d,t] - (1/Rtef[b])*bat_disch[n,s,d,t]) - (beta[b]/card(T))*EoC[n,s,d,t];

subject to estado_de_carga {n in N,  s in S, t in T, d in D, b in B: t > 1}:
    EoC[n,s,d,t] - EoC[n,s,d,t-1] = delt_t*(Rtef[b]*bat_charg[n,s,d,t] - (1/Rtef[b])*bat_disch[n,s,d,t]) - (beta[b]/card(T))*EoC[n,s,d,t];

subject to initial_sttate {n in N, d in D, s in S}:
    EoC_ini[n,s] = EoC[n,s,d,last(T)];

subject to max_energia {n in N,  s in S, t in T, d in D, b in B}:
    EoC[n,s,d,t] <= bat_module_capacity * DoD[b] * bat_inverters[n];

subject to min_energia {n in N,  s in S, t in T, d in D, b in B}:
    EoC[n,s,d,t] >= bat_module_capacity * (1-DoD[b]) * bat_inverters[n];

# PV SYSTEM allocation
subject to PV_panels_allocation {n in N}:
    PV_panels[n] <= 80*PV_allocation[n];

subject to pv_inverter_allocation {n in N}:
    PV_inverter[n] <= 20*PV_allocation[n];

subject to panels_generation {n in N, s in S, t in T}:
    PV_gen[n,s,t] = min (PV_panels[n] * pv_panel_capacity * irrad[t,s], PV_inverter[n] * pv_inverter_capacity);

# subject to PV_inverter_capacity {n in N,  s in S, t in T}:
#     PV_cut[n,s,t] <= PV_inverter[n] * pv_inverter_capacity;

subject to GD_reactive_1 {n in N, s in S, t in T}:
    - PV_gen[n,s,t] * tan(acos(theta_pf)) <= qPV_gen[n,s,t];

subject to GD_reactive_2 {n in N,  s in S, t in T}:
    qPV_gen[n,s,t] <= PV_gen[n,s,t] * tan(acos(theta_pf));


# \\\\\\\\\\\\\\\\\\\\\\\\
#  PERFORMANCE INDICATORS
# \\\\\\\\\\\\\\\\\\\\\\\\

# # Expected energy consumed
# var Ener_PD = sum {n in N, t in T, s in S, c in C} (delt_t * prob_S[s] * Pd[n]*de_curve[t,n]);

# # Pondered consumed energy
# var Ener_PD_ponderada = sum {n in N, t in T, s in S, d in D, c in C} ((delt_t * prob_S[s] * Pd[n] * de_curve[t,n])/(Pd_0[n]+ 1e-6));

# Expected energy supplied by the substation
var Ener_SE = sum {n in N, t in T, s in S, d in D} (delt_t * prob_S[s] * prob_D[d] * Pg[n,s,d,t]);

# Power loss in the system
var power_loss = sum {(l,m,n) in L, t in T, s in S, d in D} (prob_D[d] * prob_S[s] * R[l,m,n]*(P[l,m,n,s,d,t]^2 + Q[l,m,n,s,d,t]^2)/V2[n,s,d,t]);

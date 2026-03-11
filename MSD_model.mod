/*
Author: Norberto Abrante Martinez

mixed integer linear programming model for the distribution system planning 
with energy storage systems and capacitor banks allocation

*/

/*DECLARAÇÃO DE CONJUNTOS*/
set N within {0..4000}; # EDS buses/nodes
set L within {0..4000} cross N cross N; # EDS lines/branches
set S := 1..1; # set of irradiation scenarios
set D := 1..1; # set of power demand scenarios
set M := 1..12 ordered; # set of years
set B; # set of battery technologies
set T ordered; # set of time periods

param Wmax := 50; # maximum number of pieces for the current linearization
set W := 1..Wmax ; # set of linearization pieces for the current calculation

/*DECLARAÇÃO DE PARÂMETROS*/
param tag; # usado para controlar qual modelo para o calculo da corrente utilizara
param delt_t := 1; # duration of period t in hours
param prob_S {M,S}; # probability of occurrence of scenario
param prob_D {M,D}; # probability of occurrence of scenario
# current linearization calculation
param V2est {N,M,S,D,T}; # estimated voltage squared at each node
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
param curve_load {M,D,T}; # load curve
param curve_pv {M,S,T}; # load curve
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
param Pd {N}; # active power demand
param Qd {N}; # reactive power demand
param PV_allocation {N}; # csv imported param
param PV_inverter {N}; # csv imported param
param PV_panels {N}; # csv imported param
param pv_installation_cost;
param pv_inverter_cost;
param pv_panel_cost;

# EDS variables
var Pg {N,M,S,D,T}; # potência ativa fornecida pela subestação no nó i
var Qg {N,M,S,D,T}; # potência reativa fornecida pela subestação no nó i
var Qcp {N,M,S,D} >= 0; # potência reativa do capacitor no nó i
var V2 {N,M,S,D,T} >= 0; # variável que representa o quadrado de V[i]
var I2 {L,M,S,D,T}; # variável que representa o quadrado de I[i,j]
var P {L,M,S,D,T}; # fluxo de potência ativa no circuito ij
var Q {L,M,S,D,T}; # fluxo de potência ativa no circuito ij
var capacitor_allocation {N} binary; # variable that indicates if the node n allocation
var battery_allocation {N} binary; # variable that indicates if the node n allocation
var cp_modules {N} integer >= 0; # number of capacitor modules
var bat_inverters {N} integer; # number of battery inverters
var cp_operation {N,M,S,D} integer >= 0; # number of active capacitor modules
var bat_opt {N,M,S,D,T} binary; # battery operation - avoid simultaneous charge and discharge
var EoC_ini {N,M,S,D} >= 0; # battery initial state of charge
var bat_disch {N,M,S,D,T} >= 0; # battery charging power
var bat_charg {N,M,S,D,T} >= 0; # battery discharging power
var EoC {N,M,S,D,T} >= 0; # battery state of charge
var PV_gen {N,M,S,T} >= 0; # PV generation at each node
var PV_cut {N,M,S,T} >= 0; # PV curtailed power at each node
var qPV_gen {N,M,S,T}; # potência reativa fornecida pela geração distribuída no nó i
var Pp {L,M,S,D,T} >= 0; # positive part of the active power flow
var Pn {L,M,S,D,T} >= 0; # negative part of the active power flow
var Qp {L,M,S,D,T} >= 0; # positive part of the reactive power flow
var Qn {L,M,S,D,T} >= 0; # negative part of the reactive power flow
var Dp {L,M,S,D,T,W} >= 0; # linearization variable for the active power flow
var Dq {L,M,S,D,T,W} >= 0; # linearization variable for the reactive power flow

# OBJECTIVE FUNCTION
# Power loss in the system
var power_loss = 30.4 * sum {(l,m,n) in L, mes in M, t in T, s in S, d in D} (prob_D[mes,d] * prob_S[mes,s] * R[l,m,n]*(P[l,m,n,mes,s,d,t]^2 + Q[l,m,n,mes,s,d,t]^2)/V2[n,mes,s,d,t]);
var power_loss_I2 = 30.4 * sum {(l,m,n) in L, mes in M, t in T, s in S, d in D} I2[l,m,n,mes,s,d,t] * R[l,m,n] * prob_D[mes,d] * prob_S[mes,s];
var cost_loss_I2 = 30.4 * sum {(l,m,n) in L, mes in M, t in T, s in S, d in D} I2[l,m,n,mes,s,d,t] * R[l,m,n] * prob_D[mes,d] * prob_S[mes,s] * (ET[t] + SUT[t]);

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
param pv_panel_cost_total = sum {n in N} (PV_panels[n] * pv_panel_cost);
param pv_inverter_cost_total = sum {n in N} (PV_inverter[n] * pv_inverter_cost);
param pv_installation_cost_total = sum {n in N} (PV_allocation[n] * pv_installation_cost);

# total investment
var investment = investment_capacitors + investment_battery;
# subject to total_investment: investment <= 500e3;

# Expected energy supplied by the substation
var Ener_SE = 30.4 * sum {n in N, mes in M, t in T, s in S, d in D} (Pg[n,mes,s,d,t] * prob_S[mes,s] * prob_D[mes,d] * delt_t);


minimize FO1 : investment * AP + cost_loss_I2 + Ener_SE;

/* DISTRIBUTION SYSTEM RESTRICTIONS */
subject to substation_1 {m in M, d in D, s in S, t in T}:
    Pg[SE,m,s,d,t] <= S_max * 0.8;

subject to substation_2a {m in M, d in D, s in S, t in T}:
    -S_max <= Qg[SE,m,s,d,t];

subject to substation_2b {m in M, d in D, s in S, t in T}:
    Qg[SE,m,s,d,t] <= S_max;

subject to substation_3a {m in M, d in D, s in S, t in T}:
    Qg[SE,m,s,d,t] <= sqrt(2)*S_max - Pg[SE,m,s,d,t];

subject to substation_3b {m in M, d in D, s in S, t in T}:
    -sqrt(2)*S_max + Pg[SE,m,s,d,t] <= Qg[SE,m,s,d,t];

subject to balanco_potencia_ativa {n in N, mes in M, t in T, s in S, d in D}:
    sum {(l,m,n) in L} P[l,m,n,mes,s,d,t] 
    - sum {(l,n,m) in L} (P[l,n,m,mes,s,d,t] + I2[l,n,m,mes,s,d,t]*R[l,n,m]) 
    - bat_charg[n,mes,s,d,t] + bat_disch[n,mes,s,d,t] 
    + PV_gen[n,mes,s,t] 
    + Pg[n,mes,s,d,t] = Pd[n]*curve_load[mes,d,t];

subject to balanco_potencia_reativa {n in N, mes in M, t in T,s in S, d in D}:
    sum {(l,m,n) in L} Q[l,m,n,mes,s,d,t] 
    - sum {(l,n,m) in L} (Q[l,n,m,mes,s,d,t] + I2[l,n,m,mes,s,d,t]*X[l,n,m]) 
    + qPV_gen[n,mes,s,t]
    + Qcp[n,mes,s,d] + Qg[n,mes,s,d,t] = Qd[n]*curve_load[mes,d,t];

subject to queda_magnitude_tensao {(l,m,n) in L, mes in M, t in T,s in S, d in D}:
    V2[m,mes,s,d,t] - V2[n,mes,s,d,t] = 2*(P[l,m,n,mes,s,d,t] * R[l,m,n] + Q[l,m,n,mes,s,d,t] * X[l,m,n]) + I2[l,m,n,mes,s,d,t] * Z2[l,m,n];

subject to calculo_magnitude_corrente_conico {(l,m,n) in L, mes in M, t in T, s in S, d in D: tag == 2}:
    V2[n,mes,s,d,t] * I2[l,m,n,mes,s,d,t] >= P[l,m,n,mes,s,d,t]^2 + Q[l,m,n,mes,s,d,t]^2;

subject to calculo_magnitude_corrente_nl {(l,m,n) in L, mes in M, t in T, s in S, d in D: tag == 1}:
    V2[n,mes,s,d,t] * I2[l,m,n,mes,s,d,t] = P[l,m,n,mes,s,d,t]^2 + Q[l,m,n,mes,s,d,t]^2;

subject to CURRENT_FLOW_SQUARE{(l,m,n) in L, mes in M, t in T, s in S, d in D: tag == 0} :
    V2est[n,mes,s,d,t] * I2[l,m,n,mes,s,d,t] = sum{w in W} mS[w] * (Dp[l,m,n,mes,s,d,t,w] + Dq[l,m,n,mes,s,d,t,w]);
    
subject to WELTA_CURRENT_FLOW_re{(l,m,n) in L, mes in M, t in T, s in S, d in D: tag == 0} :
    Pp[l,m,n,mes,s,d,t] + Pn[l,m,n,mes,s,d,t] = sum{w in W} Dp[l,m,n,mes,s,d,t,w];

subject to WELTA_CURRENT_FLOW_im{(l,m,n) in L, mes in M, t in T, s in S, d in D: tag == 0} :
    Qp[l,m,n,mes,s,d,t] + Qn[l,m,n,mes,s,d,t] = sum{w in W} Dq[l,m,n,mes,s,d,t,w];

subject to WELTA_CURRENT_FLOW_re_soma{(l,m,n) in L, mes in M, t in T, s in S, d in D: tag == 0} :
    P[l,m,n,mes,s,d,t] = Pp[l,m,n,mes,s,d,t] - Pn[l,m,n,mes,s,d,t];

subject to WELTA_CURRENT_FLOW_im_soma{(l,m,n) in L, mes in M, t in T, s in S, d in D: tag == 0} :
    Q[l,m,n,mes,s,d,t] = Qp[l,m,n,mes,s,d,t] - Qn[l,m,n,mes,s,d,t];

subject to MAXIMUN_WELTA_CURRENT_FLOW_re {(l,m,n) in L, mes in M, w in W, t in T, s in S, d in D: tag == 0}:
    Dp[l,m,n,mes,s,d,t,w] <= WS;

subject to MAXIMUN_WELTA_CURRENT_FLOW_im {(l,m,n) in L, mes in M, w in W, t in T, s in S, d in D: tag == 0}:
    Dq[l,m,n,mes,s,d,t,w] <= WS;

# current limit
subject to corrente {(l,m,n) in L, mes in M, t in T,s in S, d in D}:
    -Imax[l,m,n]^2 <= I2[l,m,n,mes,s,d,t] <= Imax[l,m,n]^2;

# upper voltage limit
subject to tens_max {n in N, mes in M, t in T,s in S, d in D}:
    V2[n,mes,s,d,t] <= (vmax*vb)^2;

# lower voltage limit
subject to tens_min {n in N, mes in M, t in T,s in S, d in D}:
    V2[n,mes,s,d,t] >= (vmin*vb)^2;

# capacitors
subject to CapacitorAllocation {n in N}:
    capacitor_allocation[n] <= allowed_busses[n];

subject to capacitor_modules_decision {n in N}:
    cp_modules[n] <= 100*capacitor_allocation[n];

subject to capacitor_operation_decision {n in N, mes in M, s in S, d in D}:
    cp_operation[n,mes,s,d] <= cp_modules[n];

subject to reactive_power_delivered {n in N, mes in M, s in S, d in D}:
    Qcp[n,mes,s,d] = cp_size*cp_operation[n,mes,s,d];

# batteries
subject to BatteryAllocation {n in N}:
    battery_allocation[n] <= allowed_busses[n];

subject to eq_15a {n in N}:
    bat_inverters[n] <= 100*battery_allocation[n];

subject to BatteryOpt_charg {n in N, mes in M, s in S, t in T, d in D}:
    bat_charg[n,mes,s,d,t] <= 100*bat_inverter_capacity*(1 - bat_opt[n,mes,s,d,t]);

subject to BatteryOpt_disch {n in N, mes in M, s in S, t in T, d in D}:
    bat_disch[n,mes,s,d,t] <= 100*bat_inverter_capacity*bat_opt[n,mes,s,d,t];

subject to eq_23 {n in N, mes in M, s in S, t in T, d in D}:
    bat_charg[n,mes,s,d,t] <= bat_inverter_capacity*bat_inverters[n];

subject to eq_24 {n in N, mes in M, s in S, t in T, d in D}:
    bat_disch[n,mes,s,d,t] <= bat_inverter_capacity*bat_inverters[n];

subject to eq_25 {n in N, mes in M, s in S, t in T, d in D, b in B: t == 1}:
    EoC[n,mes,s,d,t] - EoC_ini[n,mes,s,d] = delt_t*(Rtef[b]*bat_charg[n,mes,s,d,t] - (1/Rtef[b])*bat_disch[n,mes,s,d,t]) - (beta[b]/card(T))*EoC[n,mes,s,d,t];

subject to estado_de_carga {n in N, mes in M, s in S, t in T, d in D, b in B: t > 1}:
    EoC[n,mes,s,d,t] - EoC[n,mes,s,d,t-1] = delt_t*(Rtef[b]*bat_charg[n,mes,s,d,t] - (1/Rtef[b])*bat_disch[n,mes,s,d,t]) - (beta[b]/card(T))*EoC[n,mes,s,d,t];

subject to initial_sttate {n in N, mes in M, d in D, s in S}:
    EoC_ini[n,mes,s,d] = EoC[n,mes,s,d,last(T)];

subject to max_energia {n in N, mes in M, s in S, t in T, d in D, b in B}:
    EoC[n,mes,s,d,t] <= bat_module_capacity * DoD[b] * bat_inverters[n];

subject to min_energia {n in N, mes in M, s in S, t in T, d in D, b in B}:
    EoC[n,mes,s,d,t] >= bat_module_capacity * (1-DoD[b]) * bat_inverters[n];

subject to PV_inveter_limit {n in N, mes in M, s in S, t in T}:
    PV_gen[n,mes,s,t] <= PV_inverter[n] * pv_inverter_capacity;

subject to PV_panels_limit {n in N, mes in M, s in S, t in T}:
    PV_gen[n,mes,s,t] <= PV_panels[n] * pv_panel_capacity * curve_pv[mes,s,t];

subject to GD_reactive_1 {n in N, mes in M, s in S, t in T}:
    - PV_gen[n,mes,s,t] * tan(acos(theta_pf)) <= qPV_gen[n,mes,s,t];

subject to GD_reactive_2 {n in N, mes in M, s in S, t in T}:
    qPV_gen[n,mes,s,t] <= PV_gen[n,mes,s,t] * tan(acos(theta_pf));

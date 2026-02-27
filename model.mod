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
# current linearization calculation
param V2est {N,Y,S,T}; # estimated voltage squared at each node
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
param de_curve {T,C}; # load curve
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
param cus_ene {T}; # energy cost
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
var Pg {N,Y,S,T}; # potência ativa fornecida pela subestação no nó i
var Qg {N,Y,S,T}; # potência reativa fornecida pela subestação no nó i
var Qcp {N,Y,S} >= 0; # potência reativa do capacitor no nó i
var V2 {N,Y,S,T} >= 0; # variável que representa o quadrado de V[i]
var I2 {L,Y,S,T} >= 0; # variável que representa o quadrado de I[i,j]
var P {L,Y,S,T}; # fluxo de potência ativa no circuito ij
var Q {L,Y,S,T}; # fluxo de potência ativa no circuito ij
var capacitor_allocation {N} binary; # variable that indicates if the node n allocation
var battery_allocation {N} binary; # variable that indicates if the node n allocation
var cp_modules {N} integer >= 0; # number of capacitor modules
var bat_inverters {N} integer; # number of battery inverters
var cp_operation {N,Y,S} integer >= 0; # number of active capacitor modules
var bat_opt {N,Y,S,T} binary; # battery operation - avoid simultaneous charge and discharge
var EoC_ini {N,Y,S} >= 0; # battery initial state of charge
var bat_disch {N,Y,S,T} >= 0; # battery charging power
var bat_charg {N,Y,S,T} >= 0; # battery discharging power
var EoC {N,Y,S,T} >= 0; # battery state of charge
var PV_allocation {N} binary; # variable that indicates if the node n allocation
var PV_inverter {N} integer >= 0;
var PV_panels {N} integer >= 0;
var PV_gen {N,Y,S,T} >= 0; # PV generation at each node
var PV_cut {N,Y,S,T} >= 0; # PV curtailed power at each node
var qPV_gen {N,Y,S,T}; # potência reativa fornecida pela geração distribuída no nó i
var Pp {L,Y,S,T} >= 0; # positive part of the active power flow
var Pn {L,Y,S,T} >= 0; # negative part of the active power flow
var Qp {L,Y,S,T} >= 0; # positive part of the reactive power flow
var Qn {L,Y,S,T} >= 0; # negative part of the reactive power flow
var Dp {L,Y,S,T,W} >= 0; # linearization variable for the active power flow
var Dq {L,Y,S,T,W} >= 0; # linearization variable for the reactive power flow
var Sd_S_influence {N,Y} >= 0; # power demand considering spatial influence
var Sd_Temp_grow {N,Y} >= 0; # variable that indicates the teporal grow of Pd
var Sd_supr {N,Y} >= 0; # variable that indicates the power demand suppressed

# OBJECTIVE FUNCTION
# social vulnerability index pre-calculated in pros_data.dat

var power_loss_I2 = sum {s in S} (prob_S[s]*sum {(l,m,n) in L, y in Y, t in T} I2[l,m,n,y,s,t]*R[l,m,n]);

minimize FO1 : power_loss_I2;

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
var investment_pv = pv_panel_cost_total + pv_inverter_cost_total + pv_installation_cost_total;
# total investment
var investment = investment_capacitors + investment_battery + investment_pv;

/* DISTRIBUTION SYSTEM RESTRICTIONS */
subject to investment_limit: investment <= 200e3;

subject to substation_1 {y in Y, s in S, t in T}:
    Pg[SE,y,s,t] <= S_max * 0.8;

subject to substation_2a {y in Y, s in S, t in T}:
    -S_max <= Qg[SE,y,s,t];

subject to substation_2b {y in Y, s in S, t in T}:
    Qg[SE,y,s,t] <= S_max;

subject to substation_3a {y in Y, s in S, t in T}:
    Qg[SE,y,s,t] <= sqrt(2)*S_max - Pg[SE,y,s,t];

subject to substation_3b {y in Y, s in S, t in T}:
    -sqrt(2)*S_max + Pg[SE,y,s,t] <= Qg[SE,y,s,t];


subject to balanco_potencia_ativa {n in N, y in Y, t in T,s in S}:
    sum {(l,m,n) in L} P[l,m,n,y,s,t] 
    - sum {(l,n,m) in L} (P[l,n,m,y,s,t] + I2[l,n,m,y,s,t]*R[l,n,m]) 
    - bat_charg[n,y,s,t] + bat_disch[n,y,s,t] 
    + PV_gen[n,y,s,t] 
    + Pg[n,y,s,t] = Pd[n]*de_curve[t,region[n]];

subject to balanco_potencia_reativa {n in N, y in Y, t in T,s in S}:
    sum {(l,m,n) in L} Q[l,m,n,y,s,t] 
    - sum {(l,n,m) in L} (Q[l,n,m,y,s,t] + I2[l,n,m,y,s,t]*X[l,n,m]) 
    + qPV_gen[n,y,s,t]
    + Qcp[n,y,s] + Qg[n,y,s,t] = Qd[n]*de_curve[t,region[n]];

subject to queda_magnitude_tensao {(l,m,n) in L, y in Y, t in T,s in S}:
    V2[m,y,s,t] - V2[n,y,s,t] = 2*(P[l,m,n,y,s,t] * R[l,m,n] + Q[l,m,n,y,s,t] * X[l,m,n]) + I2[l,m,n,y,s,t] * Z2[l,m,n];

subject to calculo_magnitude_corrente_conico {(l,m,n) in L, y in Y, t in T, s in S: tag == 2}:
    V2[n,y,s,t] * I2[l,m,n,y,s,t] >= P[l,m,n,y,s,t]^2 + Q[l,m,n,y,s,t]^2;

subject to calculo_magnitude_corrente_nl {(l,m,n) in L, y in Y, t in T, s in S: tag == 1}:
    V2[n,y,s,t] * I2[l,m,n,y,s,t] = P[l,m,n,y,s,t]^2 + Q[l,m,n,y,s,t]^2;

subject to CURRENT_FLOW_SQUARE{(l,m,n) in L, y in Y, t in T, s in S: tag == 0} :
    V2est[n,y,s,t] * I2[l,m,n,y,s,t] = sum{w in W} mS[w] * (Dp[l,m,n,y,s,t,w] + Dq[l,m,n,y,s,t,w]);
    
subject to WELTA_CURRENT_FLOW_re{(l,m,n) in L, y in Y, t in T, s in S: tag == 0} :
    Pp[l,m,n,y,s,t] + Pn[l,m,n,y,s,t] = sum{w in W} Dp[l,m,n,y,s,t,w];

subject to WELTA_CURRENT_FLOW_im{(l,m,n) in L, y in Y, t in T, s in S: tag == 0} :
    Qp[l,m,n,y,s,t] + Qn[l,m,n,y,s,t] = sum{w in W} Dq[l,m,n,y,s,t,w];

subject to WELTA_CURRENT_FLOW_re_soma{(l,m,n) in L, y in Y, t in T, s in S: tag == 0} :
    P[l,m,n,y,s,t] = Pp[l,m,n,y,s,t] - Pn[l,m,n,y,s,t];

subject to WELTA_CURRENT_FLOW_im_soma{(l,m,n) in L, y in Y, t in T, s in S: tag == 0} :
    Q[l,m,n,y,s,t] = Qp[l,m,n,y,s,t] - Qn[l,m,n,y,s,t];

subject to MAXIMUN_WELTA_CURRENT_FLOW_re {(l,m,n) in L, y in Y, w in W, t in T, s in S: tag == 0}:
    Dp[l,m,n,y,s,t,w] <= WS;

subject to MAXIMUN_WELTA_CURRENT_FLOW_im {(l,m,n) in L, y in Y, w in W, t in T, s in S: tag == 0}:
    Dq[l,m,n,y,s,t,w] <= WS;

# current limit
subject to corrente {(l,m,n) in L, y in Y, t in T,s in S}:
    -Imax[l,m,n]^2 <= I2[l,m,n,y,s,t] <= Imax[l,m,n]^2;

# upper voltage limit
subject to tens_max {n in N, y in Y, t in T,s in S}:
    V2[n,y,s,t] <= (vmax*vb)^2;

# lower voltage limit
subject to tens_min {n in N, y in Y, t in T,s in S}:
    V2[n,y,s,t] >= (vmin*vb)^2;

# capacitors
subject to CapacitorAllocation {n in N}:
    capacitor_allocation[n] <= allowed_busses[n];

subject to capacitor_modules_decision {n in N}:
    cp_modules[n] <= M*capacitor_allocation[n];

subject to capacitor_operation_decision {n in N, y in Y, s in S}:
    cp_operation[n,y,s] <= cp_modules[n];

subject to reactive_power_delivered {n in N, y in Y, s in S}:
    Qcp[n,y,s] = cp_size*cp_operation[n,y,s];

# batteries
subject to BatteryAllocation {n in N}:
    battery_allocation[n] <= allowed_busses[n];

subject to eq_15a {n in N}:
    bat_inverters[n] <= M*battery_allocation[n];

subject to BatteryOpt_charg {n in N, y in Y, s in S, t in T}:
    bat_charg[n,y,s,t] <= M*bat_inverter_capacity*(1 - bat_opt[n,y,s,t]);

subject to BatteryOpt_disch {n in N, y in Y, s in S, t in T}:
    bat_disch[n,y,s,t] <= M*bat_inverter_capacity*bat_opt[n,y,s,t];

subject to eq_23 {n in N, y in Y, s in S, t in T}:
    bat_charg[n,y,s,t] <= bat_inverter_capacity*bat_inverters[n];

subject to eq_24 {n in N, y in Y, s in S, t in T}:
    bat_disch[n,y,s,t] <= bat_inverter_capacity*bat_inverters[n];

subject to eq_25 {n in N, y in Y, s in S, t in T, b in B: t == 1}:
    EoC[n,y,s,t] - EoC_ini[n,y,s] = delt_t*(Rtef[b]*bat_charg[n,y,s,t] - (1/Rtef[b])*bat_disch[n,y,s,t]) - (beta[b]/card(T))*EoC[n,y,s,t];

subject to estado_de_carga {n in N, y in Y, s in S, t in T, b in B: t > 1}:
    EoC[n,y,s,t] - EoC[n,y,s,t-1] = delt_t*(Rtef[b]*bat_charg[n,y,s,t] - (1/Rtef[b])*bat_disch[n,y,s,t]) - (beta[b]/card(T))*EoC[n,y,s,t];

subject to initial_sttate {n in N, y in Y, s in S}:
    EoC_ini[n,y,s] = EoC[n,y,s,last(T)];

subject to max_energia {n in N, y in Y, s in S, t in T, b in B}:
    EoC[n,y,s,t] <= bat_module_capacity * DoD[b] * bat_inverters[n];

subject to min_energia {n in N, y in Y, s in S, t in T, b in B}:
    EoC[n,y,s,t] >= bat_module_capacity * (1-DoD[b]) * bat_inverters[n];

# PV SYSTEM allocation
subject to PV_panels_allocation {n in N}:
    PV_panels[n] <= 80*PV_allocation[n];

subject to pv_inverter_allocation {n in N}:
    PV_inverter[n] <= 10*PV_allocation[n];

subject to panels_generation {n in N, y in Y, s in S, t in T}:
    PV_gen[n,y,s,t] = min (PV_panels[n] * pv_panel_capacity * irrad[t,s], PV_inverter[n] * pv_inverter_capacity);

# subject to PV_inverter_capacity {n in N, y in Y, s in S, t in T}:
#     PV_cut[n,y,s,t] <= PV_inverter[n] * pv_inverter_capacity;

subject to GD_reactive_1 {n in N, y in Y, s in S, t in T}:
    - PV_gen[n,y,s,t] * tan(acos(theta_pf)) <= qPV_gen[n,y,s,t];

subject to GD_reactive_2 {n in N, y in Y, s in S, t in T}:
    qPV_gen[n,y,s,t] <= PV_gen[n,y,s,t] * tan(acos(theta_pf));


# \\\\\\\\\\\\\\\\\\\\\\\\
#  PERFORMANCE INDICATORS
# \\\\\\\\\\\\\\\\\\\\\\\\

# Expected energy consumed
var Ener_PD = sum {n in N, y in Y, t in T, s in S, c in C} (d_Y[y]* delt_t * prob_S[s] * Pd[n]*de_curve[t,region[n]]);

# Pondered consumed energy
var Ener_PD_ponderada = sum {n in N, y in Y, t in T, s in S, c in C} ((d_Y[y]* delt_t * prob_S[s] * Pd[n] * de_curve[t,region[n]])/(Pd_0[n]+ 1e-6));

# Expected energy supplied by the substation
var Ener_SE = sum {n in N, y in Y, t in T, s in S} (d_Y[y] * delt_t * prob_S[s] * Pg[n,y,s,t]);

# Power loss in the system
var power_loss = sum {s in S} (prob_S[s]*sum {(l,m,n) in L, y in Y, t in T} (R[l,m,n]*((P[l,m,n,y,s,t]^2 + Q[l,m,n,y,s,t]^2)/V2[n,y,s,t])));

# suppressed power
var Supre_power = sum {n in N, y in Y, t in T} d_Y[y] * delt_t * Sd_supr[n,y];

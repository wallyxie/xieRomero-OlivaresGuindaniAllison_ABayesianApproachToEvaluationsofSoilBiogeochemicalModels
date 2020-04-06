functions {

/*
*AWB model as depicted in Li et al., 2014
*/

/*
*Arrhenius relationship for calculating decay constants
*/

real decay(real var, real Ea, real temp, real[] x_r) {
return var * exp(-Ea / x_r[2] * (1 / temp - 1 / x_r[1]));
}

real tempMod(real var, real m_t, real temp, real[] x_r) {
return var - m_t * (temp - x_r[1]);
}

/*
*AWB ODE system described in Li et al., 2014
*C[1] is mass of SOC pool
*C[2] is mass of DOC pool
*C[3] is mass of MBC pool
*C[4] is mass of the ENZC pool
*/

real[] AWB_ODE(real t, real[] C, real[] theta, real[] x_r, int[] x_i) {

real dC_dt[4];

//parameters
real V_f; // SOC reference V_max 
real V_U_f; // DOC reference V_max
real K_f; // SOC reference K_m
real K_U_f; // DOC reference K_m
real Ea_V; // SOC V_max activation energy
real Ea_VU; // DOC V_max activation energy
real Ea_K; // SOC K_m activation energy
real Ea_KU; // DOC K_m activation energy 
real r_M; // MBC turnover rate
real r_E; // Enzyme production rate
real r_L; // Enzyme loss rate
real a_MS; // Fraction of dead MBC transferred to SOC
real E_C_f; // Reference carbon use efficiency (CUE)
real m_t; // CUE slope

/*
*Dependent variables
*/

real F_S; // SOC decomposition
real F_D; // DOC uptake
real F_M; // MBC death
real F_E; // Enzyme production
real F_L; // Enzyme loss

/*
*Match parameters to thetas
*/

V_f = theta[1]; // SOC reference V_max
V_U_f = theta[2]; // DOC reference V_max
K_f = theta[3]; // SOC reference K_m
K_U_f = theta[4]; // DOC reference K_m
Ea_V = theta[5]; // SOC V_max activation energy
Ea_VU = theta[6]; // DOC V_max activation energy
Ea_K = theta[7]; // SOC K_m activation energy
Ea_KU = theta[8]; // DOC K_m activation energy 
r_M = theta[9]; // MBC turnover rate
a_MS = theta[10]; // Fraction of dead MBC transferred to SOC
E_C_f = theta[11]; // Carbon use efficiency

r_L = x_r[9];
r_E = x_r[10];

F_S = V_f * C[4] * C[1] / (K_f + C[1]);
F_D = V_U_f * C[3] * C[2] / (K_U_f + C[2]);
F_M = r_M * C[3];
F_E = r_E * C[3];
F_L = r_L * C[4];

dC_dt[1] = x_r[3] + a_MS * F_M - F_S;
dC_dt[2] = x_r[4] + (1 - a_MS) * F_M + F_S + F_L - F_D;
dC_dt[3] = F_D * E_C_f - F_M - F_E;
dC_dt[4] = F_E - F_L;

return dC_dt;}

/*
*Function for calculating ratio of evolved flux to initial flux at each time observation point
*/

real[] model_ratios(int N_t, real t0, real[] ts, real V_ref, real V_U_ref, real Ea_V, real Ea_VU, real Ea_K, real Ea_KU, real E_C_ref, real m_t, real a_MS, real[] x_r, int[] x_i) {

real theta[12]; // ODE parameters
real C_hat[N_t,4]; // Predicted pool content
real CO2_flux_0; // Initial CO2 flux at equilibrium
real CO2_flux_hat[N_t]; // Model CO2 flux array
real CO2_flux_ratios_hat[N_t]; // Container array for ratio of evolved flux to initial flux
real C_t0[4]; // Initial conditions
real S_0; // Initial SOC
real D_0; // Initial DOC
real M_0; // Initial MIC
real E_0; // Initial ENZC
real I_S; // SOC input
real I_D; // DOC input
real r_L; // Enzyme loss rate
real r_E; // Enzyme production rate

/*
*Dependent variables
*/

real r_M;
real K_U_ref;
real K_ref;
real V_f; // V after forcing
real V_U_f;
real K_f;
real K_U_f; // K_U after forcing
real E_C_f; // E_C after forcing

I_S = x_r[3];
I_D = x_r[4];

S_0 = x_r[5];
D_0 = x_r[6];
M_0 = x_r[7];
E_0 = x_r[8];

r_L = x_r[9];
r_E = x_r[10];

C_t0[1] = S_0;
C_t0[2] = D_0;
C_t0[3] = M_0;
C_t0[4] = E_0;

r_M = (-E_C_ref * (I_D + I_S) + M_0 * r_E * (1 - E_C_ref)) / (M_0 * (E_C_ref - 1));
K_U_ref = -(D_0 * (r_M + r_E - E_C_ref * V_U_ref)) / (r_M + r_E);
K_ref = -((S_0 * (-I_S * r_E * r_L + E_C_ref * I_S * r_E * r_L - a_MS * E_C_ref * I_D * r_L * r_M - I_S * r_L * r_M + E_C_ref * I_S * r_L * r_M - a_MS * E_C_ref * I_S * r_L * r_M + E_C_ref * I_D * r_E * V_ref + E_C_ref * I_S * r_E * V_ref)) / (r_L * (-I_S * r_E + E_C_ref * I_S * r_E - a_MS * E_C_ref * I_D * r_M - I_S * r_M + E_C_ref * I_S * r_M - a_MS * E_C_ref * I_S * r_M)));

CO2_flux_0 = ((V_U_ref * C_t0[3] * C_t0[2]) / (K_U_ref + C_t0[2])) * (1 - E_C_ref);

print("S_0, D_0, M_0, E_0, V_ref, V_U_ref, K_ref, K_U_ref, Ea_V, Ea_VU, Ea_K, Ea_KU, r_E, r_L, r_M, E_C_ref, m_t, a_MS = ", S_0, ",", D_0, ",", M_0, ",", E_0, ",", V_ref, ",", V_U_ref, ",", K_ref, ",", K_U_ref, ",", Ea_V, ",", Ea_VU, ",", Ea_K, ",", Ea_KU, ",", r_E, ",", r_L, ",", r_M, ",", E_C_ref, ",", m_t, ",", a_MS);
print("CO2_flux_0 = ", CO2_flux_0);

V_f = decay(V_ref, Ea_V, x_r[1] + 3, x_r); // V after forcing calculation
V_U_f = decay(V_U_ref, Ea_VU, x_r[1] + 3, x_r); // V_U after forcing calculation
K_f = decay(K_ref, Ea_K, x_r[1] + 3, x_r); // K after forcing calculation
K_U_f = decay(K_U_ref, Ea_KU, x_r[1] + 3, x_r); // K_U after forcing calculation
E_C_f = tempMod(E_C_ref, m_t, x_r[1] + 3, x_r); // E_C after forcing calculation

print("E_C_f = ", E_C_f);

theta[1] = V_f;
theta[2] = V_U_f;
theta[3] = K_f;
theta[4] = K_U_f;
theta[5] = Ea_V;
theta[6] = Ea_VU;
theta[7] = Ea_K;
theta[8] = Ea_KU;
theta[9] = r_M;
theta[10] = a_MS;
theta[11] = E_C_f;
theta[12] = m_t;

C_hat = integrate_ode_bdf(AWB_ODE, C_t0, t0, ts, theta, x_r, x_i);

for (i in 1:N_t) {
  CO2_flux_hat[i] = ((V_U_f * C_hat[i,3] * C_hat[i,2]) / (K_U_f + C_hat[i,2])) * (1 - E_C_f);
  CO2_flux_ratios_hat[i] = CO2_flux_hat[i] / CO2_flux_0;  
}

print("CO2_flux_ratios_hat array = ", CO2_flux_ratios_hat);

return CO2_flux_ratios_hat;

}
}

data {
int<lower=1> N_t; // Number of time points to fit
int<lower=1> N_p; // Number of predictive time points
real<lower=0> t0; // Initial time
real<lower=0> ts[N_t]; // Observation times
real<lower=0> ts_p[N_p]; // Predictive observation times

vector<lower=0>[N_t] CO2_flux_ratios_vector; // Vector of data flux ratios to be fit
real<lower=0> SOC_input; // Mean hourly SOC input
real<lower=0> DOC_input; // Mean DOC input
real<lower=0> T_ref; // Reference temperature
real<lower=0> R_g; // Ideal gas constant 8.314

real<lower=0> S_0; // Steady state initial SOC pool size
real<lower=0> D_0; // Steady state initial DOC pool size
real<lower=0> M_0; // Steady state initial MIC pool size
real<lower=0> E_0; // Steady state initial ENZC pool size

real<lower=0,upper=1> r_L; // Enzyme loss rate
real<lower=0,upper=1> r_E; // Enzyme production rate
}

transformed data {
real x_r[10]; // x_r[1] is T_ref. x_r[2] is gas constant. x_r[3] is input SOC. x_r[4] is input DOC. 
int x_i[0]; // No int constants for system

x_r[1] = T_ref;
x_r[2] = R_g;
x_r[3] = SOC_input;
x_r[4] = DOC_input;

x_r[5] = S_0;
x_r[6] = D_0;
x_r[7] = M_0;
x_r[8] = E_0;

x_r[9] = r_L;
x_r[10] = r_E;
}

parameters {
real<lower=0,upper=1> V_ref; // SOC reference V_max
real<lower=0.001,upper=0.1> V_U_ref; // DOC reference V_max
real<lower=10,upper=90> Ea_V; // SOC V_max activation energy
real<lower=10,upper=90> Ea_VU; // DOC V_max activation energy
real<lower=10,upper=90> Ea_K; // SOC K_m activation energy
real<lower=10,upper=90> Ea_KU; // DOC K_m activation energy 
real<lower=0.01,upper=0.9> E_C_ref; // Carbon use efficiency
real<lower=0,upper=0.03> m_t; // CUE slope
real<lower=0.1,upper=0.9> a_MS; // Fraction of dead MBC transferred to SOC

real<lower=0> sigma; // Observation standard deviation
}

transformed parameters {
real<lower=0> CO2_flux_ratios_hat[N_t];
vector[N_t] CO2_flux_ratios_hat_vector;

CO2_flux_ratios_hat = model_ratios(N_t, t0, ts, V_ref, V_U_ref, Ea_V, Ea_VU, Ea_K, Ea_KU, E_C_ref, m_t, a_MS, x_r, x_i);
CO2_flux_ratios_hat_vector = to_vector(CO2_flux_ratios_hat);
}

model {
// informative priors
V_ref ~ normal(0.4,0.2); // SOC reference V_max
V_U_ref ~ normal(0.01,0.005); // DOC reference V_max
Ea_V ~ normal(50,25); // SOC V_max activation energy
Ea_VU ~ normal(50,25); // DOC V_max activation energy
Ea_K ~ normal(50,25); // SOC K_m activation energy
Ea_KU ~ normal(50,25); // DOC K_m activation energy 
E_C_ref ~ normal(0.4,0.2); // Reference CUE
m_t ~ normal(0.002,0.001); // CUE slope
a_MS ~ normal(0.5,0.25); // Fraction of dead MBC transferred to SOC

sigma ~ cauchy(0,1); // Residual error scale

// likelihood
CO2_flux_ratios_vector ~ normal(CO2_flux_ratios_hat_vector, sigma);
}

generated quantities{
vector[N_t] log_lik;
vector[N_t] CPOi; // Inverse CPO vector for LPML calculation
real<lower=0> CO2_flux_ratios_p[N_p];
vector[N_p] CO2_flux_ratios_new_vector; // Predictive time series vector

for (i in 1:N_t) {
  log_lik[i] = normal_lpdf(CO2_flux_ratios_vector[i] | CO2_flux_ratios_hat_vector[i], sigma);
} // Log likelihood vector for LOO package processing
CPOi = sqrt(2 * pi()) * sigma * exp(0.5 * (1 / sigma ^ 2) * square(CO2_flux_ratios_vector - CO2_flux_ratios_hat_vector));
CO2_flux_ratios_p = model_ratios(N_p, t0, ts_p, V_ref, V_U_ref, Ea_V, Ea_VU, Ea_K, Ea_KU, E_C_ref, m_t, a_MS, x_r, x_i);
// normal_rng cannot be vectorized
for (i in 1:N_p) {
  CO2_flux_ratios_new_vector[i] = normal_rng(CO2_flux_ratios_p[i], sigma);
}
}

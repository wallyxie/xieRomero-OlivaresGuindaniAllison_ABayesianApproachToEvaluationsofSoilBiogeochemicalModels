functions {

/*
*CON model as depicted in Li et al., 2014
*/

/*
*Arrhenius relationship for calculating decay constants
*/

real decay(real k_ref, real Ea, real temp, real[] x_r) {
  return k_ref * exp(-Ea / x_r[2] * (1 / temp - 1 / x_r[1]));
}

/*
*CON ODE system described in Li et al., 2014
*C[1] is mass of SOC pool
*C[2] is mass of DOC pool
*C[3] is mass of microbial pool
*/

real[] CON_ODE(real t, real[] C, real[] theta, real[] x_r, int[] x_i) {

real dC_dt[3];

//parameters
real Ea_S; // SOC activation energy
real Ea_D; // DOC activation energy
real Ea_M; // Microbial activation energy
real a_DS; // DOC to SOC transfer coefficient
real a_SD; // SOC to DOC transfer coefficient
real a_M; // Microbial to soil C transfer coefficient
real a_MS; // Fraction of dead microbes transferred to SOC
real u_M; // Microbial DOC uptake rate
real k_S_f; // k_S after forcing
real k_D_f; // k_D after forcing
real k_M_f; // k_M after forcing

/*
*Dependent variables
*/

real F_S;
real F_D;
real F_M;

/*
*Match parameters to thetas
*/

Ea_S = theta[1];
Ea_D = theta[2];
Ea_M = theta[3];
a_DS = theta[4];
a_SD = theta[5];
a_M = theta[6];
a_MS = theta[7];
k_S_f = theta[8];
k_D_f = theta[9];
k_M_f = theta[10];

u_M = x_r[8];

F_S = k_S_f * C[1];
F_D = k_D_f * C[2];
F_M = k_M_f * C[3];

dC_dt[1] = x_r[3] + a_DS * F_D + a_M * a_MS * F_M - F_S;
dC_dt[2] = x_r[4] + a_SD * F_S + a_M * (1 - a_MS) * F_M - u_M * C[2] - F_D;
dC_dt[3] = u_M * C[2] - F_M;

return dC_dt;}

/*
*Function for calculating ratio of evolved flux to initial flux at each time observation point
*/

real[] model_ratios(int N_t, real t0, real[] ts, real Ea_S, real Ea_D, real Ea_M, real a_DS, real a_SD, real a_M, real a_MS, real[] x_r, int[] x_i) {

real theta[10]; // ODE parameters
real C_hat[N_t,3]; // Predicted pool content
real CO2_flux_hat[N_t]; // Model CO2 flux array
real CO2_flux_ratios_hat[N_t]; // Container array for ratio of evolved flux to initial flux
real C_t0[3]; // Initial conditions
real CO2_flux_0; // Initial CO2 flux at equilibrium
real S_0;
real D_0;
real M_0;
real u_M;

/*
*Dependent variables
*/

real k_S_ref;
real k_D_ref;
real k_M_ref;

real k_S_f; // k_S after forcing
real k_D_f; // k_D after forcing
real k_M_f; // k_M after forcing

// print("Ea_S, Ea_D, Ea_M, a_DS, a_SD, a_M, a_MS = ", Ea_S, ",", Ea_D, ",", Ea_M,",", a_DS, ",", a_SD, ",", a_M,",", a_MS);

theta[1] = Ea_S;
theta[2] = Ea_D;
theta[3] = Ea_M;
theta[4] = a_DS;
theta[5] = a_SD;
theta[6] = a_M;
theta[7] = a_MS;

S_0 = x_r[5];
D_0 = x_r[6];
M_0 = x_r[7];
u_M = x_r[8];
k_M_ref = x_r[9];

C_t0[1] = S_0;
C_t0[2] = D_0;
C_t0[3] = M_0;

k_D_ref = (-x_r[4] - (a_SD * x_r[3]) + (D_0 * u_M) - (a_M * D_0 * u_M) + (a_M * a_MS * D_0 * u_M) - (a_M * a_MS * a_SD * D_0 * u_M)) / (D_0 * (-1 + a_DS * a_SD)); // Calculating initial k_D
k_S_ref = (x_r[3] + D_0 * (a_DS * k_D_ref + u_M * a_M * a_MS)) / x_r[5]; // Calculating initial k_S

k_S_f = decay(k_S_ref, Ea_S, x_r[1] + 3, x_r);
k_D_f = decay(k_D_ref, Ea_D, x_r[1] + 3, x_r);
k_M_f = decay(k_M_ref, Ea_M, x_r[1] + 3, x_r);

theta[8] = k_S_f;
theta[9] = k_D_f;
theta[10] = k_M_f;

// print("Calculated k_S_f, k_D_f, k_M_f = ", k_S_f, ",", k_D_f, ",", k_M_f);

CO2_flux_0 = k_S_ref * S_0 * (1 - a_SD) + k_D_ref * D_0 * (1 - a_DS) + k_M_ref * M_0 * (1 - a_M);

// print("k_S_ref, k_D_ref, k_M_ref = ", k_S_ref, ",", k_D_ref, ",", k_M_ref);

// print("CO2_flux_0 = ", CO2_flux_0);

C_hat = integrate_ode_bdf(CON_ODE, C_t0, t0, ts, theta, x_r, x_i);

for (i in 1:N_t) {
  CO2_flux_hat[i] = k_S_f * C_hat[i,1] * (1 - a_SD) + k_D_f * C_hat[i,2] * (1 - a_DS) + k_M_f * C_hat[i,3] * (1 - a_M);
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
real<lower=0> u_M; // Microbial DOC uptake rate
real<lower=0> k_M_ref; //
}

transformed data {
real x_r[9]; // x_r[1] is T_ref. x_r[2] is gas constant. x_r[3] is input SOC mean. x_r[4] is input DOC mean. x_r[5:7] are steady state values.
int x_i[0]; // No int constants for system

x_r[1] = T_ref;
x_r[2] = R_g;
x_r[3] = SOC_input;
x_r[4] = DOC_input;
x_r[5] = S_0;
x_r[6] = D_0;
x_r[7] = M_0;
x_r[8] = u_M;
x_r[9] = k_M_ref;
}

parameters {

real<lower=10,upper=90> Ea_S; // SOC activation energy
real<lower=10,upper=90> Ea_D; // DOC activation energy
real<lower=10,upper=90> Ea_M; // Microbial activation energy
real<lower=0.1,upper=0.9> a_DS; // DOC to SOC transfer coefficient
real<lower=0.1,upper=0.9> a_SD; // SOC to DOC transfer coefficient
real<lower=0.1,upper=0.9> a_M; // Microbial to soil C transfer coefficient
real<lower=0.1,upper=0.9> a_MS; // Fraction of dead microbes transferred to SOC

real<lower=0> sigma; // Observation standard deviation
}

transformed parameters {
real<lower=0> CO2_flux_ratios_hat[N_t];
vector[N_t] CO2_flux_ratios_hat_vector;

CO2_flux_ratios_hat = model_ratios(N_t, t0, ts, Ea_S, Ea_D, Ea_M, a_DS, a_SD, a_M, a_MS, x_r, x_i);
CO2_flux_ratios_hat_vector = to_vector(CO2_flux_ratios_hat);
}

model {
// informative priors
Ea_S ~ normal(50,25);
Ea_D ~ normal(50,25);
Ea_M ~ normal(50,25);
a_DS ~ normal(0.3,0.15);
a_SD ~ normal(0.3,0.15);
a_M ~ normal(0.3,0.15);
a_MS ~ normal(0.5,0.25);

sigma ~ cauchy(0,1);

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
CO2_flux_ratios_p = model_ratios(N_p, t0, ts_p, Ea_S, Ea_D, Ea_M, a_DS, a_SD, a_M, a_MS, x_r, x_i);
// normal_rng cannot be vectorized
for (i in 1:N_p) {
  CO2_flux_ratios_new_vector[i] = normal_rng(CO2_flux_ratios_p[i], sigma);
}
}

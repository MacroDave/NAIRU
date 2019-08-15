data {
  int<lower=1> T; // number of observations
  int<lower=1> J; // dimension of observations
  matrix[T,J] Y; // observations
}

parameters {
 //NAIRU Equation
	  vector[T] NAIRU; // unobserved NAIRU series
	  real<lower = 0> tau; //variance of unobserved NAIRU series

  //Inflation Equation
	  real<lower = 0, upper = 1> delta_pt; //coef on inflation expectations
	  vector[3] beta_pt; //coef on lagged inflation
	  real phi_pt; //coef on lagged ULC
	  real gamma_pt; //coef on UR gap
	  real lambda_pt; //coef on UR change
	  real alpha_pt; //coef on import prices
    real<lower = 0> eps_pt; //variance of inflation eqn error term
  
  //ULC Equation
	  real<lower = 0, upper = 1> delta_pu; //coef on inflation expectations
	  vector[3] beta_pu; //coef on lagged inflation
	  real gamma_pu; //coef on UR gap
	  real lambda_pu; //coef on UR change
	  real<lower = 0> eps_pu; //variance of ULC eqn error term 
}

model {
  //priors
  	//NAIRU
  	NAIRU[1] ~ normal(7,1);  //inital value of 'true' GDP
  	tau ~ normal(0,1);  //variance of NAIRU error term
  	 	 
  	//Inflation Equation
  	 delta_pt ~ normal(0.5,1); //coef on inflation expectations
  	 beta_pt ~ normal(0.5,1); //coef on lagged inflation
  	 phi_pt ~ normal(0.5,1); //coef on lagged ULC
  	 gamma_pt ~ normal(0,1); //coef on UR gap
  	 lambda_pt ~ normal(0,1); //coef on UR change
     alpha_pt ~ normal(0,1); //coef on import prices
     eps_pt ~ normal(0,1); //variance of inflation eqn error term
  	  
  	//ULC Equation
  	 delta_pu ~ normal(0,1); //coef on inflation expectations
  	 beta_pu ~ normal(0,1); //coef on lagged inflation
  	 gamma_pu ~ normal(0,1); //coef on UR gap
  	 lambda_pu ~ normal(0,1); //coef on UR change
  	 eps_pu ~ normal(0,1); //variance of ULC eqn error term  
  
  //likelihood
  //NAIRU
	for(t in 2:T) {
		NAIRU[t] ~ normal(NAIRU[t-1], tau);
	}
	
	//INFLATION
	for(t in 4:T) {
		Y[t,5] ~ normal(delta_pt*Y[t,6] 
						+ beta_pt[1]*Y[t-1,5] 
						+ beta_pt[2]*Y[t-2,5]
						+ beta_pt[3]*Y[t-3,5]
						+ phi_pt*Y[t-1,3]
						+ gamma_pt*(Y[t,7]-NAIRU[t])/Y[t,7]
						+ lambda_pt*(Y[t-1,7]-Y[t-2,7])/Y[t,7]
						+ alpha_pt*(Y[t-1,4]-Y[t-2,4])
						, eps_pt);
	}

	 //ULC
	 for(t in 4:T) {
		Y[t,3] ~ normal(delta_pu*Y[t,6] 
						+ beta_pu[1]*Y[t-1,3] 
						+ beta_pu[2]*Y[t-2,3]
						+ gamma_pu*(Y[t,7]-NAIRU[t])/Y[t,7]
						+ lambda_pu*(Y[t-1,7]-Y[t-2,7])/Y[t,7]
						, eps_pu);
	}
}

data {
  int<lower=1> T; // number of observations
  int<lower=1> J; // dimension of observations
  matrix[T,J] Y; // observations
}

parameters {
 //NAIRU Equation
  vector[T] NAIRU; 
  vector[1] nhat_init;
  vector[3] pthat_init;
  vector[2] puhat_init;
  
  real<lower = 0> tau; //variance of unobserved NAIRU series
  
  //Inflation Equation
	  real<lower = 0, upper = 1> delta_pt; //coef on inflation expectations
	  vector[3] beta_pt; //coef on lagged inflation
	  real gamma_pt; //coef on UR gap
	  real lambda_pt; //coef on UR change
    real<lower = 0> eps_pt; //variance of inflation eqn error term
  
  //ULC Equation
	  real<lower = 0, upper = 1> delta_pu; //coef on inflation expectations
	  vector[2] beta_pu; //coef on lagged inflation
	  real gamma_pu; //coef on UR gap
	  real lambda_pu; //coef on UR change
	  real<lower = 0> eps_pu; //variance of ULC eqn error term 

}

model {
  //priors
  	//Inflation Equation
  	 delta_pt ~ normal(0.5,1); //coef on inflation expectations
  	 beta_pt ~ normal(0.5,1); //coef on lagged inflation
  	 gamma_pt ~ normal(1,1); //coef on UR gap
  	 lambda_pt ~ normal(-0.5,1); //coef on UR change
     eps_pt ~ normal(0,1); //variance of inflation eqn error term
  	  
  	//ULC Equation
  	 delta_pu ~ normal(1,1); //coef on inflation expectations
  	 beta_pu ~ normal(0,1); //coef on lagged inflation
  	 gamma_pu ~ normal(1.5,1); //coef on UR gap
  	 lambda_pu ~ normal(-5,1); //coef on UR change
  	 eps_pu ~ normal(0,1); //variance of ULC eqn error term  

    //NAIRU
    tau ~ normal(0.1,0.2);  //variance of NAIRU error term

  nhat_init ~ normal(6, 0.2);
  pthat_init ~ normal(2.5, 1);
  puhat_init ~ normal(2.5, 1);
  
  {

  vector[T] nairu_hat;
  vector[T] pt_hat;
  vector[T] pu_hat;

  nairu_hat[1] = nhat_init[1];
  pt_hat[1:3] = pthat_init[1:3];
  pu_hat[1:2] = puhat_init[1:2];

  	for(t in 2:T) {
      nairu_hat[t] = NAIRU[t-1];
  	}
  	
  	for(t in 4:T) {
      pt_hat[t] = -gamma_pt
  		        + delta_pt*Y[t,6] 
  						+ beta_pt[1]*Y[t-1,4]
  						+ beta_pt[2]*Y[t-2,4]
  						+ beta_pt[3]*Y[t-3,4]
  						+ gamma_pt*(NAIRU[t]/Y[t,3])
						  + lambda_pt*(Y[t-1,3]-Y[t-2,3])/Y[t,3];
  	}

  	for(t in 3:T) {
      pu_hat[t] = -gamma_pu
  		        + delta_pu*Y[t,6] 
  						+ beta_pu[1]*Y[t-1,4] 
  						+ beta_pu[2]*Y[t-2,4]
  						+ gamma_pu*(NAIRU[t]/Y[t,3])
						  + lambda_pu*(Y[t-1,3]-Y[t-2,3])/Y[t,3];
  	}

    target += normal_lpdf(NAIRU | nairu_hat, tau);
    target += normal_lpdf(Y[,4] | pt_hat, eps_pt);
    target += normal_lpdf(Y[,1] | pu_hat, eps_pu);
  	
  }

}

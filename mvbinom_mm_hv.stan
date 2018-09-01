functions {
  vector interact(vector y) {
    int D = num_elements(y);
    vector[D*(D-1)/2] x;
    int k=1;
    
    for (i in 1:(D-1)) {
      for (j in (i+1):D) {
        x[k] = y[i]*y[j];
        k = k + 1;
      }
    }
    return x;
  }
  
  matrix create_Yset(int D, int D2) {
    matrix[D,D2] Yset;
    for (i in 1:D2) {
      int ind = i-1;
      
      for (j in 1:D) {
        Yset[j,i] = 2*fmod(ind,2)-1;
        ind = ind/2;
      }
    }
    return Yset;
  }
}

data {
  int N;    //total number of observations
  int M;    //number of groups
  int D;    //number of dimensions
  int D2;   //ugh wtf stan this is literally just 2^D but stan will not let that be an integer
  int n[M]; //number of observations per group
  int P;    //number of regressors
  
  matrix[D,N] Y;    //N D-dimensional data vectors
  matrix[P,M] X;    //
}

transformed data {
  vector[D] YSS[M];
  vector[D*(D-1)/2] XSS[M];
  matrix[D,D2] Yset;
  matrix[D*(D-1)/2,D2] Xset;
  matrix[D,N] Yc = 2*Y-1;
  
  {
    int start = 1;
    for (j in 1:M) {  
      XSS[j] = rep_vector(0,D*(D-1)/2);
      for (i in 1:D) YSS[j][i] = sum(Yc[i,start:(start+n[j]-1)]);
      for (i in start:(start+n[j]-1)) XSS[j] = XSS[j] + interact(Yc[:,i]);
      start = start + n[j];
    }
  }
  Yset = create_Yset(D,D2);
  for (i in 1:D2) Xset[:,i] = interact(Yset[:,i]);
  
}

parameters {
  vector[D] f1_mu;
  vector[D*(D-1)/2] f2_mu;
  vector<lower=0>[D] f1_sigma;
  vector<lower=0>[D*(D-1)/2] f2_sigma_raw;
  real<lower=0> tau_f2_sigma;
  
  matrix[P,D] f1_beta;
  
  vector[D] f1_raw[M];
  vector[D*(D-1)/2] f2_raw[M];
}

transformed parameters {
  vector[D] f1[M];
  vector[D*(D-1)/2] f2[M];
  vector[D*(D-1)/2] f2_sigma = tau_f2_sigma*f2_sigma_raw;
  
  for (i in 1:M) {
    f1[i] = f1_mu + X[,i]'*f1_beta + f1_raw[i].*f1_sigma;
    f2[i] = f2_mu + f2_raw[i].*f2_sigma;
  }
}

model {
  for (i in 1:M) {
    real logZ = log_sum_exp(f1[i]'*Yset + f2[i]'*Xset);
    target +=  YSS[i]'*f1[i] + XSS[i]'*f2[i] - n[i]*logZ;
    f1_raw[i] ~ normal(0,1);
    f2_raw[i] ~ normal(0,1);
  }
  
  f1_mu ~ normal(0,5.0);
  f2_mu ~ normal(0,1.0);
  f1_sigma ~ normal(0,1);
  tau_f2_sigma ~ normal(0,1);
  f2_sigma_raw ~ normal(0,1);
}

generated quantities {
  real lp[N];
  
  {
    int start = 1;
    for (i in 1:M) {
      int fin = n[i] + start - 1;
      real logZ = log_sum_exp(f1[i]'*Yset + f2[i]'*Xset);
      for (j in start:fin) lp[j] = Yc[:,j]'*f1[i] + interact(Yc[:,j])'*f2[i] - logZ;
      start = start + n[i];
    }
  }
}

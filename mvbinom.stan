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
  
  real calc_Z(vector f1, vector f2) {
    int D = num_elements(f1);
    int D2 = 1;
    real Z = 0;
    
    for (i in 1:D) D2 = D2*2; ##wtf stan
    
    for (i in 1:D2) {
      vector[D] y;
      int ind = i-1;
      
      for (j in 1:D) {
        y[j] = fmod(ind,2);
        ind = ind/2;
      }
      Z = Z + exp(y'*f1 + interact(y)'*f2);
    }
    return Z;
  }
  
  real multi_bernoulli_ss_lp(vector y, vector x, vector f1, vector f2, int n) {
    real Z = calc_Z(f1,f2);
    return y'*f1 + x'*f2 - n*log(Z);
  }
}

data {
  int N;    //total number of observations
  // int M;    //number of groups
  int D;    //number of dimensions
  // int n[M]; //number of observations per group
  
  matrix[D,N] Y;    //N D-dimensional data vectors
}

transformed data {
  vector[D] YSS;
  vector[D*(D-1)/2] XSS=rep_vector(0,D*(D-1)/2);
  
  for (i in 1:D) YSS[i] = sum(Y[i,:]);
  for (i in 1:N) XSS = XSS + interact(Y[:,i]);

  // {
  //   int start = 1;
  //   for (i in 1:M) {
  //     YSS[i] = Y[:,start:(start+n[i]-1)] * Y[:,start:(start+n[i]-1)]';
  //     start = start + n[i];
  //   }
  // }
}

parameters {
  vector[D] f1;
  vector[D*(D-1)/2] f2;
  // vector<lower=0>[D] f1_sigma;
  // vector<lower=0>[D2] f2_sigma;
  // 
  // vector[D] f1_raw[M];
  // vector[D2] f2_raw[M];
}

// transformed parameters {
//   vector[D] f1[M];
//   vector[D2] f2[M];
//     
//   for (i in 1:M) {
//     f1[i] = f1_mu + f1_raw[i].*f1_sigma;
//     f2[i] = f2_mu + f2_raw[i].*f2_sigma;
//   }
// }

model {
  // for (i in 1:M) {
  target += multi_bernoulli_ss_lp(YSS,XSS,f1, f2, N);
  //   f1_raw[i] ~ normal(0,1);
  //   f2_raw[i] ~ normal(0,1);
  // }
  // 
  f1 ~ normal(0,1.0);
  f2 ~ normal(0,0.5);
  // f1_sigma ~ normal(0,1); 
  // f2_sigma ~ normal(0,1);
}

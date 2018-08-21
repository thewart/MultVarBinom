functions {
  real multi_bernoulli_lpdf(matrix ss, vector f1, vector f2,  int n) {
    real fy = 0;
    real A = 1.0;
    int P = rows(ss);
    int k = 1;
    
    for (i in 1:(P-1)) {
      A = A + exp(f1[i]);
      fy = fy + f1[i]*ss[i,i];
      for (j in (i+1):P) {
        A = A + exp(f1[i] + f1[j] + f2[k]);
        fy = fy + f2[k]*ss[j,i];
        k = k + 1;
      }
    }
    A = A + exp(f1[P]);
    fy = fy + f1[P]*ss[P,P];
    
    // A = 1 + exp(f1[1]) + exp(f1[2]) + exp(f1[1] + f1[2] + f2[1]);
    // fy = ss[1,1]*f1[1] + ss[2,2]*f1[2] + ss[2,1]*f2[1];
    
    return fy - n*log(A);
  }
}

data {
  int N;    //total number of observations
  int M;    //number of groups
  int D;    //number of dimensions
  int n[M]; //number of observations per group
  
  matrix[D,N] Y;    //N D-dimensional data vectors
}

transformed data {
  matrix[D,D] YSS[M];
  int D2 = D*(D-1)/2;

  {
    int start = 1;
    for (i in 1:M) {
      YSS[i] = Y[:,start:(start+n[i]-1)] * Y[:,start:(start+n[i]-1)]';
      start = start + n[i];
    }
  }
}

parameters {
  vector[D] f1_mu;
  vector[D2] f2_mu;
  vector<lower=0>[D] f1_sigma;
  vector<lower=0>[D2] f2_sigma;
  
  vector[D] f1_raw[M];
  vector[D2] f2_raw[M];
}

transformed parameters {
  vector[D] f1[M];
  vector[D2] f2[M];
    
  for (i in 1:M) {
    f1[i] = f1_mu + f1_raw[i].*f1_sigma;
    f2[i] = f2_mu + f2_raw[i].*f2_sigma;
  }
}

model {
  for (i in 1:M) {
    YSS[i] ~ multi_bernoulli(f1[i], f2[i], n[i]);
    f1_raw[i] ~ normal(0,1);
    f2_raw[i] ~ normal(0,1);
  }
  
  f1_mu ~ normal(0,2.5);
  f2_mu ~ normal(0,2.5);
  f1_sigma ~ normal(0,1); 
  f2_sigma ~ normal(0,1);
}

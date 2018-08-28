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
        Yset[j,i] = fmod(ind,2);
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
  int D2;   //ugh wtf stan
  int n[M]; //number of observations per group
  
  matrix[D,N] Y;    //N D-dimensional data vectors
}

transformed data {
  vector[D] YSS[M];
  vector[D*(D-1)/2] XSS[M];
  matrix[D,D2] Yset;
  matrix[D*(D-1)/2,D2] Xset;
  
  {
    int start = 1;
    for (j in 1:M) {  
      XSS[j] = rep_vector(0,D*(D-1)/2);
      for (i in 1:D) YSS[j][i] = sum(Y[i,start:(start+n[j]-1)]);
      for (i in start:(start+n[j]-1)) XSS[j] = XSS[j] + interact(Y[:,i]);
      start = start + n[j];
    }
  }
  Yset = create_Yset(D,D2);
  for (i in 1:D2) Xset[:,i] = interact(Yset[:,i]);
  
  // {
    //   int start = 1;
    //   for (i in 1:M) {
      //     YSS[i] = Y[:,start:(start+n[i]-1)] * Y[:,start:(start+n[i]-1)]';
      //     start = start + n[i];
      //   }
      // }
}
    
    parameters {
    vector[D] f1[M];
    vector[D*(D-1)/2] f2[M];
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
    for (i in 1:M) {
      real logZ = log_sum_exp(f1[i]'*Yset + f2[i]'*Xset);
      target +=  YSS[i]'*f1[i] + XSS[i]'*f2[i] - n[i]*logZ;
      // f1[i] ~ normal(0,1);
      // f2[i] ~ normal(0,1);
    }
    // 
      // f1 ~ normal(0,1.0);
      // f2 ~ normal(0,0.5);
    // f1_sigma ~ normal(0,1); 
    // f2_sigma ~ normal(0,1);
    }
    
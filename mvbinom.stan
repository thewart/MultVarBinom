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
  // int M;    //number of groups
  int D;    //number of dimensions
  int D2;   //ugh wtf stan
  // int n[M]; //number of observations per group
  
  matrix[D,N] Y;    //N D-dimensional data vectors
}

transformed data {
  vector[D] YSS;
  vector[D*(D-1)/2] XSS=rep_vector(0,D*(D-1)/2);
  matrix[D,D2] Yset;
  matrix[D*(D-1)/2,D2] Xset;
  matrix[D,N] Yc = 2*Y-1;
    
  for (i in 1:D) YSS[i] = sum(Yc[i,:]);
  for (i in 1:N) XSS = XSS + interact(Yc[:,i]);
  
  Yset = create_Yset(D,D2);
  for (i in 1:D2) Xset[:,i] = interact(Yset[:,i]);
}
    
parameters {
  vector[D] f1;
  vector[D*(D-1)/2] f2;
}

model {
  real logZ = log_sum_exp(f1'*Yset + f2'*Xset);
  target +=  YSS'*f1 + XSS'*f2 - N*logZ;
  
  // f1 ~ normal(0,1.0);
  // f2 ~ normal(0,0.5);
}
    
generated quantities {
  real lp[N];
  real logZ = log_sum_exp(f1'*Yset + f2'*Xset);
  
  for (i in 1:N)
    lp[i] = Yc[:,i]'*f1 + interact(Yc[:,i])'*f2 - logZ;
}

##Make sure et/lateral are/not in log space. 
data {
 // Options
 // int<lower=0, upper=1> inc_none;
 // int<lower=0, upper=1> inc_et; // include et? 0=no, 1=yes
 // int<lower=0, upper=1> inc_lateral; // include lateral? 0=no, 1=yes
 // int<lower=0, upper=1> inc_et_lateral; // include et and lateral? 0=no, 1=yes

 // bounds on parameters
 real nInlower;
 real nInupper;
 real nOutlower;
 real nOutupper;
 real aInlower;
 real aInupper;
 real aOutlower;
 real aOutupper;
 real daInShift;
 real daOutShift;
 real nInHat;
 real nInSd;
 real aInHat;
 real aInSd;
 real nOutHat;
 real nOutSd;
 real aOutHat;
 real aOutSd;
##Change n, a, and nhat, ahat to be either real or [N]

 int N; // Sample size
 
 vector[N] sigmaIn;
 vector[N] sigmaOut;
 vector[N] qInSd;
 vector[N] qOutSd;
 vector[N] q;
 
 vector[N] da;
 vector[N] w;
 vector[N] s;
 vector[N] da2;
 vector[N] w2;
 vector[N] s2;
 vector[N] q2;
 vector[N] dv;
 vector[N] et;
 vector[N] lateral;
}

parameters {
 real < lower=nInlower, upper=nInupper > n; // mannning's n
 real < lower=aInlower, upper=aInupper > a; // Bathymetry
 real < lower=nOutlower, upper=nOutupper > nOut; // mannning's n
 real < lower=aOutlower, upper=aOutupper > aOut; // Bathymetry
 real < lower = 0 > sigma; // Error SD
 vector < lower = 0 >[N] logQ_in;
 vector < lower = 0 >[N] logQ_out;
}

transformed parameters {
  //vector[N] lhs; // LHS for Manning likelihood
  //vector[N] rhs; // log area for Manning's equation
  vector[N] lhsIn; // LHS for Manning likelihood
  vector[N] rhsIn; // log area for Manning's equation
  vector[N] lhsOut; // LHS for Manning likelihood
  vector[N] rhsOut; // log area for Manning's equation
  vector[N] lhsDV; // LHS for Manning likelihood
  vector[N] rhsDV; // log area for Manning's equation
  
  vector[N] a1;
  vector[N] a2;
  
  a1 = log(da+a);
  a2 = log(da2+aOut);
  lhsIn = (4*w)-(3*s);
  rhsIn = ((-6*n)+(10*a1))-(6*logQ_in);
  lhsOut = (4*w2)-(3*s2);
  rhsOut = ((-6*nOut)+(10*a2))-(6*logQ_out);
  // if(inc_none){
  // lhsDV = (dv/7);
  // rhsDV = exp(logQ_in)-exp(logQ_out);
  // }
  // if(inc_et){
  // lhsDV = (dv/7)+et;
  // rhsDV = exp(logQ_in)-exp(logQ_out);
  // }
  // if(inc_lateral){
  // lhsDV = (dv/7)-lateral;
  // rhsDV = exp(logQ_in)-exp(logQ_out);
  // }
  // if(inc_et_lateral){
  // lhsDV = (dv/7)-lateral+et;
  // rhsDV = exp(logQ_in)-exp(logQ_out);
  // }
  lhsDV = (dv/7)-lateral+et;
  rhsDV = exp(logQ_in)-exp(logQ_out);
}
##compare lhsDV/rhsDv vs lhs/rhs


model {
  logQ_in~normal(q,qInSd);
  logQ_out~normal(q2,qOutSd);
  n ~ normal(nInHat,nInSd);
  a+daInShift ~lognormal(aInHat, aInSd);
  nOut ~ normal(nOutHat, nOutSd);
  aOut + daOutShift~lognormal(aOutHat, aOutSd);
  
  //lhs~normal(rhs, sigma);
  lhsIn~normal(rhsIn, sigmaIn);
  lhsOut~normal(rhsOut, sigmaOut);
  lhsDV~normal(rhsDV, sigma);

}


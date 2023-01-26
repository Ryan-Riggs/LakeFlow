data {
 // Options
 // int<lower=0, upper=1> inc_none;
 // int<lower=0, upper=1> inc_et; // include et? 0=no, 1=yes
 // int<lower=0, upper=1> inc_lateral; // include lateral? 0=no, 1=yes
 // int<lower=0, upper=1> inc_et_lateral; // include et and lateral? 0=no, 1=yes

 // bounds on parameters
 real nInlower;
 real nInupper;
 real daInShift;
 real aInlower;
 real aInupper;
 
 real nInlower2;
 real nInupper2;
 real daInShift2;
 real aInlower2;
 real aInupper2;
 
 real nOutlower;
 real nOutupper;
 real aOutlower;
 real aOutupper;
 real daOutShift;

 real nInHat;
 real nInSd;
 real aInHat;
 real aInSd;
 real nIn2Hat;
 real nIn2Sd;
 real aIn2Hat;
 real aIn2Sd;
 
 real nOutHat;
 real nOutSd;
 real aOutHat;
 real aOutSd;
##Change n, a, and nhat, ahat to be either real or [N]

 int N; // Sample size
 vector[N] sigmaIn;
 vector[N] qInSd;
 vector[N] sigmaIn2;
 vector[N] qInSd2;
 vector[N] sigmaOut;
 vector[N] qOutSd;
 vector[N] q;
 vector[N] qIn2;
 
 vector[N] da;
 vector[N] w;
 vector[N] s;
 vector[N] daIn2;
 vector[N] wIn2;
 vector[N] sIn2;
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
 
 real < lower=nInlower2, upper=nInupper2 > nIn2; // mannning's n
 real < lower=aInlower2, upper=aInupper2 > aIn2; // Bathymetry

 
 real < lower=nOutlower, upper=nOutupper > nOut; // mannning's n
 real < lower=aOutlower, upper=aOutupper > aOut; // Bathymetry
 real < lower = 0 > sigma; // Error SD
 vector < lower = 0 >[N] logQ_in;
 vector < lower = 0 >[N] logQ_in2;
 vector < lower = 0 >[N] logQ_out;
}

transformed parameters {
  //vector[N] lhs; // LHS for Manning likelihood
  //vector[N] rhs; // log area for Manning's equation
  vector[N] lhsIn; // LHS for Manning likelihood
  vector[N] rhsIn; // log area for Manning's equation
  vector[N] lhsIn2; // LHS for Manning likelihood
  vector[N] rhsIn2; // log area for Manning's equation
  vector[N] lhsOut; // LHS for Manning likelihood
  vector[N] rhsOut; // log area for Manning's equation
  vector[N] lhsDV; // LHS for Manning likelihood
  vector[N] rhsDV; // log area for Manning's equation
  
  vector[N] a1;
  vector[N] a1_2;
  vector[N] a2;
  
  a1 = log(da+a);
  a1_2 = log(daIn2+aIn2);
  a2 = log(da2+aOut);
  lhsIn = (4*w)-(3*s);
  rhsIn = ((-6*n)+(10*a1))-(6*logQ_in);
  lhsIn2 = (4*wIn2)-(3*sIn2);
  rhsIn2 = ((-6*nIn2)+(10*a1_2))-(6*logQ_in2);
  
  lhsOut = (4*w2)-(3*s2);
  rhsOut = ((-6*nOut)+(10*a2))-(6*logQ_out);
  // if(inc_none){
  // lhsDV = (dv/7);
  // rhsDV = exp(logQ_in)+exp(logQ_in2)-exp(logQ_out);
  // }
  // if(inc_et){
  // lhsDV = (dv/7)+et;
  // rhsDV = exp(logQ_in)+exp(logQ_in2)-exp(logQ_out);
  // }
  // if(inc_lateral){
  // lhsDV = (dv/7)-lateral;
  // rhsDV = exp(logQ_in)+exp(logQ_in2)-exp(logQ_out);
  // }
  // if(inc_et_lateral){
  // lhsDV = (dv/7)-lateral+et;
  // rhsDV = exp(logQ_in)+exp(logQ_in2)-exp(logQ_out);
  // }
  lhsDV = (dv/7)-lateral+et;
  rhsDV = exp(logQ_in)+exp(logQ_in2)-exp(logQ_out);
}
##compare lhsDV/rhsDv vs lhs/rhs


model {
  logQ_in~normal(q,qInSd);
  logQ_in2~normal(qIn2,qInSd2);
  logQ_out~normal(q2,qOutSd);
  n ~ normal(nInHat,nInSd);
  nIn2 ~normal(nIn2Hat,nIn2Sd);
  a+daInShift ~lognormal(aInHat, aInSd);
  aIn2+daInShift2~lognormal(aIn2Hat, aIn2Sd);
  nOut ~ normal(nOutHat, nOutSd);
  aOut + daOutShift~lognormal(aOutHat, aOutSd);
  
  //lhs~normal(rhs, sigma);
  lhsIn~normal(rhsIn, sigmaIn);
  lhsIn2~normal(rhsIn2, sigmaIn2);
  lhsOut~normal(rhsOut, sigmaOut);
  lhsDV~normal(rhsDV, sigma);

}


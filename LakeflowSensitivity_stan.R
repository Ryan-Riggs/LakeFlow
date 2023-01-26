################################################################################
##LakeFlow with corrupted data. 
################################################################################
library(rstan)
library(data.table)
library(dplyr)
library(sf)
library(raster)
outpath = 'E:\\research\\RivLake\\stan_sensitivity\\'
grades=fread('E:\\research\\RivLake\\updated\\meanQ.csv')
lakeFlowCorrupted=function(page,et,q){
page = page
ryan = list()
data = readxl::read_xlsx("E:\\research\\RivLake\\updated\\Synthetic-data-11272021.xlsx", sheet=page)
data$Date = as.Date(data$date)
data$month=lubridate::month(data$Date)
data$page=page
data = left_join(data, grades,by=c('page', 'month'))
data$sum = data$V1
all=data
##Convert to 7-day sample.
toi_df = data
period = seq.Date(min(data$Date), max(data$Date), 7)
tab=list()
start=data.table(Date=period[1],dv=0,dvCorrupted=0)
for(i in 2:length(period)){
  sub=data[data$Date<=period[i]&data$Date>period[i-1],]
  tab[[i]]=data.table(Date=period[i],dv=sum(sub$`dV (m3/day)`), dvCorrupted=sum(sub$`dV corrupted (m3/day)`))
}
tab=rbindlist(tab)
tab=bind_rows(start,tab)
data=data[data$Date%in%period,]
tab=tab[tab$Date%in%data$Date,]
data$`dV (m3/day)`=tab$dv[match(data$Date,tab$Date)]
data$`dV corrupted (m3/day)`=tab$dvCorrupted[match(data$Date,tab$Date)]
data=data[!is.na(data$`dV (m3/day)`)&data$`dV (m3/day)`!=0&!is.na(data$`dV corrupted (m3/day)`)&data$`dV corrupted (m3/day)`!=0,]
########################################################################################################
##PLD
########################################################################################################
library(rgdal)
pld <- "E:\\research\\RivLake\\SWOT_Lakes\\Harmonized_SWOT_PLD_SWORD.gdb"
# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(pld)
print(fc_list)
# Read the feature class
#reaches <- readOGR(dsn=pld,layer="Harmonized_SWORD_reaches")
lakes <- st_read(dsn=pld,layer="Harmonized_PLD_lakes")
shp=list.files('E:\\research\\RivLake\\SWOT_Lakes\\LakeFlowEligible\\', full.names=TRUE,pattern='.shp')
shp = sf::st_read(shp[7])
shp = shp[shp$names=='ALLATOONA LAKE'|shp$names=='LAKE MOHAVE;COLORADO RIVER'|shp$names=='TUTTLE CREEK LAKE',]
shp = shp[!is.na(shp$names),]
lakesFilt = lakes[lakes$lake_id%in%shp$lake_id,]
lakesFilt$name=shp$names[match(lakesFilt$lake_id, shp$lake_id)]
lakesFilt$page = c(2,4,3)
tc1='74292200081'
tc2='74292200171'

lakesFilt$us_rch_id = ifelse(lakesFilt$page==4, tc1, lakesFilt$us_rch_id)
lakesFilt$us_rch_id2=tc2

########################################################################################################
##Add in geoBAM priors. 
########################################################################################################
if(page!=4){
sos = "E:\\research\\RivLake\\src\\Confluence\\sos\\sos_netcdf\\final_sos\\constrained\\na_sword_v11_SOS_priors.nc"
ryan = ncdf4::nc_open(sos)
sos_outflow = RNetCDF::open.nc(sos)
reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
index <- which(reach_ids==lakesFilt$us_rch_id[lakesFilt$page==page], arr.ind=TRUE)
gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
geoAin = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
geoNin = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
geoNInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
geoNInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
geoAInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
geoAInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
geoAinSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
geoNinSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
qHatIn <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
sigmaIn = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
qInSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
##outflow
index <- which(reach_ids==lakesFilt$ds_rch_id[lakesFilt$page==page], arr.ind=TRUE)
gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
geoAout = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
geoNout = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
geoNOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
geoNOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
geoAOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
geoAOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
geoAoutSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
geoNoutSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
qHatOut <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
qHatIn = ifelse(is.na(qHatIn), qHatOut, qHatIn)
qHatOut = ifelse(is.na(qHatOut), qHatIn, qHatOut)
sigmaOut = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
qOutSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
# qHatIn = ifelse(is.na(qHatIn), mean(data$`in Q (m3/s)`), qHatIn)
# qHatOut = ifelse(is.na(qHatOut), mean(data$`out Q (m3/s)`), qHatOut)
################################################################################
da_shift_fun=function(x) median(x) - min(x)

stan_data = list(N=nrow(data),
                 sigmaIn=rep(sigmaIn, nrow(data)),
                 sigmaOut=rep(sigmaOut, nrow(data)),
                 qInSd = rep(log(qHatIn*.1),nrow(data)),#qInSd
                 qOutSd= rep(log(qHatOut*.1),nrow(data)),#qOutSd
                 q=rep(log(qHatIn),nrow(data)),
                 sigma=rep(0.25,nrow(data)),
                 da=(data$`in dA corrupted (m2)`-min(data$`in dA corrupted (m2)`)),
                 w=log(data$`in W corrupted (m)`),
                 s=log(data$`in S corrupted (m/m)`),
                 da2=(data$`out dA corrupted (m2)`-min(data$`out dA corrupted (m2)`)),
                 w2=log(data$`out W corrupted (m)`),
                 s2=log(data$`out S corrupted (m/m)`),
                 q2=rep(log(qHatOut),nrow(data)),
                 dv=data$`dV corrupted (m3/day)`/86400,
                 et=(data$`ET (m3/s)`*et),
                 lateral=(data$sum)*q)
stan_data$da =ifelse(is.na(stan_data$da)|is.infinite(stan_data$da), 0, stan_data$da)
stan_data$da2 =ifelse(is.na(stan_data$da2)|is.infinite(stan_data$da), 0, stan_data$da2)
stan_data$s =ifelse(is.na(stan_data$s), 0, stan_data$s)
stan_data$s2 =ifelse(is.na(stan_data$s2), 0, stan_data$s2)
stan_data$lateral =ifelse(is.na(stan_data$lateral), 0, stan_data$lateral)
stan_data$et =ifelse(is.na(stan_data$et), 0, stan_data$et)
# stan_data$inc_none=none
# stan_data$inc_et=et
# stan_data$inc_lateral=q
# stan_data$inc_et_lateral=et_q
stan_data$nInlower = geoNInlower
stan_data$nInupper = geoNInupper
stan_data$aInlower=(geoAInlower)
stan_data$aInupper =(geoAInupper)
stan_data$nOutlower = geoNInlower
stan_data$nOutupper = geoNInupper
stan_data$aOutlower=(geoAOutlower)
stan_data$aOutupper = (geoAOutupper)
stan_data$daInShift = (da_shift_fun(data$`in dA corrupted (m2)`))
stan_data$daOutShift = (da_shift_fun(data$`out dA corrupted (m2)`))

stan_data$nInHat=geoNin 
stan_data$nInSd=0.25
stan_data$aInHat=geoAin
stan_data$aInSd=geoAinSD
stan_data$nOutHat=geoNout 
stan_data$nOutSd=0.25
stan_data$aOutHat=geoAout
stan_data$aOutSd=geoAoutSD
fit=stan('E:\\research\\RivLake\\src\\lakeflowStan_ifelse.stan',
         data=stan_data,
         chains=3,
         cores=5,
         iter=4000,
         control=list(stepsize=0.5,
                      adapt_delta=0.9)
)

###Comparison with true Q.
eqn1 = function(n, a, da, w, s){
  flow = (n^-1)*((a+da)^(5/3))*(w^(-2/3))*(s^0.5)
  return(flow)
}

mn = get_posterior_mean(fit)
rw=row.names(mn)
mn = data.table(mn)
mn$type=rw
mn = mn[,c('type','mean-all chains')]

model=eqn1(exp(mn$`mean-all chains`[mn$type=='n']), (mn$`mean-all chains`[mn$type=='a']), data$`in dA corrupted (m2)`,data$`in W corrupted (m)`, data$`in S corrupted (m/m)`)
modelOut=eqn1(exp(mn$`mean-all chains`[mn$type=='nOut']), (mn$`mean-all chains`[mn$type=='aOut']), data$`out dA corrupted (m2)`,data$`out W corrupted (m)`, data$`out S corrupted (m/m)`)

output_df=data.table(inflow=model,outflow=modelOut,gage_in=data$`in Q (m3/s)`, gage_out=data$`out Q (m3/s)`,
           n=exp(mn$`mean-all chains`[mn$type=='n']), a=(mn$`mean-all chains`[mn$type=='a']),
           no=exp(mn$`mean-all chains`[mn$type=='nOut']),ao=(mn$`mean-all chains`[mn$type=='aOut']),
           page=page,et=et,q=q,type='corrupt',
           geoNin=geoNin, geoAin=geoAin,geoNout=geoNout, geoAout=geoAout,
           nInTrue=all$`in n`[1], nOutTrue=all$`out n`[1],
           aInTrue=all$`in A0 (m2)`[1], aOutTrue=all$`out A0 (m2)`[1],
           date=data$Date)
#fwrite(output_df,paste0(outpath,'lake_',page, 'corrupt_','none_',none, '_q_',q,'_et_',et,'_et_q_',et_q,'.csv'))
}else{
  sos = "E:\\research\\RivLake\\src\\Confluence\\sos\\sos_netcdf\\final_sos\\constrained\\na_sword_v11_SOS_priors.nc"
  ryan = ncdf4::nc_open(sos)
  sos_outflow = RNetCDF::open.nc(sos)
  reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
  reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
  index <- which(reach_ids==lakesFilt$us_rch_id[lakesFilt$page==page], arr.ind=TRUE)
  gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
  geoAin = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
  geoNin = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
  geoNInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
  geoNInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
  geoAInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
  geoAInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
  geoAinSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
  geoNinSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
  model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
  qHatIn <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
  sigmaIn = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
  qInSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
  ##Inflow 2
  index <- which(reach_ids==lakesFilt$us_rch_id2[lakesFilt$page==page], arr.ind=TRUE)
  gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
  geoAin2 = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
  geoNin2 = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
  geoNInlower2 = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
  geoNInupper2 = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
  geoAInlower2 = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
  geoAInupper2 = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
  geoAinSD2 = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
  geoNinSD2 = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
  model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
  qHatIn2 <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
  sigmaIn2 = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
  qInSd2=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
  ##outflow
  index <- which(reach_ids==lakesFilt$ds_rch_id[lakesFilt$page==page], arr.ind=TRUE)
  gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
  geoAout = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
  geoNout = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
  geoNOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
  geoNOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
  geoAOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
  geoAOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
  geoAoutSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
  geoNoutSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
  model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
  qHatOut <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
  qHatIn = ifelse(is.na(qHatIn), qHatOut, qHatIn)
  qHatOut = ifelse(is.na(qHatOut), qHatIn, qHatOut)
  sigmaOut = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
  qOutSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
  # qHatIn = ifelse(is.na(qHatIn), mean(data$`in Q (m3/s)`), qHatIn)
  # qHatOut = ifelse(is.na(qHatOut), mean(data$`out Q (m3/s)`), qHatOut)
  ################################################################################
  da_shift_fun=function(x) median(x) - min(x)
  
  stan_data = list(N=nrow(data),
                   q=rep(log(qHatIn),nrow(data)),
                   qIn2=rep(log(qHatIn2),nrow(data)),
                   sigma=rep(0.25,nrow(data)),
                   sigmaIn=rep(sigmaIn,nrow(data)),
                   sigmaIn2=rep(sigmaIn2,nrow(data)),
                   sigmaOut=rep(sigmaOut,nrow(data)),
                   qInSd = rep(log(qHatIn*.1),nrow(data)),#qInSd
                   qInSd2 = rep(log(qHatIn2*.1),nrow(data)),#qInSd
                   qOutSd= rep(log(qHatOut*.1),nrow(data)),#qOutSd
                   da=(data$`in dA corrupted (m2)`-min(data$`in dA corrupted (m2)`)),
                   w=log(data$`in W corrupted (m)`),
                   s=log(data$`in S corrupted (m/m)`),
                   daIn2=(data$`in2 dA corrupted (m2)`-min(data$`in2 dA corrupted (m2)`)),
                   wIn2=log(data$`in2 W corrupted (m)`),
                   sIn2=log(data$`in2 S corrupted (m/m)`),
                   da2=(data$`out dA corrupted (m2)`-min(data$`out dA corrupted (m2)`)),
                   w2=log(data$`out W corrupted (m)`),
                   s2=log(data$`out S corrupted (m/m)`),
                   q2=rep(log(qHatOut),nrow(data)),
                   dv=data$`dV corrupted (m3/day)`/86400,
                   et=(data$`ET (m3/s)`*et),
                   lateral=(data$sum)*q)
  stan_data$da =ifelse(is.na(stan_data$da)|is.infinite(stan_data$da), 0, stan_data$da)
  stan_data$da2 =ifelse(is.na(stan_data$da2)|is.infinite(stan_data$da), 0, stan_data$da2)
  stan_data$s =ifelse(is.na(stan_data$s), 0, stan_data$s)
  stan_data$s2 =ifelse(is.na(stan_data$s2), 0, stan_data$s2)
  stan_data$lateral =ifelse(is.na(stan_data$lateral), 0, stan_data$lateral)
  stan_data$et =ifelse(is.na(stan_data$et), 0, stan_data$et)
  # stan_data$inc_none=none
  # stan_data$inc_et=et
  # stan_data$inc_lateral=q
  # stan_data$inc_et_lateral=et_q
  stan_data$nInlower = geoNInlower
  stan_data$nInupper = geoNInupper
  stan_data$aInlower=(geoAInlower)
  stan_data$aInupper =(geoAInupper)
  stan_data$nInlower2 = geoNInlower2
  stan_data$nInupper2 = geoNInupper2
  stan_data$aInlower2=(geoAInlower2)
  stan_data$aInupper2 =(geoAInupper2)
  stan_data$nOutlower = geoNInlower
  stan_data$nOutupper = geoNInupper
  stan_data$aOutlower=(geoAOutlower)
  stan_data$aOutupper = (geoAOutupper)
  stan_data$daInShift = (da_shift_fun(data$`in dA corrupted (m2)`))
  stan_data$daInShift2 = (da_shift_fun(data$`in2 dA corrupted (m2)`))
  stan_data$daOutShift = (da_shift_fun(data$`out dA corrupted (m2)`))
  stan_data$nInHat=geoNin
  stan_data$nInSd=0.25
  stan_data$aInHat=geoAin
  stan_data$aInSd=geoAinSD
  stan_data$nIn2Hat=geoNin2 
  stan_data$nIn2Sd=0.25
  stan_data$aIn2Hat=geoAin2
  stan_data$aIn2Sd=geoAinSD2
  stan_data$nOutHat=geoNout
  stan_data$nOutSd=0.25
  stan_data$aOutHat=geoAout
  stan_data$aOutSd=geoAoutSD
  fit=stan('E:\\research\\RivLake\\src\\lakeflowStan_ifelse_2inflows.stan',
           data=stan_data,
           chains=3,
           cores=5,
           iter=4000,
           control=list(stepsize=0.5,
                        adapt_delta=0.9)
  )
  
  ###Comparison with true Q.
  eqn1 = function(n, a, da, w, s){
    flow = (n^-1)*((a+da)^(5/3))*(w^(-2/3))*(s^0.5)
    return(flow)
  }
  
  mn = get_posterior_mean(fit)
  rw=row.names(mn)
  mn = data.table(mn)
  mn$type=rw
  mn = mn[,c('type','mean-all chains')]
  
  model=eqn1(exp(mn$`mean-all chains`[mn$type=='n']), (mn$`mean-all chains`[mn$type=='a']), data$`in dA corrupted (m2)`,data$`in W corrupted (m)`, data$`in S corrupted (m/m)`)
  modelIn2=eqn1(exp(mn$`mean-all chains`[mn$type=='nIn2']), (mn$`mean-all chains`[mn$type=='aIn2']), data$`in2 dA corrupted (m2)`,data$`in2 W corrupted (m)`, data$`in2 S corrupted (m/m)`)
  modelOut=eqn1(exp(mn$`mean-all chains`[mn$type=='nOut']), (mn$`mean-all chains`[mn$type=='aOut']), data$`out dA corrupted (m2)`,data$`out W corrupted (m)`, data$`out S corrupted (m/m)`)
  
  output_df=data.table(inflow=model,inflow2=modelIn2,outflow=modelOut,gage_in=data$`in Q (m3/s)`,gage_in2=data$`in2 Q (m3/s)`, gage_out=data$`out Q (m3/s)`,
                       n=exp(mn$`mean-all chains`[mn$type=='n']), a=(mn$`mean-all chains`[mn$type=='a']),
                       no=exp(mn$`mean-all chains`[mn$type=='nOut']),ao=(mn$`mean-all chains`[mn$type=='aOut']),
                       n2=exp(mn$`mean-all chains`[mn$type=='nIn2']),a2=(mn$`mean-all chains`[mn$type=='aIn2']),
                       page=page,et=et,q=q,type='corrupt',
                       geoNin=geoNin, geoAin=geoAin,geoNin2=geoNin2,geoAin2=geoAin2, geoNout=geoNout, geoAout=geoAout,
                       nInTrue=all$`in n`[1], nIn2True=all$`in2 n`[1], nOutTrue=all$`out n`[1],
                       aInTrue=all$`in A0 (m2)`[1], aIn2True=all$`in2 A0 (m2)`[1], aOutTrue=all$`out A0 (m2)`[1], 
                       date=data$Date)
}
fwrite(output_df,paste0(outpath,'lake_',page, 'corrupt_','_q_',q,'_et_',et,'.csv'))
}
lakeFlowCorrupted(2,0,0)
lakeFlowCorrupted(2,1,0)
lakeFlowCorrupted(2,0,1)
lakeFlowCorrupted(2,1,1)

lakeFlowCorrupted(3,0,0)
lakeFlowCorrupted(3,1,0)
lakeFlowCorrupted(3,0,1)
lakeFlowCorrupted(3,1,1)

lakeFlowCorrupted(4,0,0)
lakeFlowCorrupted(4,1,0)
lakeFlowCorrupted(4,0,1)
lakeFlowCorrupted(4,1,1)

############################################################################################################################
##LakeFlow with no measurement errors. 
############################################################################################################################
library(rstan)
library(data.table)
library(dplyr)
library(sf)
library(raster)
outpath = 'E:\\research\\RivLake\\stan_sensitivity\\'
lakeFlow=function(page,et,q){
  page = page
  ryan = list()
  data = readxl::read_xlsx("E:\\research\\RivLake\\updated\\Synthetic-data-11272021.xlsx", sheet=page)
  data$Date = as.Date(data$date)
  data$month=lubridate::month(data$Date)
  data$page=page
  data = left_join(data, grades,by=c('page', 'month'))
  data$sum = data$V1
  all=data
  ##Convert to 7-day sample.
  toi_df = data
  period = seq.Date(min(data$Date), max(data$Date), 7)
  tab=list()
  start=data.table(Date=period[1],dv=0,dvCorrupted=0)
  for(i in 2:length(period)){
    sub=data[data$Date<=period[i]&data$Date>period[i-1],]
    tab[[i]]=data.table(Date=period[i],dv=sum(sub$`dV (m3/day)`), dvCorrupted=sum(sub$`dV corrupted (m3/day)`))
  }
  tab=rbindlist(tab)
  tab=bind_rows(start,tab)
  data=data[data$Date%in%period,]
  tab=tab[tab$Date%in%data$Date,]
  data$`dV (m3/day)`=tab$dv[match(data$Date,tab$Date)]
  data$`dV corrupted (m3/day)`=tab$dvCorrupted[match(data$Date,tab$Date)]
  data=data[!is.na(data$`dV (m3/day)`)&data$`dV (m3/day)`!=0&!is.na(data$`dV corrupted (m3/day)`)&data$`dV corrupted (m3/day)`!=0,]
  ########################################################################################################
  ##PLD
  ########################################################################################################
  library(rgdal)
  pld <- "E:\\research\\RivLake\\SWOT_Lakes\\Harmonized_SWOT_PLD_SWORD.gdb"
  # List all feature classes in a file geodatabase
  subset(ogrDrivers(), grepl("GDB", name))
  fc_list <- ogrListLayers(pld)
  print(fc_list)
  # Read the feature class
  #reaches <- readOGR(dsn=pld,layer="Harmonized_SWORD_reaches")
  lakes <- st_read(dsn=pld,layer="Harmonized_PLD_lakes")
  shp=list.files('E:\\research\\RivLake\\SWOT_Lakes\\LakeFlowEligible\\', full.names=TRUE,pattern='.shp')
  shp = sf::st_read(shp[7])
  shp = shp[shp$names=='ALLATOONA LAKE'|shp$names=='LAKE MOHAVE;COLORADO RIVER'|shp$names=='TUTTLE CREEK LAKE',]
  shp = shp[!is.na(shp$names),]
  lakesFilt = lakes[lakes$lake_id%in%shp$lake_id,]
  lakesFilt$name=shp$names[match(lakesFilt$lake_id, shp$lake_id)]
  lakesFilt$page = c(2,4,3)
  tc1='74292200081'
  tc2='74292200171'
  
  lakesFilt$us_rch_id = ifelse(lakesFilt$page==4, tc1, lakesFilt$us_rch_id)
  lakesFilt$us_rch_id2=tc2
  
  ########################################################################################################
  ##Add in geoBAM priors. 
  ########################################################################################################
  if(page!=4){
    sos = "E:\\research\\RivLake\\src\\Confluence\\sos\\sos_netcdf\\final_sos\\constrained\\na_sword_v11_SOS_priors.nc"
    ryan = ncdf4::nc_open(sos)
    sos_outflow = RNetCDF::open.nc(sos)
    reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
    reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
    index <- which(reach_ids==lakesFilt$us_rch_id[lakesFilt$page==page], arr.ind=TRUE)
    gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
    geoAin = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
    geoNin = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
    geoNInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
    geoNInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
    geoAInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
    geoAInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
    geoAinSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
    geoNinSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
    model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
    qHatIn <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
    sigmaIn = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
    qInSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
    ##outflow
    index <- which(reach_ids==lakesFilt$ds_rch_id[lakesFilt$page==page], arr.ind=TRUE)
    gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
    geoAout = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
    geoNout = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
    geoNOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
    geoNOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
    geoAOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
    geoAOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
    geoAoutSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
    geoNoutSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
    model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
    qHatOut <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
    qHatIn = ifelse(is.na(qHatIn), qHatOut, qHatIn)
    qHatOut = ifelse(is.na(qHatOut), qHatIn, qHatOut)
    sigmaOut = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
    qOutSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
    # qHatIn = ifelse(is.na(qHatIn), mean(data$`in Q (m3/s)`), qHatIn)
    # qHatOut = ifelse(is.na(qHatOut), mean(data$`out Q (m3/s)`), qHatOut)
    ################################################################################
    da_shift_fun=function(x) median(x) - min(x)
    
    stan_data = list(N=nrow(data),
                     sigmaIn=rep(sigmaIn, nrow(data)),
                     sigmaOut=rep(sigmaOut, nrow(data)),
                     qInSd = rep(log(qHatIn*.1),nrow(data)),#qInSd
                     qOutSd= rep(log(qHatOut*.1),nrow(data)),#qOutSd
                     q=rep(log(qHatIn),nrow(data)),
                     sigma=rep(0.25,nrow(data)),
                     da=(data$`in dA (m2)`-min(data$`in dA (m2)`)),
                     w=log(data$`in W (m)`),
                     s=log(data$`in S (m/m)`),
                     da2=(data$`out dA (m2)`-min(data$`out dA (m2)`)),
                     w2=log(data$`out W (m)`),
                     s2=log(data$`out S (m/m)`),
                     q2=rep(log(qHatOut),nrow(data)),
                     dv=data$`dV (m3/day)`/86400,
                     et=(data$`ET (m3/s)`*et),
                     lateral=(data$sum*q))
    stan_data$da =ifelse(is.na(stan_data$da)|is.infinite(stan_data$da), 0, stan_data$da)
    stan_data$da2 =ifelse(is.na(stan_data$da2)|is.infinite(stan_data$da), 0, stan_data$da2)
    stan_data$s =ifelse(is.na(stan_data$s), 0, stan_data$s)
    stan_data$s2 =ifelse(is.na(stan_data$s2), 0, stan_data$s2)
    stan_data$lateral =ifelse(is.na(stan_data$lateral), 0, stan_data$lateral)
    stan_data$et =ifelse(is.na(stan_data$et), 0, stan_data$et)
    # stan_data$inc_none=none
    # stan_data$inc_et=et
    # stan_data$inc_lateral=q
    # stan_data$inc_et_lateral=et_q
    stan_data$nInlower = geoNInlower
    stan_data$nInupper = geoNInupper
    stan_data$aInlower=(geoAInlower)
    stan_data$aInupper =(geoAInupper)
    
    stan_data$nOutlower = geoNInlower
    stan_data$nOutupper = geoNInupper
    stan_data$aOutlower=(geoAOutlower)
    stan_data$aOutupper = (geoAOutupper)
    stan_data$daInShift = (da_shift_fun(data$`in dA (m2)`))
    stan_data$daOutShift = (da_shift_fun(data$`out dA (m2)`))

    stan_data$nInHat=geoNin 
    stan_data$nInSd=0.25
    stan_data$aInHat=geoAin
    stan_data$aInSd=geoAinSD
    stan_data$nOutHat=geoNout 
    stan_data$nOutSd=0.25
    stan_data$aOutHat=geoAout
    stan_data$aOutSd=geoAoutSD
    fit=stan('E:\\research\\RivLake\\src\\lakeflowStan_ifelse.stan',
             data=stan_data,
             chains=3,
             cores=5,
             iter=4000,
             control=list(stepsize=0.5,
                          adapt_delta=0.9)
    )
    
    ###Comparison with true Q.
    eqn1 = function(n, a, da, w, s){
      flow = (n^-1)*((a+da)^(5/3))*(w^(-2/3))*(s^0.5)
      return(flow)
    }
    
    mn = get_posterior_mean(fit)
    rw=row.names(mn)
    mn = data.table(mn)
    mn$type=rw
    mn = mn[,c('type','mean-all chains')]
    
    model=eqn1(exp(mn$`mean-all chains`[mn$type=='n']), (mn$`mean-all chains`[mn$type=='a']), data$`in dA (m2)`,data$`in W (m)`, data$`in S (m/m)`)
    modelOut=eqn1(exp(mn$`mean-all chains`[mn$type=='nOut']), (mn$`mean-all chains`[mn$type=='aOut']), data$`out dA (m2)`,data$`out W (m)`, data$`out S (m/m)`)
    
    output_df=data.table(inflow=model,outflow=modelOut,gage_in=data$`in Q (m3/s)`, gage_out=data$`out Q (m3/s)`,
                         n=exp(mn$`mean-all chains`[mn$type=='n']), a=(mn$`mean-all chains`[mn$type=='a']),
                         no=exp(mn$`mean-all chains`[mn$type=='nOut']),ao=(mn$`mean-all chains`[mn$type=='aOut']),
                         page=page,et=et,q=q,type='synthetic',
                         geoNin=geoNin, geoAin=geoAin,geoNout=geoNout, geoAout=geoAout,
                         nInTrue=all$`in n`[1], nOutTrue=all$`out n`[1],
                         aInTrue=all$`in A0 (m2)`[1], aOutTrue=all$`out A0 (m2)`[1],
                         date=data$Date)
    #fwrite(output_df,paste0(outpath,'lake_',page, 'synthetic_','none_',none, '_q_',q,'_et_',et,'_et_q_',et_q,'.csv'))
  }else{
    sos = "E:\\research\\RivLake\\src\\Confluence\\sos\\sos_netcdf\\final_sos\\constrained\\na_sword_v11_SOS_priors.nc"
    ryan = ncdf4::nc_open(sos)
    sos_outflow = RNetCDF::open.nc(sos)
    reach_grp <- RNetCDF::grp.inq.nc(sos_outflow, "reaches")$self
    reach_ids <- RNetCDF::var.get.nc(reach_grp, "reach_id")
    index <- which(reach_ids==lakesFilt$us_rch_id[lakesFilt$page==page], arr.ind=TRUE)
    gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
    geoAin = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
    geoNin = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
    geoNInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
    geoNInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
    geoAInlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
    geoAInupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
    geoAinSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
    geoNinSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
    model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
    qHatIn <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
    sigmaIn = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
    qInSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
    ##Inflow 2
    index <- which(reach_ids==lakesFilt$us_rch_id2[lakesFilt$page==page], arr.ind=TRUE)
    gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
    geoAin2 = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
    geoNin2 = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
    geoNInlower2 = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
    geoNInupper2 = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
    geoAInlower2 = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
    geoAInupper2 = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
    geoAinSD2 = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
    geoNinSD2 = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
    model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
    qHatIn2 <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
    sigmaIn2 = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
    qInSd2=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
    ##outflow
    index <- which(reach_ids==lakesFilt$ds_rch_id[lakesFilt$page==page], arr.ind=TRUE)
    gbp_grp <- RNetCDF::grp.inq.nc(sos_outflow, "gbpriors/reach")$self
    geoAout = RNetCDF::var.get.nc(gbp_grp, "logA0_hat")[index]
    geoNout = RNetCDF::var.get.nc(gbp_grp, "logn_hat")[index]
    geoNOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_logn")[index]
    geoNOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_logn")[index]
    geoAOutlower = RNetCDF::var.get.nc(gbp_grp, "lowerbound_A0")[index]
    geoAOutupper = RNetCDF::var.get.nc(gbp_grp, "upperbound_A0")[index]
    geoAoutSD = RNetCDF::var.get.nc(gbp_grp, "logA0_sd")[index]
    geoNoutSD = RNetCDF::var.get.nc(gbp_grp, "logn_sd")[index]
    model_grp <- RNetCDF::grp.inq.nc(sos_outflow, "model")$self
    qHatOut <- RNetCDF::var.get.nc(model_grp, "mean_q")[index]
    qHatIn = ifelse(is.na(qHatIn), qHatOut, qHatIn)
    qHatOut = ifelse(is.na(qHatOut), qHatIn, qHatOut)
    sigmaOut = RNetCDF::var.get.nc(gbp_grp, "sigma_man")[index]
    qOutSd=RNetCDF::var.get.nc(gbp_grp, "logQ_sd")[index]
    # qHatIn = ifelse(is.na(qHatIn), mean(data$`in Q (m3/s)`), qHatIn)
    # qHatOut = ifelse(is.na(qHatOut), mean(data$`out Q (m3/s)`), qHatOut)
    ################################################################################
    da_shift_fun=function(x) median(x) - min(x)
    
    stan_data = list(N=nrow(data),
                     q=rep(log(qHatIn),nrow(data)),
                     qIn2=rep(log(qHatIn2),nrow(data)),
                     sigma=rep(0.25,nrow(data)),
                     sigmaIn=rep(sigmaIn,nrow(data)),
                     sigmaIn2=rep(sigmaIn2,nrow(data)),
                     sigmaOut=rep(sigmaOut,nrow(data)),
                     qInSd = rep(log(qHatIn*.1),nrow(data)),#qInSd
                     qInSd2 = rep(log(qHatIn2*.1),nrow(data)),#qInSd
                     qOutSd= rep(log(qHatOut*.1),nrow(data)),#qOutSd
                     da=(data$`in dA (m2)`-min(data$`in dA (m2)`)),
                     w=log(data$`in W (m)`),
                     s=log(data$`in S (m/m)`),
                     daIn2=(data$`in2 dA (m2)`-min(data$`in2 dA (m2)`)),
                     wIn2=log(data$`in2 W (m)`),
                     sIn2=log(data$`in2 S (m/m)`),
                     da2=(data$`out dA (m2)`-min(data$`out dA (m2)`)),
                     w2=log(data$`out W (m)`),
                     s2=log(data$`out S (m/m)`),
                     q2=rep(log(qHatOut),nrow(data)),
                     dv=data$`dV (m3/day)`/86400,
                     et=(data$`ET (m3/s)`*et),
                     lateral=(data$sum*q))
    stan_data$da =ifelse(is.na(stan_data$da)|is.infinite(stan_data$da), 0, stan_data$da)
    stan_data$da2 =ifelse(is.na(stan_data$da2)|is.infinite(stan_data$da), 0, stan_data$da2)
    stan_data$s =ifelse(is.na(stan_data$s), 0, stan_data$s)
    stan_data$s2 =ifelse(is.na(stan_data$s2), 0, stan_data$s2)
    stan_data$lateral =ifelse(is.na(stan_data$lateral), 0, stan_data$lateral)
    stan_data$et =ifelse(is.na(stan_data$et), 0, stan_data$et)
    # stan_data$inc_none=none
    # stan_data$inc_et=et
    # stan_data$inc_lateral=q
    # stan_data$inc_et_lateral=et_q
    stan_data$nInlower = geoNInlower
    stan_data$nInupper = geoNInupper
    stan_data$aInlower=(geoAInlower)
    stan_data$aInupper =(geoAInupper)
    stan_data$nInlower2 = geoNInlower2
    stan_data$nInupper2 = geoNInupper2
    stan_data$aInlower2=(geoAInlower2)
    stan_data$aInupper2 =(geoAInupper2)
    stan_data$nOutlower = geoNInlower
    stan_data$nOutupper = geoNInupper
    stan_data$aOutlower=(geoAOutlower)
    stan_data$aOutupper = (geoAOutupper)
    stan_data$daInShift = (da_shift_fun(data$`in dA (m2)`))
    stan_data$daInShift2 = (da_shift_fun(data$`in2 dA (m2)`))
    stan_data$daOutShift = (da_shift_fun(data$`out dA (m2)`))
    
    stan_data$nInHat=geoNin
    stan_data$nInSd=0.25
    stan_data$aInHat=geoAin
    stan_data$aInSd=geoAinSD
    stan_data$nIn2Hat=geoNin2 
    stan_data$nIn2Sd=0.25
    stan_data$aIn2Hat=geoAin2
    stan_data$aIn2Sd=geoAinSD2
    stan_data$nOutHat=geoNout
    stan_data$nOutSd=0.25
    stan_data$aOutHat=geoAout
    stan_data$aOutSd=geoAoutSD
    fit=stan('E:\\research\\RivLake\\src\\lakeflowStan_ifelse_2inflows.stan',
             data=stan_data,
             chains=3,
             cores=5,
             iter=4000,
             control=list(stepsize=0.5,
                          adapt_delta=0.9)
    )
    
    ###Comparison with true Q.
    eqn1 = function(n, a, da, w, s){
      flow = (n^-1)*((a+da)^(5/3))*(w^(-2/3))*(s^0.5)
      return(flow)
    }
    
    mn = get_posterior_mean(fit)
    rw=row.names(mn)
    mn = data.table(mn)
    mn$type=rw
    mn = mn[,c('type','mean-all chains')]
    
    model=eqn1(exp(mn$`mean-all chains`[mn$type=='n']), (mn$`mean-all chains`[mn$type=='a']), data$`in dA (m2)`,data$`in W (m)`, data$`in S (m/m)`)
    modelIn2=eqn1(exp(mn$`mean-all chains`[mn$type=='nIn2']), (mn$`mean-all chains`[mn$type=='aIn2']), data$`in2 dA (m2)`,data$`in2 W (m)`, data$`in2 S (m/m)`)
    modelOut=eqn1(exp(mn$`mean-all chains`[mn$type=='nOut']), (mn$`mean-all chains`[mn$type=='aOut']), data$`out dA (m2)`,data$`out W (m)`, data$`out S (m/m)`)
    
    output_df=data.table(inflow=model,inflow2=modelIn2,outflow=modelOut,gage_in=data$`in Q (m3/s)`,gage_in2=data$`in2 Q (m3/s)`, gage_out=data$`out Q (m3/s)`,
                         n=exp(mn$`mean-all chains`[mn$type=='n']), a=(mn$`mean-all chains`[mn$type=='a']),
                         no=exp(mn$`mean-all chains`[mn$type=='nOut']),ao=(mn$`mean-all chains`[mn$type=='aOut']),
                         n2=exp(mn$`mean-all chains`[mn$type=='nIn2']),a2=(mn$`mean-all chains`[mn$type=='aIn2']),
                         page=page,et=et,q=q,type='synthetic',
                         geoNin=geoNin, geoAin=geoAin,geoNin2=geoNin2,geoAin2=geoAin2, geoNout=geoNout, geoAout=geoAout,
                         nInTrue=all$`in n`[1], nIn2True=all$`in2 n`[1], nOutTrue=all$`out n`[1],
                         aInTrue=all$`in A0 (m2)`[1], aIn2True=all$`in2 A0 (m2)`[1], aOutTrue=all$`out A0 (m2)`[1],
                         date=data$Date)
  }
  fwrite(output_df,paste0(outpath,'lake_',page, 'synthetic_','_q_',q,'_et_',et,'.csv'))
}
lakeFlow(2,0,0)
lakeFlow(2,1,0)
lakeFlow(2,0,1)
lakeFlow(2,1,1)

lakeFlow(3,0,0)
lakeFlow(3,1,0)
lakeFlow(3,0,1)
lakeFlow(3,1,1)

lakeFlow(4,0,0)
lakeFlow(4,1,0)
lakeFlow(4,0,1)
lakeFlow(4,1,1)



library(data.table)
library(dplyr)
library(ggplot2)
library(hydroGOF)
source("E:/research/2019_08_30_rivObs/git/src/Error_stats_functions.R")

################################################################################
##Plots for LakeFlow manuscript. 
################################################################################
files=list.files('E:\\research\\RivLake\\stan_sensitivity\\', full.names=TRUE)
#files = files[-grep('daily', files)]
comb=lapply(files, fread)
jnd = rbindlist(comb, fill=TRUE)
jnd$lake=jnd$page
jnd$model = NA
jnd$model=ifelse(jnd$et==1, 'E', jnd$model)
jnd$model=ifelse(jnd$q==1, 'Q', jnd$model)
jnd$model=ifelse(jnd$q==1&jnd$et==1, 'EQ', jnd$model)
jnd$model=ifelse(is.na(jnd$model), 'none', jnd$model)
jnd$type=ifelse(jnd$type=='corrupt', 'Corrupted', jnd$type)
jnd$type=ifelse(jnd$type=='synthetic', 'Synthetic', jnd$type)


jnd$lake=ifelse(jnd$lake==3, 'Mohave', jnd$lake)
jnd$lake=ifelse(jnd$lake==2, 'Allatoona', jnd$lake)
jnd$lake=ifelse(jnd$lake==4, 'Tuttle Creek', jnd$lake)
secondInflow = jnd[!is.na(jnd$inflow2),]
secondInflow$lake = 'Tuttle Creek Inflow 2'
secondInflow$gage_in = secondInflow$gage_in2
secondInflow$inflow=secondInflow$inflow2
jnd = bind_rows(jnd, secondInflow)

jnd$type=factor(jnd$type, levels=c('Synthetic', 'Corrupted'))
jnd$lake=factor(jnd$lake, levels=c('Mohave', 'Tuttle Creek','Tuttle Creek Inflow 2', 'Allatoona'))
jnd$model=factor(jnd$model, levels=c('none', 'E', 'Q', 'EQ'))


inflow=jnd[,c('inflow','gage_in', 'type','lake','model', 'date')]
inflowTall=melt(inflow, measurement.vars=c('inflow'),id.vars=c('gage_in', 'type','lake','model', 'date'))
inflowTall$lakeFlow=paste0(inflowTall$lake, ' inflow')
inflowTall$gage=inflowTall$gage_in
inflowTall$lf=inflowTall$value
outflow=jnd[,c('outflow', 'gage_out', 'type','lake', 'model','date')]
outflowTall=melt(outflow, measurement.vars=c('outflow'),id.vars=c('gage_out', 'type','lake','model','date'))
outflowTall$lakeFlow=paste0(outflowTall$lake, ' outflow')
outflowTall$gage=outflowTall$gage_out
outflowTall$lf=outflowTall$value
combined=bind_rows(inflowTall,outflowTall)
cunique=combined[!duplicated(combined[,c('lakeFlow','model','type','lake', 'gage', 'date')]),]
cunique=cunique[cunique$lakeFlow!='Tuttle Creek Inflow 2 outflow',]
cunique$lakeFlow=factor(cunique$lakeFlow,levels=c('Allatoona inflow','Allatoona outflow',
                                          'Mohave inflow','Mohave outflow',
                                          'Tuttle Creek inflow', 'Tuttle Creek outflow','Tuttle Creek Inflow 2 inflow'))

plot(cunique[cunique$lakeFlow=='Tuttle Creek inflow']$lf, 
     cunique[cunique$lakeFlow=='Tuttle Creek Inflow 2 inflow']$lf)
abline(0,1)
##############################################################################################################
##Q scatter plot and performance metrics. 
##############################################################################################################
modelPlot=ggplot(data=cunique)+
  # geom_rect(data = background,aes(fill = model),alpha=0.25,xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf)+
  # scale_fill_manual(values=c('white','#C48672','#5694DD','#BD69B8'))+
  scale_y_log10()+
  scale_x_log10()+
  coord_cartesian(xlim=c(1,1000),ylim=c(1,1000))+
  geom_abline(slope=1,intercept = 0,lty=2,size=0.5)+
  geom_point(aes(x=gage,y=lf,color=lakeFlow),size=0.5,alpha=0.2)+
  #scale_shape_manual(values=c(1, 2,3))+
  #geom_text(data=cunique[,rmse(lf,gage),by=list(model,type)],aes(x=7, y=900,label=paste0('RMSE=',sprintf("%.2f", round(V1,2)))),size=3)+
  geom_text(data=cunique[,mean(abs(lf-gage),na.rm=TRUE),by=list(model,type)],aes(x=7, y=900,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),size=3)+
  geom_text(data=cunique[,mean((lf-gage),na.rm=TRUE),by=list(model,type)],aes(x=7, y=500,label=paste0('Bias=',sprintf("%.2f", round(V1,2)))),size=3)+
  geom_text(data=cunique[,sqrt(mean((lf-gage)^2)),by=list(model,type)],aes(x=7, y=200,label=paste0('RMSE=',sprintf("%.2f", round(V1,2)))),size=3)+
  facet_wrap(~type+model, ncol=4,labeller = label_wrap_gen(multi_line=FALSE))+
  scale_color_brewer(palette='Dark2')+
  #scale_color_manual(values=c(dk2[4],dk2[5],dk2[6]))+
  #scale_color_manual(values=c('purple','yellow','orange'))+
  ylab('LakeFlow Discharge (cms)')+
  xlab('Gauge Discharge (cms)')+
  theme_bw()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        axis.title=element_text(size=12, color="black"),
        axis.text=element_text(size=10, color="black"),
        legend.text=element_text(size=9, color="black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=12, color="black"),
        aspect.ratio = 1,panel.grid = element_line(size=0.25),
        panel.spacing = unit(0.5, "lines"),plot.margin = margin(5.5,7.25,5.5,5.5))
modelPlot
# ggsave("E:\\research\\RivLake\\Figures\\modelScatterPlotWide.pdf",modelPlot,
#        dpi=1000,units="in",width=6.5,height=8)

bar=cunique
bar=bar[!is.na(lf),list(nse=NSE(lf,gage)*100,rbias=rBias(lf,gage),nrmse=NRMSE(lf,gage)),by=list(model,type,lakeFlow)]
bar=melt(bar,measurement.vars=c('nse','rbias','nrmse'))
bar[,median(value),by=list(model,type,variable)]

bar[variable=='nse',list(quantile(value,0.25)),by=list(model,type)]


source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

barPlot=ggplot(data=bar)+
  # geom_rect(data = background,aes(fill = model),alpha=0.25,xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf)+
  # scale_fill_manual(values=c('white','#C48672','#5694DD','#BD69B8'))+  
  #stat_ecdf(aes(y=value, color=variable))+
  geom_flat_violin(aes(x=variable,y=value), fill='grey90',color='grey90')+
  geom_boxplot(aes(x=variable,y=value),fill='lightblue',outlier.shape = NA,width=0.3,lwd=0.25)+
  coord_cartesian(ylim=c(-50,125))+
  # geom_point(aes(x=gage,y=lf,color=lakeFlow,alpha=0.9),size=0.25)+
  # geom_text(data=c[,rmse(lf,gage),by=list(model,type)],aes(x=7, y=900,label=paste0('RMSE=',sprintf("%.2f", round(V1,2)))),size=3)+
  # geom_text(data=c[,mean(abs(lf-gage),na.rm=TRUE),by=list(model,type)],aes(x=7, y=300,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),size=3)+
  facet_wrap(~type+model, ncol=4,labeller = label_wrap_gen(multi_line=FALSE))+
  #scale_fill_brewer(palette='Dark2')+
  ylab('Percent')+
  xlab('')+
  theme_bw()+
  theme(legend.position = 'none',
        legend.title = element_blank(),
        axis.title=element_text(size=12, color="black"),
        axis.text=element_text(size=10, color="black"),
        legend.text=element_text(size=9, color="black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=12, color="black"),
        aspect.ratio = 1,panel.grid = element_line(size=0.25),
        panel.spacing = unit(0.5, "lines"),plot.margin = margin(5.5,7.25,5.5,5.5))
barPlot
library(ggpubr)
mdlBarPlot=ggarrange(modelPlot,barPlot,nrow=2,common.legend = TRUE,labels='auto') 
mdlBarPlot
ggsave("E:\\research\\RivLake\\Figures\\mdlBarPlotWide.pdf",mdlBarPlot,
       dpi=1000,units="in",width=6.5,height=8)
####################################################################################################################
##Hydrographs. 
####################################################################################################################
library(scales)
allatoona = readxl::read_xlsx("E:\\research\\RivLake\\updated\\Synthetic-data-11272021.xlsx", sheet=2)
allInf = data.frame(Date=allatoona$date,gage= allatoona$`in Q (m3/s)`,lakeFlow='Allatoona inflow')
allOut = data.frame(Date=allatoona$date,gage= allatoona$`out Q (m3/s)`,lakeFlow='Allatoona outflow')
mohave = readxl::read_xlsx("E:\\research\\RivLake\\updated\\Synthetic-data-11272021.xlsx", sheet=3)
moInf = data.frame(Date=mohave$date,gage= mohave$`in Q (m3/s)`,lakeFlow='Mohave inflow')
moOut = data.frame(Date=mohave$date,gage= mohave$`out Q (m3/s)`,lakeFlow='Mohave outflow')
tc=readxl::read_xlsx("E:\\research\\RivLake\\updated\\Synthetic-data-11272021.xlsx", sheet=4)
tcInf = data.frame(Date=tc$date,gage= tc$`in Q (m3/s)`,lakeFlow='Tuttle Creek inflow')
tcOut = data.frame(Date=tc$date,gage= tc$`out Q (m3/s)`,lakeFlow='Tuttle Creek outflow')
tcInf2 = data.frame(Date=tc$date,gage= tc$`in2 Q (m3/s)`,lakeFlow='Tuttle Creek Inflow 2 inflow')
gageData=bind_rows(allInf, allOut,moInf, moOut,tcInf,tcOut,tcInf2)
gageData$Date = as.Date(gageData$Date)


hyPlot2=ggplot(gageData)+#c[!is.na(gage)][model=='none'])+
  facet_wrap(~factor(lakeFlow,levels=c('Allatoona inflow','Allatoona outflow',
                                       'Mohave inflow','Mohave outflow',
                                       'Tuttle Creek inflow', 'Tuttle Creek outflow','Tuttle Creek Inflow 2 inflow'
  )), scales='free',ncol=2,labeller = label_wrap_gen(multi_line=FALSE))+
  geom_line(aes(x=Date,y=gage),size=0.5,col='grey50')+
  geom_line(data=cunique[!is.na(lf)&type=='Corrupted'][model=='none'|model=='EQ'],
            aes(x=date,y=lf,color=model),alpha=0.5,size=0.5)+
  #geom_point(data=c[c$model=='EQ'&!is.na(lf)&type=='corrupted',],aes(x=Date, lf),col='blue4',size=0.25)+
  scale_color_manual(values=c('blue4','red3'))+
  scale_x_date(breaks = scales::pretty_breaks(n = 3)) +
  #geom_text(data=hyText,aes(x=V3, y=Inf,label=paste0('NSE=',sprintf("%.2f", round(NSE,2)), ', ','rBias(%)=',sprintf("%.0f", round(rBias,2)),', ','NRMSE(%)=',sprintf("%.0f", round(nrmse,2)))),size=3,hjust=0,vjust=1,nudge_x = (365*0.5),nudge_y=-500)+
  ylab('Discharge (cms)')+
  xlab(NULL)+
  #geom_text(data=jnd[jnd$type=='synthetic',],aes(x=median(jnd$Date),max(jnd$`out Q (m3/s)`), label=paste0('NSE=',round(nseOut,2))))+
  theme_bw()+
  theme(legend.position = 'top',
        legend.title=element_blank(),
        axis.title=element_text(size=12, color="black"),
        axis.text=element_text(size=10, color="black"),
        legend.text=element_text(size=11, color="black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=12, color="black"))
hyPlot2
ggsave("E:\\research\\RivLake\\Figures\\normHydro.pdf",hyPlot2,
       dpi=1000,units="in",width=6.5,height=8)
###########################################################################################################
##FLP plots.
###########################################################################################################
names=c('n','a','nTrue','aTrue','geoN','geoA','type','lake','model')
inflow=jnd[,c('n','a','nInTrue','aInTrue','geoNin','geoAin','type','lake','model')]
colnames(inflow)=names
inflowTall=melt(inflow, measurement.vars=c('n','a'),id.vars=c('nTrue','aTrue','geoN','geoA','type','lake','model'))
inflowTall$lakeFlow=paste0(inflowTall$lake, ' inflow')
inflowTall = inflowTall[inflowTall$lake!='Tuttle Creek Inflow 2']


inflow2=jnd[,c('n2','a2','nIn2True','aIn2True','geoNin2','geoAin2','type','lake','model')]
colnames(inflow2)=names
inflow2=inflow2[complete.cases(inflow2),]
inflow2 = inflow2[inflow2$lake=='Tuttle Creek Inflow 2',]
inflow2Tall=melt(inflow2, measurement.vars=c('n','a'),id.vars=c('nTrue','aTrue','geoN','geoA','type','lake','model'))
inflow2Tall$lakeFlow=paste0(inflow2Tall$lake, ' inflow')


outflow=jnd[,c('no','ao','nOutTrue','aOutTrue','geoNout','geoAout','type','lake','model')]
colnames(outflow)=names
outflowTall=melt(outflow, measurement.vars=c(c('n','a')),id.vars=c('nTrue','aTrue','geoN','geoA','type','lake','model'))
outflowTall$lakeFlow=paste0(outflowTall$lake, ' outflow')
combined=bind_rows(inflowTall,inflow2Tall,outflowTall)
combined$gage = combined$nTrue
combined$gage = ifelse(combined$variable=='a', combined$aTrue, combined$gage)
combined$geo=combined$geoN
combined$geo=ifelse(combined$variable=='a', combined$geoA, combined$geo)
combined$geo = exp(combined$geo)

cunique=combined[!duplicated(combined[,c('lakeFlow','model','type','variable')]),]
cunique=cunique[cunique$lakeFlow!='Tuttle Creek Inflow 2 outflow',]
cunique$lakeFlow=factor(cunique$lakeFlow,levels=c('Allatoona inflow','Allatoona outflow',
                                                  'Mohave inflow','Mohave outflow',
                                                  'Tuttle Creek inflow', 'Tuttle Creek outflow','Tuttle Creek Inflow 2 inflow'))


aPlot=ggplot(data=cunique[cunique$variable=='a'])+
  coord_cartesian(xlim=c(0.3,300),ylim=c(0.3,300))+
  scale_y_log10()+
  scale_x_log10()+
  geom_abline(slope=1,intercept = 0,lty=2,size=0.5)+
  geom_point(aes(x=gage,y=value,color=lakeFlow),size=1,alpha=0.5)+
  scale_shape_manual(values=c(1, 2,3))+
  geom_text(data=cunique[variable=='a',mean(abs(value-gage),na.rm=TRUE),by=list(model,type)],aes(x=3, y=300,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),size=3)+
  geom_text(data=cunique[variable=='a',mean(abs(geo-gage),na.rm=TRUE),by=list(model,type)],aes(x=3, y=150,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),col='red',size=3)+
  facet_wrap(~type+model, ncol=4,labeller = label_wrap_gen(multi_line=FALSE))+
  scale_color_brewer(palette='Dark2')+
  ylab('LakeFlow Bathymetry (m2)')+
  xlab('True Bathymetry (m2)')+
  theme_bw()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        axis.title=element_text(size=12, color="black"),
        axis.text=element_text(size=10, color="black"),
        legend.text=element_text(size=9, color="black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=12, color="black"),
        aspect.ratio = 1,panel.grid = element_line(size=0.25),
        panel.spacing = unit(0.5, "lines"),plot.margin = margin(5.5,7.25,5.5,5.5))
aPlot

ndf = cunique[cunique$variable=='n']
ndf$gage=(ndf$gage)
ndf$geo=(ndf$geo)
ndf$value=(ndf$value)
nPlot=ggplot(data=ndf)+
  # geom_rect(data = background,aes(fill = model),alpha=0.25,xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf)+
  # scale_fill_manual(values=c('white','#C48672','#5694DD','#BD69B8'))+
  geom_density(aes(value, ..scaled..),
               fill = 'grey20',
               alpha = 0.25)+
  scale_x_log10()+
  geom_vline(data=ndf,aes(xintercept=gage),lty=1,size=0.5,col='blue')+
  geom_vline(data=ndf,aes(xintercept=geo),lty=2,size=0.5,col='red')+
  # geom_text(data=c[,rmse(log(flp$n_in),log(flp$n_true)),by=list(model,type)],aes(x=0.02, y=4,label=paste0('RMSE=',sprintf("%.2f", round(V1,2)))),size=3)+
  # geom_text(data=c[,mean(abs(flp$n_in-flp$n_true),na.rm=TRUE),by=list(model,type)],aes(x=0.02, y=3,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),size=3)+
  # geom_text(data=c[,mean(abs(flp$n_in-flp$n_geoBam),na.rm=TRUE),by=list(model,type)],aes(x=0.02, y=3,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),size=3)+
  facet_wrap(~type+model, ncol=4,labeller = label_wrap_gen(multi_line=FALSE))+
  scale_color_brewer(palette='Dark2')+
  ylab('Density')+
  xlab('LakeFlow Mannings n')+
  theme_bw()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        axis.title=element_text(size=12, color="black"),
        axis.text=element_text(size=10, color="black"),
        legend.text=element_text(size=9, color="black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=12, color="black"),
        aspect.ratio = 1,panel.grid = element_line(size=0.25),
        panel.spacing = unit(0.5, "lines"),plot.margin = margin(5.5,7.25,5.5,5.5))
nPlot

cPlot=ggarrange(aPlot,nPlot,nrow=2,common.legend = TRUE,labels='auto') 
cPlot
ggsave("E:\\research\\RivLake\\Figures\\combPlotWide.png",cPlot,
       dpi=1000,units="in",width=6.5,height=8)






ndf = cunique[cunique$variable=='n']
ndf$gage=log(ndf$gage)
ndf$geo=log(ndf$geo)
ndf$value=log(ndf$value)
nPlot=ggplot(data=ndf)+
  coord_cartesian(xlim=c(-3.7,-3),ylim=c(-3.7,-3))+
  #scale_y_log10()+
  #scale_x_log10()+
  geom_abline(slope=1,intercept = 0,lty=2,size=0.5)+
  geom_point(aes(x=gage,y=value,color=lakeFlow),size=1,alpha=0.5)+
  scale_shape_manual(values=c(1, 2,3))+
  geom_text(data=ndf[,mean(abs(value-gage),na.rm=TRUE),by=list(model,type)],aes(x=-3.5, y=-3.2,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),size=3)+
  geom_text(data=ndf[,mean(abs(geo-gage),na.rm=TRUE),by=list(model,type)],aes(x=-3.5, y=-3.3,label=paste0('MAE=',sprintf("%.2f", round(V1,2)))),col='red',size=3)+
  facet_wrap(~type+model, ncol=4,labeller = label_wrap_gen(multi_line=FALSE))+
  scale_color_brewer(palette='Dark2')+
  ylab('log n')+
  xlab('log n')+
  theme_bw()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        axis.title=element_text(size=12, color="black"),
        axis.text=element_text(size=10, color="black"),
        legend.text=element_text(size=9, color="black"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=12, color="black"),
        aspect.ratio = 1,panel.grid = element_line(size=0.25),
        panel.spacing = unit(0.5, "lines"),plot.margin = margin(5.5,7.25,5.5,5.5))
nPlot
cPlot=ggarrange(aPlot,nPlot,nrow=2,common.legend = TRUE,labels='auto') 
cPlot
ggsave("E:\\research\\RivLake\\Figures\\combPlotWide.pdf",cPlot,
       dpi=1000,units="in",width=6.5,height=8)


##MAE comparison. 
ryan=cunique[,mean(abs(value-gage),na.rm=TRUE),by=list(model,type, variable)]
ryan = ryan[ryan$variable=='a']
ryan$geo=222.23
pc=function(old,new){
  ((new-old)/old)*100
}
ryan$pc=pc(ryan$geo, ryan$V1)
mean(ryan$pc)



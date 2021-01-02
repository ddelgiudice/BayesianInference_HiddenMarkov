######################################
# Cal & Unc Ana for iid, bias descr, multipliers with 2 rain scen; 3 ev for calibration
# Dario Del Giudice - Eawag - ETHZ - Carnegie Institution
# Cite as: Del Giudice et al. (2016), Describing the catchment-averaged precipitation as a stochastic process improves parameter and input estimation, WRR
# last modified: Feb 2016
######################################

rm(list=ls(all=TRUE)); graphics.off()

path.comp<- 'C:/Users/Dario/Desktop/CodeData_SIP' # change path
setwd(path.comp)
hodie   <- Sys.Date(); hodie <- format(hodie, format="%Y%b%d")
print.pdf =F

####################################
# Select your rainfall scenario and error model

RaSc <- "malaP"  # Sc2
# RaSc <- "bonaP" # Sc1

# .............
# EM   <- "iid" 
# EM   <- "Bias" 
# EM   <- "mult"
EM   <- "SIP"

# choose output transformation
transf  <- "ls" 

# name output variable
var      = "Q" # our output is flow 

# selected rain events 
cal.ev <- c(2,6,9)

# define iterations in calibration 
iter.inf   = 50 # just to get started; it should be then increased.

# define fix parameters
if(EM == "iid"||EM == "mult") {par.fix   <- c(Del.Max=0,sd.B_Q=0, ks_Q=0, corrlen=0, Delta=0)
}else           {par.fix   <- c(Del.Max=20, ks_Q=0, Delta=0)} # considering bias previous parameters are calibrated


#####################################
# load packages and define functions
source("Functions.R")   # functions to run calibration
if (!require("deSolve")) install.packages("deSolve") ; library(deSolve) #ODE solver
if (!require("kimisc")) install.packages("kimisc") ;library(kimisc) #convert time with hms.to.seconds(x)
if (!require("adaptMCMC")) install.packages("adaptMCMC") ;library("adaptMCMC") # adaptive MCMC sampling
if (!require("GenSA")) install.packages("GenSA") ;library("GenSA") # http://cas.stat.ucla.edu/page_attachments/0000/0079/Mullen_nov27_2012.pdf
if (!require("tmvtnorm")) install.packages("tmvtnorm") ;library(tmvtnorm) #truncated multivariate normal distribution
source("http://www.phaget4.org/R/myImagePlot.R") 



#####################################
# import precipitation (from pluvio 2,3, and from STP) and discharge data
load("DataParam.RData")

# manipulate data to put them in the right format
inp.7  = inp.7[which(inp.7[,1]==71),] # pluviometer from MS in Fluntern (7km distant)
t.7.or = paste(paste(inp.7[,2],inp.7[,3],inp.7[,4], sep="-"),paste(inp.7[,5],inp.7[,6], sep=":"), sep=" ")
t.7.or = as.POSIXlt(t.7.or, tz="UTC")
inp.7  = data.frame(date.time=t.7.or,rain=inp.7[,ncol(inp.7)])

inp.data[,1]<-as.POSIXct(inp.data[,1],tz="UTC")
inp.2[,1]   <-as.POSIXct(inp.2[,1],tz="UTC")   # Resolution (0.001 mm/min) 
inp.3[,1]   <-as.POSIXct(inp.3[,1],tz="UTC") # Resolution (0.025 mm/min) 
out.data[,1]<-as.POSIXct(out.data[,1],format="%d.%m.%Y %H:%M",tz="UTC")


##############################################################################
# This part rearranges the data and defines time vectors for the 3 calibration events
##############################################################################
## Period I - Calibration

# select event 
ev.nr       = cal.ev[1] 
dt.mod      = 4/60 # same as flow meas [hr]
ind.per.cal = seq(from=which(inp.data[,1]==selected_events1[ev.nr,2]-3600),to=which(inp.data[,1]==selected_events1[ev.nr,3]+3600*2))

da.ti.cal   = inp.data[ind.per.cal,1] # UTC time of the data (every min)
out.tabl    = out.data[match(da.ti.cal,out.data[,1]),] # l/s (with 1 min res)
P3.tabl     = inp.3[match(da.ti.cal,inp.3[,1]),] # 1 min data from P3 (with same index as P1)
P2.tabl.1   = inp.2[match(da.ti.cal,inp.2[,1]),]; P2.tabl.1[,2] = apply(cbind(P3.tabl[,2],P2.tabl.1[,2]), 1, mean)
P7.tabl     = inp.7.2m[match(da.ti.cal,inp.7.2m[,1]),]
P7.tabl.raw = inp.7[match(da.ti.cal,inp.7[,1]),]
P7.tabl.raw = P7.tabl.raw[complete.cases(P7.tabl.raw),] 
inp.tabl    = inp.data[ind.per.cal,]; inp.tabl[,1]=out.tabl[,1]; inp.tabl[,2]=NA
inp.tabl.FL = inp.tabl

# calculate the mean rain rate over last 4 mins
for (i in 1:nrow(out.tabl))
{
  if(!is.na(out.tabl[i,1])) 
  {
    inp.tabl[i,2]     <- mean(c(P2.tabl.1[i,2],P2.tabl.1[max((i-1),1),2],
                                P2.tabl.1[min(i+1,nrow(out.tabl)),2],P2.tabl.1[min(i+2,nrow(out.tabl)),2]))
    inp.tabl.FL[i,2]  <- mean(c(P7.tabl[i,2],P7.tabl[max((i-1),1),2],
                                P7.tabl[min(i+1,nrow(out.tabl)),2],P7.tabl[min(i+2,nrow(out.tabl)),2]))
  }  
}

out.ta.th   = out.tabl[complete.cases(out.tabl),] # selected out table without Nas

if(RaSc == "bonaP") {inp.ta.th   = inp.tabl[complete.cases(inp.tabl),]  # eliminate NA-rows within the input
                     inp.1min.1  = P2.tabl.1[,2]
                     inp.meas.1  = P2.tabl.1[,2]
                     pluv.tabl   = P2.tabl.1
                     dt.inp      = as.numeric(P2.tabl.1[2,1]-P2.tabl.1[1,1])/60
}else               {inp.ta.th   = inp.tabl.FL[complete.cases(inp.tabl.FL),] 
                     inp.1min.1  = P7.tabl[,2]
                     inp.meas.1  = P7.tabl.raw[,2]/10 # important! these data are not mm/min but mm/10min
                     pluv.tabl   = P7.tabl.raw
                     dt.inp      = 10/60} # change when I will start using these data

out.meas.1    = out.ta.th[,11]
d.t.cal.th    = out.ta.th[,1] # now thinned at 4 min
inp4mod.1     <- inp.ta.th[,2]*60/1000 # mm/min to m/hr
imber.vel     <- inp.ta.th[,2]/60      # mm/min to mm/s
imber_i        = c(rep(0,par.fix["Del.Max"]),imber.vel) # shifted to account for delay betw inp-out


# set time & create a calibration layout
# ---------------------------------------
origo     = hms.to.seconds(format(d.t.cal.th[1], format="%H:%M:%S",tz=""))/3600#hr
ult.cal   = (length(d.t.cal.th)-1)*dt.mod+origo #hr -- end of event 1
t.grid.Ome.1= format(seq(origo,ult.cal,dt.mod), nsmall=6)


# for input labelling with time
quom.prin.ra <- hms.to.seconds(format(pluv.tabl[1,1], format="%H:%M:%S",tz=""))/3600
ult.grid.1.ra  <- (nrow(pluv.tabl)-1)*dt.inp+quom.prin.ra #hr 
t.grid.Ime.1   = format(seq(quom.prin.ra, ult.grid.1.ra,dt.inp), nsmall=6) #vector with elapsed hr from simulation beginning


# for input labelling with finest time (1min)
t.grid.1mi.1= format(seq(hms.to.seconds(format(P2.tabl.1[1,1], format="%H:%M:%S",tz=""))/3600,
                                    hms.to.seconds(format(P2.tabl.1[nrow(P2.tabl.1),1], format="%H:%M:%S",tz=""))/3600,
                                    1/60),nsmall=6) 


# check if output measurements come first than rain; if yes, remove the first value
if(as.numeric(t.grid.Ome.1[1])<as.numeric(t.grid.Ime.1[1]))
{
  out.meas.1 = out.meas.1[2:length(t.grid.Ome.1)] # eliminate 1st output measurement coming before rains begins
  inp4mod.1  = inp4mod.1[2:length(t.grid.Ome.1)] # reduce the length of the "artificially" aggregated rain
  t.grid.Ome.1     = t.grid.Ome.1[2:length(t.grid.Ome.1)]
  out.meas.1.br = out.meas.1
  inp4mod.1.br  = inp4mod.1 # reduce the length of the "artificially" aggregated rain
  t.grid.Ome.1.br     = t.grid.Ome.1
}

L.obs.1   = paste("Q",t.grid.Ome.1,sep="_") 
ind.comp.xi.1 = match(t.grid.Ime.1,t.grid.1mi.1,5)
names(out.meas.1)  = L.obs.1
names(inp.1min.1)  = t.grid.1mi.1 # meas rain [mm/min] every min [for out lkh test]
names(inp.meas.1)  = t.grid.Ime.1 # meas rain [mm/min] every min or 10 (dep on pluv) [for inference]
names(inp4mod.1)   = t.grid.Ome.1 # aggregated rain [mm/min] every 4min [for mod test]

##############################################################################
# Period II - Calibration
ev.nr   = cal.ev[2] 
ind.per.cal = seq(from=which(inp.data[,1]==selected_events1[ev.nr,2]-3600),to=which(inp.data[,1]==selected_events1[ev.nr,3]+3600*2))

da.ti.cal   = inp.data[ind.per.cal,1] # UTC time of the data (every min)
out.tabl    = out.data[match(da.ti.cal,out.data[,1]),] # l/s (with 1 min res)
P3.tabl     = inp.3[match(da.ti.cal,inp.3[,1]),] # 1 min data from P3 (with same index as P1)
P2.tabl.2   = inp.2[match(da.ti.cal,inp.2[,1]),]; P2.tabl.2[,2] = apply(cbind(P3.tabl[,2],P2.tabl.2[,2]), 1, mean)
P7.tabl     = inp.7.2m[match(da.ti.cal,inp.7.2m[,1]),]
P7.tabl.raw = inp.7[match(da.ti.cal,inp.7[,1]),]
P7.tabl.raw = P7.tabl.raw[complete.cases(P7.tabl.raw),] 
inp.tabl    = inp.data[ind.per.cal,]; inp.tabl[,1]=out.tabl[,1]; inp.tabl[,2]=NA
inp.tabl.FL = inp.tabl

# calculate the mean rain rate over last 4 mins
for (i in 1:nrow(out.tabl))
{
  if(!is.na(out.tabl[i,1])) 
  {
    inp.tabl[i,2]     <- mean(c(P2.tabl.2[i,2],P2.tabl.2[max((i-1),1),2],
                                P2.tabl.2[min(i+1,nrow(out.tabl)),2],P2.tabl.2[min(i+2,nrow(out.tabl)),2]))
    inp.tabl.FL[i,2]  <- mean(c(P7.tabl[i,2],P7.tabl[max((i-1),1),2],
                                P7.tabl[min(i+1,nrow(out.tabl)),2],P7.tabl[min(i+2,nrow(out.tabl)),2]))
  }  
}


out.ta.th   = out.tabl[complete.cases(out.tabl),] # selected out table without Nas
if(RaSc == "bonaP") {inp.ta.th   = inp.tabl[complete.cases(inp.tabl),]  # eliminate NA-rows within the input
                     inp.1min.2  = P2.tabl.2[,2]
                     inp.meas.2  = P2.tabl.2[,2]
                     pluv.tabl   = P2.tabl.2
}else               {inp.ta.th   = inp.tabl.FL[complete.cases(inp.tabl.FL),] 
                     inp.1min.2  = P7.tabl[,2]
                     inp.meas.2  = P7.tabl.raw[,2]/10 # important! these data are not mm/min but mm/10min
                     pluv.tabl   = P7.tabl.raw
}
out.meas.2  = out.ta.th[,11]

d.t.ca.th.2 = out.ta.th[,1] # now thinned at 4 min
inp.cal.2   = inp.ta.th[,2] # better averaged inp [mm/min] every 4 minutes (probably useless)
inp4mod.2   = inp.cal.2*60/1000# mm/min to m/hr
imber.vel   = inp.cal.2/60    # mm/min to mm/s
imber_ii    = c(rep(0,par.fix["Del.Max"]),imber.vel) # shifted to account for delay betw inp-out


# set time & create a calibration layout
# ---------------------------------------
D.dies    <- (as.numeric(as.Date(d.t.ca.th.2[1])-as.Date(d.t.cal.th[1])))*24
# for output labelling with time
D.hodie  <-  hms.to.seconds(format(d.t.ca.th.2[1], format="%H:%M:%S",tz=""))/3600
quom.prin <- D.dies + D.hodie # hr elapsed from the midn preceding I ev to beg of this one
ult.grid.Ome.2  <- (length(d.t.ca.th.2)-1)*dt.mod+quom.prin #hr 
t.grid.Ome.2 = format(seq(quom.prin,ult.grid.Ome.2,dt.mod), nsmall=6)

# for input labelling with time (1 or 10 min)
D.hodie  <-  hms.to.seconds(format(pluv.tabl[1,1], format="%H:%M:%S",tz=""))/3600
quom.prin.ra <- D.dies + D.hodie # hr elapsed from the midn preceding I ev to beg of this one
ult.grid.2.ra  <- (nrow(pluv.tabl)-1)*dt.inp+quom.prin.ra #hr 
t.grid.Ime.2 =   format(seq(quom.prin.ra,ult.grid.2.ra,dt.inp), nsmall=6)

# for input labelling with finest time (1min)
D.hodie  <-  hms.to.seconds(format(P2.tabl.2[1,1], format="%H:%M:%S",tz=""))/3600
quom.prin.ra <- D.dies + D.hodie # hr elapsed from the midn preceding I ev to beg of this one
ult.grid.2.ra  <- (nrow(P2.tabl.2)-1)*1/60+quom.prin.ra #hr 
t.grid.1mi.2 =   format(seq(quom.prin.ra,ult.grid.2.ra,1/60), nsmall=6)

# check if output measurements come first than rain; if yes, remove the first value
if(as.numeric(t.grid.Ome.2[1])<as.numeric(t.grid.Ime.2[1]))
{
  out.meas.2 = out.meas.2[2:length(t.grid.Ome.2)]
  inp4mod.2  = inp4mod.2[2:length(t.grid.Ome.2)] # reduce the length of the "artificially" aggregated rain
  t.grid.Ome.2     = t.grid.Ome.2[2:length(t.grid.Ome.2)]
  out.meas.2.br = out.meas.2
  inp4mod.2.br  = inp4mod.2 # reduce the length of the "artificially" aggregated rain
  t.grid.Ome.2.br     = t.grid.Ome.2
}

L.obs.2   = paste("Q",t.grid.Ome.2,sep="_") 
ind.comp.xi.2 = match(t.grid.Ime.2,t.grid.1mi.2)
names(out.meas.2)  = L.obs.2
names(inp.1min.2)  = t.grid.1mi.2 # meas rain [mm/min] every min [for out lkh test]
names(inp.meas.2)  = t.grid.Ime.2 # meas rain [mm/min] every min or 10 (dep on pluv) [for inference]
names(inp4mod.2)   = t.grid.Ome.2 # aggregated rain [mm/min] every 4min [for mod test]

##############################################################################
# Period III - Calibration

ev.nr   = cal.ev[3] 
ind.per.cal = seq(from=which(inp.data[,1]==selected_events1[ev.nr,2]-3600),to=which(inp.data[,1]==selected_events1[ev.nr,3]+3600*2))

da.ti.cal   = inp.data[ind.per.cal,1] # UTC time of the data (every min)
out.tabl    = out.data[match(da.ti.cal,out.data[,1]),] # l/s (with 1 min res)
P3.tabl     = inp.3[match(da.ti.cal,inp.3[,1]),] # 1 min data from P3 (with same index as P1)
P2.tabl.3   = inp.2[match(da.ti.cal,inp.2[,1]),]; P2.tabl.3[,2] = apply(cbind(P3.tabl[,2],P2.tabl.3[,2]), 1, mean)
P7.tabl     = inp.7.2m[match(da.ti.cal,inp.7.2m[,1]),]
P7.tabl.raw = inp.7[match(da.ti.cal,inp.7[,1]),]
P7.tabl.raw = P7.tabl.raw[complete.cases(P7.tabl.raw),] 
inp.tabl    = inp.data[ind.per.cal,]; inp.tabl[,1]=out.tabl[,1]; inp.tabl[,2]=NA
inp.tabl.FL = inp.tabl

# calculate the mean rain rate over last 4 mins
for (i in 1:nrow(out.tabl))
{
  if(!is.na(out.tabl[i,1])) 
  {
    inp.tabl[i,2]     <- mean(c(P2.tabl.3[i,2],P2.tabl.3[max((i-1),1),2],
                                P2.tabl.3[min(i+1,nrow(out.tabl)),2],P2.tabl.3[min(i+2,nrow(out.tabl)),2]))
    inp.tabl.FL[i,2]  <- mean(c(P7.tabl[i,2],P7.tabl[max((i-1),1),2],
                                P7.tabl[min(i+1,nrow(out.tabl)),2],P7.tabl[min(i+2,nrow(out.tabl)),2]))
  }  
}

out.ta.th   = out.tabl[complete.cases(out.tabl),] # selected out table without Nas
if(RaSc == "bonaP") {inp.ta.th   = inp.tabl[complete.cases(inp.tabl),]  # eliminate NA-rows within the input
                     inp.1min.3  = P2.tabl.3[,2]
                     inp.meas.3  = P2.tabl.3[,2]
                     pluv.tabl   = P2.tabl.3
}else               {inp.ta.th   = inp.tabl.FL[complete.cases(inp.tabl.FL),] 
                     inp.1min.3  = P7.tabl[,2]
                     inp.meas.3  = P7.tabl.raw[,2]/10 # important! these data are not mm/min but mm/10min
                     pluv.tabl   = P7.tabl.raw
}
out.meas.3  = out.ta.th[,11]

d.t.ca.th.3 = out.ta.th[,1] # now thinned at 4 min
inp.cal.3   = inp.ta.th[,2] # better averaged inp [mm/min] every 4 minutes (probably useless)
inp4mod.3   = inp.cal.3*60/1000# mm/min to m/hr
imber.vel   = inp.cal.3/60    # mm/min to mm/s
imber_iii   = c(rep(0,par.fix["Del.Max"]),imber.vel) # shifted to account for delay betw inp-out

# set time & create a calibration layout
# ---------------------------------------
# D.dies: hrs elapsed from the midnight following the first/reference event to the midnight previous to this event
# D.fin.in: hrs remaining in the day of the first event after its end + hrs of from midnight to the beginning of this event 
# ult.cal: hrs from the 00:00 of the day of the beginning of the I event to the end of the I event
D.dies    <- (as.numeric(as.Date(d.t.ca.th.3[1])-as.Date(d.t.cal.th[1])))*24
# for output labelling with time
D.hodie  <-  hms.to.seconds(format(d.t.ca.th.3[1], format="%H:%M:%S",tz=""))/3600
quom.prin <- D.dies + D.hodie # hr elapsed from the midn preceding I ev to beg of this one
ult.grid.Ome.3  <- (length(d.t.ca.th.3)-1)*dt.mod+quom.prin #hr 
t.grid.Ome.3 = format(seq(quom.prin,ult.grid.Ome.3,dt.mod), nsmall=6)

# for input labelling with time (1 or 10 min)
D.hodie  <-  hms.to.seconds(format(pluv.tabl[1,1], format="%H:%M:%S",tz=""))/3600
quom.prin.ra <- D.dies + D.hodie # hr elapsed from the midn preceding I ev to beg of this one
ult.grid.3.ra  <- (nrow(pluv.tabl)-1)*dt.inp+quom.prin.ra #hr 
t.grid.Ime.3 = format(seq(quom.prin.ra,ult.grid.3.ra,dt.inp), nsmall=6)

# for input labelling with finest time (1min)
D.hodie  <-  hms.to.seconds(format(P2.tabl.3[1,1], format="%H:%M:%S",tz=""))/3600
quom.prin.ra <- D.dies + D.hodie # hr elapsed from the midn preceding I ev to beg of this one
ult.grid.3.ra  <- (nrow(P2.tabl.3)-1)*1/60+quom.prin.ra #hr 
t.grid.1mi.3 =   format(seq(quom.prin.ra,ult.grid.3.ra,1/60), nsmall=6)

# check if output measurements come first than rain; if yes, remove the first value
if(as.numeric(t.grid.Ome.3[1])<as.numeric(t.grid.Ime.3[1]))
{
  out.meas.3 = out.meas.3[2:length(t.grid.Ome.3)]
  inp4mod.3  = inp4mod.3[2:length(t.grid.Ome.3)] # reduce the length of the "artificially" aggregated rain
  t.grid.Ome.3     = t.grid.Ome.3[2:length(t.grid.Ome.3)]
  out.meas.3.br = out.meas.3
  inp4mod.3.br  = inp4mod.3 # reduce the length of the "artificially" aggregated rain
  t.grid.Ome.3.br     = t.grid.Ome.3
}

L.obs.3   = paste("Q",t.grid.Ome.3,sep="_") 
ind.comp.xi.3 = match(t.grid.Ime.3,t.grid.1mi.3)
names(out.meas.3)  = L.obs.3
names(inp.1min.3)  = t.grid.1mi.3 # meas rain [mm/min] every min [for out lkh test]
names(inp.meas.3)  = t.grid.Ime.3 # meas rain [mm/min] every min or 10 (dep on pluv) [for inference]
names(inp4mod.3)   = t.grid.Ome.3 # aggregated rain [mm/min] every 4min [for mod test]

##############################################################################
# Data preparation complete - calibration preparation starts
##############################################################################


model.function <- lin_1box_ww # simulator

# Define parameters
# ----------------------------

par <- par.opt.ini[-which(names(par.opt.ini)=="sd.B_Q")] # initial precalibrated parameters

ref.out   <- c(y.ref_Q = 50) # high output [l/s] (used for transformation)
par.init  <- c(par) 

if(EM == "mult") par.init  <- c(par, beta1=1, beta2=1, beta3=1, sdv_b=0.1)
if(EM == "Bias") par.init  <- c(par, par.opt.ini["sd.B_Q"], corrlen=7*dt.mod)
if(EM == "SIP"){ # inp transf parameters + correlation of the OU process 
  avg.inp.par <- c(corrl_inp=tau.opt.dim)
  sd.inp.par <- c(corrl_inp=tau.sd.dim)
}

pri_sd         <- par.init/10 # 10% CoV for most parameters
pri_sd["nco1"] <- 1.860135e-02 # all values from a previous dry weather calibration
pri_sd["co2"]  <- 1.049271e-02
pri_sd["bf"]   <- 4.830275e-02
pri_sd["nsi1"] <- 9.400194e-02
pri_sd["nsi2"] <- 1.865796e-02
pri_sd["k"]    <- par.init["k"]*.2
pri_sd["sd.Eps_Q"]    <- par.init["sd.Eps_Q"]*.2
lim.sx.ks_Q    <- NA

if (transf=="BC")  par.tr = c(l1=1, l2=1)
if (transf=="ls") par.tr = c(alpha=25, beta=50) # chosen option

# Tranform error model hyperparameters
# --------------------------------------

if(is.na(par.tr["alpha"])) 
{  # BC or no transformation
  if (!is.na(pri_sd["corrlen"]))   par.init[paste("sd.B",var, sep="_")]  <- par.init[paste("sd.B",var, sep="_")]*sysanal.boxcox.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                                                                                      par.tr["l1"],par.tr["l2"])
  par.init[paste("sd.Eps",var, sep="_")]<- par.init[paste("sd.Eps",var, sep="_")]*sysanal.boxcox.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                                                       par.tr["l1"],par.tr["l2"])
  if (!is.na(par.init[paste("ks",var, sep="_")])) 
  { # we have input dependence via an heteroskedastik EM
    
    lim.sx.ks_Q = as.numeric(sysanal.boxcox(0,par.tr["l1"],par.tr["l2"])/max(imber.vel))
    
    pri_sd["ks_Q"]= (-lim.sx.ks_Q+(sysanal.boxcox(par.init["ks_Q"],par.tr["l1"],par.tr["l2"])/max(imber.vel)))
    
    par.init["ks_Q"] = 0.1*pri_sd["ks_Q"]} # start from a small value
  
}  else  { # log-sinh
  if (!is.na(pri_sd["corrlen"]))  par.init[paste("sd.B",var, sep="_")]  <- par.init[paste("sd.B",var, sep="_")]*sysanal.logsinh.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                                                                                      par.tr["alpha"],par.tr["beta"])
  par.init[paste("sd.Eps",var, sep="_")]<- par.init[paste("sd.Eps",var, sep="_")]*sysanal.logsinh.deriv(ref.out[paste("y.ref",var, sep="_")],
                                                                                                        par.tr["alpha"],par.tr["beta"])
  if (!is.na(par.init[paste("ks",var, sep="_")]))  
    
  {
    lim.sx.ks_Q = as.numeric(sysanal.logsinh(0,par.tr["alpha"],par.tr["beta"])/max(imber.vel))
    
    pri_sd["ks_Q"]= (-lim.sx.ks_Q+(sysanal.logsinh(par.init["ks_Q"],par.tr["alpha"],par.tr["beta"])/max(imber.vel)))
    
    par.init["ks_Q"] = 0.1*pri_sd["ks_Q"]
    
  } 
} 


# define prior parameter distribution
# --------------------------------------

prior.pbdis<- list(A       =c("Lognormal", par.init["A"],pri_sd["A"]),
                   bf      =c("Lognormal", par.init["bf"],pri_sd["bf"]),
                   nco1    =c("Lognormal", par.init["nco1"],abs(pri_sd["nco1"])),
                   co2     =c("Lognormal", par.init["co2"],pri_sd["co2"]),
                   k       =c("Lognormal", par.init["k"],pri_sd["k"] ),
                   nsi1    =c("Lognormal", par.init["nsi1"],abs(pri_sd["nsi1"])),
                   nsi2    =c("Lognormal", par.init["nsi2"],abs(pri_sd["nsi2"])),
                   sdv_b   =c("Lognormal", par.init["sdv_b"],.2*par.init["sdv_b"]),
                   beta1    =c("Lognormal", par.init["beta1"], "sdv_b"),
                   beta2    =c("Lognormal", par.init["beta2"], "sdv_b"),
                   beta3    =c("Lognormal", par.init["beta3"], "sdv_b"),
                   corrlen =c("Lognormal", par.init["corrlen"],pri_sd["corrlen"]), # [hr]; same units as layout
                   sd.Eps_Q=c("Lognormal", par.init["sd.Eps_Q"],.1*par.init["sd.Eps_Q"]),#param of the observ err model
                   sd.B_Q  =c("NormalTrunc", 0,par.init["sd.B_Q"],0, 1e+08 ),
                   ks_Q    =c("NormalTrunc", lim.sx.ks_Q, pri_sd["ks_Q"],lim.sx.ks_Q, 1e+08 ), # not used here, kept fix
                   Delta   =c("Exponential", par.init["Delta"]) # not used here, kept fix
)

#define objective function (logposterior)

likel.fun  <- sysanal.loglikeli.bias.inp # Eq. 1
if(EM == "mult") likel.fun  <- sysanal.loglikeli.mult

logposterior.3 <- function(par) # for classical error models (i.e. not SIP)
{ names(par)=names(par.init)
  out=sysanal.logposterior.3(par,
                             model         = model.function,
                             L.1           = L.obs.1,
                             y.1           = out.meas.1,
                             L.2           = L.obs.2,
                             y.2           = out.meas.2,
                             L.3           = L.obs.3,
                             y.3           = out.meas.3,
                             prior.dist    = "indep",
                             prior.def     = prior.pbdis,
                             loglikeli     = likel.fun,
                             par.fix       = par.fix,
                             par.tr        = par.tr,
                             Var.Bs        = sysanal.Var.Bs, # output covariance
                             inp.1         = imber_i, # in case one uses an heteroskedastik lkh
                             Inp.1         = inp4mod.1,
                             inp.2         = imber_ii, 
                             Inp.2         = inp4mod.2,
                             inp.3         = imber_iii, 
                             Inp.3         = inp4mod.3,
                             dt            = dt.mod,
                             sd.Eps        = sysanal.sd.Eps.L)
  return(out)}


# # >> cov mat of the jump distribution of the out parameters
out.cand.cov <- diag((pri_sd/2)^2,length(par.init)) # intial guess (use a better one if you have it)


########################################################
# FUNCTIONS REQUIRED FOR SIP INFERENCE

if(EM == "SIP")
{

# >> prior density of output parameters evaluated at the proposed samples
out.logpri <- function(out.par) # f(th,psi_y) [used in Eq. 15]
{ out=sysanal.calcpdf_mv(z       = out.par,
                         dist    = "indep",
                         distdef = prior.pbdis)
  return(out)} 


# >> definition of the prior distribution of the input parameter
inp.pri.par <- list(corrl_inp       =c("Lognormal", avg.inp.par["corrl_inp"],sd.inp.par["corrl_inp"])) # psi_x


# >> pdf the true input (alias: density of the probabilistic model of the true input)  [used in Eq. 16]
inp.pro.den <- function(t.grid.cal,inp.par,sim.inp)
{ out=dou(t   = t.grid.cal,
          tau = inp.par,
          y   = sim.inp) # modeled stochastic rainrate (in the space of the rainfall potential)
  return(out)}   


# >> inp transf funct from OU to ppt [Eq. 11]
inp.trans <- function(OU.traj)
{ out=tr.norm2ppt(OU.traj,par=par.tr.opt)
  return(out)}   


# >> inverse inp transf funct: from rain [mm/min] to potential rain (only for rainrate>0) [used as described after Eq. 17]
rain.trans <- function(RI,par=par.tr.opt)
{
  RI.cr = tr.norm2ppt(x=par.tr.opt["x0"],par=par.tr.opt) #mm/min
  b2 <- par["x0"] - ( par["alpha1"]*par["a1"]/(par["alpha2"]*par["a2"])*(par["x0"]-par["b1"])^(par["alpha1"]-1) )^(1/(par["alpha2"]-1))
  c2 <- par["a1"]*(par["x0"]-par["b1"])^par["alpha1"] - par["a2"]*(par["x0"]-b2)^par["alpha2"]
  xi <- rep(NA, length(RI))
  for (i in 1:length(RI))
  {if(RI[i]>RI.cr) xi[i]  <- b2 + ((RI[i]-c2)/par["a2"])^(1/par["alpha2"])
   if(RI[i]<RI.cr && RI[i]>0) xi[i]  <- par["b1"] + ((RI[i]-0)/par["a1"])^(1/par["alpha1"])}
  return(xi)
}


# >> Probabilistic hydrological model feeded with true input (alias: output likelihood) [used in Eq. 15]
out.loglikh <- function(out.par, L, y.obs, sim.inp,t.grid.cal, displ.res=F) # f(yo|th,psi_y,x)
{ names(out.par)=names(par.init)
  
  # use the right time resolution for the hydrological model (res[inp]=<res[out data])
  mod.tim.ind <- match(sysanal.decode(L)$val,as.numeric(t.grid.cal)) # match output and input time
  mod.inp     <- sim.inp[mod.tim.ind]   
out=sysanal.loglikeli.bias.inp(par       = out.par,
                                 Inp           = mod.inp, # input of the model (rainfall) [mm/hr]
                                 model         = model.function,
                                 L             = L,
                                 y.obs         = y.obs,
                                 par.fix       = par.fix,
                                 par.tr        = par.tr,
                                 Var.Bs        = sysanal.Var.Bs,
                                 dt            = dt.mod, # model time step
                                 sd.Eps        = sysanal.sd.Eps.L)
  if (displ.res == T) {
    plot(t.grid.cal[mod.tim.ind],y.obs, ylim = c(0,100))  
    lines(t.grid.cal[mod.tim.ind],model.function(par=out.par,L=L, Inp=mod.inp, dt=dt.mod)  )
    legend("topright", legend=paste("Lkh =", out))
  }

  return(out)} 

# >> prior parameters of the input error model [used in Eq. 13]
inp.ob.pa.ini <- c(var.xi  = .4, beta.xi = 100000000000) # beta.xi is actually not used [Fig. S1 - S2]

# >> cov mat of the jump distribution of the inp parameters [used in Eq. 16]
inp.cand.cov <- diag((inp.ob.pa.ini*.1)^2,length(inp.ob.pa.ini))
colnames(inp.cand.cov) = names(inp.ob.pa.ini);rownames(inp.cand.cov) = names(inp.ob.pa.ini)

# >> prior density of input parameters evaluated at the proposed samples [used in Eq. 16]
inp.logpri <- function(inp.par) # fpsi_x(psi_x)
{ out=sysanal.calcpdf_mv(z       = inp.par,
                         dist    = "indep",
                         distdef = list(var.xi  =c("Lognormal", inp.ob.pa.ini["var.xi"],inp.ob.pa.ini["var.xi"]*.5),
                                        beta.xi =c("Lognormal", inp.ob.pa.ini["beta.xi"],inp.ob.pa.ini["beta.xi"]*.5)))
  return(out)}   





##############################################################################
# Preparation complete - Calibration can start
##############################################################################


interv.pts = c(100,50) # defines the initial and final average number of time points in every subinterval for inp sampling
thin=1

ptm <- proc.time() 
Post.Par.Inp <-   mcmc.sip.3e(out.par.ini = par.init,
                               inp.par.fix = c(avg.inp.par,thresh=as.numeric(par.tr.opt["b1"])),
                               inp.par.ini = inp.ob.pa.ini, #   var.xi  
                               t.grid.1mi.1  = t.grid.1mi.1,
                               t.grid.Ime.1  = t.grid.Ime.1,
                               t.grid.1mi.2  = t.grid.1mi.2,
                               t.grid.Ime.2  = t.grid.Ime.2,
                               t.grid.1mi.3  = t.grid.1mi.3,
                               t.grid.Ime.3  = t.grid.Ime.3,
                               inp.cand.cov= inp.cand.cov,
                               inp.logpri  = inp.logpri,
                               L.1  = L.obs.1, 
                               y.1  = out.meas.1,
                               L.2  = L.obs.2, 
                               y.2  = out.meas.2,
                               L.3  = L.obs.3, 
                               y.3 = out.meas.3,
                               interv.pts  = interv.pts,
                               adapt_par =c(100,500,0.5,0.5),
                               ppt2norm.tra  = rain.trans,
                               out.cand.cov= out.cand.cov,
                               out.loglikh = out.loglikh,
                               out.logpri  = out.logpri,
                               meas.rain.1 = inp.meas.1,#[mm/min]
                               meas.rain.2 = inp.meas.2,#[mm/min]
                               meas.rain.3 = inp.meas.3,#[mm/min]
                               sampsize    = iter.inf,
                              traj.Inp.pr.1= NA, # if you have a good OU realization to start, use it here
                              traj.Inp.pr.2= NA,
                              traj.Inp.pr.3= NA,
                              thin = thin)
comp.time = proc.time() - ptm
comp.time

# Plot results

Post.obs.inp.par = Post.Par.Inp$post.inp.par[nrow(Post.Par.Inp$post.inp.par),]
dim.sam = nrow(Post.Par.Inp$post.inp.1)
len.ts.1= ncol(Post.Par.Inp$post.inp.1)
len.ts.2= ncol(Post.Par.Inp$post.inp.2)
len.ts.3= ncol(Post.Par.Inp$post.inp.3)
real2plo= min(dim.sam/2,200)

.pardefault <- par() 
# Figure 1
# visualize inference results
ind.sav =2:iter.inf/thin
Post.par.samp = cbind(Post.Par.Inp$post.out.par[ind.sav,-2],var.xi = Post.Par.Inp$post.inp.par[ind.sav,1], 
                      Post.Par.Inp$out.log.dens[ind.sav,], 
                      xi.1At30pct = Post.Par.Inp$post.inp.1[ind.sav,len.ts.1/3], xi.1At50pct = Post.Par.Inp$post.inp.1[ind.sav,len.ts.1/2],
                      xi.2At30pct = Post.Par.Inp$post.inp.2[ind.sav,len.ts.2/3], xi.2At50pct = Post.Par.Inp$post.inp.2[ind.sav,len.ts.2/2],
                      xi.3At30pct = Post.Par.Inp$post.inp.3[ind.sav,len.ts.1/3], xi.3At50pct = Post.Par.Inp$post.inp.3[ind.sav,len.ts.3/2])
sysanal.plot.chains(Post.par.samp)


# FIGURE 2.1 (just for first calibration period)
# posterior output predictions without measurement noise
# visually check the inferred rain (true space + meas rain) [mm/min]


par(mfrow=c(2,1), mar=c(0,1.5,0.1,0), oma= c(2,0,0.0,0), mgp = c(2.5, .5, 0))
for (l in seq(from = dim.sam/2, to = dim.sam, length.out = real2plo)) 
{ plot(t.grid.1mi.1,inp.trans(Post.Par.Inp$post.inp.1[l,]), col=gray(.81), xaxt = "n", type="l", ylim=c(0,max(P2.tabl.2[,2])*1.01))
  par(new=T)}
# abline(v=t.grid.1mi.2[max.II+35], lty="dotted", col="darkorchid")
points(t.grid.Ime.1,inp.meas.1, cex=.7, pch="|")

R2.inp <- 1-(sum((inp.meas.1-inp.trans(Post.Par.Inp$post.inp.1[dim.sam*.915,ind.comp.xi.1]))^2, na.rm = T))/sum((inp.meas.1-mean(inp.meas.1, na.rm = T))^2, na.rm = T)
legend("right", legend=paste("Typic Inp R2 =", R2.inp), bty="n")
par(new=F)

for (l in seq(from = dim.sam/2, to = dim.sam, length.out = real2plo)) 
{ sim.inp = inp.trans(Post.Par.Inp$post.inp.1[l,])
  out.par = Post.Par.Inp$post.out.par[l,]
  L       = L.obs.1
  mod.tim.ind <- match(sysanal.decode(L)$val,as.numeric(t.grid.1mi.1))
  mod.inp     <- sim.inp[mod.tim.ind]*60/1000 
  plot(t.grid.Ome.1,model.function(par=out.par,L=L, Inp=mod.inp, dt=dt.mod), col=gray(.81), type="l", ylim=c(0,90),xlim=as.numeric(c(t.grid.1mi.1[1],t.grid.1mi.1[length(t.grid.1mi.1)])))
  par(new=T)}
out = out.loglikh(Post.Par.Inp$post.out.par[dim.sam*.915,],sim.inp=sim.inp*60/1000 ,
                  t.grid.cal=as.numeric(t.grid.1mi.1),L=L, y.obs=out.meas.1,  displ.res=F)+ out.logpri(Post.Par.Inp$post.out.par[dim.sam*.915,])
legend("top", legend=paste("Typic Out Lkh =", out), bty="n")
NS.opt <- 1-(sum((out.meas.1-model.function(par=Post.Par.Inp$post.out.par[dim.sam*.915,],L=L, Inp=mod.inp, dt=dt.mod))^2, na.rm = T))/sum((out.meas.1-mean(out.meas.1, na.rm = T))^2, na.rm = T)
legend("right", legend=paste("NS =", NS.opt), bty="n")
points(t.grid.Ome.1,out.meas.1, cex=.5)


# FIGURE 3
# visually check the inference results
par(mfrow=c(3,1), mar=c(0,1.5,0.1,0), oma= c(2,0,0.0,0), mgp = c(2.5, .5, 0))

plot(Post.Par.Inp$inp.log.dens.1, cex=.2, xaxt = "n") 
legend("bottom", legend=paste("max inp dens.1 =", max(Post.Par.Inp$inp.log.dens.1)), bty="n")
plot(Post.Par.Inp$inp.log.dens.2, cex=.2, xaxt = "n") 
legend("bottom", legend=paste("max inp dens.2 =", max(Post.Par.Inp$inp.log.dens.2)), bty="n")
legend("center", legend=paste("inference cost =", round(comp.time[3]/3600,digits=2),"[hr]"), bty="n", text.col="blue") 
plot(Post.Par.Inp$inp.log.dens.3, cex=.2) 
legend("bottom", legend=paste("max inp dens.3 =", max(Post.Par.Inp$inp.log.dens.3)), bty="n")


} else {# run previous part only if EM == SIP; RUN this with OTHER ERROR MODELS

  # Posterior exploration with previous error models
  
jump.cov = out.cand.cov
  
  if(EM == "mult") {
    par.optim  <- c(par.init, beta1=1, beta2=1, beta3=1, sdv_b=0.1)
    iid.cov    <- jump.cov 
    jump.cov   <- diag((pri_sd/2)^2,length(par.init))
    jump.cov[1:nrow(iid.cov),1:ncol(iid.cov)] <- iid.cov
    colnames(jump.cov) <- names(par.init);rownames(jump.cov) <- names(par.init) 
  }
  

  if(EM == "Bias")  {
    par.init  <- c(par.init,sd.B_Q=as.numeric(par.opt.ini["sd.Eps_Q"]*2/3), corrlen=7*dt.mod)
    iid.cov    <- jump.cov 
    jump.cov   <- diag((pri_sd/2)^2,length(par.init))
    jump.cov[1:nrow(iid.cov),1:ncol(iid.cov)] <- iid.cov
    colnames(jump.cov) <- names(par.init);rownames(jump.cov) <- names(par.init) 
  
  }


ptm <- proc.time()
  MHa <- Metro_Hastings_0815(li_func    = logposterior.3 , 
                             pars       = par.init,
                             prop_sigma = jump.cov,
                             iterations = iter.inf, 
                             burn_in    = 1,
                             adapt_par  = c(100,100,0.5,.49),
                             quiet      = F)
comp.time = proc.time() - ptm
comp.time 
  

sampl.propa= min(iter.inf/4, 2001)

# Eliminate adaptive burn in phase
MCMC.propa = MHa$trace[seq(max(1,nrow(MHa$trace)-2*sampl.propa),nrow(MHa$trace),2),] # thinned every two elements

  
dim.sam = nrow(MHa$trace)
len.ts.1= length(inp4mod.1)
len.ts.2= length(inp4mod.2)
len.ts.3= length(inp4mod.3)
ParScen = paste("EM=", EM,"_RaSc", RaSc,sep="")
real2plo= min(dim.sam/2,200)

.pardefault <- par() 

# FIGURE 1.1
# posterior output predictions without measurement noise
# visually check the inferred rain (true space + meas rain) [mm/min]


par(mfrow=c(2,1), mar=c(0,1.5,0.1,0), oma= c(2,0,0.0,0), mgp = c(2.5, .5, 0))
plot(t.grid.1mi.1,rep(NA,length(t.grid.1mi.1)), col=gray(.81), xaxt = "n", type="l", ylim=c(0,max(P2.tabl.2[,2])*1.01))
points(t.grid.Ime.1,inp.meas.1, cex=.7, pch="|", xaxt = "n", ylim=c(0,max(P2.tabl.2[,2])*1.01))

for (l in seq(from = dim.sam/2, to = dim.sam, length.out = real2plo)) 
{ out.par = MHa$trace[l,]
  L       = L.obs.1
  mod.inp = inp4mod.1 
  plot(t.grid.Ome.1,model.function(par=out.par,L=L, Inp=mod.inp, dt=dt.mod), col=gray(.81), type="l", ylim=c(0,90),xlim=as.numeric(c(t.grid.1mi.1[1],t.grid.1mi.1[length(t.grid.1mi.1)])))
  par(new=T)}
NS.opt <- 1-(sum((out.meas.1-model.function(par=MHa$trace[dim.sam*.915,],L=L, Inp=mod.inp, dt=dt.mod))^2, na.rm = T))/sum((out.meas.1-mean(out.meas.1, na.rm = T))^2, na.rm = T)
legend("right", legend=paste("NS =", NS.opt), bty="n")
points(t.grid.Ome.1,out.meas.1, cex=.5)



# FIGURE 2
# visually check the inference results
Post.par.samp = cbind(MHa$trace,log.pdf=MHa$unnor.log.post)
sysanal.plot.chains(Post.par.samp)

}

# calibration completed



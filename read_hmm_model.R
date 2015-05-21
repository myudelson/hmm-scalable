
read_hmm_model <- function(filename) {
  raw <- readLines(filename)
  #
  # read header
  #
  solver_structure <- as.numeric(substr(strsplit(raw[1],"\t")[[1]][[2]],1,1))
  nPIg <- 0
  nAg  <- 0
  nBg  <- 0
  nK <- as.numeric(strsplit(raw[2],"\t")[[1]][[2]])
  nG <- as.numeric(strsplit(raw[3],"\t")[[1]][[2]])
  nS <- as.numeric(strsplit(raw[4],"\t")[[1]][[2]])
  nO <- as.numeric(strsplit(raw[5],"\t")[[1]][[2]])
  
#   STRUCTURE_PIgk    = 4,  // 4 - PI by skll&group, A,B by skill
#   STRUCTURE_PIAgk   = 5,  // 5 - PI, A by skll&group, B by skill
#   STRUCTURE_Agk     = 6,  // 6 - A by skll&group, PI,B by skill
#   STRUCTURE_PIABgk  = 7,  // 5 - PI, A, B by skll&group
  if( solver_structure %in% c(4,5,7) )
    nPIg <- nG
  if( solver_structure %in% c(5,6,7) )
    nAg <- nG
  if( solver_structure %in% c(7) )
    nBg <- nG
  
  line_no_nulls <- 6 + nPIg + nAg + nBg # last line before skills/students
  model <- data.frame(kc=character(nK), unit=character(nK), pLo=numeric(nK), pF=numeric(nK), pT=numeric(nK), pS=numeric(nK), pG=numeric(nK))
  model$kc <- as.character(model$kc)
  model$unit <- as.character(model$unit)
  for( i in 1:nK ) { # for all KCs
      ln <- line_no_nulls + (i-1) * 4 + 1
      model$kc[i]   <- as.character(strsplit(raw[ln],"\t")[[1]][[2]])
      model$unit[i] <-          strsplit(model$kc[i],"__")[[1]][[1]]
      model$pLo[i]  <- as.numeric(strsplit(raw[ln+1],"\t")[[1]][[2]])
      model$pF[i]   <- as.numeric(strsplit(raw[ln+2],"\t")[[1]][[3]])
      model$pT[i]   <- as.numeric(strsplit(raw[ln+2],"\t")[[1]][[4]])
      model$pS[i]   <- as.numeric(strsplit(raw[ln+3],"\t")[[1]][[3]])
      model$pG[i]   <- as.numeric(strsplit(raw[ln+3],"\t")[[1]][[4]])
  } # for all KCs
  return ( model )
}

four_densities <- function( data, a_legend, a_main, label, colors) {
  xr <- range(density(data[[1]])$x)
  yr <- range(density(data[[1]])$y)
  if( xr[1]>range(density(data[[2]])$x)[1] )
    xr[1] <- range(density(data[[2]])$x)[1]
  if( xr[1]>range(density(data[[3]])$x)[1] )
    xr[1] <- range(density(data[[3]])$x)[1]
  if( xr[1]>range(density(data[[4]])$x)[1] )
    xr[1] <- range(density(data[[4]])$x)[1]
  if( yr[1]>range(density(data[[2]])$y)[1] )
    yr[1] <- range(density(data[[2]])$y)[1]
  if( yr[1]>range(density(data[[3]])$y)[1] )
    yr[1] <- range(density(data[[3]])$y)[1]
  if( yr[1]>range(density(data[[4]])$y)[1] )
    yr[1] <- range(density(data[[4]])$y)[1]
  if( xr[2]<range(density(data[[2]])$x)[2] )
    xr[2] <- range(density(data[[2]])$x)[2]
  if( xr[2]<range(density(data[[3]])$x)[2] )
    xr[2] <- range(density(data[[3]])$x)[2]
  if( xr[2]<range(density(data[[4]])$x)[2] )
    xr[2] <- range(density(data[[4]])$x)[2]
  if( yr[2]<range(density(data[[2]])$y)[2] )
    yr[2] <- range(density(data[[2]])$y)[2]
  if( yr[2]<range(density(data[[3]])$y)[2] )
    yr[2] <- range(density(data[[3]])$y)[2]
  if( yr[2]<range(density(data[[4]])$y)[2] )
    yr[2] <- range(density(data[[4]])$y)[2]
  
  
  plot(density(data[[1]]),xlim=xr,ylim=yr,main=a_main,xlab=label, col=colors[1])
  par(new=T)
  plot(density(data[[2]]),xlim=xr,ylim=yr,col=colors[2],main="",xlab="",ylab="")
  par(new=T)
  plot(density(data[[3]]),xlim=xr,ylim=yr,col=colors[3],main="",xlab="",ylab="")
  par(new=T)
  plot(density(data[[4]]),xlim=xr,ylim=yr,col=colors[4],main="",xlab="",ylab="")
  legend("topright",legend=a_legend,col=colors,lty=c(1,1,1,1))
}

model_a89uskts1.2  <- read_hmm_model("a89_model_uskts_multi1.2.txt")
model_a89cli  <- read_hmm_model("2008_CT_model_a89shrt.txt")
model_a89irtfx <- read_hmm_model("a89_model_uskts_multi1.2pLoIRTfx.txt")
model_a89irtft <- read_hmm_model("a89_model_uskts_multi1.2pLoIRT.txt")

model_UOPXuskts1.2 <- read_hmm_model("uopx_model208_uskts_multi1.2.txt")
model_UOPXcli <- read_hmm_model("uopx_CLI_model_shrt208.txt")
model_UOPXirtfx <- read_hmm_model("uopx_model208_uskts_multi1.2pLoLLfx.txt")
model_UOPXirtft <- read_hmm_model("uopx_model208_uskts_multi1.2pLoLL.txt")


tier = "208"
if(TRUE) {
mod_uopx_cli <- read_hmm_model("monet208/monet208_usktsCLI_model_shrt.txt")
mod_uopx_fit <- read_hmm_model("monet208/monet208_uskts_model_multi1.2CLI.txt")
mod_uopx_lls <- read_hmm_model("monet208/monet208_uskts_model_multi1.2CLIpLoLL.txt")
mod_uopx_llf <- read_hmm_model("monet208/monet208_uskts_model_multi1.2CLIpLoLLfx.txt")

label = paste("UOPX",tier)
layout(matrix(c(1,2,3,4),2,2,byrow=T))
four_densities( list(mod_uopx_cli$pLo, mod_uopx_fit$pLo, mod_uopx_llf$pLo, mod_uopx_lls$pLo), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                label,"pLo",
                c("black","red","blue","magenta"))
four_densities( list(mod_uopx_cli$pT, mod_uopx_fit$pT, mod_uopx_llf$pT, mod_uopx_lls$pT), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                label,"pT",
                c("black","red","blue","magenta"))
four_densities( list(mod_uopx_cli$pS, mod_uopx_fit$pS, mod_uopx_llf$pS, mod_uopx_lls$pS), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                label,"pS",
                c("black","red","blue","magenta"))
four_densities( list(mod_uopx_cli$pG, mod_uopx_fit$pG, mod_uopx_llf$pG, mod_uopx_lls$pG), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                label,"pG",
                c("black","red","blue","magenta"))
}

plot(mod_uopx_cli$pLo,mod_uopx_fit$pLo)
plot(mod_uopx_cli$pT,mod_uopx_fit$pT)
plot(mod_uopx_cli$pS,mod_uopx_fit$pS)
plot(mod_uopx_cli$pG,mod_uopx_fit$pG)

# UOPX
# shipped RMSE=0.401979 Acc=0.789994
# fit     RMSE=0.369913 Acc=0.819893
# irt fx  RMSE=0.372148 Acc=0.818399
# irt fit RMSE=0.369913 Acc=0.819893

# A89
# shipped RMSE=0.384395 Acc=0.823886
# fit     RMSE=0.359026 Acc=0.829359
# irt fx  RMSE=0.360302 Acc=0.828155
# irt fit RMSE=0.358839 Acc=0.829406



plot(density(model_a89uskts1.2$pLo),col="red")
par(new=T)
xr <- range(density(model_a89uskts1.2$pLo)$x)
yr <- range(density(model_a89uskts1.2$pLo)$y)
plot(density(model_a89cli$pLo),xlim=xr,ylim=yr)

# UoPX
layout(matrix(c(1,2,3,4),2,2,byrow=T))
four_densities( list(model_UOPXcli$pLo, model_UOPXuskts1.2$pLo, model_UOPXirtfx$pLo, model_UOPXirtft$pLo), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "UoPX","pLo",
                c("black","red","blue","magenta"))
four_densities( list(model_UOPXcli$pT, model_UOPXuskts1.2$pT, model_UOPXirtfx$pT, model_UOPXirtft$pT), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "UoPX","pT",
                c("black","red","blue","magenta"))
four_densities( list(model_UOPXcli$pS, model_UOPXuskts1.2$pS, model_UOPXirtfx$pS, model_UOPXirtft$pS), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "UoPX","pS",
                c("black","red","blue","magenta"))
four_densities( list(model_UOPXcli$pG, model_UOPXuskts1.2$pG, model_UOPXirtfx$pG, model_UOPXirtft$pG), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "UoPX","pG",
                c("black","red","blue","magenta"))

# A89
layout(matrix(c(1,2,3,4),2,2,byrow=T))
four_densities( list(model_a89cli$pLo, model_a89uskts1.2$pLo, model_a89irtfx$pLo, model_a89irtft$pLo), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "KDD A89","pLo",
                c("black","red","blue","magenta"))
four_densities( list(model_a89cli$pT, model_a89uskts1.2$pT, model_a89irtfx$pT, model_a89irtft$pT), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "KDD A89","pT",
                c("black","red","blue","magenta"))
four_densities( list(model_a89cli$pS, model_a89uskts1.2$pS, model_a89irtfx$pS, model_a89irtft$pS), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "KDD A89","pS",
                c("black","red","blue","magenta"))
four_densities( list(model_a89cli$pG, model_a89uskts1.2$pG, model_a89irtfx$pG, model_a89irtft$pG), 
                c("shipped","fit","irt pLo fixed","irt pLo fit"),
                "KDD A89","pG",
                c("black","red","blue","magenta"))




# akc1.2 <- read.table("a891.2kcs.txt")
# ukc1.2 <- read.table("uopx1.2kcs.txt")
# akc7.2 <- read.table("a897.2kcs.txt")
# ukc7.2 <- read.table("uopx7.2kcs.txt")
# 
# model_a89uskts1.2$rmse  <-akc1.2[,4]
# model_a89uskts1.2$acc   <-akc1.2[,6]
# model_UOPXuskts1.2$rmse <-ukc1.2[,4]
# model_UOPXuskts1.2$acc  <-ukc1.2[,6]
# 
# model_a89uskts1.2$rmse  <-akc7.2[,4]
# model_a89uskts1.2$acc   <-akc7.2[,6]
# model_UOPXuskts1.2$rmse <-ukc7.2[,4]
# model_UOPXuskts1.2$acc  <-ukc7.2[,6]
# 
# a89_uopx_merge <- merge(model_a89uskts1.2, model_UOPXuskts1.2, by=c("kc"), suffixes=c(".A",".U"))
# plot(a89_uopx_merge$rmse.A,a89_uopx_merge$rmse.U)
# abline(lm(a89_uopx_merge$rmse.U~a89_uopx_merge$rmse.A))
# cor.test(a89_uopx_merge$rmse.A,a89_uopx_merge$rmse.U) # 0.6845783 for 1.2, 0.6824376  for 7.2
# 
# plot(a89_uopx_merge$acc.A,a89_uopx_merge$acc.U)
# abline(lm(a89_uopx_merge$acc.U~a89_uopx_merge$acc.A))
# cor.test(a89_uopx_merge$acc.A,a89_uopx_merge$acc.U) # 0.7197323 for 1.2, 0.6983421 for 7.2


stderr <- function(x) {
  return ( sd(x)/sqrt(length(x)) )
}

agAmu <- aggregate(cbind(pLo,pT,pS,pG)~unit,data=model_a89uskts1.2,  FUN=mean)
agUmu <- aggregate(cbind(pLo,pT,pS,pG)~unit,data=model_UOPXuskts1.2, FUN=mean)
agAse <- aggregate(cbind(pLo,pT,pS,pG)~unit,data=model_a89uskts1.2,  FUN=stderr)
agUse <- aggregate(cbind(pLo,pT,pS,pG)~unit,data=model_UOPXuskts1.2, FUN=stderr)
agA <- merge(agAmu, agAse,by=c("unit"),suffixes=c("mu","se"))
agU <- merge(agUmu, agUse,by=c("unit"),suffixes=c("mu","se"))


agAU<-merge(agA,agU,by=c("unit"), suffixes=c(".A",".U"))

plot(agAU$pLo.U,type="l",col="blue")
par(new=T)
lines(agAU$pLo.A,type="l")

agAU_pLo <- agAU[,c("unit","pLomu.A","pLose.A","pLomu.U","pLose.U")]
agAU_pLo$sig <- 0
for(i in 1:dim(agAU_pLo)[1]) {
  if( (agAU_pLo$pLomu.A+agAU_pLo$pLose.A*1.96)<=(agAU_pLo$pLomu.U-agAU_pLo$pLose.U*1.96) ||
      (agAU_pLo$pLomu.A-agAU_pLo$pLose.A*1.96)>=(agAU_pLo$pLomu.U+agAU_pLo$pLose.U*1.96)  )
    agAU_pLo$sig[i] <- 1
}


agAU_pT <- agAU[,c("unit","pTmu.A","pTse.A","pTmu.U","pTse.U")]
agAU_pT$sig <- 0
for(i in 1:dim(agAU_pT)[1]) {
  if( (agAU_pT$pTmu.A+agAU_pT$pTse.A*1.96)<=(agAU_pT$pTmu.U-agAU_pT$pTse.U*1.96) ||
        (agAU_pT$pTmu.A-agAU_pT$pTse.A*1.96)>=(agAU_pT$pTmu.U+agAU_pT$pTse.U*1.96)  )
    agAU_pT$sig[i] <- 1
}

agAU_pS <- agAU[,c("unit","pSmu.A","pSse.A","pSmu.U","pSse.U")]
agAU_pS$sig <- 0
for(i in 1:dim(agAU_pS)[1]) {
  if( (agAU_pS$pSmu.A+agAU_pS$pSse.A*1.96)<=(agAU_pS$pSmu.U-agAU_pS$pSse.U*1.96) ||
        (agAU_pS$pSmu.A-agAU_pS$pSse.A*1.96)>=(agAU_pS$pSmu.U+agAU_pS$pSse.U*1.96)  )
    agAU_pS$sig[i] <- 1
}

agAU_pG <- agAU[,c("unit","pGmu.A","pGse.A","pGmu.U","pGse.U")]
agAU_pG$sig <- 0
for(i in 1:dim(agAU_pG)[1]) {
  if( (agAU_pG$pGmu.A+agAU_pG$pGse.A*1.96)<=(agAU_pG$pGmu.U-agAU_pG$pGse.U*1.96) ||
        (agAU_pG$pGmu.A-agAU_pG$pGse.A*1.96)>=(agAU_pG$pGmu.U+agAU_pG$pGse.U*1.96)  )
    agAU_pG$sig[i] <- 1
}

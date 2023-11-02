#'
#'@import dplyr
#'@import ggplot2
#'@import devtools
#'@import tidyr
#'@import readr
#'@import stringr
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
NULL

#' 1. Calculate spr, spr0, ypr by age ---
#' 
#' @export
cal.SPRetc3 <- function(age,waa,maa,M,F,waa.mid,maa.mid,init.num=1,plus=TRUE){
  if(min(age)==1){
    res <- data.frame(age=age,waa=waa,maa=maa,M=M,F=F,waa.mid=waa.mid,maa.mid=maa.mid)
  }
  else{
    if(min(age)>1){  # assume constant M for age < min(age)
      tmp <- min(age)-1
      res <- data.frame(age=1:max(age),waa=c(rep(0,tmp),waa),maa=c(rep(0,tmp),maa),
                        maa.mid=c(rep(0,tmp),maa.mid),waa.mid=c(rep(0,tmp),waa.mid),
                        M=c(rep(M[1],tmp),M),F=c(rep(F[1],tmp),F))
    }
    else{
      cat("age is less than 1!!")
    }
  }
  
  #res$waa.mid <- (res$waa+c(res$waa[-1],rev(res$waa)[1]))/2
  #res$maa.mid <- (res$maa+c(res$maa[-1],rev(res$maa)[1]))/2 
  
  res$remain.n <- res$remain.n0 <-  res$fished.n <- 0 
  
  if(isTRUE(plus)){ # for stocks with plus-group
  res$remain.n <- res$remain.n0 <-  res$fished.n <- 0 
  for(i in 1:nrow(res)-1){
      res$remain.n[i] <- exp(-sum(M[1:i])-sum(F[1:i])) 
      res$remain.n0[i] <- exp(-sum(M[1:i]))
  }
      res$remain.n[nrow(res)] <- exp(-sum(M[1:nrow(res)])-sum(F[1:nrow(res)])) * (1-exp(-((Inf-0)-(nrow(res)-2))*(M[nrow(res)]+F[nrow(res)])))/(1-exp(-M[nrow(res)]-F[nrow(res)])) 
	  
      res$remain.n0[nrow(res)] <- exp(-sum(M[1:nrow(res)])) * (1-exp(-((Inf-0)-(nrow(res)-2))*(M[nrow(res)]+0)))/(1-exp(-M[nrow(res)]-0))
  
  for(i in 2:nrow(res)-1){
      res$fished.n[i] <- #F[i]/(F[i]+M[i]) * (1-exp(-F[i]-M[i])) * res$remain.n[i-1]
        F[i]/(F[i]+M[i]) * (1-exp(-F[i]-M[i])) * exp(-sum(M[1:(i-1)])-sum(F[1:(i-1)]))
  }
      res$fished.n[nrow(res)] <- F[nrow(res)]/(F[nrow(res)]+M[nrow(res)]) * (1-exp(-F[nrow(res)]-M[nrow(res)])) *exp(-sum(M[1:nrow(res)])-sum(F[1:nrow(res)])) * (1-exp(-((Inf-0)-(nrow(res)-2))*(M[nrow(res)]+F[nrow(res)])))/(1-exp(-M[nrow(res)]-F[nrow(res)]))
  }
  else{# for stocks with NO plus-group
  for(i in 1:nrow(res)){
    res$remain.n[i] <- exp(-sum(M[1:i])-sum(F[1:i])) # survival number at age 生残率
    res$remain.n0[i] <- exp(-sum(M[1:i]))
  }
  # scaling of remain.n 
  #res$remain.n <- res$remain.n/res$remain.n[2]
  #res$remain.n0 <- res$remain.n0/res$remain.n0[2]
  
  for(i in 2:nrow(res)){
    res$fished.n[i] <- #F[i]/(F[i]+M[i]) * (1-exp(-F[i]-M[i])) * res$remain.n[i-1]
      F[i]/(F[i]+M[i]) * (1-exp(-F[i]-M[i])) * exp(-sum(M[1:(i-1)])-sum(F[1:(i-1)]))
  }
  }
  
  res$ssb.w <- res$remain.n * res$maa.mid * res$waa.mid # spr by age
  res$ssb0.w <- res$remain.n0 * res$maa.mid * res$waa.mid # spr0 by age
  res$fished.w <- res$fished.n * res$waa # ypr by age  
  
  return(res)
}

#---------------------------------------------------------
#' 2. Calculate spr, spr0, ypr, yield for the stock ---
#' 
#' @export
calYPR.simple2 <- function(waa,maa,M,F,Fmulti,age,waa.mid,maa.mid,SR.coef=1,is.ricker=F,steepness=1,plus){
  ypr <- spr <- spr0 <- rep(0,length(Fmulti))
  for(i in 1:length(Fmulti)){
    res <- cal.SPRetc3(waa=waa,maa=maa,M=M,F=F*Fmulti[i],age=age,waa.mid=waa.mid,maa.mid=maa.mid,plus=plus) 
    ypr[i] <- sum(res$fished.w,na.rm=T) # ypr for stock: equation 2 in paper
    spr[i] <- sum(res$ssb.w,na.rm=T)# spr for stock: equation 6.1 in paper 
    spr0[i] <- sum(res$ssb0.w,na.rm=T)# spr0 (no fishing) for stock 
  }
  
  # Beverton-holt
  if(!is.ricker){
    SperS0 <- 1-1/SR.coef[1]*(1-spr/spr0) # equation 5.1 in paper
    SperS0 <- ifelse(SperS0<0,NA,SperS0)
    RperR0 <- SperS0*(1-SR.coef[1]*(1-SperS0))^{-1}
  }
  else{
    # Ricker
    SperS0 <- 1-1/SR.coef[1]*log(spr0/spr) # equation 5.2 in paper
    SperS0 <- ifelse(SperS0<0,NA,SperS0)  
    RperR0 <- SperS0*exp(SR.coef[1]*(1-SperS0)) # equation 3.2 in paper
  }
  
  ab.yield <- ypr*RperR0 # equation 1 in paper
  
  return(list(ypr=ypr,spr=spr,spr0=spr0,RperR0=RperR0,
              ab.yield=ab.yield,SperS0=SperS0))
}

#----------------------------------------------------------------------
#' @export
my.normalize <- function(x,FUN2="mean",...){
  tmpfunc <- get(FUN2)
  sweep(x,2,apply(x,2,tmpfunc,...),FUN="/")
}

#'3. Further interface to get result of Clark (1991)
#' 
#' @export
get.clark.data2 <- function(bpara, SR.coef=c(0.750,0.875,0.938,1.386,2.079,2.773),is.ricker=c(F,F,F,T,T,T),f.vec=seq(from=0,to=1,by=0.01),is.graph=F,steepness){
  mat.tmp <- mat.relativeS <- RperR0 <- perspr <- matrix(0,length(f.vec),length(SR.coef))
  one_Abr <- rep(0,length(SR.coef))
  
  for(i in 1:length(SR.coef)){
    res <- calYPR.simple2(waa=bpara$waa,maa=bpara$maa,M=bpara$M,F=bpara$F,age=bpara$Age,
                          waa.mid=bpara$waa.mid,maa.mid=bpara$maa.mid,
                          Fmulti=f.vec,SR.coef=SR.coef[i],is.ricker=is.ricker[i],steepness=steepness[i],plus=bpara$plus[1])
    mat.tmp[,i] <- res$ab.yield
    mat.relativeS[,i] <- res$SperS0
    RperR0[,i] <- res$RperR0
    one_Abr[i]<-steepness[i]
    perspr[,i]<-res$spr/res$spr0 # note: this does not depend on SR.coef
  }
  
  f.max<-f.vec[res$ypr==max(res$ypr)]
  
  # calc MMY
  x0 <- res$spr/res$spr0 
  
  x <- my.normalize(mat.tmp,FUN2="max",na.rm=T) #mat.tmp is abs.yield for each S-R, x is normalized by MSY (equation 7 in paper)
  
  x <- ifelse(is.na(x),0,x) #replace NA with 0
  
  x <- apply(x,1,min) #choose the minimum x among the SR.coef for each F
  
  spr.mmy <- x0[x==max(x)]  #%SPRmmy
  f.mmy <- f.vec[spr.mmy==x0]
  
  if(is.graph==T){
    par(mfrow=c(3,2)) 
    # YPR (panel A)
	
    plot(f.vec,res$ypr,type="l",xlab="F multiplier",ylab="YPR",cex.lab=2, cex.axis=2) 
    points(f.max,max(res$ypr),pch=1)
    legend("bottomright",pch=1,legend=paste("Fmax=",f.vec[res$ypr==max(res$ypr)]),cex=1.5)
	
    lty.tmp <- as.numeric(is.ricker)+1
    
    # F vs. Absolute yield (panel B)
    matplot(f.vec,mat.tmp,type="l",lty=lty.tmp,col=1,xlab="F multiplier",
            ylab="Absolute yield",cex.lab=2, cex.axis=2) 
    
    # Relative spawning biomass vs relative yield (panel C)
    matplot(mat.relativeS,my.normalize(mat.tmp,FUN2="max",na.rm=T),col=1,lty=lty.tmp,
            type="l",xlab="Relative spawning biomass",ylab="Relative yield",cex.lab=2, cex.axis=2) 
    
    # Relative spawning biomass per recruit vs relative yield (panel D)
    matplot(x0,my.normalize(mat.tmp,FUN2="max",na.rm=T),col=1,lty=lty.tmp,
            type="l",xlab="SPR/SPR0",ylab="Relative yield",cex.lab=2, cex.axis=2) 
    points(x0,x,type="b")
    lines(c(x0[x==max(x)],x0[x==max(x)]),c(0,max(x)),col="gray")
    text(spr.mmy,max(x)*0.2,paste("SPR/SPR0 at MMY =",round(spr.mmy,2)),pos=4,cex=2)
    
    # Relative spawning biomass per recruit & F (panel E)
    plot(f.vec,x0,type="l",xlab="F multiplier",ylab="SPR/SPR0",cex.lab=2, cex.axis=2)
    points(f.mmy,spr.mmy)
    abline(v=f.mmy,col="gray")
    text(f.mmy,spr.mmy*1.2,paste("F_MMY =",round(f.mmy,2)),pos=4,cex=2)
    lines(c(0,f.mmy),c(spr.mmy,spr.mmy),col="gray",lty=2)
    text(0,spr.mmy*0.7,paste("SPR_MMY =\n",round(spr.mmy,2)),pos=4,cex=2)    
    
    # Relative yiled vs Fishing mortality
    matplot(f.vec,my.normalize(mat.tmp,FUN2="max",na.rm=T),col=1,lty=lty.tmp,
            type="l",xlab="F multiplier",ylab="Relative yield",cex.lab=2, cex.axis=2) 
    abline(v=f.mmy,col="gray")
    text(f.mmy,0.1,paste("F_MMY =",round(f.mmy,2)),pos=4,cex=2)
    abline(v=mean(bpara$M),col="gray")
    text(mean(bpara$M),0.2,paste("F=M (mean)",round(mean(bpara$M),3)),pos=4,cex=2)
    abline(v=f.max,col="gray")
    text(f.max,0.3,paste("Fmax=",round(mean(f.max),2)),pos=4,cex=2)
    
    jpgfile<-paste(bpara$Stock,".jpg",sep="")
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 1200, 
             units = "px", 
             pointsize = 18,
             quality = 100,
             res=80,
             antialias="cleartype")
    
    dev.off()
  }
  
  # calculate %SPR_Fmax
  F_x0 <- as.data.frame(cbind(f.vec,x0))
  colnames(F_x0) <- c("F","spr_spr0")
  
  f.max <- signif(f.max,digits=2) 
  
  as.data.frame(Fmax_x0 <- F_x0 %>% filter(F==f.max))
  spr.max <- Fmax_x0$spr_spr0 #　%SPR_Fmax
  
  # calculate F and %SPR for each SR when MSY
  F_Ymax <- SPR_Ymax <- Ymax <- mmyY <- mmy_MSY_Yloss <- sango_MSY_Yloss <- yonjyu_MSY_Yloss <-rep(0,length(SR.coef))
  
  for(i in 1:length(SR.coef)){
    
    F_Y<-as.data.frame(cbind(f.vec,x0,mat.tmp[,i])) #x0 <- res$spr/res$spr0
    maxF_Yvec<- F_Y %>% filter(V3==max(V3, na.rm=TRUE))
    
    
    F_Ymax[i]<-maxF_Yvec$f.vec # Fmsy
    SPR_Ymax[i]<-maxF_Yvec$x0 # %SPRmsy
    Ymax[i]<-maxF_Yvec$V3 # MSY
    mmyF_Yvec<- F_Y %>% filter(f.vec==f.mmy) # MMY
    mmyY[i] <- mmyF_Yvec$V3
    mmy_MSY_Yloss[i] <- mmyF_Yvec$V3/maxF_Yvec$V3 # Ymmy/MSY
    
    F_Y$x0<-round((F_Y$x0),2)
    sango_MSY_Ylossvec <- F_Y %>% mutate(diffs=abs(F_Y$x0-0.35)) %>% filter(diffs==min(diffs)) # 35%SPR_Yield
    yonjyu_MSY_Ylossvec <- F_Y %>% mutate(diffs=abs(F_Y$x0-0.40)) %>% filter(diffs==min(diffs)) # 40%SPR_Yield
    sango_MSY_Yloss[i] <-sango_MSY_Ylossvec$V3/maxF_Yvec$V3 # 35%SPR_Yield/MSY
    yonjyu_MSY_Yloss[i] <-yonjyu_MSY_Ylossvec$V3/maxF_Yvec$V3 # 40%SPR_Yield/MSY
  }
  
  #invisible(list(c(f.max,f.mmy,spr.mmy,spr.max),cbind(f.vec,x,x0),res,mat.tmp,mat.relativeS,mat.relativeS,RperR0,one_Abr,tameshi))
  
  c(f.max,spr.mmy*100,f.mmy,F_Ymax,SPR_Ymax*100,mmy_MSY_Yloss,sango_MSY_Yloss,yonjyu_MSY_Yloss,steepness)
}
#----------------------------------------------------------------------------------------------------------
#' @export
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}






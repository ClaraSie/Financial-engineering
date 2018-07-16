###歐式################################################################################
data2022<-read.table("C:\\Users\\Clara\\Desktop\\2022.csv",col.names = "stock",sep = ",")
data2022before<-read.table("C:\\Users\\Clara\\Desktop\\2022before.csv",col.names = "beforestock",header = TRUE,sep = ",")

M<-100000  #面額
R<-M/6 #轉換比率
r<-1.035/100  #無風險利率
T<-c();returnRate<-c() 
for(i in 0:29){returnRate[i+1]<-(data2022before$beforestock[i+2]-data2022before$beforestock[i+1])/data2022before$beforestock[i+1]}
#sig<-0.17
sig<-sd(returnRate,na.rm = FALSE)*sqrt(252) #年化sig
EurCB<-c()
for(j in 1:30){ #轉換時間第67天
  T<-(252*3-66-j+1)/252
  d1<-(log(data2022$stock[j+66]*R/M)+(T*(r+0.5*sig^2)))/(sig*sqrt(T))
  d2<-d1-sig*sqrt(T)
  C0<-data2022$stock[j+66]*pnorm(d1)-M/R*exp(-r*T)*pnorm(d2)
  EurCB[j]<-M*exp(-r*T)+R*C0
}

###樹狀(美式)#########################################################################

k<-6 #轉換價格
TT<-(252*3-66)/252 #年數
N<-4 #期數
pu<-0.5
pd<-0.5
dt<-TT/N 
u<-exp((r-0.5*sig^2)*dt+sig*sqrt(dt))
d<-exp((r-0.5*sig^2)*dt-sig*sqrt(dt))
disc<-exp(-r*dt)

count<-N+1
St<-diag(c(rep(0,each=N+1)))
C<-diag(c(rep(0,each=N+2)))
TreeCB<-c()
for(x in 1:30){
  S<-data2022$stock[x+66]
  for(i in count:1){ 
    for(j in 1:i){ 
      St[j,i]<-S*(u^(i-j)*d^(j-1))
      C[j,i]<-max((St[j,i]-k),(disc*(pu*C[j,i+1]+pd*C[j+1,i+1])))
    }
  }
  TreeCB[x]<-M*exp(-r*TT)+R*C[1,1]
}

###最小平方蒙地卡羅法(美式)###########################################################

m<-5000 #模擬次數
nudt<-(r-0.5*sig^2)*dt
sigsdt<-sig*sqrt(dt)
MC_CB<-matrix(0,nrow = 30)
for(x in 1:30){
  lnS<-log(data2022$stock[66+x])
  for(i in 1:m){   #模擬m次
    lnSt1<-lnS
    lnSt2<-lnS
    for(j in 1:N){ #切N期
      error<-matrix(rnorm(m*N,mean = 0, sd=1),ncol = N)  #誤差項矩陣
      lnSt1<-lnSt1+nudt+(sigsdt*error)
      lnSt2<-lnSt2+nudt+(sigsdt*(-error))
    }
    St1<-matrix(exp(lnSt1),ncol = N,byrow = TRUE) #找出原來股價
    St2<-matrix(exp(lnSt2),ncol = N,byrow = TRUE)
    MC_St<-rbind(St1,St2)
    MC_C<-MC_St-k
    MC_C[which(MC_C<0)]<-0 #每日選擇權價格
  }

  payoff_sum<-0;conti_test<-c()
  MC_C1<-matrix(0,nrow = m,ncol = N)
  for(i in 3:1){
    X1<-MC_St[1:m,i]; Y<-MC_C[1:m,i+1]
    x1<-X1[-which(MC_C[1:m,i]==0)]     #St
    y<-Y[-which(MC_C[1:m,i]==0)]*disc  #C
    x2<-x1^2
    
    redata<-data.frame(y,x1,x2)
    re<-lm(formula = y~x1+x2, data = redata) #迴歸
    model<-data.frame(re$coefficients)       #x1x2
    conti_test<-c(model$re.coefficients[1]+model$re.coefficients[2]*X1+model$re.coefficients[3]*X1^2) 
    exer_test<-MC_C[1:m,i]
    contrast<-cbind(exer_test,conti_test)
    contrast[which(contrast[1:m,1]<contrast[1:m,2])]=0
    MC_C1[1:m,N]<-MC_C[1:m,N]     #time N
    MC_C1[1:m,i]<-contrast[1:m,1]
    MC_C1[which(MC_C1[1:m,i]!=0),(i+1):N]=0
  }
  for(i in 1:N){
    w<-which(MC_C1[1:m,i]!=0)
    payoff<-exp(-r*dt*(i-1))*MC_C[w,i]
    payoff_sum<-payoff_sum+sum(payoff)
  }
  MC_C0<-payoff_sum/m
  MC_CB[x,1]<-M*exp(-r*TT)+R*MC_C0  #蒙地卡羅法CB
}
Price<-data.frame(EurCB/1000,TreeCB/1000,MC_CB/1000)
write.table(Price,"C:\\Users\\Clara\\Desktop\\CBPrice.csv",sep = ",")  #輸出csv檔


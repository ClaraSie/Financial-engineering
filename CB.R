###�ڦ�################################################################################
data2022<-read.table("C:\\Users\\Clara\\Desktop\\2022.csv",col.names = "stock",sep = ",")
data2022before<-read.table("C:\\Users\\Clara\\Desktop\\2022before.csv",col.names = "beforestock",header = TRUE,sep = ",")

M<-100000  #���B
R<-M/6 #�ഫ��v
r<-1.035/100  #�L���I�Q�v
T<-c();returnRate<-c() 
for(i in 0:29){returnRate[i+1]<-(data2022before$beforestock[i+2]-data2022before$beforestock[i+1])/data2022before$beforestock[i+1]}
#sig<-0.17
sig<-sd(returnRate,na.rm = FALSE)*sqrt(252) #�~��sig
EurCB<-c()
for(j in 1:30){ #�ഫ�ɶ���67��
  T<-(252*3-66-j+1)/252
  d1<-(log(data2022$stock[j+66]*R/M)+(T*(r+0.5*sig^2)))/(sig*sqrt(T))
  d2<-d1-sig*sqrt(T)
  C0<-data2022$stock[j+66]*pnorm(d1)-M/R*exp(-r*T)*pnorm(d2)
  EurCB[j]<-M*exp(-r*T)+R*C0
}

###��(����)#########################################################################

k<-6 #�ഫ����
TT<-(252*3-66)/252 #�~��
N<-4 #����
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

###�̤p����X�a�dù�k(����)###########################################################

m<-5000 #��������
nudt<-(r-0.5*sig^2)*dt
sigsdt<-sig*sqrt(dt)
MC_CB<-matrix(0,nrow = 30)
for(x in 1:30){
  lnS<-log(data2022$stock[66+x])
  for(i in 1:m){   #����m��
    lnSt1<-lnS
    lnSt2<-lnS
    for(j in 1:N){ #��N��
      error<-matrix(rnorm(m*N,mean = 0, sd=1),ncol = N)  #�~�t���x�}
      lnSt1<-lnSt1+nudt+(sigsdt*error)
      lnSt2<-lnSt2+nudt+(sigsdt*(-error))
    }
    St1<-matrix(exp(lnSt1),ncol = N,byrow = TRUE) #��X��Ӫѻ�
    St2<-matrix(exp(lnSt2),ncol = N,byrow = TRUE)
    MC_St<-rbind(St1,St2)
    MC_C<-MC_St-k
    MC_C[which(MC_C<0)]<-0 #�C�����v����
  }

  payoff_sum<-0;conti_test<-c()
  MC_C1<-matrix(0,nrow = m,ncol = N)
  for(i in 3:1){
    X1<-MC_St[1:m,i]; Y<-MC_C[1:m,i+1]
    x1<-X1[-which(MC_C[1:m,i]==0)]     #St
    y<-Y[-which(MC_C[1:m,i]==0)]*disc  #C
    x2<-x1^2
    
    redata<-data.frame(y,x1,x2)
    re<-lm(formula = y~x1+x2, data = redata) #�j�k
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
  MC_CB[x,1]<-M*exp(-r*TT)+R*MC_C0  #�X�a�dù�kCB
}
Price<-data.frame(EurCB/1000,TreeCB/1000,MC_CB/1000)
write.table(Price,"C:\\Users\\Clara\\Desktop\\CBPrice.csv",sep = ",")  #��Xcsv��

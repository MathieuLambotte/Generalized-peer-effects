# This code is a companion of the paper of Lambotte, Mathieu, A Generalized Model of Peer Effects for Binary Outcomes: 
# Strategic Complementarity and Pressure to Conform (December 11, 2023). Available at SSRN: 
# http://dx.doi.org/10.2139/ssrn.4661240 
# The paper and the code have NOT been peer-reviewed yet. As the paper is currently under review,
# please do not share this code without informing me.

#Setting ####
#Some packages are not useful
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr,haven,fastDummies,margins,reshape,data.table,ggplot2)
packages<-c("dplyr","haven","fastDummies","margins","reshape","data.table","ggplot2")
lapply(packages, require, character.only = TRUE)
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#Load data ([Y,X])
X<-read.csv("fake_X.csv")
Y<-read.csv("fake_Y.csv")

#Load network matrix A
A<-read.csv("fake_A.csv")

#create the row normalized matrix G
G<-A/rowSums(A)
G<-ifelse(is.nan(G),0,G)

#if several subnetworks, create a spliting variable to associate each household in the data to a subnetwork
#df$network indicates to which network each household belongs
networks<-split(seq_len(nrow(df)),df$network)
#if X does not contain network dummy variables, create them
fastDummies::dummy_cols(df,remove_selected_columns=TRUE)

#create dummy for isolated or not
X$not_isolated<-ifelse(rowSums(A)>0,1,0)
X<-as.matrix(X)
Y<-as.matrix(Y)
# names of the individual variables, without the dummy not_isolated and the dummy for network fixed effects
k<-colnames(X)[!c(grepl("networks", colnames(X)))]
k<-k[-which(k=="not_isolated")]

#GXAY (local-aggregate peer effect) ####
dataset <- as.data.frame(cbind(Y, X[,!colnames(X) %in%c("not_isolated")],G%*%X[,k],(A%*%Y)))
colnames(dataset) <- c("y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"Ay")

#starting values are parameters from a model without endogenous peer effects and a zero for the endogenous peer effect
b1<-c(coef(glm(y~.-Ay,dataset,family=binomial(link='logit'))),0)

colnames(Y)<-"Y"
Mg1<-A%*%Y
w1<-A
w2<-G
RE_Y<-matrix(0,nrow(X),1)
flag=1
count=1
SSE_b=matrix(0,100,1)
RE=matrix(0,100,1)
B<-matrix(0,100,1)
#Iterative MLE
while(flag==1){
  for(i in 1:length(networks)){
    Xt<-X[networks[[i]],]
    Yt<-Y[networks[[i]]]
    w1t<-A[networks[[i]],networks[[i]]]
    w2t<-G[networks[[i]],networks[[i]]]
    Mg1t<-Mg1[networks[[i]]]
    mr<-nrow(Xt)
    euclid=1
    x0=matrix(1,mr,1)/4
    endog=matrix(0,mr,1)
    
    while (euclid>0.00000001){
      xx_new=logit2prob(Xt[,!colnames(Xt)%in%c("not_isolated")]%*%matrix(b1[2:(length(networks)+length(k))])+w2t%*%Xt[,k]%*%b1[(length(networks)+length(k)+1):(length(networks)+(2*length(k)))]+w1t%*%x0%*%b1[(length(networks)+(2*length(k))+1)]+b1[1]*matrix(1,nrow=mr,ncol=1)) #equation (8)
      
      euclid=max(sum((xx_new-x0)*(xx_new-x0)))
      x0=xx_new
    }
    
    endog=xx_new
    RE_Y[networks[[i]]]<-endog
    Mg1[networks[[i]]]<-w1t%*%endog
  }
  RE[count]=mean(RE_Y)
  datafold<-dataset
  dataf<-setNames(data.frame(cbind(Y,X[,!colnames(X)%in%c("not_isolated")],G%*%X[,k],Mg1)),c("Y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"A_Y"))
  formulaf<-as.formula(paste("Y", paste("~",paste("1+",paste(paste0(colnames(X[,!colnames(X) %in%c("not_isolated")]),sep="",collapse = "+"),paste(paste0(paste0("G_",colnames(X[,k])),sep="",collapse = "+")), "A_Y",sep="+")))))
  
  bnew_model<-glm(formulaf,data=dataf,family=binomial(link='logit'))
  B[count,]<-bnew_model$coefficients[c("A_Y")]
  SSE=t(bnew_model$coef-b1)%*%(bnew_model$coef-b1)
  if (is.na(SSE)){
    flag=0
  }  else if (SSE<0.00000001){
    flag=0
  }  else if (count>1000){
    flag=0
  } else {flag=1}
  
  if (flag==1){
    SSE_b[count]=SSE
    b1=(bnew_model$coefficients-b1)/(ceiling(count/2))+b1
    count=count+1
  }
}

print(RE[count])
bnew_model_GXAY<-bnew_model
naive_ME_GXAY<-margins(bnew_model_GXAY)
round(logLik(bnew_model_GXAY),2)
round(BIC(bnew_model_GXAY),2)
round(AIC(bnew_model_GXAY),2)

#SE
y<-Y
p<-predict(bnew_model_GXAY,type="response")
kappa<-bnew_model_GXAY$coefficients
graf<-matrix(0,length(kappa),length(kappa))
for(z in 1:length(networks)){
  #temp1 is the direct effect partial p_i/ partial kappa
  temp1<-cbind(matrix(1,nrow=nrow(X[networks[[z]],]),ncol=1),X[networks[[z]],!colnames(X) %in%c("not_isolated")],G[networks[[z]],networks[[z]]]%*%X[networks[[z]],k],Mg1[networks[[z]]])
  #temp2 is indirect effect through strategic complementarity
  temp2<-kappa['A_Y']*(p[networks[[z]]]*(1-p[networks[[z]]]))*A[networks[[z]],networks[[z]]]%*%temp1
  
  #temp4 is partial of loglikelihood / partial kappa
  temp4<-as.numeric(y[networks[[z]]]-p[networks[[z]]])*as.matrix(temp1+temp2)
  
  graf<-graf+t(temp4)%*%temp4*solve(length(networks)/(length(networks)-1))[,1]
}

#cluster adjustment cameron et miller Robust Inference with Clustered Data p8

msscore = graf/(nrow(X)) #(length(networks)*(nrow(X)-1))/((length(networks)-1)*(nrow(X)-(60)))
se=sqrt(diag(solve(msscore)/nrow(X)))
result_GXAY=cbind(kappa, se, 2*(1-pnorm(abs(kappa/se)))) #bnew$par/se=t-ratio, 2*1-pnorm(t)<0.05 iff t>1.96
write.csv(result_GXAY,"Results GXAY.csv")

result_GXAY=coeftemp(bnew_model,vcov=solve(msscore)/nrow(X))

#marginal effects
#ME Average ME for theta, beta and psi can be obtained with margins(bnew_model_AYGY)
#but ME for gamma are direct (from margins()) + indirect
ME_theta<-marginal_effects(bnew_model_GXAY)[(length(k)+length(networks)):(2*length(k)+length(networks)-1)]
ME_psi<-marginal_effects(bnew_model_GXAY)[(2*length(k)+length(networks))]
ME_gamma_direct_all<-marginal_effects(bnew_model_GXAY)[1:(length(k)+length(networks)-1)]
ME_gamma_direct<-marginal_effects(bnew_model_GXAY)[c(1:length(k))]

temp1<-matrix(0,nrow(X),length(k))
temp2<-matrix(0,nrow(X),length(k))
temp3<-matrix(0,nrow(X),length(k))
for(z in 1:length(networks)){
  for(i in 1:nrow(A[networks[[z]],networks[[z]]])){ 
    temp1[networks[[z]],]<-as.matrix(diag(length(networks[[z]]))[,i]*ME_gamma_direct[networks[[z]],]+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],])
    #\mathbb{1}_{l=i}(\hat{\gamma}_k)+  g_{li}\hat{\theta}_k
    temp2[networks[[z]],]<-as.matrix(p[networks[[z]]]*(1-p[networks[[z]]])*(kappa['A_Y']*A[networks[[z]],networks[[z]]]%*%temp1[networks[[z]],])+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],]) 
    #[p^*_j (1-p^*_j)(\hat{\psi}A_j \times temp1+\hat{\beta}G_j \times temp1]+g_{ji}\hat{\theta}_k  => big accolade 
    temp3[networks[[z]],][i,]<-((kappa['A_Y']*p[networks[[z]]][i]*(1-p[networks[[z]]][i])*A[networks[[z]],networks[[z]]])%*%temp2[networks[[z]],])[i,]
    #\biggl(\hat{\psi}p^*_i(1-p^*_i)A_i+ \hat{\beta}p^*_i(1-p^*_i)G_i\biggr) \times big accolade (temp2)
  }
} 
ME_gamma_indirect<-data.frame(temp3)
colnames(ME_gamma_indirect)<-paste0(colnames(X[,k]),"_indirect")

ME_gamma_total<-ME_gamma_indirect+ME_gamma_direct
colnames(ME_gamma_total)<-paste0(colnames(X[,k]),"_total")

ME_all_GXAY<-c(colMeans(ME_gamma_direct_all),colMeans(ME_gamma_total),colMeans(ME_theta),colMeans(ME_psi))
write.csv(ME_all_GXAY,"Marginal Effects GXAY.csv")

#GXGY (local-average peer effect) ####
dataset <- as.data.frame(cbind(Y, X[,!colnames(X) %in%c("not_isolated")],G%*%X[,k],X[,colnames(X) %in%c("not_isolated")]*(G%*%Y-0.5)))
colnames(dataset) <- c("y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"Gy")
#starting values are parameters from a model without endogenous peer effects and zeros for endogenous peer effects
b1<-c(coef(glm(y~.-Gy,dataset,family=binomial(link='logit'))),0)


colnames(Y)<-"Y"
Mg2<-X[,colnames(X) %in%c("not_isolated")]*(G%*%Y-0.5)
w2<-G
RE_Y<-matrix(0,nrow(X),1)
flag=1
count=1
SSE_b=matrix(0,100,1)
RE=matrix(0,100,1)
B<-matrix(0,100,2)
#Iterative MLE
while(flag==1){
  for(i in 1:length(networks)){
    Xt<-X[networks[[i]],]
    Yt<-Y[networks[[i]]]
    w2t<-G[networks[[i]],networks[[i]]]
    Mg2t<-Mg2[networks[[i]]]
    mr<-nrow(Xt)
    euclid=1
    x0=matrix(1,mr,1)/4
    endog=matrix(0,mr,1)
    
    while (euclid>0.00000001){
      xx_new=logit2prob(Xt[,!colnames(Xt) %in%c("not_isolated")]%*%matrix(b1[2:(length(networks)+length(k))])+w2t%*%Xt[,k]%*%b1[(length(networks)+length(k)+1):(length(networks)+(2*length(k)))]+Xt[,colnames(Xt) %in%c("not_isolated")]*(w2t%*%x0-0.5)%*%b1[(length(networks)+(2*length(k))+1)]+b1[1]*matrix(1,nrow=mr,ncol=1)) #equation (8)
      
      euclid=max(sum((xx_new-x0)*(xx_new-x0)))
      x0=xx_new
    }
    
    endog=xx_new
    RE_Y[networks[[i]]]<-endog
    Mg2[networks[[i]]]<-(w2t%*%endog-0.5) #update Mg
  }
  RE[count]=mean(RE_Y)
  datafold<-dataset
  dataf<-setNames(data.frame(cbind(Y,X[,!colnames(X) %in%c("not_isolated")],G%*%X[,k],X[,colnames(X) %in%c("not_isolated")]*Mg2)),c("Y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"G_Y"))
  formulaf<-as.formula(paste("Y", paste("~",paste("1+",paste(paste0(colnames(X[,!colnames(X) %in%c("not_isolated")]),sep="",collapse = "+"),paste(paste0(paste0("G_",colnames(X[,k])),sep="",collapse = "+")), "G_Y",sep="+")))))
  
  bnew_model<-glm(formulaf,data=dataf,family=binomial(link='logit'))
  B[count,]<-bnew_model$coefficients[c("A_Y","G_Y")]
  SSE=t(bnew_model$coef-b1)%*%(bnew_model$coef-b1)
  if (is.na(SSE)){
    flag=0
  }  else if (SSE<0.00000001){
    flag=0
  }  else if (count>1000){
    flag=0
  } else {flag=1}
  
  if (flag==1){
    SSE_b[count]=SSE
    b1=(bnew_model$coefficients-b1)/(ceiling(count/2))+b1
    count=count+1
  }
}

print(RE[count])
bnew_model_GXGY<-bnew_model
naive_ME_GXGY<-margins(bnew_model_GXGY)
round(logLik(bnew_model_GXGY),2)
round(BIC(bnew_model_GXGY),2)
round(AIC(bnew_model_GXGY),2)

#SE

y<-Y
p<-predict(bnew_model_GXGY,type="response")
kappa<-bnew_model_GXGY$coefficients
graf<-matrix(0,length(kappa),length(kappa))
for(z in 1:length(networks)){
  #temp1 is the direct effect partial p_i/ partial kappa
  temp1<-cbind(matrix(1,nrow=nrow(X[networks[[z]],]),ncol=1),X[networks[[z]],!colnames(X) %in%c("not_isolated")],G[networks[[z]],networks[[z]]]%*%X[networks[[z]],k],X[networks[[z]],colnames(X) %in%c("not_isolated")]*Mg2[networks[[z]]])
  #temp3 is indirect effect through conformity
  temp3<-kappa['G_Y']*(p[networks[[z]]]*(1-p[networks[[z]]]))*G[networks[[z]],networks[[z]]]%*%temp1
  #temp4 is partial of loglikelihood / partial kappa
  temp4<-as.numeric(y[networks[[z]]]-p[networks[[z]]])*as.matrix(temp1+temp3)
  
  graf<-graf+t(temp4)%*%temp4*solve(length(networks)/(length(networks)-1))[,1]
}

#cluster adjustment cameron et miller Robust Inference with Clustered Data p8

msscore = graf/(nrow(X)) *(length(networks)*(nrow(X)-1))/((length(networks)-1)*(nrow(X)-(60)))
se=sqrt(diag(solve(graf)/nrow(X)))
result_GXGY=cbind(kappa, se, 2*(1-pnorm(abs(kappa/se)))) #bnew$par/se=t-ratio, 2*1-pnorm(t)<0.05 iff t>1.96
write.csv(result_GXGY,"Results GXGY.csv")

#marginal effects
#ME Average ME for theta, beta and psi can be obtained with margins(bnew_model_AYGY)
#but ME for gamma are direct (from margins()) + indirect
ME_theta<-marginal_effects(bnew_model_GXAY)[(length(k)+length(networks)):(2*length(k)+length(networks)-1)]
ME_beta<-marginal_effects(bnew_model_GXAY)[(2*length(k)+length(networks))]
ME_gamma_direct_all<-marginal_effects(bnew_model_GXAY)[1:(length(k)+length(networks)-1)]
ME_gamma_direct<-marginal_effects(bnew_model_GXAY)[c(1:length(k))]

temp1<-matrix(0,nrow(X),length(k))
temp2<-matrix(0,nrow(X),length(k))
temp3<-matrix(0,nrow(X),length(k))
for(z in 1:length(networks)){
  for(i in 1:nrow(A[networks[[z]],networks[[z]]])){ 
    temp1[networks[[z]],]<-as.matrix(diag(length(networks[[z]]))[,i]*ME_gamma_direct[networks[[z]],]+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],])
    #\mathbb{1}_{l=i}(\hat{\gamma}_k)+  g_{li}\hat{\theta}_k
    temp2[networks[[z]],]<-as.matrix(p[networks[[z]]]*(1-p[networks[[z]]])*(kappa['G_Y']*A[networks[[z]],networks[[z]]]%*%temp1[networks[[z]],])+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],]) 
    #[p^*_j (1-p^*_j)(\hat{\psi}A_j \times temp1+\hat{\beta}G_j \times temp1]+g_{ji}\hat{\theta}_k  => big accolade 
    temp3[networks[[z]],][i,]<-((kappa['G_Y']*p[networks[[z]]][i]*(1-p[networks[[z]]][i])*G[networks[[z]],networks[[z]]])%*%temp2[networks[[z]],])[i,]
    #\biggl(\hat{\psi}p^*_i(1-p^*_i)A_i+ \hat{\beta}p^*_i(1-p^*_i)G_i\biggr) \times big accolade (temp2)
  }
} 
ME_gamma_indirect<-data.frame(temp3)
colnames(ME_gamma_indirect)<-paste0(colnames(X[,k]),"_indirect")

ME_gamma_total<-ME_gamma_indirect+ME_gamma_direct
colnames(ME_gamma_total)<-paste0(colnames(X[,k]),"_total")

ME_all_GXGY<-c(colMeans(ME_gamma_direct_all),colMeans(ME_gamma_total),colMeans(ME_theta),colMeans(ME_beta))
write.csv(ME_all_GXGY,"Marginal Effects GXGY.csv")

#GXGYAY (generalized peer effect) ####
dataset <- as.data.frame(cbind(Y, X[,!colnames(X) %in%c("not_isolated")],G%*%X[,k],A%*%Y,X[,colnames(X) %in%c("not_isolated")]*(G%*%Y-0.5)))
colnames(dataset) <- c("y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"Ay","Gy")
#starting values are parameters from a model without endogenous peer effects and zeros for endogenous peer effects
b1<-c(coef(glm(y~.-Ay-Gy,dataset,family=binomial(link='logit'))),0,0)

colnames(Y)<-"Y"
Mg1<-A%*%Y
Mg2<-X[,colnames(X) %in%c("not_isolated")]*(G%*%Y-0.5)
w1<-A
w2<-G
RE_Y<-matrix(0,nrow(X),1)
flag=1
count=1
SSE_b=matrix(0,100,1)
RE=matrix(0,100,1)
B<-matrix(0,100,2)
#Iterative MLE
while(flag==1){
  for(i in 1:length(networks)){
    Xt<-X[networks[[i]],]
    Yt<-Y[networks[[i]]]
    w1t<-A[networks[[i]],networks[[i]]]
    w2t<-G[networks[[i]],networks[[i]]]
    Mg1t<-Mg1[networks[[i]]]
    Mg2t<-Mg2[networks[[i]]]
    mr<-nrow(Xt)
    euclid=1
    x0=matrix(1,mr,1)/4
    endog=matrix(0,mr,1)
    
    while (euclid>0.00000001){
      xx_new=logit2prob(Xt[,!colnames(Xt) %in%c("not_isolated")]%*%matrix(b1[2:(length(networks)+length(k))])+w2t%*%Xt[,k]%*%b1[(length(networks)+length(k)+1):(length(networks)+(2*length(k)))]+w1t%*%x0%*%b1[(length(networks)+(2*length(k))+1)]+Xt[,colnames(Xt) %in%c("not_isolated")]*(w2t%*%x0-0.5)%*%b1[(length(networks)+(2*length(k))+2)]+b1[1]*matrix(1,nrow=mr,ncol=1)) #equation (8)
      
      euclid=max(sum((xx_new-x0)*(xx_new-x0)))
      x0=xx_new
    }
    
    endog=xx_new
    RE_Y[networks[[i]]]<-endog
    Mg1[networks[[i]]]<-w1t%*%endog
    Mg2[networks[[i]]]<-(w2t%*%endog-0.5) #update Mg
  }
  RE[count]=mean(RE_Y)
  datafold<-dataset
  dataf<-setNames(data.frame(cbind(Y,X[,!colnames(X) %in%c("not_isolated")],G%*%X[,k],Mg1,X[,colnames(X) %in%c("not_isolated")]*Mg2)),c("Y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"A_Y","G_Y"))
  formulaf<-as.formula(paste("Y", paste("~",paste("1+",paste(paste0(colnames(X[,!colnames(X) %in%c("not_isolated")]),sep="",collapse = "+"),paste(paste0(paste0("G_",colnames(X[,k])),sep="",collapse = "+")), "A_Y+G_Y",sep="+")))))
  
  bnew_model<-glm(formulaf,data=dataf,family=binomial(link='logit'))
  B[count,]<-bnew_model$coefficients[c("A_Y","G_Y")]
  SSE=t(bnew_model$coef-b1)%*%(bnew_model$coef-b1)
  if (is.na(SSE)){
    flag=0
  }  else if (SSE<0.00000001){
    flag=0
  }  else if (count>1000){
    flag=0
  } else {flag=1}
  
  if (flag==1){
    SSE_b[count]=SSE
    b1=(bnew_model$coefficients-b1)/(ceiling(count/2))+b1
    count=count+1
  }
}


print(RE[count])
bnew_model_GXAYGY<-bnew_model
naive_ME_GXAYGY<-margins(bnew_model_GXAYGY)
round(logLik(bnew_model_GXAYGY),2)
round(BIC(bnew_model_GXAYGY),3)
round(AIC(bnew_model_GXAYGY),3)

#graph RE
colnames(B)<-c('expression(psi)','expression(beta)')
dataB<-as.data.frame(B[c(1:count),])
dataB$iter <- c(1:nrow(dataB))
meltdataB<-reshape2::melt(dataB,id="iter")
line_colors <- c("royalblue", "forestgreen")
Paraplot<-ggplot(meltdataB,aes(x=iter,y=value,colour=variable,group=variable)) + geom_line()+ ylab("Estimated value")+
  scale_x_log10()+xlab("Iteration (log scale)")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_colour_manual(name = '',values = line_colors,
                      labels = expression(hat(psi),hat(beta),Delta~beta))+
  theme(legend.text=element_text(size=12))+theme(aspect.ratio=1/2) 
ggsave("Para evo.png", plot = Paraplot, width = 6, height = 4, dpi = 600)

dataRE<-as.data.frame(RE[c(2:count)])
colnames(dataRE)<-"Average rational \n expectations equilibrium"
dataRE$iter <- c(1:nrow(dataRE))
meltdataRE<-reshape2::melt(dataRE,id="iter")
REplot<-ggplot(meltdataRE,aes(x=iter,y=value,colour=variable,group=variable)) + geom_line()+ ylab("P")+
  scale_x_log10()+xlab("Iteration (log scale)")+theme_bw()+scale_color_manual(name="",values='darkred')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.text=element_text(size=12))+theme(aspect.ratio=1/2) 
ggsave("Re evo.png", plot = REplot, width = 6, height = 4, dpi = 600)

#SE

y<-Y
p<-predict(bnew_model_GXAYGY,type="response")
kappa<-bnew_model_GXAYGY$coefficients
graf<-matrix(0,length(kappa),length(kappa))
for(z in 1:length(networks)){
  #temp1 is the direct effect partial p_i/ partial kappa
  temp1<-cbind(matrix(1,nrow=nrow(X[networks[[z]],]),ncol=1),X[networks[[z]],!colnames(X) %in%c("not_isolated")],G[networks[[z]],networks[[z]]]%*%X[networks[[z]],k],Mg1[networks[[z]]],X[networks[[z]],colnames(X) %in%c("not_isolated")]*Mg2[networks[[z]]])
  #temp2 is indirect effect through strategic complementarity
  temp2<-kappa['A_Y']*(p[networks[[z]]]*(1-p[networks[[z]]]))*A[networks[[z]],networks[[z]]]%*%temp1
  #temp3 is indirect effect through conformity
  temp3<-kappa['G_Y']*(p[networks[[z]]]*(1-p[networks[[z]]]))*G[networks[[z]],networks[[z]]]%*%temp1
  #temp4 is partial of loglikelihood / partial kappa
  temp4<-as.numeric(y[networks[[z]]]-p[networks[[z]]])*as.matrix(temp1+temp2+temp3)
  
  graf<-graf+t(temp4)%*%temp4*solve(length(networks)/(length(networks)-1))[,1]
}

#cluster adjustment cameron et miller Robust Inference with Clustered Data p8

msscore = graf/(nrow(X)) #(length(networks)*(nrow(X)-1))/((length(networks)-1)*(nrow(X)-(60)))
se=sqrt(diag(solve(msscore)/nrow(X)))
result_GXAYGY=cbind(kappa, se, 2*(1-pnorm(abs(kappa/se)))) #bnew$par/se=t-ratio, 2*1-pnorm(t)<0.05 iff t>1.96
write.csv(result_GXAYGY,"Results GXAYGY.csv")

#marginal effects
#ME Average ME for theta, beta and psi can be obtained with margins(bnew_model_AYGY)
#but ME for gamma are direct (from margins()) + indirect
ME_theta<-marginal_effects(bnew_model_GXAY)[(length(k)+length(networks)):(2*length(k)+length(networks)-1)]
ME_psi<-marginal_effects(bnew_model_GXAY)[(2*length(k)+length(networks))]
ME_beta<-marginal_effects(bnew_model_GXAY)[(2*length(k)+length(networks)+1)]
ME_gamma_direct_all<-marginal_effects(bnew_model_GXAY)[1:(length(k)+length(networks)-1)]
ME_gamma_direct<-marginal_effects(bnew_model_GXAY)[c(1:length(k))]

temp1<-matrix(0,nrow(X),length(k))
temp2<-matrix(0,nrow(X),length(k))
temp3<-matrix(0,nrow(X),length(k))
for(z in 1:length(networks)){
  for(i in 1:nrow(A[networks[[z]],networks[[z]]])){ 
    temp1[networks[[z]],]<-as.matrix(diag(length(networks[[z]]))[,i]*ME_gamma_direct[networks[[z]],]+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],])
    #\mathbb{1}_{l=i}(\hat{\gamma}_k)+  g_{li}\hat{\theta}_k
    temp2[networks[[z]],]<-as.matrix(p[networks[[z]]]*(1-p[networks[[z]]])*(kappa['A_Y']*A[networks[[z]],networks[[z]]]%*%temp1[networks[[z]],]+
                                                                              kappa['G_Y']*A[networks[[z]],networks[[z]]]%*%temp1[networks[[z]],])+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],]) 
    #[p^*_j (1-p^*_j)(\hat{\psi}A_j \times temp1+\hat{\beta}G_j \times temp1]+g_{ji}\hat{\theta}_k  => big accolade 
    temp3[networks[[z]],][i,]<-((kappa['A_Y']*p[networks[[z]]][i]*(1-p[networks[[z]]][i])*A[networks[[z]],networks[[z]]]+
                                   kappa['G_Y']*p[networks[[z]]][i]*(1-p[networks[[z]]][i])*G[networks[[z]],networks[[z]]])%*%temp2[networks[[z]],])[i,]
    #\biggl(\hat{\psi}p^*_i(1-p^*_i)A_i+ \hat{\beta}p^*_i(1-p^*_i)G_i\biggr) \times big accolade (temp2)
  }
} 
ME_gamma_indirect<-data.frame(temp3)
colnames(ME_gamma_indirect)<-paste0(colnames(X[,k]),"_indirect")

ME_gamma_total<-ME_gamma_indirect+ME_gamma_direct
colnames(ME_gamma_total)<-paste0(colnames(X[,k]),"_total")

ME_all_GXAYGY<-c(colMeans(ME_gamma_direct_all),colMeans(ME_gamma_total),colMeans(ME_theta),colMeans(ME_psi),colMeans(ME_beta))
write.csv(ME_all_GXAYGY,"Marginal Effects GXAYGY.csv")
#GXGY model of Lee et al. (2014) (includes isolated households and only a local average component) ####
dataset <- as.data.frame(cbind(Y, X[,!colnames(X) %in%c("not_isolated")],G%*%X[,k],G%*%Y))
colnames(dataset) <- c("y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"Gy")
#starting values are parameters from a model without endogenous peer effects and zeros for endogenous peer effects
b1<-c(coef(glm(y~.-Gy,dataset,family=binomial(link='logit'))),0)

colnames(Y)<-"Y"
Mg2<-(G%*%Y)
w2<-G
RE_Y<-matrix(0,nrow(X),1)
flag=1
count=1
SSE_b=matrix(0,100,1)
RE=matrix(0,100,1)
B<-matrix(0,100,1)
#Iterative MLE
while(flag==1){
  for(i in 1:length(networks)){
    Xt<-X[networks[[i]],]
    Yt<-Y[networks[[i]]]
    w2t<-G[networks[[i]],networks[[i]]]
    Mg1t<-Mg1[networks[[i]]]
    Mg2t<-Mg2[networks[[i]]]
    mr<-nrow(Xt)
    euclid=1
    x0=matrix(1,mr,1)/4
    endog=matrix(0,mr,1)
    while (euclid>0.00000001){
      xx_new=logit2prob(Xt[,!colnames(Xt) %in%c("not_isolated")]%*%matrix(b1[2:(length(networks)+length(k))])+
                          w2t%*%Xt[,k]%*%b1[(length(networks)+length(k)+1):(length(networks)+(2*length(k)))]+
                          (w2t%*%x0)%*%b1[(length(networks)+(2*length(k))+1)]+b1[1]*matrix(1,nrow=mr,ncol=1)) #equation (8)
      euclid=max(sum((xx_new-x0)*(xx_new-x0)))
      x0=xx_new
    }
    
    endog=xx_new
    RE_Y[networks[[i]]]<-endog
    Mg2[networks[[i]]]<-(w2t%*%endog) #update Mg
  }
  RE[count]=mean(RE_Y)
  datafold<-dataset
  dataf<-setNames(data.frame(cbind(Y,X[,!colnames(X) %in%c("not_isolated")],G%*%X[,k],Mg2)),c("Y",colnames(X[,!colnames(X) %in%c("not_isolated")]),paste0("G_",colnames(X[,k])),"G_Y"))
  formulaf<-as.formula(paste("Y", paste("~",paste("1+",paste(paste0(colnames(X[,!colnames(X) %in%c("not_isolated")]),sep="",collapse = "+"),paste(paste0(paste0("G_",colnames(X[,k])),sep="",collapse = "+")), "G_Y",sep="+")))))
  
  bnew_model<-glm(formulaf,data=dataf,family=binomial(link='logit'))
  B[count,]<-bnew_model$coefficients[c("G_Y")]
  SSE=t(bnew_model$coef-b1)%*%(bnew_model$coef-b1)
  if (is.na(SSE)){
    flag=0
  }  else if (SSE<0.00000001){
    flag=0
  }  else if (count>1000){
    flag=0
  } else {flag=1}
  
  if (flag==1){
    SSE_b[count]=SSE
    b1=(bnew_model$coefficients-b1)/(ceiling(count/2))+b1
    count=count+1
  }
}



print(RE[count])
bnew_model_GXGYLLL<-bnew_model
naive_ME_GXGYLLL<-margins(bnew_model_GXGYLLL)
round(logLik(bnew_model_GXGYLLL),2)
round(BIC(bnew_model_GXGYLLL),3)
round(AIC(bnew_model_GXGYLLL),3)

#SE

y<-Y
p<-predict(bnew_model_GXGYLLL,type="response")
kappa<-bnew_model_GXGYLLL$coefficients
graf<-matrix(0,length(kappa),length(kappa))
for(z in 1:length(networks)){
  #temp1 is the direct effect partial p_i/ partial kappa
  temp1<-cbind(matrix(1,nrow=nrow(X[networks[[z]],]),ncol=1),X[networks[[z]],!colnames(X) %in%c("not_isolated")],
               G[networks[[z]],networks[[z]]]%*%X[networks[[z]],k],Mg2[networks[[z]]])
  #temp3 is indirect effect through conformity
  temp3<-kappa['G_Y']*(p[networks[[z]]]*(1-p[networks[[z]]]))*G[networks[[z]],networks[[z]]]%*%temp1
  #temp4 is partial of loglikelihood / partial kappa
  temp4<-as.numeric(y[networks[[z]]]-p[networks[[z]]])*as.matrix(temp1+temp3)
  
  graf<-graf+t(temp4)%*%temp4*solve(length(networks)/(length(networks)-1))[,1]
}

#cluster adjustment cameron et miller Robust Inference with Clustered Data p8

msscore = graf/(nrow(X)) #(length(networks)*(nrow(X)-1))/((length(networks)-1)*(nrow(X)-(60)))
se=sqrt(diag(solve(msscore)/nrow(X)))
result_GXGYLLL=cbind(kappa, se, 2*(1-pnorm(abs(kappa/se)))) #bnew$par/se=t-ratio, 2*1-pnorm(t)<0.05 iff t>1.96
write.csv(result_GXGYLLL,"Results GXGYLLL.csv")

#marginal effects
#ME Average ME for theta, beta and psi can be obtained with margins(bnew_model_AYGY)
#but ME for gamma are direct (from margins()) + indirect
ME_theta<-marginal_effects(bnew_model_GXAY)[(length(k)+length(networks)):(2*length(k)+length(networks)-1)]
ME_beta<-marginal_effects(bnew_model_GXAY)[(2*length(k)+length(networks))]
ME_gamma_direct_all<-marginal_effects(bnew_model_GXAY)[1:(length(k)+length(networks)-1)]

ME_gamma_direct<-marginal_effects(bnew_model_GXAY)[c(1:length(k))]
temp1<-matrix(0,nrow(X),length(k))
temp2<-matrix(0,nrow(X),length(k))
temp3<-matrix(0,nrow(X),length(k))
for(z in 1:length(networks)){
  for(i in 1:nrow(A[networks[[z]],networks[[z]]])){ 
    temp1[networks[[z]],]<-as.matrix(diag(length(networks[[z]]))[,i]*ME_gamma_direct[networks[[z]],]+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],])
    #\mathbb{1}_{l=i}(\hat{\gamma}_k)+  g_{li}\hat{\theta}_k
    temp2[networks[[z]],]<-as.matrix(p[networks[[z]]]*(1-p[networks[[z]]])*(kappa['G_Y']*A[networks[[z]],networks[[z]]]%*%temp1[networks[[z]],])+G[networks[[z]],networks[[z]]][,i]*ME_theta[networks[[z]],]) 
    #[p^*_j (1-p^*_j)(\hat{\psi}A_j \times temp1+\hat{\beta}G_j \times temp1]+g_{ji}\hat{\theta}_k  => big accolade 
    temp3[networks[[z]],][i,]<-((kappa['G_Y']*p[networks[[z]]][i]*(1-p[networks[[z]]][i])*G[networks[[z]],networks[[z]]])%*%temp2[networks[[z]],])[i,]
    #\biggl(\hat{\psi}p^*_i(1-p^*_i)A_i+ \hat{\beta}p^*_i(1-p^*_i)G_i\biggr) \times big accolade (temp2)
  }
} 
ME_gamma_indirect<-data.frame(temp3)
colnames(ME_gamma_indirect)<-paste0(colnames(X[,k]),"_indirect")

ME_gamma_total<-ME_gamma_indirect+ME_gamma_direct
colnames(ME_gamma_total)<-paste0(colnames(X[,k]),"_total")

ME_all_GXGY<-c(colMeans(ME_gamma_direct_all),colMeans(ME_gamma_total),colMeans(ME_theta),colMeans(ME_beta))
write.csv(ME_all_GXGY,"Marginal Effects GXAYGYLLL.csv")




############# ROC functions #############
ROC<-function(res){  
  # 1st column is class has 0 and 1 only 
  # 2nd colum is their scores 
  ord<-order(res[,2],decreasing=T)
  score<-res[ord,2]
  class<-res[ord,1]
  
  temp1<-unique(score)
  n2<-length(temp1)
  n<-length(class)
  class0<-which(class==0)
  class1<-which(class==1)
  n1<-length(class1)
  n0<-length(class0)
  
  Sen<-rep(0,(n2+1)) #Sensitivity
  Spe<-rep(1,(n2+1)) #Specificity
  for (i in 1:n2){
    tmp1<-which(score>=temp1[i])
    tmp2<-setdiff(1:n,tmp1)
    Sen[(i+1)]<-length(intersect(tmp1,class1))/n1
    Spe[(i+1)]<-length(intersect(tmp2,class0))/n0   
  }
  out<-data.frame(Sen=Sen,Spe=Spe)
  out
}

ROC.score<-function(Sen,Spe){
  n<-length(Sen)-1
  tmp1<-1-Spe
  tmp2<-diff(tmp1,lag=1)
  tmp3<-rep(0,n)
  for (i in 1:n){ tmp3[i]<-(Sen[i]+Sen[(i+1)])/2 }
  out<-tmp3%*%tmp2
  out
}


##import data
temp1 <- read.csv("/home/alex/Documents/Stat project/Grant application/unimelb_training.csv",header = TRUE)
id    <- temp1[,1]
y     <- temp1[,2]
temp1 <- temp1[,-c(1,2)]


############# Distribution Analysis######
#cbind(colnames(temp1),seq(1:ncol(temp1)))

# 8708 rows, 252 observations
dim(temp1)

# colnames and dependent variable is Grant.Status
names(temp1) 

# distribution of response variable(almost balanced)
#    0    1 
# 4716 3992 
table(y)

############# Data cleaning #############

### Start date
# chron library has days/years function
library(chron)
date0 <- as.Date(c("01/01/04"),"%d/%m/%y")
date1 <- as.Date(temp1[,4],"%d/%m/%y")    
#date  <- as.numeric(date1-date0)%%365

month <- months(date1)
day <- days(date1)
year <- as.numeric(years(date1)) - as.numeric(years(date0))

daytable <- table(day, y) # check contigency table of day and y
monthtable <- table(month, y)
yeartable <- table(year, y)

chisq.test(daytable) # independent test
chisq.test(monthtable)
chisq.test(yeartable)

### Sponsor code
spon.code1 <- as.character(temp1[,1])
spon.code1[spon.code1==""]<-"999" #Set the missing value to sponsor code 999
spon.cate  <- unique(spon.code1)
spon.le    <- length(spon.cate)   #Number of sponsors
spon.numb  <- rep(0,spon.le)      #Number of applications for each sponsor
spon.numbs <- rep(0,spon.le)      #Number of success applications for each sponsor
spon.rate  <- rep(0,spon.le)      #Success rate for each sponsor
for (i in 1:spon.le){
  tmp1         <- which(spon.code1==spon.cate[i])
  spon.numb[i] <- length(tmp1)
  spon.numbs[i] <- length(which(y[tmp1]==1))
  spon.rate[i]  <- spon.numbs[i]/spon.numb[i]
}
spon.chis  <- rep(0,spon.le) #Find the p-value of chi-square test for each sponsor

for (i in 1:spon.le){
  if (spon.numb[i]<=2){
    spon.chis[i] <- NA
  }else {
    tmp1 <- which(spon.code1==spon.cate[i])
    tmp2 <- which(y[tmp1]==0)
    if ((length(tmp2)==0)||length(tmp2)==length(tmp1)){ 
      spon.chis[i] <- 0 
    } else{
      spon.chis[i] <- chisq.test(table(spon.code1[tmp1],y[tmp1]))$p.value
    }
  }
}

#Big sponsors
spon.cate1 <- spon.cate[which(spon.numb>100)] 

#Small and favor to success (collapse to one group)
spon.cate2 <- spon.cate[which((spon.numb<=100)&(spon.rate>.5)&(spon.chis<0.2))] # why compare to 0.2?

#Small and favor to failure (collapse to one group)
spon.cate3 <- spon.cate[which((spon.numb<=100)&(spon.rate<.5)&(spon.chis<0.2))]

#Small and no obvious preference or tiny sponsors (collapse to one group)
spon.cate4 <- spon.cate[which((spon.numb<=100)&(spon.chis>=0.2)|is.na(spon.chis))]

n <- nrow(temp1)
spon.code2 <- rep(NA,n)
for (i in 1:n){
  tmp1 <- which(spon.cate1==spon.code1[i])
  tmp2 <- length(which(spon.cate2==spon.code1[i]))
  tmp3 <- length(which(spon.cate3==spon.code1[i]))
  tmp4 <- length(which(spon.cate4==spon.code1[i]))
  if (length(tmp1)==1){
    spon.code2[i] <- spon.code1[i]
  } else{
    if (tmp2==1){
      spon.code2[i] <- "AAA"
    } else{
      if (tmp3==1){
        spon.code2[i] <- "BBB"
      } else{
        if (tmp4==1){spon.code2[i] <- "CCC"}
      }
    }
  }
}
spon.code <- as.factor(spon.code2)

### Grant code
grant.code1 <- as.character(temp1[,2])
grant.code1[grant.code1==""] <- "999" #999 as missing value
granttable <- table(grant.code1,y)
grantmargintable <- margin.table(granttable,1)
grantproptable <- prop.table(granttable,1)
cbind(granttable,total = grantmargintable, prop = grantproptable[,1])

#  0    1 total      prop
#10A 2478 1475  3953 0.6268657
#10B  251  196   447 0.5615213
#20A   60   68   128 0.4687500
#20C  202  205   407 0.4963145
#30A    0    3     3 0.0000000#
#30B  894  413  1307 0.6840092
#30C  150  208   358 0.4189944
#30D  126   52   178 0.7078652
#30E    6    5    11 0.5454545
#30F    1    0     1 1.0000000#
#30G   48   12    60 0.8000000
#40C    0    6     6 0.0000000#
#50A  337  600   937 0.3596585
#999  163  749   912 0.1787281

grant.code1[grant.code1=="30A"] <- "999"
grant.code1[grant.code1=="40C"] <- "999"
grant.code1[grant.code1=="30E"] <- "10B"  #different from Dr. Li  
grant.code1[grant.code1=="30F"] <- "30G"  #different from Dr. Li 
grant.code <- as.factor(grant.code1)

### Contract value
cvalue1 <- as.character(temp1[,3])
cvalue1[cvalue1==""] <- "999" #999 as missing value
contracttable <- table(cvalue1,y)

contractmargintable <- margin.table(contracttable,1)
contractproptable <- prop.table(contracttable,1)
cbind(contracttable,total = contractmargintable, prop = contractproptable[,1])
#  0    1 total      prop
#999 2913  650  3563 0.8175695
#A    875 1601  2476 0.3533926
#B    332  326   658 0.5045593
#C    166  284   450 0.3688889
#D    121  318   439 0.2756264
#E     69  244   313 0.2204473
#F     58  208   266 0.2180451
#G    103  294   397 0.2594458
#H     37   41    78 0.4743590#
#I     14    7    21 0.6666667
#J     23    3    26 0.8846154
#K      0    6     6 0.0000000#
#L      0    2     2 0.0000000#
#M      1    1     2 0.5000000#
#O      1    1     2 0.5000000#
#P      2    0     2 1.0000000
#Q      1    6     7 0.1428571#

cvalue1[which((cvalue1=="H ") | (cvalue1 == "M ") | (cvalue1 == "O "))]<-"B "   
cvalue1[which((cvalue1=="I ")|(cvalue1=="J ")|(cvalue1=="P "))]<-"999"
cvalue1[which((cvalue1=="K ")|(cvalue1=="L ")|(cvalue1=="Q "))]<-"F "
cvalue <- as.factor(cvalue1)

#####Number of success and failure####
nsuccess <- temp1[,34]
median <- summary(nsuccess)[3]
nsuccessind <- rep(0,length(nsuccess))
nsuccessind[which(is.na(nsuccess))] = 1
nsuccess[which(is.na(nsuccess))] = median

nfail    <- temp1[,35]
median <- summary(nfail)[3]
nfailind <- rep(0,length(nfail))
nfailind[which(is.na(nfail))] = 1
nfail[which(is.na(nfail))] = median

######### Has PhD#####
has.PhD <- temp1$With.PHD.1
levels(has.PhD)
has.PhDtable <- table(has.PhD,y)
has.PhDtable
chisq.test(has.PhDtable)

######## Country ####
country <- temp1$Country.of.Birth.1
levels(country)
countrytable <- table(country,y)
countrytable
chisq.test(countrytable)

####### Has.ID ###
has.ID <- temp1$Person.ID.1
has.ID[which(!is.na(has.ID))] = 1 # is.na retruns a logical vector
has.ID[which(is.na(has.ID))] = 0 # have to handle non missing value first
has.IDtable <- table(has.ID,y)
has.IDtable
chisq.test(has.IDtable)

###### Role ###
role <- temp1$Role.1
levels(role)
roletable <- table(role,y)
role.margin.table <- margin.table(roletable,1)
role.prop.table <- prop.table(roletable,1)
cbind(roletable,sum = role.margin.table, prop0 = role.prop.table[,1])
chisq.test(roletable)


role[which(role == "")] = "DELEGATED_RESEARCHER"  # Merge certain levels
role[which(role == "EXTERNAL_ADVISOR")] = "STUD_CHIEF_INVESTIGATOR"
role[which(role == "HONVISIT")] = "STUD_CHIEF_INVESTIGATOR"

# Department No.
dept <- as.factor(temp1$Dept.No..1)

###### Papers#########
Astar <- temp1$A..1
summary(Astar)
# range(Astar, na.rm = TRUE)
table(Astar,y)
median <- summary(Astar)[3]
Astarind <- rep(0,length(Astar))
Astarind[which(is.na(Astar))] = 1
Astar[which(is.na(Astar))] = median

A <- temp1$A.1
summary(A)
table(A,y)
median <- summary(A)[3]
Aind <- rep(0,length(A))
Aind[which(is.na(A))] = 1
A[which(is.na(A))] = median

B <- temp1$B.1
summary(B)
table(B,y)
median <- summary(B)[3]
Bind <- rep(0,length(B))
Bind[which(is.na(B))] = 1
B[which(is.na(B))] = median

C <- temp1$C.1
summary(C)
table(C,y)
median <- summary(C)[3]
Cind <- rep(0,length(C))
Cind[which(is.na(C))] = 1
C[which(is.na(C))] = median

A..col <- which(colnames(temp1) == "A..1" | colnames(temp1) == "A..2" |
                  colnames(temp1) == "A..3" | colnames(temp1) == "A..4" | 
                  colnames(temp1) == "A..5" | colnames(temp1) == "A..6" | 
                  colnames(temp1) == "A..7" | colnames(temp1) == "A..8" | 
                  colnames(temp1) == "A..9" | colnames(temp1) == "A..10" | 
                  colnames(temp1) == "A..11" | colnames(temp1) == "A..12" | 
                  colnames(temp1) == "A..13" | colnames(temp1) == "A..14" | 
                  colnames(temp1) == "A..15")
A.. <- temp1[,A..col]
maxAstar <- apply(A..,1,function(x) max(x,na.rm = TRUE))
table(maxAstar,y)

median <- summary(maxAstar)[3]
maxAstarind <- rep(0,length(maxAstar))
maxAstarind[which(maxAstar == -Inf)] = 1
maxAstar[which(maxAstar == -Inf)] = median


C.col <- which(colnames(temp1) == "C.1" | colnames(temp1) == "C.2" |
                 colnames(temp1) == "C.3" | colnames(temp1) == "C.4" | 
                 colnames(temp1) == "C.5" | colnames(temp1) == "C.6" | 
                 colnames(temp1) == "C.7" | colnames(temp1) == "C.8" | 
                 colnames(temp1) == "C.9" | colnames(temp1) == "C.10" | 
                 colnames(temp1) == "C.11" | colnames(temp1) == "C.12" | 
                 colnames(temp1) == "C.13" | colnames(temp1) == "C.14" | 
                 colnames(temp1) == "C.15")
C. <- temp1[,C.col]
maxC <- apply(C.,1,function(x) max(x,na.rm = TRUE))
table(maxC,y)

median <- summary(maxC)[3]
maxCind <- rep(0,length(maxC))
maxCind[which(maxC == -Inf)] = 1
maxC[which(maxC == -Inf)] = median

####Num of people involved ####
personcol <- which(colnames(temp1) == "Role.1" | colnames(temp1) == "Role.2" |
                     colnames(temp1) == "Role.3" | colnames(temp1) == "Role.4" | 
                     colnames(temp1) == "Role.5" | colnames(temp1) == "Role.6" | 
                     colnames(temp1) == "Role.7" | colnames(temp1) == "Role.8" | 
                     colnames(temp1) == "Role.9" | colnames(temp1) == "Role.10" | 
                     colnames(temp1) == "Role.11" | colnames(temp1) == "Role.12" | 
                     colnames(temp1) == "Role.13" | colnames(temp1) == "Role.14" | 
                     colnames(temp1) == "Role.15")
persons <- temp1[,personcol]
numPeople <- rowSums(persons != "", na.rm = TRUE)
numPeopletable <- table(numPeople,y)
numPeopletable
chisq.test(numPeopletable)



data1 <- data.frame(y=y,day = day, month = month, year = year,sponsor=spon.code,grant=grant.code,
                    cvalue=cvalue,nsuccess=nsuccess,nsuccessind = nsuccessind, nfail=nfail,nfailind = nfailind,
                    has.PhD = has.PhD, country = country, has.ID = has.ID, role = role, Astar = Astar,
                    Astarind = Astarind, A = A, Aind = Aind, B = B, Bind = Bind, C = C, Cind = Cind,
                    maxAstar = maxAstar, maxAstarind = maxAstarind, maxC = maxC, maxCind = maxCind,
                    numPeople = numPeople)

########## randomForest ###############
data1$y <- as.factor(data1$y)
library(randomForest)
library(caret)
library(pROC)

#tuning mtry
#tuneRF(x = data1[tra,2:ncol(data1)], y = data1[tra,1],trace = TRUE, ntreeTry = 1000, stepFactor = 1.5, improve = 0.0000001, plot = TRUE, doBest = FALSE,)
#mtry  OOBError
#3.OOB    3 0.1085714
#4.OOB    4 0.1057143
#5.OOB    5 0.1078571
#7.OOB    7 0.1082857

set.seed(1) #set the random seed, so that you can repeat the same results
ind <- sample(1:nrow(data1),size=nrow(data1),replace=F)
tra <- ind[1:7000]
val <- ind[7001:nrow(data1)]

gbmGrid <- expand.grid(mtry = c(3,4,5,6))
nrow(gbmGrid)

#fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
ptm <- proc.time() 
set.seed(2)
rfFit2 <- train(x = data1[tra,2:ncol(data1)], y = data1[tra,1],method = "rf", tuneGrid = gbmGrid)
(proc.time()-ptm)/60 # 19 mins
rfFit2
#plot(rfFit2,cex.axis = 5)
plot(rfFit2,cex.axis = 5, xlab = "Number of Randomly Selected Predictors")

ptm <- proc.time() 
set.seed(1)
newdata1_A.rf <- randomForest (x = data1[tra,2:ncol(data1)], y = data1[tra,1], mtry = 4, importance = TRUE)
(proc.time()-ptm)/60
print(newdata1_A.rf) # OOB: 10.73%
#plot(newdata1_A.rf, main = "Random Forest Error Rate VS Number of Trees")

#plot(newdata1_A.rf, xlab = "Trees", main = NULL)
par(mar = c(3.7,5,0.7,2), mgp = c(2.5,0.8,0))
matplot(1:newdata1_A.rf$ntree, newdata1_A.rf$err.rate, type = "l", xlab = "Trees", ylab = "Error")
legend("topright", inset = 0.1, c("Error Rate", "5% CI Upper Bound", "5% CI Lower Bound"), lty = c(1,2,2), col = c("black", "red","green"),
       bty = "n")

pred1 <- predict(newdata1_A.rf,data1[val,2:ncol(data1)],type = "prob")[,2]
pred2<- as.numeric(pred1>=0.5)
true <- data1[val,1]
temp3 <- cbind.data.frame(pred2,true)
error <- (nrow(temp3)-length(which(temp3[,1]==temp3[,2])))/nrow(temp3)
error #when mtry = 4, error is 0.1018735;
#when mtry = 5, error is 0.09894614

#temp<-ROC(cbind.data.frame(data1[val,1],pred1)) #factor 0/1 converts to 1/2 in cbind
#as.vector(ROC.score(temp$Sen,temp$Spe)) #0.9409988
rf.ROC <- roc(data1[val,1],pred1)
par(las = 1,mgp = c(3,1,0),mar = c(5.1,5.1,4.1,2.1), oma = c(1,2,0,0))
dev.off()
plot(rf.ROC, print.auc = FALSE, col = "red", print.thres.adj=c(1,-0.5), print.auc.adj=c(-1.5,5),
     cex = 0.7, ylab = NA)
##used to adjust the distance between axis labels and tick mark labels
mtext(side = 2, "Sensitivity",las = 0, line = 2.3)

#plot(rf.ROC,asp = NA, pin = c(2,2))

#variable importance
#par(las = 1,mgp = c(2,3,0),mar = c(5,6,2,2))
#varImpPlot(newdata1_A.rf, sort = TRUE, main = NULL, type = 1,cex = 0.6)
imp <- importance(newdata1_A.rf, type = 1)
imp.order <- imp[order(imp[,1], decreasing = TRUE),]
#par(las = 1,mgp = c(3,1,0),mar = c(5.1,4.1,4.1,2.1))
#dotchart(imp.order, labels = names(imp.order), xlab = "Relative Importance", cex = 0.4)
library(Hmisc)
par(mar = c(3,2,0.5,2), mgp = c(2.5,0.8,0))
dotchart2(imp.order, labels = names(imp.order), xlab = "Relative Importance of RF", 
          width.factor = 3, dotsize =2, pch = 21, cex = 0.75)

##GBM
library(gbm)
set.seed(1) #set the random seed, so that you can repeat the same results
ind <- sample(1:nrow(data1),size=nrow(data1),replace=F)
tra <- ind[1:7000]
val <- ind[7001:nrow(data1)]

####tuning parameters ####
data1$y <- as.factor(data1$y)
library(caret)
library(e1071)
gbmGrid <- expand.grid(interaction.depth = c(6,8,10), n.trees = (15:20)*50,shrinkage = c(0.02,0.03))
nrow(gbmGrid)

fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
ptm <- proc.time() 
set.seed(2)
gbmFit2 <- train(x = data1[tra,2:ncol(data1)], y = data1[tra,1],method = "gbm", verbose = FALSE,
                 tuneGrid = gbmGrid)
(proc.time()-ptm)/60
gbmFit2

#The final values used for the model were n.trees
#= 200, interaction.depth = 9 and shrinkage = 0.01.

#interaction.depth = seq(1,10,1), n.trees = (1:20)*50,shrinkage = seq(0,0.01,0.001)
#The final values used for the model were n.trees
#= 1000, interaction.depth = 10 and shrinkage = 0.01. 

#The final values used for the model were n.trees
#= 900, interaction.depth = 8 and shrinkage = 0.02.

#The final values used for the model were n.trees
#= 1000, interaction.depth = 10 and shrinkage = 0.02.
#########
nu  <- 0.001
D <- 3

data1 <- data.frame(y=y,day = day, month = month, year = year,sponsor=spon.code,grant=grant.code,
                    cvalue=cvalue,nsuccess=nsuccess,nsuccessind = nsuccessind, nfail=nfail,nfailind = nfailind,
                    has.PhD = has.PhD, country = country, has.ID = has.ID, role = role, Astar = Astar,
                    Astarind = Astarind, A = A, Aind = Aind, B = B, Bind = Bind, C = C, Cind = Cind,
                    maxAstar = maxAstar, maxAstarind = maxAstarind, maxC = maxC, maxCind = maxCind,
                    numPeople = numPeople)
set.seed(1)
fit.gbm<-gbm.fit(data1[tra,2:ncol(data1)], data1[tra,1], n.tree=500, interaction.depth=D, 
                 shrinkage=nu, distribution="bernoulli", verbose=FALSE)   # How to choose these parameters
best.iter<-gbm.perf(fit.gbm, plot.it=FALSE, method="OOB")
while(fit.gbm$n.trees-best.iter<50){
  fit.gbm<-gbm.more(fit.gbm, 100)           # do another 50 iterations
  best.iter<-gbm.perf(fit.gbm,plot.it=FALSE,method="OOB")
}
best.iter

pred1<- predict.gbm(fit.gbm,data1[val,2:ncol(data1)],n.trees = best.iter,type="response")
pred2<- as.numeric(pred1>=0.5)
tab  <- table(data1[val,1],pred2)
tab
(tab[1,2]+tab[2,1])/sum(tab) #0.1282201

gbm.ROC <- roc(data1[val,1],pred1)
plot(gbm.ROC, add = TRUE, print.auc = FALSE, col = "blue", 
     lty = 2, cex = 0.7 )

legend("bottomright", inset = 0.05, c("RF", "GBM"), lty = c(1,2), col = c("red", "blue"),
       bty = "n")


# variable importance
#par(las = 1,mgp = c(1.5,0.3,0),mar = c(5,7,2,2), cex= 0.75) #mgp: adjust margin of axis, label 
imp2 <- summary(fit.gbm, n.trees = best.iter, plotit = FALSE)
imp2.order <- imp2[order(imp2[,2], decreasing = TRUE),]
#par(las = 1,mgp = c(3,1,0),mar = c(5.1,4.1,4.1,2.1))
library(Hmisc)
#par(mar = c(4,0,2,4), mgp = c(2.5,0.3,0)) 
par(mar = c(3,2,0.5,2), mgp = c(2.5,0.8,0))
dotchart2(imp2.order$rel.inf, labels = rownames(imp2.order), xlab = "Relative Importance of GMB", 
          width.factor = 3, dotsize =2, pch = 21, cex = 0.75)

# Cvalue
cvaluetable <- table(data1$y,data1$cvalue)
#par(las = 1,mgp = c(2.5,0.3,0),mar = c(5,7,2,2), cex= 0.75)
par(mar = c(3.7,6,0.5,2), mgp = c(2.5,0.8,0))
barplot(cvaluetable, main= NULL,
        xlab="Contract Value", ylab = "Counts", col=c("darkblue","red"),
        legend = rownames(cvaluetable), beside=TRUE)


nsuccesstable <- table(data1$y,data1$nsuccess)
barplot(nsuccesstable, main=NULL,
        xlab="Number of Success History", ylab = "Counts", col=c("darkblue","red"),
        legend = rownames(nsuccesstable), beside=TRUE)

sponsortable <- table(data1$y,data1$sponsor)
barplot(sponsortable, main=NULL,
        xlab="Sponsor", ylab = "Counts", col=c("darkblue","red"),
        legend = rownames(sponsortable), beside=TRUE)

#rookie
temp <- cbind.data.frame(nsuccess,nfail,y)
temp1 <- temp[which(temp$nfail == 0) ,]
temp2 <- temp1[which(temp1$nsuccess == 0) ,]
rookietable <- table(temp2$y)
barplot(rookietable, main=NULL,
        xlab="Rookie", ylab = "Counts", col = NA)

par(mar = c(3.7,6,1.5,2), mgp = c(2.5,0.8,0))
nfailtable <- table(data1$y,data1$nfail)
barplot(nfailtable, main=NULL,
        xlab="Number of Fail History", ylab = "Counts", col=c("darkblue","red"),
        legend = rownames(nfailtable), beside=TRUE)

par(mar = c(3.7,5,0.5,2), mgp = c(2.5,0.8,0))
monthtable <- table(data1$y,factor(data1$month,levels=month.name))
#legend(xjust = 1.5)
#par(las = 1,mgp = c(2,0.3,0),mar = c(5,7,2,2), cex= 0.75)
barplot(monthtable, main= NULL,
        xlab="Month", ylab = "Counts", col=c("darkblue","red"),
        legend = rownames(monthtable), beside=TRUE)

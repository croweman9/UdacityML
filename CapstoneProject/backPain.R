library(HSAUR3)
library(rpart)
library(party)
library(gclus) 
library(partykit)
library(adabag)
library(randomForest)
library(MASS)

# Read in the data
pain<-read.table(url("http://mathsci.ucd.ie/~brendan/data/Physio.txt"))

pain=na.omit(pain)

is.factor(pain$assigned.labels)

pain=data.frame(pain)

####  Classify the data

fit.r <- rpart(pain$assigned.labels~.,data=pain)
plot(as.party(fit.r))

plot(fit.r)
text(fit.r,use.n=TRUE,cex=0.5,xpd=TRUE,col="red")

pred <- predict(fit.r,type="class")

tab <- table(pain$assigned.labels,pred)
tab
# Work out the accuracy
sum(diag(tab))/sum(tab)


#Manual Bagging
set.seed(1234)
res=matrix(NA,100,2)

iterlim <- 100
for (iter in 1:iterlim)
{
  N <- nrow(pain)
  indtrain <- sample(1:N,replace = TRUE)
  indtrain <- sort(indtrain)
  indtest <- setdiff(1:N,indtrain)


  fit.r2 <- rpart(assigned.labels~.,data=pain,subset=indtrain)
  pred2 <- predict(fit.r2,type="class",newdata=pain)
  
  # Look at table for the test data only (rows=truth, cols=prediction)
  tab2 <- table(pain$assigned.labels[indtest],pred2[indtest])
  tab3 <- table(pain$assigned.labels[indtrain],pred2[indtrain])


  res[iter,1]=sum(diag(tab2))/sum(tab2)
  res[iter,2]=sum(diag(tab3))/sum(tab3)

  }

res
colnames(res)<-c("test","train")
apply(res,2,summary)







###### Check against a multinomial model
set.seed(245)
res2<-matrix(NA,100,4)

# Start simulation to look at this 
iterlim <- 100
for (iter in 1:iterlim)
{
  # Sample 50% of the data as training data
  # Sample 25% of the data as validation 
  # Let the remaining 25% data be test data
  
  N <- nrow(pain)
  indtrain <- sample(1:N,size=0.50*N,replace=FALSE)
  indtrain <- sort(indtrain)
  indvalid <- sample(setdiff(1:N,indtrain),size=0.25*N)
  indvalid <- sort(indvalid)
  indtest <- setdiff(1:N,union(indtrain,indvalid))
  
  fit.r4 <- rpart(pain$assigned.labels~.,data=pain,subset=indtrain)
  
  # Fit a logistic regression to the training data only too
  mult1 <- multinom(pain$assigned.labels~., data=pain,subset=indtrain)
  fit.m=stepAIC(mult1,scope = ~ ., direction = "forward")
  
  # Classify for ALL of the observations
  pred.r4 <- predict(fit.r4,type="class",newdata=pain)
  pred.m <- predict(fit.m,type="class",newdata=pain)
  
  # Look at table for the validation data only (rows=truth, cols=prediction)
  tab.r4 <- table(pain$assigned.labels[indvalid],pred.r4[indvalid])
  tab.m <- table(pain$assigned.labels[indvalid],pred.m[indvalid])
  
  # Work out the accuracy
  acc.r4 <- sum(diag(tab.r4))/sum(tab.r4)
  acc.m <- sum(diag(tab.m))/sum(tab.m)
  
  # Store the results
  res2[iter,1] <- acc.r4
  res2[iter,2] <- acc.m
  
  # Look at the method that did best on the validation data 
  # when applied to the test data
  if (acc.r4>acc.m)
  {
    tab5 <- table(pain$assigned.labels[indtest],pred.r4[indtest])
    acc <- sum(diag(tab5))/sum(tab5)
    res2[iter,3] <- 1
    res2[iter,4] <- acc
  }else
  {
    tab5 <- table(pain$assigned.labels[indtest],pred.m[indtest])
    acc <- sum(diag(tab5))/sum(tab5)
    res2[iter,3] <- 2
    res2[iter,4] <- acc
  }
  
} 

# Check out the error rate summary statistics.
colnames(res2)<-c("valid.r","valid.m","chosen","test")
apply(res2,2,summary)
table(res2[,3])






# A priori (all im looking for is the combos that these appear in)

pain<-read.table(url("http://mathsci.ucd.ie/~brendan/data/Physio.txt"))
is.factor(pain2)
pain=data.frame(pain)
pain2=pain[,-37]
pain2=1*(pain2==1)

N=nrow(pain)
indices <- 1:N

fit <- apriori(pain2,parameter=list(confidence=0.9,support=0.6))
fit
painmat <- as.matrix(pain2)

#Get a table of the outcome for each symptoms
apply(painmat,2,table)

disabilitycount <- apply(painmat,1,sum)
hist(disabilitycount)
summary(disabilitycount)

pain2=as.matrix(pain2)




pain2=pain[,-37]


pain2 <- na.omit(pain2)
pain2 = sort(pain2)
pain2=as.matrix(pain2)

pain3 <- as(pain2,"transactions")

x=apply(pain2,2,table)
x=x[2,]
x=x/425
x

plot(x[2,],type = h)

symptomcount <- apply(pain2,1,sum)
hist(symptomcount)
summary(symptomcount)


fit <- apriori(pain2,parameter=list(confidence=0.9,support=0.4
                                        ,minlen=6,maxlen=8))
fit <- sort(fit,by="lift")
inspect(fit)



############## Clustering ###############


pain2<-pain[,-37]
pain2=na.omit(pain2)

#Plot the data. Do you see clusters?
pairs(pain2)
# Run k-means (with K=3) and store the resulting fit in fitkm.
pain2<-1*(pain2==1)    #Making it a matrix

fitkm<-kmeans(pain2,centers=3)
fitkm
# Find optimal number clusters
N=10  
storemat=matrix(NA,10,2)
for (k in 1:N) 
{
  fitkm<-kmeans(pain2,centers=k,nstart=10)
  storemat[k,1]=fitkm$tot.withinss
  storemat[k,2]=k
}
storemat

plot(storemat[,2],storemat[,1],type = "b",main="Optimal Clusters",xlab="Number clusters"
     ,ylab="Within SS")




# Silhouettes

d=dist(pain2, method = "binary")^2
for (i in 2:N) {
  fitkmtest=kmeans(pain2,centers = i,nstart = 10)
  sil = silhouette(fitkmtest$cluster,d)
  plot(sil)  
}

#Average Silhouette Values
vec=c(0.56,.48,.4,.34,.28,.26,.24,.23,.23)
vec2=c(2,3,4,5,6,7,8,9,10)

plot(vec2,vec,type = "b",main = "Avg. Sil Vs. No. Clusters"
     ,xlab = "Number Clusters",ylab = "Avg. Silhouette")


# Check cross tab between k-means and k-medoid
fitkm<-kmeans(pain2,centers=3,nstart=10)
fitkmed<-pam(d,k=3)

tab<-table(fitkm$cluster,fitkmed$clustering)
tab

classAgreement(tab)

#Check if similar to back pain disorders, [37]

tabfinal=table(fitkm$cluster,pain[,37])
tabfinal

acc.km=sum(diag(tabfinal))/sum(tabfinal)
acc.km
classAgreement(tabfinal)


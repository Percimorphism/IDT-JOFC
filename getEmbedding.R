#This is the first part of the IDT*JOFC part, i.e. JOFC
#IDT * JOFC
#We use the neuron data in Carey and youngser's science paper

source("fastsmacof.R")
source("idt-function.R")

#the following outcome is "Znorm"
print(load("Znorm-single-37780x632-f8.Rbin"))

head(Znorm)
dim(Znorm)


#data sets D1 and D2
D1=Znorm[,80:158]
D2=Znorm[,554:632]

samp=sample(1:37780, 3000, replace=FALSE)

D1=D1[samp,]
D2=D2[samp,]

#we compute the dissimilarity matrix
Delta1=as.matrix(dist(D1, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
Delta2=as.matrix(dist(D2, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))

Delta1=Delta1/norm(Delta1, "F")
Delta2=Delta2/norm(Delta2, "F")


Delta=list(Delta1, Delta2)

embedding=fast.smacof(0.9, Delta, ndim=500)

save(embedding, file="embeddingF=2,8N=3000D=500.Rbin")



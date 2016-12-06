#This is the second half of the idt*jofc, i.e. IDT. 
#IDT * JOFC
#We use the neuron data in Carey and youngser's science paper

source("fastsmacof.R")
source("idt-function.R")


print(load("embeddingF=2,8N=1000D=500.Rbin"))

typeof(embedding)

embedconflist=embedding$conf
#print(embedconflist)


###############################################################
#it returns a list of matrices (4 in the iris dataset)

#stack the final configs on top of each other
embedconf=data.frame()
for (i in embedconflist){
  #print(i)
  embedconf=rbind(embedconf, i)
}

rownames(embedconf) = NULL
embedconf=as.matrix(embedconf)
#print((config))
#print(rownames(config))

#PLEASE CENTER THE DATA FIRST! This is important!
embedconf=scale(embedconf, center = TRUE)
attr(embedconf,"scaled:center")<-NULL 
attr(embedconf,"scaled:scale")<-NULL 
##############################################################



embedconf=as.data.frame(embedconf)
print("")
print("the dim of the embedding conf is")
print(dim(embedconf))


embedconf.pca=prcomp(embedconf, retx=TRUE, center=TRUE)

dim=getElbows(embedconf.pca$sdev)[3]
print("")
print("dim we selected in scree plot is:")
print(dim)

idtdata=embedconf.pca$x #[,1:dim]

resultlabel=doIdt(idtdata, maxdepth=5, minnum=50)
#plottree(idtdata, maxdepth=5, minnum=50 )

print(resultlabel[1:1000])

print("  seperating line  ")

print(resultlabel[1001:2000])


adjustedRandIndex(resultlabel[1:1000], resultlabel[1001:2000])

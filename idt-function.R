source("getElbows.R")
library(mclust)
library(data.tree)


################################################################################
#we create a function for splitting a given dataset after pca, output should be a list of
#new dataset, each is a subset of the ORIGINAL dataset
#there is a minnum of elements in the dataset, if under that threshhold, we will
#not continue doing pca on the feat
splitnode = function(feat, minnum){
  #we plot everything into a pdf file
  pdf("plots.pdf")
  
  #put it as a dataframe
  as.data.frame(feat)
  
  #check if the feat set is larger than minnum
  if (nrow(feat) <= minnum){
    splittedfeat=list()
    newlabel=rep(1, nrow(feat))
    feat_with_label=cbind(feat, newlabel)
    feat_with_label=as.data.frame(feat_with_label)
    splittedfeat[[1]]=feat_with_label
    return(splittedfeat)
  }
  
  
  else{
    
    # first we do PCA on the data, note that feat.pca$x is the data after pca (transformed data)
    #centerfeat=scale(feat, center=TRUE)
    #attr(centerfeat,"scaled:center")<-NULL 
    #attr(centerfeat,"scaled:scale")<-NULL 
    ########################################################################
    feat.pca=prcomp(feat, retx=TRUE, center=TRUE, scale=TRUE)
    
    #we select the 2nd elbow as the dimension for clustering, also get the scree plot
    if (!is.na(getElbows(feat.pca$sdev)[2])){
      dim=getElbows(feat.pca$sdev)[2]
    }
    else{
      dim=getElbows(feat.pca$sdev)[1]
    }
    #print(dim)
    
    #print(feat.pca$x)
    #we do clustering on the data feat.pca$x[,1:dim]
    feat_clusted = Mclust(feat.pca$x[,1:dim], G=1:2, modelNames=c("EEE","EEV","VEV","VVV"))
    
    #plot the clusted data
    #plot(feat_clusted, what="classification")
    
    #print out the labels after clustering
    #print(feat_clusted$classification)
    
    #reattached the labels to the feat_pca, prepare for break-up by label
    feat.pca_newdatawithlabel=cbind(feat.pca$x, feat_clusted$classification)
    
    #make sure the resulted data is a dataframe
    feat.pca_newdatawithlabel=as.data.frame(feat.pca_newdatawithlabel)
    
    #we split the data accorading to the lables
    feat2_withlabel = split(feat.pca_newdatawithlabel , f = feat.pca_newdatawithlabel[,ncol(feat.pca_newdatawithlabel)])
    
    dev.off()
    
    #remember here we return a list of subset of original dataset
    #not the data after pca
    splittedfeat=list()
    labelindex=1
    for (i in feat2_withlabel){
      rownames=rownames(i)
      newlabel=rep(labelindex, nrow(feat[rownames,]))
      subcluster_with_label=cbind(feat[rownames,], newlabel)
      subcluster_with_label=as.data.frame(subcluster_with_label)
      labelindex=labelindex+1
      #print(feat[rownames,])
      #print(subcluster_with_label)
      splittedfeat[[length(splittedfeat)+1]]=subcluster_with_label
    }
    #return a list of subclusters
    return(splittedfeat)
  }
}

#####################################################################################
#we grow tree, where mynode is an node object
getNode=function(mynode, minnum){
  #print(mynode$data)
  #we split the node (clustering on the data in the node)
  #print(minnum)
  children=splitnode(mynode$data, minnum)
  
  #we make them into a list of node objects
  child=list()
  index=1
  for (i in children){
    #print(typeof(i))
    #############################################################
    #note the syntax here!!! It took me a while to fiture out!!!#
    #############################################################
    #n=ncol(i)
    child[[index]]=Node$new(index, data=i[, (-ncol(i))], labels=i[, ncol(i)])
    #print(child[[index]]$data)
    index=index+1
  }
  return(child)
}
#print(getNode(idtTree))

#####################################################################################
#This function recursively grow the tree
#treedepth is where we keep track of the depth of the recursion
buildIdtTree=function(mynode, maxdepth, minnum, treedepth=0){
  #we returned a list of nodes that need to be added as children of mynode
  #note here we return a list of nodes!
  #print(minnum)
  children=getNode(mynode, minnum)
  #for(i in children){
  #  print(i$labels)
  #}
  
  childnode=list()
  index=1
  for(i in children){
    if(treedepth < maxdepth){
      childnode[[index]]=mynode$AddChildNode(i)
      #print(childnode[[index]]$data)
      #print(childnode[[index]]$labels)
      #note here we keep track of the depth of the recursion
      buildIdtTree(childnode[[index]], maxdepth, minnum, treedepth+1)
      index=index+1
    }
  }
  return(mynode)
}

###############################################################################################
#we do a complete idt Algorithm
doIdt=function(feat, maxdepth, minnum){
  #print(minnum)
  #we build out the tree, from the node
  idtTree=Node$new("idtTree", data=as.data.frame(feat))
  #the resulted tree after idt
  result=buildIdtTree(idtTree, maxdepth, minnum)
  #print("we are going to plot the tree")
  #plot(result)
  #print("we printed the tree")
  
  #the really stupid way of making the labels different
  #Anyone can think of a clever way?
  #get the list of all leaf nodes
  finalclusters=Traverse(result, filterFun = isLeaf)
  #creat a list of all labels
  finalclusterlabels=c()
  labelindex=1
  for (i in finalclusters){
    #print(rownames(i$data))
    #i$labels=rep(labelindex, length(i$labels))
    for (j in rownames(i$data)){
      finalclusterlabels[as.numeric(j)]=labelindex
    }
    labelindex=labelindex+1
  }
  print("the labels after we do idt is:")
  print(finalclusterlabels)
  return(finalclusterlabels)
  
}

plottree=function(feat, maxdepth, minnum){
  idtTree=Node$new("idtTree", data=as.data.frame(feat))
  #the resulted tree after idt
  result=buildIdtTree(idtTree, maxdepth, minnum)
  
  plot(result)
  
}





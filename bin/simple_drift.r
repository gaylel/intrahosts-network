pop=matrix(rep(0,100000),ncol=100)

size.chain=100
n.pop.peak=10000
mutation.rate=.0001
effective.BN=100

pop=matrix(rep(0,size.chain*n.pop.peak),ncol=size.chain)


#init.pop[1,]=1

evolve.init<-function(init)
{
  new.pop=t(matrix(rep(t(init),n.pop.peak/effective.BN),nrow=size.chain))
  
  result=t(apply(new.pop,1,function(x)
    {
    n.mut=rbinom(1,size.chain,mutation.rate)
    if(n.mut>0)
    {
      mut=sample(1:size.chain,n.mut)
      x[mut]=(1-x[mut])^2
    }
    return(x)
  }))
  
  #sum(result!=new.pop)
  
  return(result)
}

consensus<-function(pop)
{
  apply(pop,2,median)  
}

summary_pop<-function(pop)
{
  print("Consensus:")
  con=consensus(pop)
  print(con)
  
  minor=colSums(t(t(pop)-consensus(pop)))
  pos.var=which(minor!=0)
  
  for(i in pos.var)
  {
    print(paste(abs(minor[i])*100/n.pop.peak,"% variant in pos",i))
  }
}

list_unique_seq<-function(pop)
{
  res<-NULL
  
  #extract the consensus
  con<-consensus(pop)
  diff=rowSums(abs(t(t(pop)-con)))
  
  pop_it=pop[order(diff),]
  diff=sort(diff)
 
  n_seq=sum(diff==0)
  res<-list(list(seq=con,n=n_seq,dist=0))  
  
  remaining_seq=length(pop_it[,1])-n_seq
  if(remaining_seq>0)
    pop_it=pop_it[-(which(diff==0)),]
    
  #iterates for minority sequences
  while(remaining_seq>0)
  {
    if(!is.vector(pop_it))
    {
      diff=rowSums(abs(t(t(pop_it)-pop_it[1,]))) 
    
      n_seq=sum(diff==0)
      p_res<-list(seq=pop_it[1,],n=n_seq,dist=sum(abs(con-pop_it[1,])))
      res<-append(res,list(p_res))
    
      remaining_seq=length(pop_it[,1])-n_seq
      if(remaining_seq>0)
        pop_it=pop_it[-(which(diff==0)),]
    } else
    {
      res<-append(res,list(list(seq=pop_it,n=1,dist=sum(abs(con-pop_it)))))
      remaining_seq=0
    }  
  }
  
  return(res)
}  
  
n_chain=100
chain=rep(0,n_chain)
for(i in 1:n_chain)
{
  init.pop=pop[sample(1:n.pop.peak,effective.BN),]  
  pop=evolve.init(init.pop)
  #print(table(rowSums(pop))*100/n.pop.peak)
  chain[i]=sum(rowSums(pop)==0)
}

uni_seq<-list_unique_seq(pop)
n_seq=length(uni_seq)
dist_mat=matrix(rep(0,n_seq^2),ncol=n_seq)
for(i in 2:n_seq)
  for(j in 1:(i-1))
  {
    dist_mat[i,j]=sum(abs(uni_seq[[i]]$seq-uni_seq[[j]]$seq))
    dist_mat[j,i]=dist_mat[i,j]
  }

M<-mst(dist_mat)
require(PEcAn.all)
logger.setQuitOnSevere(FALSE)
settings <- read.settings("gr.settings.xml")
td <- get.trait.data(settings$pfts,settings$run$dbfiles,settings$database,TRUE)

## rescale trait data
trait.file = file.path(settings$pfts$pft$outdir, 'trait.data.Rdata')

load(trait.file)
for(i in 1:length(trait.data)){
  trait.data[[i]]$mean = trait.data[[i]]$mean/100
  trait.data[[i]]$stat = trait.data[[i]]$stat/100
}
save(trait.data,file=trait.file)

##PEcAn - get posterior priors
run.meta.analysis(td, settings$meta.analysis$iter, settings$run$dbfiles, settings$database)
load(file.path(settings$pfts$pft$outdir,"trait.mcmc.Rdata"))
load(file.path(settings$pfts$pft$outdir,"post.distns.Rdata"))

#########################################

c <- read.csv("cost.csv")
cost <- c$CO2Produced
NC <- length(cost) # #Components

## Convert gCO2 to gC
## gCO2*(12gC/44gCO2)
cost = cost*(12/44)


leafvariables = c('l_carbohydrates','l_lignin','l_lipids','l_minerals','l_organicacids','l_protein')
stemvariables = c('s_carbohydrates','s_lignin','s_lipids','s_minerals','s_organicacids','s_protein')
rootvariables = c('r_carbohydrates','r_lignin','r_lipids','r_minerals','r_organicacids','r_protein')

variables=matrix(c(leafvariables,stemvariables,rootvariables),NC,3)
NV=length(variables) #Number of variables

##########################################################
##FUNCTION FOR ANALYSIS
## function input: trait data and construction costs
## function outopt: list of 4 matrices of growth respiration values for uninformative prior, leaf, root, and stem
###########################################################

getdistribution <- function(trait.mcmc,post.distns,cost,variables) {
  NC=length(cost)
  NV=length(variables)

  #R = Rl + Rs + Rr
  
  #leaf
  #Rl = kl*Gl
  #kl = cost (g C produced) * pcompl (percent composition of leaf components)

  ## calc mean and sd from meta-analysis
  mean = matrix(NA,NC,3)
  var = matrix(NA,NC,3)
  
  for(i in 1:NV){
    if(variables[i] %in% names(trait.mcmc)){
      y = as.matrix((trait.mcmc[[variables[i]]]))[,"beta.o"]
      mean[i]= mean(y)
      var[i]= var(y)
    } else {
      ## use the prior
      row = which(rownames(post.distns) == variables[i])
      if(length(row)==1 & post.distns$distn[row] == 'beta'){
        x = post.distns[row,]
        mean[i] = x$parama/(x$parama+x$paramb)
        var[i]  = (x$parama*x$paramb)/((x$parama+x$paramb)^2*(x$parama+x$paramb+1)^2)    
      }
    }
  }

  ## moment matching to est. alpha
             
              # USING DIRICHLET:
              # mean[i]=a[i]/a0
              # var[i]=a[i]*(a0-a[i])/((a0)^2*(a0+1))
                # a = matrix(NA,NC,3)
                # for(i in 1:length(variables)){
                  # a[i]=mean[i]*(((mean[i]-mean[i]^2)/var[i])-1)
                  # }
  # USING BETA
    # E[x]=M=a/(a+B)
      # B=a(1-M)/M
    # Var[X]=aB/[(a+B)^2(a+B+1)]
  a=B=matrix(NA,NC,3)
  for(i in 1:NV) {
    a[i]=(1-mean[i])*mean[i]^2/var[i]-mean[i]
    B[i]= a[i]*(1-mean[i])/mean[i]
  }
  
  ########## functions to rescale percent composition to sum to 1 #############
  
  NewP.oldDoesntWork <- function(k,p,a,b){
    # calculate current quantile
    q0 = pbeta(p,a,b)
    qm = pbeta(a/(a+b),a,b)
  
    # adjust by k
    qnew = qm + k*(q0-qm)
    qnew[qnew<0] = 0
    qnew[qnew>1] = 1
  
    # convert back to p
    pnew = qbeta(qnew,a,b)
    return(pnew)  
  }

  NewP <- function(k,p,a,b){
    # calculate current quantile
    q0 = pbeta(p,a,b)

    # calc SD equivalent of current quantile
    sd0 = qnorm(q0)
  
    # adjust by k
    sd.new = sd0 + k
  
    # calc new quantile
    q.new = pnorm(sd.new)
    
    # convert back to p
    pnew = qbeta(q.new,a,b)
    return(pnew)  
  }

  SumToOneFactor <- function(k,p,a,b){
    pnew = NewP(k,p,a,b)
    # assess sum to 1
    return((sum(pnew)-1)^2)
  }

  N = 5000 # Iterations
  ## l=leaf; s=stem; r=root; nd=assuming no parameter data
  G=Gl=Gs=Gr=matrix(1,N,1)
  Rl=Rs=Rr=Rnd=matrix(NA,N,1)
  pcompl=pcomps=pcompr=pcompnd=matrix(NA,N,NC) #storage for % composition
  kl=ks=kr=knd=matrix(NA,N,1) #cost*%composition

  # get percent composition using alpha and beta
  for(i in 1:N){
    # rdirichlet(1,c(,1,1,1,1,1))
    # pcompl[i,]=rdirichlet(1,c(a[,1]))
    # pcomps[i,]=rdirichlet(1,c(a[,2]))
    # pcompr[i,]=rdirichlet(1,c(a[,3]))
    for (j in 1:NC) {
      pcompnd[i,j]=rbeta(1,1,5)  
      pcompl[i,j]=rbeta(1,a[j,1],B[j,1]) 
      pcomps[i,j]=rbeta(1,a[j,2],B[j,2]) 
      pcompr[i,j]=rbeta(1,a[j,3],B[j,3])
    }
    ## Rescale pcomp output so sums to 1
    kopt = optimize(SumToOneFactor,c(-10,10),p=pcompnd[i,],a=1,b=6)
    popt = NewP(kopt$minimum,p=pcompnd[i,],a=1,b=6)
  
    koptl = optimize(SumToOneFactor,c(-10,10),p=pcompl[i,],a=a[,1],b=B[,1])
    poptl = NewP(koptl$minimum,p=pcompl[i,],a=a[,1],b=B[,1])
  
    kopts = optimize(SumToOneFactor,c(-10,10),p=pcomps[i,],a=a[,2],b=B[,2])
    popts = NewP(kopts$minimum,p=pcomps[i,],a=a[,2],b=B[,2])
  
    koptr = optimize(SumToOneFactor,c(-10,10),p=pcompr[i,],a=a[,3],b=B[,3])
    poptr = NewP(koptr$minimum,p=pcompr[i,],a=a[,3],b=B[,3])
  
    knd[i,]=sum(cost*popt)
    kl[i,]=sum(cost*poptl)
    ks[i,]=sum(cost*popts)
    kr[i,]=sum(cost*poptr)
  
    if(i %% 1000 == 0) print(i)

  }

  # Calculate growth respiration for leaf, stem, and root
  Rnd=knd*G  ## UNINFORMATIVE PRIOR; no percent composition data
  Rl=kl*Gl
  Rs=ks*Gs
  Rr=kr*Gr

  R<- list("Rnd"=Rnd,"Rl"=Rl,"Rs"=Rs,"Rr"=Rr,"var"=var)
  return(R)
} #end of function
##########################################################################

R.allplants <- getdistribution(trait.mcmc,post.distns,cost,variables)

########## Create Plot of Distributions ##################################

cols = 1:4
dRnd = density(R.allplants$Rnd)
plot(density(R.allplants$Rl),xlim=range(dRnd$x),col=cols[2])
lines(dRnd,col=cols[1])
lines(density(R.allplants$Rs),col=cols[3])
lines(density(R.allplants$Rr),col=cols[4])
legend("topright",legend=c("Null","Leaf","Stem","Root"),col=cols,lwd=2)


########### Variance Decomposition ####################################################
## sum(Pcomp^2*Var(cost) + sum(cost^2*Var(Pcomp))  ## no variance in construction costs

vd = matrix(NA,NC,3)

for (i in 1:NC){
  vd[i,1]=cost[i]^2*var(pcompl[,i])
  vd[i,2]=cost[i]^2*var(pcomps[,i])
  vd[i,3]=cost[i]^2*var(pcompr[,i])
}

## alternative that doesn't have sum to 1 constraints
for (i in 1:NC){
  vd[i,1]=cost[i]^2*R.allplants$var[i,1]
  vd[i,2]=cost[i]^2*R.allplants$var[i,2]
  vd[i,3]=cost[i]^2*R.allplants$var[i,3]
}

colnames(vd) <- c("leaf","stem","root")
rownames(vd) <- c("carb","lignin","lipid","mineral","OA","protein")

totvar <- apply(vd,2,sum)
t(vd)/totvar  ##  % variance

totsd <- apply(sqrt(vd),2,sum)
t(sqrt(vd))/totsd *100 ##  % sd


##########################################################
  
## Build covariates table

ctable=matrix(NA,0,NV)
colnames(ctable) <-variables

for(i in 1:length(variables)){
  ##Find variable in trait data
  if(variables[i]%in%names(trait.data)){
    tr=which(names(trait.data)==variables[i])
    ##Create unique ID for trait
    v=paste(trait.data[[tr]]$specie_id,trait.data[[tr]]$site_id,sep="#")
    for(j in 1:length(v)){
      ##if ID is already has a row in the table
      if(v[j]%in%rownames(ctable)){
        rownumber=which(rownames(ctable)==v[j])
        if(is.na(ctable[rownumber,i])) {
          ctable[rownumber,i]=trait.data[[i]]$mean[j]
        } else  {
        ####But if space in table already full
        ##average current and new value
          ctable[rownumber,i]==mean(c(ctable[rownumber,i],trait.data[[i]]$mean[j]))
        }
      } else{
        ##if ID is new
        newrow=matrix(NA,1,length(variables))
        rownames(newrow)=v[j]
        newrow[,i]=trait.data[[i]]$mean[j]
        ctable=rbind(ctable,newrow)
      }
    }
  }
}

## fit missing data model to estimate NAs
MissingData = "
model{
  for(i in 1:n){
    x[i,] ~ dmnorm(mu,tau)
  }
  mu ~ dmnorm(m0,t0)
  tau ~ dwish(R,k)
  x[2,4] <- xmis 
  xmis ~ dnorm(0.2,1)
}
"
w = ncol(ctable)
data <- list(x = ctable,n=nrow(ctable),m0=rep(1/6,w),t0 = diag(1,w),R = diag(1e-6,w),k=w)

#test
#w = 4
#data <- list(x = ctable[1:2,1:w],n=2,m0=rep(1/6,w),t0 = diag(1,w),R = diag(1e-6,w),k=w)


j.model = jags.model(file=textConnection(MissingData),
                     data = data,
                     n.chains=1,
                     n.adapt=10,
                     inits = list(xmis = 0.1))


logit <- function(p){
  log(p/(1-p))
}

ilogit <- function(x){
  exp(x)/(1+exp(x))
}


Z = logit(ctable)
m = nrow(ctable)

###set up covariates
Zorig <- as.matrix(Z)
ncov <- ncol(as.matrix(Z))
#find Zobs
Zobs <- apply(Z,2,mean,na.rm=TRUE)
## HACK##
if(is.nan(Zobs[10])) Zobs[10] = Zobs[4]
if(is.nan(Zobs[11])) Zobs[11] = Zobs[5]
if(is.nan(Zobs[13])) Zobs[13] = Zobs[1]
if(is.nan(Zobs[16])) Zobs[16] = Zobs[4]
if(is.nan(Zobs[17])) Zobs[17] = Zobs[5]
if(is.nan(Zobs[18])) Zobs[18] = Zobs[12]

n.Z <- nrow(Z)
for(i in 1:ncov){
  Z[is.na(Zorig[,i]),i] <- Zobs[i]
}

## initial guess
Z.init = ilogit(Z)

#priors for Zmis

#mean mu
muZ.ic <- Zobs
mu.Z0 <- rep(logit(1/6),ncol(Z))  #post-normalization
M.Z0 <- diag(rep(10,ncol(Z)))
IM.Z0 <- solve(M.Z0)
#cov V
V.Z.ic <- diag(cov(Z,use="pairwise.complete.obs"))
x.Z <- ncov + 2
V.Z0.all <-  M.Z0*x.Z
V.Z0 <- diag(V.Z0.all)
IV.Z0 <- solve(V.Z0.all)
mu.Z <- mu.Z0
V.Z <- V.Z0.all 
IV.Z <- solve(V.Z)

library(MCMCpack)
library(mvtnorm)

## set storage
start = 1
ngibbs = 500
muZgibbs <- matrix(0,nrow=ngibbs,ncol=ncov)
VZgibbs <- matrix(0,nrow=ngibbs,ncol=ncov*(ncov+1)/2)
Zgibbs <- Z*0

#gibbs loop
btimes <- 0
for(g in start:ngibbs){
  print(g)

  ##missing Z's - mean
  bigv <- try(solve(n.Z*IV.Z + IM.Z0))
  if(is.numeric(bigv)){
    smallv <- apply(Z %*% IV.Z,2,sum) + IM.Z0 %*% mu.Z0
    mu.Z <- rmvnorm(1,bigv %*% smallv,bigv)   
  }
  muZgibbs[g,] <- mu.Z

  
##missing Z's - Variance
u <- 0
for(i in 1:m){ u <- u + crossprod(Z[i,]-mu.Z) }
V.Z.orig <- V.Z
IV.Z.orig <- IV.Z
V.Z <- riwish(x.Z + n.Z, V.Z0.all + u)
IV.Z <- try(solve(V.Z))
if(!is.numeric(IV.Z)){
  IV.Z <- IV.Z.orig
  V.Z <- V.Z.orig
}
VZgibbs[g,] <- vech(V.Z)

##missing Z's - draw missing values
for(i in 1:m){
  for(j in 1:ncov){
    if(is.na(Zorig[i,j])){
      bigv <- 1/IV.Z[j,j]
      smallv <- mu.Z[j]*IV.Z[j,j]
      zcols <- 1:ncov; zcols <- zcols[zcols != j]
      for(k in zcols){
        smallv <- smallv + (Z[i,k] - mu.Z[k])*IV.Z[k,j]
      }
      Z[i,j] <- rnorm(1,bigv * smallv, sqrt(bigv))
    }  	
  }
}
  Zgibbs = Zgibbs + Z
  
  if(g %% 500 == 0){ save.image("GR.RData")}
} #end Z.fillmissing
Zbar = ilogit(Zgibbs/g)

sum(is.na(Zbar))
cbind(apply(Zbar,2,mean),apply(Z.init,2,mean),
apply(ctable,2,mean,na.rm=TRUE))
pdf("muZgibb.pdf")
plot(as.mcmc(ilogit(muZgibbs)))
dev.off()

##################################################################
## PCA & Cluster Analysis
data.leaf <- Z.init[,1:6]
data.stem <- Z.init[,7:12]
data.root <- Z.init[,13:18]

## cluster analysis on raw leaf data
cluster.leaf <- kmeans(data.leaf,2)
plot(Z.init[,2],Z.init[,3])
plot(Z.init[,2],Z.init[,3],col=cluster.leaf$cluster)

cluster.stem <- kmeans(data.stem,2)
cluster.root <- kmeans(data.root,2)

## cluster analysis on leaf data weighted by construction costs
cluster.leaf.cost <- kmeans(t(t(data.leaf)*cost),2)
cluster.stem.cost <- kmeans(t(t(data.stem)*cost),2)
cluster.root.cost <- kmeans(t(t(data.root)*cost),2)
plot(Z.init[,2],Z.init[,3],col=cluster.leaf.cost$cluster)

## principal component analysis on raw leaf data
pca.leaf <- prcomp(data.leaf,retx=TRUE)
pca.leaf$sdev/sum(pca.leaf$sdev)*100

plot(pca.leaf)
plot(pca.leaf$x[,1],pca.leaf$x[,2])

## principal component analysis on leaf data weighed by construction costs
pca.leaf.cost <- prcomp(data.leaf,scale=cost,retx=TRUE)


##
cluster.pca.leaf <- kmeans(t(t(pca.leaf$x)*pca.leaf$sdev^2),2)
plot(pca.leaf$x[,1],pca.leaf$x[,2],col=cluster.pca.leaf$cluster)

cluster.pca.leaf <- kmeans(t(t(pca.leaf$x)*pca.leaf$sdev^2),2)

phenol = rbinom(nrow(pca.leaf$x),1,0.5) ## replace this with real data
## phenol.char = c("E","D")

plot(pca.leaf$x[,1],pca.leaf$x[,2],col=cluster.pca.leaf$cluster)


library(MASS)
library(vegan)
library("RPostgreSQL")

dbparms <- list(driver="PostgreSQL" , user = "bety", dbname = "bety", password = "bety")
con     <- db.open(dbparms)

## species category gymnosperm?
input = db.query(paste('SELECT "id","scientificname","commonname","Category","GrowthForm" FROM species'),con)
#input = db.query(paste("SELECT * FROM species"),con)

## vectors of categoires corresponding to data
categories = growthform = speciesnames = commonname = id = vector()

for (i in 1:length(input$id)) {
  k = grep(input$id[i], rownames(data.leaf), fixed=TRUE) ######## rownames(data.leaf) can't be redefined (should be id#site)
  categories[k]=input$Category[i]
  growthform[k]=input$GrowthForm[i]
  speciesnames[k]=input$scientificname[i]
  commonname[k]=input$commonname[i]
  id[k]=input$id[i]
}




#################### Fill missing categories ##################################
growthform[1]="Single Crown"####????
growthform[2]="Bunch" ####??????
  
#American Beech
categories[3]="Dicot"
growthform[3]="Single Stem"
commonname[3]="American beech"

growthform[9]="Single Stem"
growthform[17]="Single Stem"

categories[18]="Dicot"
growthform[18]="Single Crown" ###CHECK
commonname[18]="Yellow Alpine Pasqueflower"

growthform[21]="Single Stem"
growthform[22:26]="Bunch"

categories[29]="Dicot"
growthform[29]="Single Stem"
commonname[29]="Oak"

growthform[30]="Single Stem"
growthform[35]="Single Crown" ###???

woody <- vector()
for (i in 1:length(growthform)) {
  if (growthform[i]=="Single Stem") {
    woody[i]=TRUE
  } else if (growthform[i]=="Multiple Stem") {
    woody[i]=TRUE
  } else if (growthform[i]=="Single Crown") {
    woody[i]=TRUE
  } else if (growthform[i]=="Rhizomatous") {
    woody[i]=TRUE
  } else if (growthform[i]=="Bunch") {
    woody[i]=FALSE
  } else {
    woody[i]=NA
  }
}

sel.vect = which(!is.na(categories))


rownames(data.leaf)=rownames(data.stem)=rownames(data.root)=speciesnames #################### renaming rownames(data.leaf) for pca labels

characteristics = cbind(categories[sel.vect]=="Monocot",woody[sel.vect]==TRUE)
colnames(characteristics)=c("Monocot","Woody")

## fit species charactaristics to compositional pca
cluster.leaf.cost <- kmeans(t(t(data.leaf)*cost),2)

pca.leaf.cost <- prcomp(data.leaf[sel.vect,],scale=cost,retx=TRUE)
#ef.leaf <- envfit(pca.leaf.cost,as.factor(categories[sel.vect]),na.rm=TRUE)
ef.leaf <- envfit(pca.leaf.cost,characteristics,na.rm = TRUE)
biplot(pca.leaf.cost,cex=0.6,col=cluster.pca.leaf$cluster)
plot (ef.leaf,cex=0.5)


cluster.stem.cost <- kmeans(t(t(data.stem)*cost),2)

pca.stem.cost <- prcomp(data.stem[sel.vect,],scale=cost,retx=TRUE)
ef.stem <- envfit(pca.stem.cost,characteristics,na.rm=TRUE)
biplot(pca.stem.cost,cex=0.8)
plot(ef.stem,cex=0.5)


cluster.root.cost <- kmeans(t(t(data.root)*cost),2)

pca.root.cost <- prcomp(data.root[sel.vect,],scale=cost,retx=TRUE)
ef.root <- envfit(pca.root.cost,characteristics,na.rm=TRUE)
biplot(pca.root.cost,cex=0.8)
plot(ef.root,cex=0.5)

##################################
## Split into 2 distributions
##################################

####### Monocot vs. Dicot ########

## query species for char.
j=which(categories=="Monocot")

trait.data.mono = list()
trait.data.dicot = list()
for (i in 1:length(trait.data)) { 
  sel.mono = which(trait.data[[i]]$specie_id %in% id[j])
  trait.data.mono[[i]] = trait.data[[i]][sel.mono,]
  trait.data.dicot[[i]] = trait.data[[i]][-sel.mono,]
  #  for(l in 1:length(trait.data[[i]]$specie_id)) {
  #    if(trait.data[[i]]$specie_id[l]%in%id[j]) {
  #      #m=as.matrix(trait.data[[i]])
  #    }
  #  }
  #j = which(trait.data[[i]]$specie_id%in%id[j]) 
}  
  names(trait.data.mono) = names(trait.data)
  names(trait.data.dicot) = names(trait.data)

## split trait.data
#trait.data.foo1 = trait.data[k]
#trat.data.foo2 = trait.data[!k]

td.mono = td
td.mono$pft$outdir = "/home/carya/pecan/pft/gr.mono"

td.dicot = td
td.dicot$pft$outdir = "/home/carya/pecan/pft/gr.dicot/"


## save
#save(trait.data.foo1,file=foo1)
#save(trait.data.foo2,file=foo2)
save(trait.data.mono,file=file.path(td.mono$pft$outdir, 'trait.data.Rdata'))
save(trait.data.dicot,file=file.path(td.dicot$pft$outdir, 'trait.data.Rdata'))


##PEcAn - get posterior priors
#run.meta.analysis()
run.meta.analysis(td.mono, settings$meta.analysis$iter, settings$run$dbfiles, settings$database)
load(file.path(settings$pfts$pft$outdir,"trait.mcmc.Rdata"))
load(file.path(settings$pfts$pft$outdir,"post.distns.Rdata"))

R.mono = getdistribution(trait.mcmc,post.distns,cost,variables)

run.meta.analysis(td.dicot, settings$meta.analysis$iter, settings$run$dbfiles, settings$database)
load(file.path(settings$pfts$pft$outdir,"trait.mcmc.Rdata"))
load(file.path(settings$pfts$pft$outdir,"post.distns.Rdata"))

R.dicot = getdistribution(trait.mcmc,post.distns,cost,variables)

#### Probability distributions different?
ks.test(R.mono$Rl,R.dicot$Rl)
ks.test(R.mono$Rs,R.dicot$Rs)
ks.test(R.mono$Rr,R.dicot$Rr)

cols = 1:4
dR.monond = density(R.mono$Rnd)
plot(density(R.mono$Rl),xlim=range(dR.monond$x),col=cols[2])
lines(dR.monond,col=cols[1])
lines(density(R.mono$Rs),col=cols[3])
lines(density(R.mono$Rr),col=cols[4])
lines(density(R.dicot$Rnd),col=cols[1],lty=2)
lines(density(R.dicot$Rl),col=cols[2],lty=2)
lines(density(R.dicot$Rs),col=cols[3],lty=2)
lines(density(R.dicot$Rr),col=cols[4],lty=2)
legend("topright",legend=c("Null","Leaf","Stem","Root","Monocot","Dicot"),col=c(cols,1,1),lwd=2,lty=c(1,1,1,1,1,2))



###### Woody vs Nonwoody #######

## query species for char.
j=which(woody==TRUE)

trait.data.woody = list()
trait.data.nonwoody = list()
for (i in 1:length(trait.data)) { 
  sel.woody = which(trait.data[[i]]$specie_id %in% id[j])
  trait.data.woody[[i]] = trait.data[[i]][sel.woody,]
  trait.data.nonwoody[[i]] = trait.data[[i]][-sel.woody,]
}
names(trait.data.woody) = names(trait.data)
names(trait.data.nonwoody) = names(trait.data)


td.woody = td
td.woody$pft$outdir = "/home/carya/pecan/pft/gr.woody/"

td.nonwoody = td
td.nonwoody$pft$outdir = "/home/carya/pecan/pft/gr.nonwoody/"

## save
save(trait.data.woody,file=file.path(td.woody$pft$outdir, 'trait.data.Rdata'))
save(trait.data.nonwoody,file=file.path(td.nonwoody$pft$outdir, 'trait.data.Rdata'))

#run.meta.analysis()
run.meta.analysis(td.woody, settings$meta.analysis$iter, settings$run$dbfiles, settings$database)
load(file.path(settings$pfts$pft$outdir,"trait.mcmc.Rdata"))
load(file.path(settings$pfts$pft$outdir,"post.distns.Rdata"))

R.woody = getdistribution(trait.mcmc,post.distns,cost,variables)

run.meta.analysis(td.nonwoody, settings$meta.analysis$iter, settings$run$dbfiles, settings$database)
load(file.path(settings$pfts$pft$outdir,"trait.mcmc.Rdata"))
load(file.path(settings$pfts$pft$outdir,"post.distns.Rdata"))

R.nonwoody = getdistribution(trait.mcmc,post.distns,cost,variables)

#### Test probability distributions

ks.test(R.woody$Rl,R.nonwoody$Rl)
ks.test(R.woody$Rs,R.nonwoody$Rs)
ks.test(R.woody$Rr,R.nonwoody$Rr)

cols = 1:4
dR.woodynd = density(R.woody$Rnd)
plot(density(R.woody$Rl),xlim=range(dR.woodynd$x),col=cols[2])
lines(dR.woodynd,col=cols[1])
lines(density(R.woody$Rs),col=cols[3])
lines(density(R.woody$Rr),col=cols[4])
lines(density(R.nonwoody$Rnd),col=cols[1],lty=2)
lines(density(R.nonwoody$Rl),col=cols[2],lty=2)
lines(density(R.nonwoody$Rs),col=cols[3],lty=2)
lines(density(R.nonwoody$Rr),col=cols[4],lty=2)
legend("topright",legend=c("Null","Leaf","Stem","Root","Woody","Nonwoody"),col=c(cols,1,1),lwd=2,lty=c(1,1,1,1,1,2))

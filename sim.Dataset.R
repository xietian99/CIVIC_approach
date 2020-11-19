library(dplyr)
# Simulated Data Set 1
##---- For binary outcome
pi = c(0.2, 0.8)
tau.matrix = matrix(c(0.5, 0.1, 0.7, 0.1, 0.3, 0.1 ), nrow=2)

rownames(tau.matrix) = c("Class 1: 0.2", "Class 2: 0.8")
colnames(tau.matrix) = c("HA", "HAI", "int.NA")
n = 200 # 
set.seed(8090)
sim.data = matrix(rep(0, 4*n), ncol=4, nrow=n)
colnames(sim.data) = c("int.ha", "int.hai", "int.na","h1")

sim.class = sample(c(1,2), n, replace = T, prob = pi)

for (i in 1:n){
  cls = sim.class[i]
  for (j in 1:3) {
    sim.data[i,j] = rbinom(n=1, size = 1, prob = tau.matrix[cls,j])
  }
  cls.Prob =sim.data[i,1]+sim.data[i,2]
  if( cls.Prob == 0){ 
    sim.data[i,4] = rbinom(n=1, size = 1, prob = 0.9)
  } else if (cls.Prob == 2) {
    sim.data[i,4] = rbinom(n=1, size = 1, prob = 0.2)
  } else {
    if (sim.data[i,1] == 0) { 
      sim.data[i,4] = rbinom(n=1, size = 1, prob = 0.7)
    } else { 
      sim.data[i,4] = rbinom(n=1, size = 1, prob = 0.5) }
  }
}
##-----

#---- Simulate Dataset 2
##---For muti-level
pi = c(0.2, 0.8)
ParTable = matrix( c(1:12,
                   rep(0:2, each = 4),
                   rep(rep(0:1, each = 2), 3),
                   rep(0:1, 6), 
                   rep(c(0.9, 0.7, 0.3, 0.25, 0.15, 0.1), each = 2),
                   rep(0, 12),
                   rep(0, 12)
                   ),
                   ncol = 7
                   )
colnames(ParTable) = c("Group", "HAI", "HA", "NA", "Pr(H1=1)", "MeanH1", "Weight")
sim.class = sample(c(1,2), n, replace = T, prob = pi)
l11 = list(c(0.3,0.4,0.3))
l12 = list(c(0.5,0.5))
l13 = list(c(0.7,0.3))
l21 = list(c(0.8,0.15,0.05))
l22 = list(c(0.9,0.1))
l23 = list(c(0.9,0.1))
tau.matrix = matrix(c(l11,l12,l13,
                      l21,l22,l23), byrow = T, ncol =3)
rownames(tau.matrix) = c("Class 1: 0.2", "Class 2: 0.8")
colnames(tau.matrix) = c("HAI(0,1,2)", "HA(0,1)", "int.NA(0,1)")


n = 100000
set.seed(2763)
sim.data = matrix(rep(0, 4*n), ncol=4, nrow=n)
colnames(sim.data) = c("int.hai", "int.ha", "int.na","h1")
MarkerLevel = list(c(0:2), c(0:1), c(0,1))
sim.class = sample(c(1,2), n, replace = T, prob = pi)
for (i in 1:n){
  cls = sim.class[i]
  for (j in 1:3){
    sim.data[i,j] = sample(MarkerLevel[[j]], 
                           size = 1, 
                           prob = unlist(tau.matrix[cls, j]))
  }
  grp = 4*sim.data[i,1] + 2*sim.data[i,2] + 1*sim.data[i,3] + 1
  sim.data[i,4] = rbinom(n=1, size=1, prob = ParTable[grp,"Pr(H1=1)"])
}

ans = as.data.frame(sim.data) %>% 
  group_by(int.hai, int.ha, int.na) %>%  
  summarise(meanH1= mean(h1), weight = n()/n)

# For convient computing, bacause previous code we use dataset (ha, hai, na)
sim.data[,c(1,2,3,4)] = sim.data[,c(2,1,3,4)]
colnames(sim.data) = c("int.ha", "int.hai", "int.na","h1")
int2cls = function(x){return(as.factor(x))}
sim.data2 = as.data.frame(sim.data) %>% mutate_at(c("int.ha", "int.hai", "int.na"), int2cls) 

ParTable[,"Weight"] = ans$weight
ParTable[,"MeanH1"] = ans$meanH1
ParTable1 = ParTable
save.image("../SaveData/n1500.sim2.RData")

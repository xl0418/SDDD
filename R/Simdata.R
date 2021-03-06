# Arguments introduction
# n: number of locations in the model. Multiple locations model is already incorporated.
# parsN: a vector of number of species on each location at initial stage. In total, the number of species
#        in the system should be 2.
# For example: if set n = 2, parsN can be set as parsN = c(1,1) which means 1 species for each location
#              parsN = c(2,0) which means 2 species on location 1 and 0 species on location 2.
#              In this case, we can set migration rate to 0 to approach our spatial model to non-spatial model.
# age: crown age. As the simulation is conditioning on the survival of two ancestor specis, the crown age is fixed.
# pars: a vector of parameters. (lambda, mu, K'( or K)).
#              lambda: intrisic speciation rate
#              mu: extinction rate
#              K': ecological limit when argument K_fix = 0
#              K : carrying capacity when argument K_fix = 1
# seed: seeds for 100 simulation for this specific parameter setting.
# lambda_allo0: intrisic allopatric speciation rate.
# M0: intrisic migration rate.
# rn: the number of simulations for this specific parameter setting. It should be less than the maximum of seed.
# num: folder number in which the data is saving.
# K_fix: fix carrying capacity (K_fix = 1) or fix ecological limit (K_fix = 0)
# dirname: the directory in which you want to save the data. Please create a folder named by dirname under your main directory (getwd() can show the main directory).
#          Then the subfolders will be created automatically and named in this format "2Modeltestgroup15age112". "15" indicates the crown age. "112" indicates the group number.
#          The default setting for number of the groups is:
#          number: 101 - 109  Scenario 1 for extinction 0 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 111 - 119  Scenario 1 for extinction 0.1 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 121 - 129  Scenario 1 for extinction 0.2 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 131 - 139  Scenario 1 for extinction 0.4 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 141 - 149  Scenario 2 for extinction 0 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 151 - 159  Scenario 2 for extinction 0.1 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 161 - 169  Scenario 2 for extinction 0.2 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 171 - 179  Scenario 2 for extinction 0.4 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 181 - 189  Scenario 3 for extinction 0 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 191 - 199  Scenario 3 for extinction 0.1 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 201 - 209  Scenario 3 for extinction 0.2 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 211 - 219  Scenario 3 for extinction 0.4 varying migration rate from 0,0.05,0.1,0.15,0.3,0.5,1,5,1000
#          number: 221 - 224  Scenario 4 for migration 0 varying extinction rate from 0,0.1,0.2,0.4
#          number: 231 - 234  Scenario 5 for migration 0 varying extinction rate from 0,0.1,0.2,0.4


SimMdata <- function(n, parsN ,age=20,pars=c(0.8,0.3,20), seed = c(1:100), lambda_allo0=0.2, M0=0,rn=100,num,K_fix = 1,dirname){
  mainDir = dirname
  subDir = paste(n,"Modeltestgroup",age,"age",num,sep = "")
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  for (i in 1:rn){
    seed_fun = seed[i]
    result = sddsim(n = n, parsN=parsN,age=age,pars = pars, seed_fun = seed_fun,  lambda_allo0=lambda_allo0, M0=M0,K_fix = K_fix)
    file = paste(mainDir,"/",n,"Modeltestgroup",age,"age",num,"/out",i,"sim.Rdata",sep = "")
    save(result,file = file)
  }
}

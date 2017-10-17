
# Code for simulating all scenarios.
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

# x: if x = y, the package simulation starts.
simpackage = function(dirname){
for(base in c(100,140,180)){
  if(base == 100) pars = c(0.8,0,20)   #Scenario 1
  else if(base == 140) pars = c(0.8,0,40)    #Scenario 2
  else if(base == 180) pars = c(0.8,0,20,40)   #Scenario 3

  for(i in (base + 1):(base + 9)){
    M0 = c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000)
    SimMdata(n = 2, pars = pars ,parsN = c(1,1), age = 15, M0= M0[i-base],K_fix = 0,num = i,dirname = dirname)

  }
  pars[2] = pars[2]+0.1

  for(i in (base + 11):(base + 19)){
    M0 = c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000)

    SimMdata(n = 2, pars = pars ,parsN = c(1,1), age = 15, M0= M0[i-base-10],K_fix = 0,num = i,dirname = dirname)

  }

  pars[2] = pars[2]+0.1

  for(i in (base + 21):(base + 29)){
    M0 = c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000)

    SimMdata(n = 2, pars = pars ,parsN = c(1,1), age = 15, M0= M0[i-base-20],K_fix = 0,num = i,dirname = dirname)

  }
  pars[2] = pars[2]+0.2

  for(i in (base + 31):(base + 39)){
    M0 = c(0,0.05,0.1,0.15,0.3,0.5,1,5,1000)

    SimMdata(n = 2, pars = pars ,parsN = c(1,1), age = 15, M0= M0[i-base-30],K_fix = 0,num = i,dirname = dirname)

  }

}


mu = c(0,0.1,0.2,0.4)

# Scenario 4
for(i in 221:224){
  pars = c(0.8,0,20)
  pars[2] = mu[i-220]
  SimMdata(n = 2, pars = pars ,parsN = c(2,0), age = 15, M0= 0,K_fix = 0,num = i,dirname = dirname)

}
# Scenario 5
for(i in 231:234){
  pars = c(0.8,0,40)
  pars[2] = mu[i-230]
  SimMdata(n = 2, pars = pars ,parsN = c(2,0), age = 15, M0= 0,K_fix = 0,num = i,dirname = dirname)

}
  }


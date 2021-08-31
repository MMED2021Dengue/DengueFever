MinAgeSample <- 4
MaxAgeSample <- 11
Pop <- 10000
MaxAge <- 99

FoIs <- seq(0, .1, length = 1001)
TotalDiseases <- rep(0, 1001)
ObservedSeroprevalence <- rep(0, 1001)

for (run in 1:1001){
  FoI <- FoIs[run]
  
  ImmuneHistory <- matrix(NA, MaxAge, 5)
  ImmuneHistory[1, 1] <- Pop
  ImmuneHistory[1, 2:5] <- 0
  
  ProbInfection <- c(1 - (1 - FoI)^4,
                     1 - (1 - FoI)^3,
                     1 - (1 - FoI)^2,
                     FoI,
                     0)
  
  for (age in 2:MaxAge){
    ImmuneHistory[age, 1] <- ImmuneHistory[age - 1, 1] * (1 - ProbInfection[1])
    for (infection in 1:4){
      ImmuneHistory[age, infection + 1] <- ImmuneHistory[age - 1, infection] * ProbInfection[infection] + 
                                              ImmuneHistory[age - 1, infection + 1] * (1 - ProbInfection[infection + 1])
    }
  }
  # ImmuneHistory
  # rowSums(ImmuneHistory)
  
  SeroPrevalence <- 1 - ImmuneHistory[,1] / Pop
  ObservedSeroprevalence[run] <- mean(SeroPrevalence[MinAgeSample:MaxAgeSample])
  
  InfectionCount <- matrix(0, MaxAge, 4)
  for (age in 2:MaxAge){
    for (infection in 1:4){
      InfectionCount[age, infection] <- ImmuneHistory[age - 1, infection] * ProbInfection[infection]
    }
  }
  
  # InfectionCount
  # Other option
  # ProbOfDisease <- c(0, .05, 0, 0)
  ProbOfDisease <- c(0.005, 0.05, 0.0005, 0.0005)
  
  DiseaseByInfectionNumber <- rep(0,4)
  
  for (infection_num in 1:4){
    DiseaseByInfectionNumber[infection_num] <- sum(InfectionCount[,infection_num]) * ProbOfDisease[infection_num]
  }
  TotalDiseases[run] <- sum(DiseaseByInfectionNumber)

}
par(mfrow=c(1,2))
plot(FoIs, TotalDiseases, type='l')
plot(FoIs, ObservedSeroprevalence, type='l')


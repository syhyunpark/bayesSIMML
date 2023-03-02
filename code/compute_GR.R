# load("C:/Users/Danni/OneDrive - NYU Langone Health/BayesSIMML/results_testlinearcode/v7brms/all_24_scenarios_brms.rda")
load("all_24_scenarios_brms.rda")  # read-in the saved data 

results <- unlist2d(results.aggregated2, idcols = "replicate",DT = TRUE)#the first columns: scenario's id;the second column:simulation's id for each secnario
pri_scen <- function(x){
  paste0(x["replicate"],". ","n=",x["n"],","," p=",x["p"],","," g.choice=",x["g.choice"],",",
         " m.choice=",x["m.choice"])
}
results$scenario <- apply(results,1,pri_scen)
class(results)
# save(results,file="./all_24_scenarios_table.rda")
# # load("./all_24_scenarios_table.rda")

ans <- results[,
               .(beta.psrf.1.03=mean(GR.beta.psrf.1.03),
                 gamma.psrf.1.03=round(mean(GR.gamma.psrf.1.03),3),
                 m.psrf.1.03=mean(GR.m.psrf.1.03),
                 
                 beta.psrf.1.05=mean(GR.beta.psrf.1.05),
                 gamma.psrf.1.05=round(mean(GR.gamma.psrf.1.05),3),
                 m.psrf.1.05=mean(GR.m.psrf.1.05),
                 
                 beta.psrf.1.07=mean(GR.beta.psrf.1.07),
                 gamma.psrf.1.07=round(mean(GR.gamma.psrf.1.07)),
                 m.psrf.1.07=mean(GR.m.psrf.1.07)),
               keyby = .(m.choice, g.choice, p,n)]

ans
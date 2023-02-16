#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
require(rms)
library(fst)
library(gridExtra)
library(grid)

#INPUT
args = commandArgs()



print(args)                                                                                                

# logistic_and_quantile_compare.R $task.cpus ${ENH_list} ${cohort_fam} \
#     ${TS_ENH_GWAS_compartment_originalOR_summary} ${TS_ENH_GWAS_compartment_originalOR_best}\
#     ${TS_ENH_GWAS_compartment_OR_by_measure1_summary} ${TS_ENH_GWAS_compartment_OR_by_measure1_best}\
#     ${TS_ENH_GWAS_compartment_OR_by_measure2_summary} ${TS_ENH_GWAS_compartment_OR_by_measure2_best}\
#     ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
#     ${merged_GWAS_summary} ${merged_GWAS_best}\
#     ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  \
#     ${original_LOO_GWAS_summary} ${original_LOO_GWAS_prsice} ${original_LOO_GWAS_best}\
#     ${modif_name_1} ${modif_name_2}
nthreads = as.numeric(args[8])
#set max CPU processes
setDTthreads(nthreads)
threads_fst(nr_of_threads = round(nthreads/3*2))

(ENH_list = args[9])
(diagnosis = fread(args[10], header=F, col.names = c("FID", "IID", "IIDf", "IIDm", "sex", "dx" )) %>%
    dplyr::select("FID", "IID", "dx"))

TS_ENH_GWAS_compartment_originalOR_summary = args[11]
TS_ENH_GWAS_compartment_originalOR_best = 
  fread(args[12], select=c("FID", "IID", "PRS")) %>% 
  dplyr::rename(TS_ENH_GWAS_compartment_originalOR_best_PRS = PRS)

TS_ENH_GWAS_compartment_OR_by_measure1_summary = args[13]
TS_ENH_GWAS_compartment_OR_by_measure1_best = 
  fread(args[14], select=c("FID", "IID", "PRS")) %>% 
  dplyr::rename(TS_ENH_GWAS_compartment_OR_by_measure1_best_PRS = PRS)

TS_ENH_GWAS_compartment_OR_by_measure2_summary = args[15]
TS_ENH_GWAS_compartment_OR_by_measure2_best = 
  fread(args[16], select=c("FID", "IID", "PRS")) %>% 
  dplyr::rename(TS_ENH_GWAS_compartment_OR_by_measure2_best_PRS = PRS)

residual_GWAS_compartment_summary = args[17]
residual_GWAS_compartment_best = 
  fread(args[18], select=c("FID", "IID", "PRS"))  %>% 
  dplyr::rename(residual_GWAS_compartment_best_PRS = PRS)

merged_GWAS_summary = args[19]
merged_GWAS_best = 
  fread(args[20], select=c("FID", "IID", "PRS")) %>% 
  dplyr::rename(merged_GWAS_best_PRS = PRS)

# ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  
TS_ENH_GWAS_compartment_originalOR_prsice = args[21]
TS_ENH_GWAS_compartment_OR_by_measure1_prsice = args[22]
TS_ENH_GWAS_compartment_OR_by_measure2_prsice = args[23]
residual_GWAS_compartment_prsice = args[24]
merged_GWAS_prsice = args[25]

original_LOO_GWAS_summary = args[26]
original_LOO_GWAS_prsice = args[27]
(original_LOO_GWAS_best = fread(args[28], select=c("FID", "IID", "PRS"))  %>% 
  dplyr::rename(original_LOO_GWAS_best_PRS = PRS))

modif_name_1 = args[29]
modif_name_2 = args[30]
threshold = args[31]
condition_name = args[32] #SCZ or HCM

#set input variables
number_quantiles = 5

# pop_prev = population prevalence
pop_prev = ifelse(test = condition_name == "SCZ", yes = 0.01, no = 0.007)



#OUTPUT
analysis_output_txt = paste0(ENH_list, "_",threshold, "_", modif_name_1,"_", modif_name_2,"_", Sys.Date(),"_model_fit.txt")
model_fit_plot      = paste0(ENH_list, "_",threshold, "_", modif_name_1,"_", modif_name_2,"_", Sys.Date(),"_MODEL_FIT_PLOT.pdf")
PRS_double_QUANTILE_PLOT  = paste0(ENH_list, "_",threshold, "_", modif_name_1,"_", modif_name_2,"_", Sys.Date(),"_PRS_double_QUANTILE_PLOT.pdf")
# PRS_comparison_figure_path = paste0(ENH_list, "_",threshold, "_", modif_name_1,"_", modif_name_2,"_", Sys.Date(), "_PRS_comparison_plot.pdf")
CoD_per_SNP_plot_scaled= paste0(ENH_list, "_",threshold, "_", modif_name_1,"_", modif_name_2,"_", Sys.Date(), "_scaled_CoD_per_snp_plot.pdf")
CoD_ALL_plots = paste0(ENH_list, "_",threshold, "_", modif_name_1,"_", modif_name_2,"_", Sys.Date(), ".pdf")




r_color <- colors()


#SUMMARY TABLE
#import thresholds and SNP N for each summary
(summary_table = 
    rbind(
      "TS_ENH_GWAS_compartment_originalOR_summary" = 
        data.frame(fread(TS_ENH_GWAS_compartment_originalOR_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
      "TS_ENH_GWAS_compartment_OR_by_measure1_summary" = 
        data.frame(fread(TS_ENH_GWAS_compartment_OR_by_measure1_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
      "TS_ENH_GWAS_compartment_OR_by_measure2_summary" = 
        data.frame(fread(TS_ENH_GWAS_compartment_OR_by_measure2_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
      "residual_GWAS_compartment_summary"=
        data.frame(fread(residual_GWAS_compartment_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
      "original_LOO_GWAS_summary"=
        data.frame(fread(original_LOO_GWAS_summary, select=c("Threshold", "PRS.R2","Num_SNP")))#"PRS.R2.adj",
    ) %>%  rownames_to_column(var = "compartment"))



#BEST TABLE
#create total PRS score
(BEST_PRS_score_per_UKBB_participant <- TS_ENH_GWAS_compartment_originalOR_best %>%
    left_join(TS_ENH_GWAS_compartment_OR_by_measure1_best) %>%
    left_join(TS_ENH_GWAS_compartment_OR_by_measure2_best) %>%
    left_join(residual_GWAS_compartment_best) %>%
    # left_join(merged_GWAS_best) %>%
    left_join(original_LOO_GWAS_best) %>%
    left_join(diagnosis) %>%
    mutate(dx=factor(dx), IID=factor(IID)) %>%
    select(-FID) %>%
    remove_missing() #%>% head(n=50000)
)


(scaled_BEST_PRS_score_per_UKBB_participant <- BEST_PRS_score_per_UKBB_participant)
# fwrite(BEST_PRS_score_per_UKBB_participant,"BEST_PRS_score_per_UKBB_participant.txt")
scaled_BEST_PRS_score_per_UKBB_participant[,c(2:6)] <-  data.frame(scale(BEST_PRS_score_per_UKBB_participant[,c(2:6)], center = T, scale = T))+10
# head(scaled_BEST_PRS_score_per_UKBB_participant)
# fwrite(scaled_BEST_PRS_score_per_UKBB_participant,"scaled_BEST_PRS_score_per_UKBB_participant.txt")

(scaled_BEST_PRS_score_per_UKBB_participant<-
    scaled_BEST_PRS_score_per_UKBB_participant %>% 
    # mutate(weight_total_PRS_best = 
    #          (TS_ENH_GWAS_compartment_originalOR_best_PRS * summary_table[summary_table$compartment=="TS_ENH_GWAS_compartment_originalOR_summary",]$Num_SNP  /summary_table[summary_table$compartment=="merged_GWAS_summary",]$Num_SNP) + 
    #          (residual_GWAS_compartment_best_PRS     * summary_table[summary_table$compartment=="residual_GWAS_compartment_summary",]$Num_SNP      /summary_table[summary_table$compartment=="merged_GWAS_summary",]$Num_SNP)) %>%
    remove_missing() %>% 
    #generate quantiles
    mutate(original_GWAS_q = 
             factor(ntile(original_LOO_GWAS_best_PRS, n = number_quantiles))) %>% 
    # mutate(merged_GWAS_q = 
    #          factor(ntile(merged_GWAS_best_PRS, n = number_quantiles))) %>% 
    mutate(residual_GWAS_compartment_q = 
             factor(ntile(residual_GWAS_compartment_best_PRS, n = number_quantiles))) %>% 
    mutate(TS_ENH_compartment_originalOR_q = 
             factor(ntile(TS_ENH_GWAS_compartment_originalOR_best_PRS, number_quantiles))) %>% 
    # mutate(weight_total_q = 
    #          factor(ntile(weight_total_PRS_best, number_quantiles)))%>% 
    mutate(TS_ENH_compartment_OR_by_measure1_q = 
             factor(ntile(TS_ENH_GWAS_compartment_OR_by_measure1_best_PRS, number_quantiles)))%>% 
    mutate(TS_ENH_compartment_OR_by_measure2_q = 
             factor(ntile(TS_ENH_GWAS_compartment_OR_by_measure2_best_PRS, number_quantiles)))
)





# str(scaled_BEST_PRS_score_per_UKBB_participant)
# scaled_BEST_PRS_score_per_UKBB_participant[rowSums(is.na(scaled_BEST_PRS_score_per_UKBB_participant)) > 0,]


#### measures of model fit 
## https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.21614
# nt = total number of the sample
(nt = NROW(diagnosis))
# ncase = number of cases
(ncase = NROW(diagnosis[diagnosis$dx==2,]))
# ncont = number of controls
(ncont = NROW(diagnosis[diagnosis$dx==1,]))
# case_prev_in_sample = proportion of cases in the case-control samples
(case_prev_in_sample = ncase/nt)
# thd = the threshold on the normal distribution which truncates the proportion of disease prevalence
(thd = -qnorm(pop_prev,0,1))
(zv = dnorm(thd)) #z (normal density)
(mv = zv/pop_prev) #mean liability for case
(mv2 = -mv*pop_prev/(1-pop_prev)) #mean liability for controls

#R2 on the observed scale
theta = mv*(case_prev_in_sample-pop_prev)/(1-pop_prev)*(mv*(case_prev_in_sample-pop_prev)/(1-pop_prev)-thd) #θ in equation 15
cv = pop_prev*(1-pop_prev)/zv^2*pop_prev*(1-pop_prev)/(case_prev_in_sample*(1-case_prev_in_sample)) #C inequation 15


# Start writing to an output file
(CoD_per_SNP = data.frame())
summary_table

# sink(analysis_output_txt)
# ## original_GWAS
lrm(dx ~ original_LOO_GWAS_best_PRS, 
                                      data = scaled_BEST_PRS_score_per_UKBB_participant)
#logistic model
(pmv = glm(dx ~ original_LOO_GWAS_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant,family = binomial(probit))) #probit model
summary(pmv)
# R2 on the liability scale using the transformation
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(original_GWAS_logistic_model_R2 = R2O*cv/(1+R2O*theta*cv))
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  c("0",original_GWAS_logistic_model_R2,
    summary_table[summary_table$compartment=="original_LOO_GWAS_summary",]$Num_SNP, NA)
))


# (merged_GWAS_logistic_model <- lrm(dx ~ merged_GWAS_best_PRS, 
#                                    data = scaled_BEST_PRS_score_per_UKBB_participant))
# #logistic model
# # (lmv = glm(dx ~ merged_GWAS_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(logit)))
# pmv = glm(dx ~ merged_GWAS_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant,family = binomial(probit)) #probit model
# # R2 on the liability scale using the transformation
# (R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
# #R2 on the observed scale
# (merged_GWAS_logistic_model_R2 = R2O*cv/(1+R2O*theta*cv))
# (CoD_per_SNP = rbind(
#   CoD_per_SNP,
#   c("X",merged_GWAS_logistic_model_R2,
#     summary_table[summary_table$compartment=="merged_GWAS_summary",]$Num_SNP, NA)
# ))

##residual
lrm(dx ~ residual_GWAS_compartment_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant)
#logistic model
pmv = glm(dx ~ residual_GWAS_compartment_best_PRS, 
           data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(probit)) #logistic model
summary(pmv)
# R2 on the liability scale using the transformation
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(residual_GWAS_compart_logistic_model_R2 = R2O*cv/(1+R2O*theta*cv))
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  c("1",residual_GWAS_compart_logistic_model_R2,summary_table[summary_table$compartment=="residual_GWAS_compartment_summary",]$Num_SNP, NA)
))

# TS ENH original OR
lrm(dx ~ TS_ENH_GWAS_compartment_originalOR_best_PRS, data = BEST_PRS_score_per_UKBB_participant)
#logistic model
(pmv = glm(dx ~ TS_ENH_GWAS_compartment_originalOR_best_PRS, 
           data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(probit))) #logistic model
summary(pmv)
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(TS_ENH_originalOR_compart_logistic_model_R2 = R2O*cv/(1+R2O*theta*cv))
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  c("2",TS_ENH_originalOR_compart_logistic_model_R2,summary_table[summary_table$compartment=="TS_ENH_GWAS_compartment_originalOR_summary",]$Num_SNP, NA)
))

# TS ENH  OR by measure 1
lrm(dx ~ TS_ENH_GWAS_compartment_OR_by_measure1_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant)
#logistic model
(pmv = glm(dx ~ TS_ENH_GWAS_compartment_OR_by_measure1_best_PRS, 
           data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(probit))) #logistic model
summary(pmv)
# R2 on the liability scale using the transformation
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(TS_ENH_OR_by_measure1_compart_logistic_model_R2 = R2O*cv/(1+R2O*theta*cv))
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  c("2b",TS_ENH_OR_by_measure1_compart_logistic_model_R2,
    summary_table[summary_table$compartment=="TS_ENH_GWAS_compartment_OR_by_measure1_summary",]$Num_SNP, NA)
))

# TS ENH  OR by measure 2
lrm(dx ~ TS_ENH_GWAS_compartment_OR_by_measure2_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant)
#logistic model
(pmv = glm(dx ~ TS_ENH_GWAS_compartment_OR_by_measure2_best_PRS, 
           data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(probit))) #logistic model
# R2 on the liability scale using the transformation
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(TS_ENH_OR_by_measure2_compart_logistic_model_R2 = R2O*cv/(1+R2O*theta*cv))
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  c("2c",TS_ENH_OR_by_measure2_compart_logistic_model_R2,
    summary_table[summary_table$compartment=="TS_ENH_GWAS_compartment_OR_by_measure2_summary",]$Num_SNP, NA)
))


##simple additive model
lrm(dx ~ residual_GWAS_compartment_best_PRS + TS_ENH_GWAS_compartment_originalOR_best_PRS, 
    data = scaled_BEST_PRS_score_per_UKBB_participant)
#logistic model
pmv = glm(dx ~ residual_GWAS_compartment_best_PRS + 
             TS_ENH_GWAS_compartment_originalOR_best_PRS, 
           data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(probit)) #logistic model
summary(pmv)
# R2 on the liability scale using the transformation
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(residual_GWAS_plus_TS_ENH_originalOR_logistic_model_R2 = R2O*cv/(1+R2O*theta*cv))


##logistic full factorial design 
lrm(dx ~ residual_GWAS_compartment_best_PRS*TS_ENH_GWAS_compartment_originalOR_best_PRS, 
                                                                 data = scaled_BEST_PRS_score_per_UKBB_participant)
#logistic model
pmv = glm(dx ~ residual_GWAS_compartment_best_PRS*TS_ENH_GWAS_compartment_originalOR_best_PRS, 
           data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(probit))
# R2 on the liability scale using the transformation
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(logistic_full_factorial_design_model_R2 = R2O*cv/(1+R2O*theta*cv))


#logistic full factorial design including interactions and non-linear interactions
pmv = glm(dx ~ residual_GWAS_compartment_best_PRS*TS_ENH_GWAS_compartment_originalOR_best_PRS +
             I(residual_GWAS_compartment_best_PRS^2) + I(TS_ENH_GWAS_compartment_originalOR_best_PRS^2), 
           data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial(probit))
summary(pmv)
# R2 on the liability scale using the transformation
(R2O = var(pmv$fitted.values)/(ncase/nt*ncont/nt))
#R2 on the observed scale
(logistic_full_factorial_design_nonlinear_interactions_model_R2 = R2O*cv/(1+R2O*theta*cv))

# Stop writing to the file
# sink()

# CoD PER SNP df
colnames(CoD_per_SNP)=c("partition","CoD","Num_SNP","CoD_per_SNP")
CoD_per_SNP[c(2:4)]<-sapply(CoD_per_SNP[c(2:4)],as.numeric)
CoD_per_SNP$CoD_per_SNP = (CoD_per_SNP$CoD / CoD_per_SNP$Num_SNP)*10^5
CoD_per_SNP

addline_format <- function(x,...){
  gsub('\\s|__','\n',x)
}

#create DF for plotting
(df_plot<- data.frame(
  partition=c(factor(c("0",
                       #"1",
                       "1","2","2b","2c","3","3b",
                       "3c"
                       ))),
  partition_name= factor(c("0-original_GWAS",
                           #"1-merged_GWAS",
                           "1-residual_GWAS", 
                           "2-TS_ENH original_OR",
                           paste0("2b-TS_ENH_OR* ",modif_name_1),
                           paste0("2c-TS_ENH_OR* ",modif_name_2),
                           "3-residual_GWAS+ TS_ENH_original_OR",
                           "3b-residual_GWAS*TS_ENH_original_OR (FFD)",
                           "3c-residual_GWAS*TS_ENH_original_OR FFD+quadratic_factors"
  )),
  R2=c(original_GWAS_logistic_model_R2,
       #merged_GWAS_logistic_model_R2, #1
       residual_GWAS_compart_logistic_model_R2,#2 
       TS_ENH_originalOR_compart_logistic_model_R2, #3
       TS_ENH_OR_by_measure1_compart_logistic_model_R2,
       TS_ENH_OR_by_measure2_compart_logistic_model_R2,
       residual_GWAS_plus_TS_ENH_originalOR_logistic_model_R2, #4
       logistic_full_factorial_design_model_R2,logistic_full_factorial_design_nonlinear_interactions_model_R2
  ),
            colour=c("navyblue",
                     "green",
                     "green",
                     "green",
                     "green",
                     "darkgreen",
                     "darkgreen",                   "darkgreen"
  )) %>% left_join(CoD_per_SNP, by="partition") %>% dplyr::select(-CoD)  %>% 
    mutate_at(c("Num_SNP"), ~replace_na(.,-1)) %>%
    mutate(xlabel=factor(paste0(addline_format(partition_name), "\nN_SNP ", Num_SNP)))
)


p1 <- ggplot(data = df_plot, aes(
  y=reorder(xlabel, desc(xlabel)), 
  x=R2, 
  label=round(R2,4))
) +  
  geom_segment(aes(yend = reorder(xlabel, desc(xlabel))), xend = 0, colour = "grey50") +                      
  geom_point(color=df_plot$colour, size=4) + 
  ggrepel::geom_label_repel(size = rel(3), fill = "azure", col="black",
                            hjust = 1, nudge_y = -0.2, point.padding = NA, box.padding = 0.5)+ 
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .15))) + 
  scale_y_discrete(expand = expansion(add = .7)) + #.6 is default - add a little spacing to the sides
  xlab("") +  ylab("")+
  theme_bw() +
  theme(
    axis.text.y = element_text(lineheight = 0.8, angle = 25, size = rel(1)),
    plot.margin = margin(t = 0, r = 0, b = 1, l = 0, "cm"),
    panel.grid.major.y = element_blank()
  )


f1<-arrangeGrob(textGrob("A)", just = "left",
                         gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                textGrob(paste("Coefficients of determination"), 
                         gp = gpar(fontsize = 16, fontface = "bold", col="darkgreen")), 
                #textGrob("diagnosis ~ PRS, probit link function \nProportion of the total variance explained by the genetic factor on the liability scale, \ncorrected for ascertainment, as per Lee et al 2012", gp = gpar(fontsize = 10)), 
                p1,
                layout_matrix=rbind(c(1,2),
                                    c(3,3)),
                widths = c(0.05, 1), heights = c(0.05, 1))


p2 <-ggplot(data = df_plot[!is.na(df_plot$Num_SNP),], 
            aes(
              y=reorder(xlabel, desc(xlabel)),
              x=CoD_per_SNP, 
              label=round(CoD_per_SNP,3))) +  
  geom_segment(aes(yend = reorder(xlabel, desc(xlabel))), xend = 0, colour = "grey50") +
  geom_point(color=df_plot$colour, size=4) + 
  ggrepel::geom_label_repel(size = rel(3), fill = "azure", col="black",
                            hjust = 1, nudge_y = -0.2, point.padding = NA, box.padding = 0.5)+
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) + 
  xlab("") +  ylab("")+
  theme_bw() +
  theme(
    axis.text.y = element_text(lineheight = 0.8, angle = 25, size = rel(1)),
    plot.margin = margin(t = 0, r = 1, b = 1, l = 0, "cm"),
    panel.grid.major.y = element_blank()
  )


f2<-arrangeGrob(textGrob("B)", just = "left",
                         gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                textGrob(paste("Coefficients of determination per SNP * 10^5"), 
                         gp = gpar(fontsize = 16, fontface = "bold", col="darkblue")), 
                #textGrob("diagnosis ~ PRS, probit link function \nProportion of the total variance explained by the genetic factor on the liability scale, \ncorrected for ascertainment, as per Lee et al 2012", gp = gpar(fontsize = 10)), 
                p2,
                layout_matrix=rbind(c(1,2),
                                    c(3,3)),
                widths = c(0.05, 1), heights = c(0.05, 1))


#save all plots
ggsave(filename = model_fit_plot, f1,  width = 8, height = 7)
ggsave(filename = paste0(ENH_list, "_",threshold, "_", modif_name_1,"_", modif_name_2,"_", Sys.Date(),"_MODEL_FIT_PLOT_per_snp.pdf"), 
        f2,  width = 8, height = 7)



#scaled plot
scale1 <- function(x) scale(x, center = F)[,1]
p3 <- df_plot%>% 
  mutate(scaled_CoD = R2, scaled_N = Num_SNP, scaled_CoD_per_SNP = CoD_per_SNP) %>% 
  mutate_at(vars(contains("scaled")), scale1) %>% drop_na() %>% 
  dplyr::select("partition_name","scaled_CoD","scaled_N","scaled_CoD_per_SNP") %>% 
  pivot_longer(!partition_name) %>% 
  ggplot(aes(x=value, fill=name, y=addline_format(partition_name))) + 
  geom_col(position = "dodge") +
  scale_fill_brewer(name="",palette="Set1")+#coord_flip()+
  theme_bw()+
  xlab("") +  ylab("")+
  theme(
    axis.text.y = element_text(lineheight = 0.8, angle = 25, size = rel(1)),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, "cm"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )


f3<-arrangeGrob(textGrob("C)", just = "left",
                         gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                textGrob(paste("Relative number of SNPs, total CoD, and CoD per SNP"), 
                         gp = gpar(fontsize = 16, fontface = "bold", col="darkblue")), 
                #textGrob("diagnosis ~ PRS, probit link function \nProportion of the total variance explained by the genetic factor on the liability scale, \ncorrected for ascertainment, as per Lee et al 2012", gp = gpar(fontsize = 10)), 
                p3,
                layout_matrix=rbind(c(1,2),
                                    c(3,3)),
                widths = c(0.05, 1), heights = c(0.05, 1))








# OR based plots
# double quantile plot for interactions

#https://cran.r-project.org/web/packages/samplesizeCMH/vignettes/samplesizeCMH-introduction.html
#https://stats.stackexchange.com/questions/593123/can-i-add-up-ors-for-specific-predictors/593130#593130
#https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704-ep713_confounding-em/BS704-EP713_Confounding-EM7.html
scaled_BEST_PRS_score_per_UKBB_participant

# CALCULATE ORS
#merged_GWAS_q_OR
summary(logistic<-glm(formula = dx ~ original_GWAS_q, 
                      data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial, na.action = "na.omit"))
(original_GWAS_q_OR<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
colnames(original_GWAS_q_OR) <- c("quantile", "OR", "LCI", "UCI")
original_GWAS_q_OR

(ORs <- original_GWAS_q_OR)
ORs[1,]<-list("1",1,1,1)
ORs[,1]<-list(1:nrow(ORs))
ORs$original_OR_quant <- paste("all_thresh:",summary_table[summary_table$compartment=="original_LOO_GWAS_summary",c("Threshold")],"_N=", summary_table[summary_table$compartment=="original_LOO_GWAS_summary",c("Num_SNP")])
ORs


ORs_original_OR <- ORs
for  (i in 1:number_quantiles) {
  print(i)
  logistic<-glm(formula = dx ~ original_GWAS_q, 
                data = scaled_BEST_PRS_score_per_UKBB_participant[scaled_BEST_PRS_score_per_UKBB_participant$TS_ENH_compartment_originalOR_q==i,], 
                family = binomial, na.action = "na.omit")
  (OR<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
  OR[1,]<-list("1",1,1,1)
  OR[,1]<-list(1:nrow(OR))
  colnames(OR) <- c("quantile", "OR", "LCI", "UCI")
  OR$original_OR_quant <- paste0("TS_ENH_q_",i)
  
  ORs_original_OR<-rbind(ORs_original_OR,OR)
}
ORs_original_OR$comp=paste0("2) TS_ENH: Original OR",#"__thresh_",
                            #summary_table[summary_table$compartment=="TS_ENH_GWAS_compartment_originalOR_summary",c("Threshold")],
                            "\nSNP N=",NROW(scaled_BEST_PRS_score_per_UKBB_participant[scaled_BEST_PRS_score_per_UKBB_participant$TS_ENH_compartment_originalOR_q==i,]))
ORs_original_OR


ORs_OR_by_ES <- ORs
for  (i in 1:number_quantiles) {
  print(i)
  logistic<-glm(formula = dx ~ original_GWAS_q, 
                data = scaled_BEST_PRS_score_per_UKBB_participant[scaled_BEST_PRS_score_per_UKBB_participant$TS_ENH_compartment_OR_by_measure1_q==i,], 
                family = binomial, na.action = "na.omit")
  (OR<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
  OR[1,]<-list("1",1,1,1)
  OR[,1]<-list(1:nrow(OR))
  colnames(OR) <- c("quantile", "OR", "LCI", "UCI")
  OR$original_OR_quant <- paste0("TS_ENH_q_",i)
  
  
  ORs_OR_by_ES<-rbind(ORs_OR_by_ES,OR)
}
ORs_OR_by_ES$comp=paste0("2b) TS_ENH OR * ",modif_name_1,#"__thresh_",
                         #summary_table[summary_table$compartment=="TS_ENH_GWAS_compartment_OR_by_measure1_summary",c("Threshold")],
                         "\nSNP N=",NROW(scaled_BEST_PRS_score_per_UKBB_participant[scaled_BEST_PRS_score_per_UKBB_participant$TS_ENH_compartment_OR_by_measure1_q==i,]))
ORs_OR_by_ES


ORs_OR_by_exp <- ORs
for  (i in 1:number_quantiles) {
  print(i)
  logistic<-glm(formula = dx ~ original_GWAS_q, 
                data = scaled_BEST_PRS_score_per_UKBB_participant[scaled_BEST_PRS_score_per_UKBB_participant$TS_ENH_compartment_OR_by_measure2_q==i,], 
                family = binomial, na.action = "na.omit")
  (OR<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
  OR[1,]<-list("1",1,1,1)
  OR[,1]<-list(1:nrow(OR))
  colnames(OR) <- c("quantile", "OR", "LCI", "UCI")
  OR$original_OR_quant <- paste0("TS_ENH_q_",i)
  
  ORs_OR_by_exp<-rbind(ORs_OR_by_exp,OR)
}
ORs_OR_by_exp$comp=paste0("2c) TS_ENH OR * ",modif_name_2,#"__thresh_",
                          #summary_table[summary_table$compartment=="TS_ENH_GWAS_compartment_OR_by_measure2_summary",c("Threshold")],
                          "\nSNP N=",NROW(scaled_BEST_PRS_score_per_UKBB_participant[scaled_BEST_PRS_score_per_UKBB_participant$TS_ENH_compartment_OR_by_measure2_q==i,]))
ORs_OR_by_exp

(all_ORs<-rbind(
  ORs_original_OR,ORs_OR_by_ES,ORs_OR_by_exp) %>% 
    mutate(
      comp=factor(comp),
      comp=forcats::fct_relevel(comp)
    ))
table(all_ORs$comp)


# pdf(file = PRS_double_QUANTILE_PLOT, width = 11, height = 7)
p4 = ggplot(data = all_ORs , aes(y= OR, ymin = LCI, ymax=UCI, 
                                 x=factor(quantile), colour=original_OR_quant, group=original_OR_quant)) + 
  facet_wrap(facets = vars((comp)))+
  scale_colour_manual(name="ENH compartment quantile", values = c("tomato","#ccece6", "#99d8c9", "#41ae76","#006d2c", "#00441b",r_color))+
  geom_pointrange(position = position_dodge(width = 0.3))  + 
  ylab(paste0("OR for ",condition_name))+   xlab('Original PRS quantile')+
  # labs(title =  paste("Participant distribution by HCM OR by original PGC GWAS quantile\nand further by", ENH_list, "quantile"))+ 
  theme_bw() +
  theme(
    strip.text.x = element_text(size = rel(1.3)),
    axis.text = element_text(size = rel(1)),
    axis.title = element_text(size = rel(1.5)),
    plot.margin = margin(t = 0, r = 1, b = 0, l = 1, "cm"),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

# dev.off()
f4<-arrangeGrob(textGrob("D)", just = "left",
                         gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                textGrob(paste("Participant distribution by ",condition_name," OR by original PGC GWAS quantile\nand further by", ENH_list, "quantile"), 
                         gp = gpar(fontsize = 16, fontface = "bold",col="maroon")), 
                p4,
                layout_matrix=rbind(c(1,2),
                                    c(3,3)),
                widths = c(0.05, 1), heights = c(0.1, 1))


ggsave(filename = CoD_per_SNP_plot_scaled, arrangeGrob(f3, f4, ncol = 2),  width = 17, height = 7)


## save overall plot, adding threshold, date, 
ggsave(filename = paste0("CoD_ALL_plots_",CoD_ALL_plots), 
       arrangeGrob(
         textGrob(ENH_list,
                  gp = gpar(fontsize = 18, fontface = "bold",col="navyblue")),
         f1, f2, f3, f4, 
         textGrob(paste("PRSice --score avg - C+T threshold:",threshold,"plot created on:",Sys.Date()),
                  gp = gpar(fontsize = 10, col="maroon")),
         layout_matrix=rbind(c(1,1),
                             c(2,3),
                             c(4,5),
                             c(6,6)),
         heights = c(0.08, 0.5, 0.5, 0.04)
       ),  
       width = 17, height = 14)

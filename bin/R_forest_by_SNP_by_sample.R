#!/usr/bin/env Rscript

library(tidyverse)
library(metafor)
library(data.table)

#INPUT
args = commandArgs()

# print(args)

meta_analytic_results_filelist = readLines(args[8])


(SNPlist = c("rs1888998","rs198817","rs397553","rs34380086","rs55682899","rs759168","rs198817","rs11717383","rs2213806"))#))#"rs10921330","rs2491405","rs1447514","rs2894249","rs436031","rs6017135","rs4962696","rs10901813","rs12117347","rs1811852","rs9531841","rs11934877","rs7327579","rs2834230","rs10759370","rs13008496","rs12510182","rs11053603","rs4821101","rs3129716","rs72717438","rs140540991","rs2067098","rs385303","rs35341285","22:45820268_CG_C","rs7015469","rs12462897","rs4912782","5:64796837_CT_C","rs13131819","rs12993075","rs2972553","rs322353","rs4426626","rs8118848","rs12548201","rs6737950","rs558721608","rs1713073","rs1682809",,"rs6960438","rs4971212","rs12130986","rs6656884","rs926734","rs2787350","rs8070990","rs2651077","rs2617865","rs2143594","rs6935354","rs552180","rs1554457","rs881375","rs12475055","rs59237912","rs237897","rs2216258","rs1870073","rs4652564","rs2151119","rs2169501","rs2301220","rs73136048","rs2809342","rs73368107","rs17826073","rs11651521","rs368186","rs7972413","rs2413047","rs71311393","rs6765857","rs9838091","rs9268832","rs1480380","rs72694957","rs10514971","rs242557","rs12741781","rs10411755","rs62652294","rs6038300","rs75219778","rs11617366","rs10926978","rs196588","rs59513960","rs7657658","rs7148171","rs77080586","rs62086904","rs7032218","rs9616819"


df <- data.frame(cohort=character(),
                MARKERNAME=character(), 
                 OR=numeric(), OR_95L=numeric(), OR_95U=numeric(), 
                 N=numeric(), 
                 EA=character()) 

for (SNP in SNPlist) {
    print(SNP)
    for (i in 1:length(meta_analytic_results_filelist)) {
        # print(i)
        # print(meta_analytic_results_filelist[[i]])
        
        # #MARKERNAME	EA	NEA	OR	OR_95L	OR_95U	N	EAF	STRAND
        cohort = #str_split(meta_analytic_results_filelist[[i]], "_",simplify = T)[1]
                str_match(meta_analytic_results_filelist[[i]],"REC_ALL_BRAIN_EPs_([\\w]+)_eur_associations_GWAMA_format")[2]
        # pop = str_split(meta_analytic_results_filelist[[i]], "_",simplify = T)[2]
        
        print(d <- fread(file=meta_analytic_results_filelist[[i]], select = c("MARKERNAME","OR","OR_95L","OR_95U","N","EA")) %>%
            dplyr::filter(MARKERNAME==SNP))
        # check if df not empty (no snp in dataset, 0 rows in d)
        if(!is.null(d)){
            # print("not null")
            if(dim(as.data.frame(d))[1]>0){
                df = rbind(
                    df,
                    cbind(cohort,  d[1,])
                )
            }
        }
    }
}
#remove rows with missing ORs and make ORs and CIs numeric
(df<- df %>% mutate(across(OR:OR_95U, as.numeric)))
#remove rows with Inf or NA
# df = df[complete.cases(df[c("OR","OR_95L","OR_95U")]),]
df<-df[is.finite(rowSums(df[,c("OR","OR_95L","OR_95U")])),]
#also remove very large values (v uncertain estimates)
df<-df[rowSums(df[,c("OR","OR_95L","OR_95U")])<1000,]

fwrite( file= paste0(Sys.Date(),"_meta_results_per_study_per_SNP.csv"), x=df )


#check that every SNP has same EA across cohorts and remove cohorts with only 1 dataset
(to_include = df %>% mutate(rs_allele=paste0(MARKERNAME,"_",EA)) %>%
    group_by(rs_allele) %>% count(EA) %>% filter(n>1) %>% ungroup())
(df = df %>% mutate(rs_allele=paste0(MARKERNAME,"_",EA)) %>%
    filter(rs_allele %in% unique(to_include$rs_allele)))

#metafor
for(snippo in unique(df$rs_allele)){
    # library(dplyr)
    print(snippo)
    
    print(df_meta <- df %>% 
        #filter by SNP
        dplyr::filter(rs_allele==snippo) %>%
        #yi is log(OR)
        mutate(yi= log(OR), slab=cohort) %>% 
        #5) The confidence interval bounds can be converted into the standard errors with: 
        #https://www.metafor-project.org/doku.php/tips:assembling_data_or
        mutate(sei= (log(OR_95U) - log(OR_95L))/(2*1.96)) %>%
        # 6) Any missing values for the vi variable (the sampling variances) can now be replaced with: 
        mutate(vi = sei^2, zi=NULL, sei = NULL)%>%
        mutate(across(c(OR:OR_95U, yi, vi), as.numeric) ) 
    )
    

    print(summary(res <- rma(yi, vi, method="EE", slab=slab, data=df_meta)))
    # print(predict(res, transf=exp, digits=2))
    # saveRDS(res, paste0(snippo,"res.rds"))

    png(filename=paste0(snippo,"_",Sys.Date(),"_by_cohort.png"), width = 2500, height = 3000, res=200)
    forest(res, atransf=exp, slab=slab, order="obs",
        at=log(c(.05, .25, 1, 4)), alim=log(c(0.01,10)), xlim=c(-5,5),
        refline=0, cex=.9, header=TRUE,
        #title
        main=paste(snippo, "recessive association with schizophrenia in PGC cohorts + UKBB"))
    dev.off()
    
}
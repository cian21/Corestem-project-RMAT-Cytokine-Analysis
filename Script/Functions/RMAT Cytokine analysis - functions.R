# Corestem project - RMAT Cytokine Analysis.proj
# 20210216 Sehwan Chun at Corestem, Inc.
# 2.1. cleaning process for RMAT Cytokine analysis

#### 1. Library Loading ####
packs = c("data.table", "ggpubr", "leaps")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. cpm calculation ####
reg_file_cleaning = function(reg_file, impute = TRUE, month){
    reg_file = subset(reg_file, is.na(TGFB1_V5) == FALSE)
    row.names(reg_file) = 1:nrow(reg_file)
    
    reg_file = reg_file[,-c(35:40)] #TNF-A, IL-1B REMOVED
    reg_file$X1M_slope = trunc((reg_file$X1M - reg_file$X0M.baseline.) / 1 * 10^3) / 10^3
    reg_file$X2M_slope = trunc((reg_file$X2M - reg_file$X0M.baseline.) / 2 * 10^3) / 10^3
    reg_file$X3M_slope = trunc((reg_file$X3M - reg_file$X0M.baseline.) / 3 * 10^3) / 10^3
    reg_file$X4M_slope = trunc((reg_file$X4M - reg_file$X0M.baseline.) / 4 * 10^3) / 10^3
    if(impute == TRUE){
        reg_file$X6M_slope = trunc((as.numeric(reg_file$X6M_Imp) - reg_file$X0M.baseline.) / 6 * 10^3) / 10^3
    }else{
        reg_file$X6M_slope = trunc((as.numeric(reg_file$X6M) - reg_file$X0M.baseline.) / 6 * 10^3) / 10^3
    }
    reg_file$MINUS3M_slope = trunc((reg_file$X0M.baseline. - reg_file$minus_3M) / 3 * 10^3) / 10^3
    reg_file$progression_rate_change_6m = reg_file$X6M_slope - reg_file$MINUS3M_slope
    reg_file$progression_rate_change_4m = reg_file$X4M_slope - reg_file$MINUS3M_slope
    
    reg_file$GR6 = ifelse(reg_file$X6M_slope >= (reg_file$MINUS3M_slope * 0.5), "GR","PR")
    reg_file$GR4 = ifelse(reg_file$X4M_slope >= (reg_file$MINUS3M_slope * 0.5), "GR","PR")
    
    #reg_file = reg_file[reg_file$GR == "GR",]
    assign("reg_file", reg_file, envir = globalenv())
    
    reg_x_values_v5 = reg_file[,c(8,seq(20,37,3))]
    assign("reg_x_values_v5", reg_x_values_v5, envir = globalenv())
    
    if(month == 4){
        reg_x_values_diff_percent = reg_file[,seq(22,37,3)] / reg_file[,seq(20,37,3)] * 100
        reg_x_values_diff_percent = cbind(reg_file["progression_rate_change_4m"], reg_x_values_diff_percent)
        assign("reg_x_values_diff_percent", reg_x_values_diff_percent, envir = globalenv())
        
        GR4_reg = subset(reg_file, GR4 == "GR")
        assign("GR4_reg", GR4_reg, envir = globalenv())
        
        GR4_reg_x_values_v5 = GR4_reg[,c(8,seq(20,37,3))]
        assign("GR4_reg_x_values_v5", GR4_reg_x_values_v5, envir = globalenv())
        
        GR4_reg_x_values_diff_percent = GR4_reg[,seq(22,37,3)] / GR4_reg[,seq(20,37,3)] * 100
        GR4_reg_x_values_diff_percent = cbind(GR4_reg["progression_rate_change_4m"], GR4_reg_x_values_diff_percent)
        assign("GR4_reg_x_values_diff_percent", GR4_reg_x_values_diff_percent, envir = globalenv())
    }else if(month == 6){
        reg_x_values_diff_percent = reg_file[,seq(22,37,3)] / reg_file[,seq(20,37,3)] * 100
        reg_x_values_diff_percent = cbind(reg_file["progression_rate_change_6m"], reg_x_values_diff_percent)
        assign("reg_x_values_diff_percent", reg_x_values_diff_percent, envir = globalenv())
        
        GR6_reg = subset(reg_file, GR6 == "GR")
        assign("GR6_reg", GR6_reg, envir = globalenv())
        
        GR6_reg_x_values_v5 = GR6_reg[,c(8,seq(20,37,3))]
        assign("GR6_reg_x_values_v5", GR6_reg_x_values_v5, envir = globalenv())
        
        GR6_reg_x_values_diff_percent = GR6_reg[,seq(22,37,3)] / GR6_reg[,seq(20,37,3)] * 100
        GR6_reg_x_values_diff_percent = cbind(GR6_reg["progression_rate_change_6m"], GR6_reg_x_values_diff_percent)
        assign("GR6_reg_x_values_diff_percent", GR6_reg_x_values_diff_percent, envir = globalenv())
    }
}

    


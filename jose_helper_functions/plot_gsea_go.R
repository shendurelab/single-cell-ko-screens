plot_gsea_go<-function(gsa_res, 
                       mode="distinct", 
                       top_g_gene_sets=15, 
                       q_thresh=0.05, 
                       fill_color="black", 
                       whitelist=NULL, 
                       plot_q_vals=FALSE)
{
  # pvl_mat <- get_gene_set_p_vals(detection_df, go_gs, alternative=alternative)
  resTab <- GSAsummaryTable(gsa_res)
  #resTab <- as.data.frame(gsa_res$resTab)
  #colnames(resTab) <- c("p_value", "q_value", "inc_in_set", "excl_in_set", "incl_not_in_set", "excl_not_in_set")
  #resTab$gene_set_name <- row.names(gsa_res$resTab)
  # pvl_bp_go <- arrange(pvl_bp_go, desc(-log10(q_value)))
  # pvl_bp_go <- transform(pvl_bp_go, gene_set_name=reorder(gene_set_name, -log10(q_value)) ) 
  
  res_df <- resTab[,c("Name", "Genes (tot)", "Genes (up)", "Genes (down)")]
  
  if (mode == "distinct"){
    resTab <- resTab[with(resTab, order(-log10(resTab[,"p (dist.dir.up)"]))), ]
    res_df$stat_up <- resTab[,"Stat (dist.dir.up)"]
    res_df$stat_dn <- -resTab[,"Stat (dist.dir.dn)"]
    res_df$p_up <- resTab[,"p (dist.dir.up)"]
    res_df$p_dn <- resTab[,"p (dist.dir.dn)"]
    res_df$q_up <- resTab[,"p adj (dist.dir.up)"]
    res_df$q_dn <- resTab[,"p adj (dist.dir.dn)"]
  }
  else if (mode == "mixed"){
    resTab <- resTab[with(resTab, order(-log10(resTab[,"p (mix.dir.up)"]))), ]
    res_df$stat_up <- resTab[,"Stat (mix.dir.up)"]
    res_df$stat_dn <- -resTab[,"Stat (mix.dir.dn)"]
    res_df$p_up <- resTab[,"p (mix.dir.up)"]
    res_df$p_dn <- resTab[,"p (mix.dir.dn)"]
    res_df$q_up <- resTab[,"p adj (mix.dir.up)"]
    res_df$q_dn <- resTab[,"p adj (mix.dir.dn)"]
  }
  
  
  #if (is.null(whitelist)){
  to_draw_df_dn <- arrange(subset(res_df, q_dn < 1), desc(q_dn))
  to_draw_df_dn <- to_draw_df_dn[,c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat_dn", "p_dn", "q_dn")]
  colnames(to_draw_df_dn) <- c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat", "p_val", "q_val")
  to_draw_df_dn$q_val <- log10(to_draw_df_dn$q_val + 1e-3)
  
  to_draw_df_up <- arrange(subset(res_df, q_up < 1), desc(q_up))
  to_draw_df_up <- to_draw_df_up[,c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat_up", "p_up", "q_up")]
  colnames(to_draw_df_up) <- c("Name", "Genes (tot)", "Genes (up)", "Genes (down)", "stat", "p_val", "q_val")
  to_draw_df_up$q_val <- -log10(to_draw_df_up$q_val + 1e-3)
  
  res_df <- rbind(to_draw_df_up, to_draw_df_dn)
  #print (res_df)
  
  # }else{
  #   print(subset(pvl_bp_go, gene_set_name %in% whitelist))
  #   res_df <- subset(resTab, gene_set_name %in% whitelist)
  # }
  
  res_df <- arrange(res_df, desc(q_val))
  res_df <- transform(res_df, Name=reorder(Name, q_val) ) 
  
  if (is.null(whitelist) == FALSE){
    g <- qplot(Name, q_val, data=res_df, size=I(0.65), color = Name %in% whitelist) 
    
    up_gene_sets <- subset(res_df, stat > 0 & Name %in% whitelist)
    if (nrow(up_gene_sets) > 0){
      g <- g + geom_text(aes(Name, q_val, label=Name, hjust="right", check_overlap=TRUE), data=up_gene_sets)
    }
    down_gene_sets <- subset(res_df, stat < 0 & Name %in% whitelist)
    if (nrow(down_gene_sets) > 0){
      g <- g + geom_text(aes(Name, q_val, label=Name, hjust="left", check_overlap=TRUE), data=down_gene_sets)
    }
    
    print (subset(res_df, Name %in% whitelist))
  }
  else{
    g <- qplot(Name, q_val, data=res_df)
  }
  g <- g + theme(axis.text.x = element_blank()) + 
    theme(axis.ticks.x = element_blank()) 
  g <- g + geom_hline(yintercept=-log10(q_thresh)) 
  g <- g + geom_hline(yintercept=log10(q_thresh)) 
  g <- g + scale_color_manual(values=c("black", "red"))
  g <- g + theme(legend.position="none") 
  #g <- g + theme(plot.background = element_blank()) 
  g <- g + theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA))
  #g <- g + theme(panel.border = element_blank(), axis.line = element_line(size=0.2)) 
  
  g
}
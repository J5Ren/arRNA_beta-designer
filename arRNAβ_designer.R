library(optparse)
library(pheatmap)
library(ggplot2)
library(ggrastr)
library(egg)
library(dplyr)
parser <- OptionParser()
parser <- add_option(parser, c("-s", "--seq"), help="sequence")
parser <- add_option(parser, c("-t", "--taget"), help="target site (1 base)")
RNA <- toupper(parse_args(parser)$seq) #input: mRNA sequence
RNA=strsplit(RNA,'')[[1]]
RNA[RNA%in%"T"]="U"
RNA=paste0(RNA,collapse = "")
site <- as.numeric(parse_args(parser)$taget) #input: target A position

# 1. Potential Bystander Site ---------------------------------------------
All_A <- as.numeric(gregexpr("A", RNA)[[1]])
RNA_len <- nchar(RNA)
All_A <- setdiff(All_A,c(1:4,(RNA_len-19):RNA_len)) #Left 4, right 20 excluded
pos=c()
for(i in 1:length(All_A)){
  if(substr(RNA,All_A[i]-1,All_A[i]-1)%in%"G"){
    pos=c(pos)
  }else{
    pos=c(pos,i)
  }
}
All_A=All_A[pos]#GA excluded
All_A=setdiff(All_A,site)
df <- data.frame(letters = unlist(strsplit(RNA, "")), Position = 1:RNA_len,
                 bystander=rep('no',RNA_len))
df[All_A,'bystander']='yes'
pdf('1.potential bystander site.pdf',width = RNA_len/3)
ggplot(df, aes(x = Position, y = 1)) +
  geom_text(aes(label = letters, color = factor(bystander)), size = 10) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(color = "potential bystander sites")+
  ylab("")+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

All_A_minwindow=All_A[which((All_A>=site-4)&(All_A<=site+20))]

# 151nt arRNA candidate ---------------------------------------------------
start <- max(site + 20 - 150,
             1)
end <- min(site - 4,
           RNA_len - 150)
start_range=start:end
df=data.frame(id=1:(end-start+1),
              bystander_num=0)
for(i in 1:nrow(df)){
  df[i,2]=length(All_A[All_A%in%c(start_range[i]:(start_range[i]+150))])
}
candidate=df[df$bystander_num%in%min(df$bystander_num),][1,1]
candidate_window=start_range[candidate]:(start_range[candidate]+150)
candidate_bystander=All_A[All_A%in%candidate_window]
All_A=candidate_bystander
All_A_minwindow=All_A_minwindow[All_A_minwindow%in%candidate_window]
# candidate out bulge -------------------------------------------------------------------
Eliminate_bystander_out <- function(bystander_position) {
  bulge3_start_range <- (bystander_position - 6):(bystander_position + 20)
  bulge4_start_range <- (bystander_position - 7):(bystander_position + 20)
  bulge_df <- data.frame(
    start = c(bulge3_start_range, bulge4_start_range),
    end = c(bulge3_start_range + 2, bulge4_start_range + 3)
  )
  bulge_df <- bulge_df[bulge_df$start >= min(candidate_window) & bulge_df$end <= max(candidate_window), ]
  bulge_df <- bulge_df[order(bulge_df$start, -bulge_df$end), ]
  keep <- logical(nrow(bulge_df))
  last_end <- -Inf
  for (i in seq_len(nrow(bulge_df))) {
    if (bulge_df$end[i] > last_end) {
      keep[i] <- TRUE
      last_end <- bulge_df$end[i]
    }
  }
  res=bulge_df[keep, ]
  res$id=paste(res$start,res$end,sep='-')
  res=res[!duplicated(res$id),]
  res=res[,-3]
  return(res)
}
bulge_out=list()
All_A=setdiff(All_A,All_A_minwindow)
for(pos in c(All_A,site)){
  bulge_out[[as.character(pos)]]=Eliminate_bystander_out(pos)
}
combined_bulge_out <- bind_rows(bulge_out, .id = "source")
combined_bulge_out$affect_target=0
combined_bulge_affect_target=combined_bulge_out[combined_bulge_out$source%in%as.character(site),]
affect_target_site_out <- unique(unlist(
  lapply(1:nrow(combined_bulge_affect_target), function(i) {
    seq(combined_bulge_affect_target[i, "start"], combined_bulge_affect_target[i, "end"])
  })
))
for(i in 1:nrow(combined_bulge_out)){
  if(any(affect_target_site >= combined_bulge_out[i,'start'] & affect_target_site <= combined_bulge_out[i,'end'])){
    combined_bulge_out[i,'affect_target']=1
  }
}
combined_bulge_out=combined_bulge_out[combined_bulge_out$affect_target%in%0,]
# candidate in bulge ------------------------------------------------------
Eliminate_bystander_in <- function(bystander_position) {
  bulge_start_range <- (bystander_position - 1):(bystander_position + 15)
  bulge_df <- data.frame(
    start = bulge_start_range,
    end = bulge_start_range
  )
  bulge_df <- bulge_df[bulge_df$start >= min(candidate_window) & bulge_df$end <= max(candidate_window), ]
  bulge_df <- bulge_df[order(bulge_df$start, -bulge_df$end), ]
  keep <- logical(nrow(bulge_df))
  last_end <- -Inf
  for (i in seq_len(nrow(bulge_df))) {
    if (bulge_df$end[i] > last_end) {
      keep[i] <- TRUE
      last_end <- bulge_df$end[i]
    }
  }
  res=bulge_df[keep, ]
  res$id=paste(res$start,res$end,sep='-')
  res=res[!duplicated(res$id),]
  res=res[,-3]
  return(res)
}
bulge_in=list()
for(pos in c(All_A_minwindow,site)){
  bulge_in[[as.character(pos)]]=Eliminate_bystander_in(pos)
}
combined_bulge_in <- bind_rows(bulge_in, .id = "source")
combined_bulge_in$affect_target=0
combined_bulge_in[combined_bulge_in$start%in%combined_bulge_in[combined_bulge_in$source%in%as.character(site),'start'],'affect_target']=1
combined_bulge_in=combined_bulge_in[combined_bulge_in$affect_target%in%0,]
# bulge_selection ---------------------------------------------------------
combined_bulge <- rbind(combined_bulge_in,combined_bulge_out)
row_counts <- combined_bulge %>%
  group_by(across(-source)) %>%  
  summarise(
    count = n(),  
    sources = paste(unique(source), collapse = ","),  
    .groups = "drop"
  ) %>%
  arrange(desc(count))
set_list <- lapply(row_counts$sources, function(row) {
  unlist(strsplit(row, ","))  # 按逗号分割并转为向量
})
names(set_list)=1:nrow(row_counts)
find_minimum_cover <- function(set_list, target_set) {
  selected_sets <- list()
  selected_names <- character()
  remaining_chars <- target_set
  while (length(remaining_chars) > 0) {
    coverage <- sapply(set_list, function(s) length(intersect(s, remaining_chars)))
    best_index <- which.max(coverage)
    best_name <- names(set_list)[best_index]
    selected_sets <- append(selected_sets, list(set_list[[best_index]]))
    selected_names <- c(selected_names, best_name)
    remaining_chars <- setdiff(remaining_chars, set_list[[best_index]])
    set_list <- set_list[-best_index]
  }
  return(selected_names)
}


bulge_final=row_counts[find_minimum_cover(set_list, c(All_A,All_A_minwindow)),]


# arRNA generate ---------------------------------------------------------
bulge_final_arRNA=bulge_final
bulge_final_arRNA$start=as.numeric(RNA_len-bulge_final$end+1)
bulge_final_arRNA$end=as.numeric(RNA_len-bulge_final$start+1)
site=RNA_len-as.numeric(site)+1
bulge_final_arRNA=bulge_final_arRNA[order(bulge_final_arRNA$start),]
MyRevComRNASeq=function(rna_sequence){
  rna_sequence=toupper(rna_sequence)
  complement <- chartr("ACGU", "UGCA", rna_sequence)
  reverse_complement <- paste(rev(strsplit(complement, "")[[1]]), 
                              collapse = "")
  return(reverse_complement)
}
arRNA=MyRevComRNASeq(RNA)
arRNA=strsplit(arRNA,'')[[1]]
arRNA[site]="C"
candidate_window=RNA_len-candidate_window+1
arRNA=arRNA[candidate_window[order(candidate_window)]]
bulge_sites=c()
for(i in 1:nrow(bulge_final_arRNA)){
  bulge_sites=c(bulge_sites,as.numeric(bulge_final_arRNA[i,1]):as.numeric(bulge_final_arRNA[i,2]))
}
#arRNA=paste0(arRNA[setdiff(1:length(arRNA),bulge_sites)],collapse = '') #deletion
mismatch=function(base){
  if(base%in%"A"){
    return("U")
  }else if(base%in%"U"){
    return("A")
  }else if(base%in%"C"){
    return("G")
  }else{
    return("C")
  }
}
arRNA[bulge_sites]=sapply(arRNA[bulge_sites],mismatch)
arRNA=paste0(arRNA,collapse = '') #mismatch

options <- commandArgs()

IPS <- options[which(grepl("-IPS", options))+1]
nog <- options[which(grepl("-eggnog", options))+1]
arg <- options[which(grepl("-deeparg", options))+1]

tmhmmfile <- options[which(grepl("-hmm", options))+1]

sigpfile <- options[which(grepl("-sigp", options))+1]

labelling_file <- options[which(grepl("-clusts", options))+1]

libloc <- "/home/projects/group-a/bin/miniconda3/lib/R/library"

library(data.table, lib.loc = libloc)
library(Biostrings, lib.loc = libloc)

#Interpro has multiple annots per prot on account of the different DBs... needs to have either multiple cols or filters.
interproscan <- read.csv(IPS, sep = "\t", header = F)

deeparg <- read.csv(arg, sep = "\t", header = T, stringsAsFactors = F)


eggnog <- read.csv(nog, sep = "\t", skip = 3, stringsAsFactors = F)
eggnog <- eggnog[eggnog$seed_ortholog_evalue <= 10e-5,]

tmhmm <- read.csv(tmhmmfile, sep = "\t", stringsAsFactors = F, header = F)

signalP <- read.csv(sigpfile, sep = "\t", skip = 1, header = F, stringsAsFactors = F)

original_seqs <- read.csv(labelling_file, sep = "\t", header = F, stringsAsFactors = F)



representatives <- as.character(original_seqs$V9[original_seqs$V10 == "*"])
representatives <- unlist(strsplit(as.character(representatives), split=" "))
representatives <- representatives[c(T,F,F)]

seq_names <- unlist(strsplit(as.character(original_seqs$V9), split=" "))
seq_names <- seq_names[c(T,F,F)]

translator <- data.frame(seq = as.character(original_seqs$V9), represent = as.character(original_seqs$V10), stringsAsFactors = F)
translator$represent[translator$represent == "*"] <- translator$seq[translator$represent == "*"]

translator$proper_name <- unlist(strsplit(translator$represent, split = " "))[c(T,F,F)]
translator_proper_seq <- seq_names

translator$sample <-  unlist(strsplit(seq_names, split=";"))[c(T,F)]
translator$samp_ID <- unlist(strsplit(seq_names, split=";"))[c(F,T)]


eggnog_output <- eggnog[match(translator$proper_name, eggnog$X.query_name),]
eggnog_output$X.query_name <- translator$seq
eggnog_output$sample <- translator$sample
eggnog_output$seq_in_samp <- translator$samp_ID
eggnog_output <- eggnog_output[!is.na(eggnog_output$seed_eggNOG_ortholog),]

eggnog_output_uniq <- unique(eggnog_output$X.query_name)
eggnog_output <- eggnog_output[match(eggnog_output_uniq, eggnog_output$X.query_name),]


interpro_output <- interproscan[match(translator$proper_name, interproscan$V1),]
interpro_output$V1 <- translator$seq
interpro_output$sample <- translator$sample
interpro_output$seq_in_samp <- translator$samp_ID
interpro_output <- interpro_output[!is.na(interpro_output$V3),]

interpro_output_uniq <- unique(interpro_output$V1)
interpro_output <- interpro_output[match(interpro_output_uniq, interpro_output$V1),]


deeparg_output <- deeparg[match(translator$proper_name, deeparg$read_id),]

deeparg_output$read_id <- translator$seq
deeparg_output$sample <- translator$sample
deeparg_output$seq_in_samp <- translator$samp_ID
deeparg_output <- deeparg_output[!is.na(deeparg_output$query.start),]


deeparg_output_uniq <- unique(deeparg_output$read_id)
deeparg_output <- deeparg_output[match(deeparg_output_uniq, deeparg_output$read_id),]


signalP_output <- signalP[match(translator$proper_name, signalP$V1),]
signalP_output$V1 <- translator$seq
signalP_output$sample <- translator$sample
signalP_output$seq_in_samp <- translator$samp_ID
signalP_output <- signalP_output[!is.na(signalP_output$V2),]

signalP_output_uniq <- unique(signalP_output$V1)
signalP_output <- signalP_output[match(signalP_output_uniq, signalP_output$V1),]

hmm_output <- tmhmm[match(translator$proper_name, tmhmm$V1),]

hmm_output <- lapply(unique(translator$proper_name), function(x){
  return(tmhmm[tmhmm$V1 == x,])
})




names(hmm_output) = unique(translator$proper_name)


bigger_hmm_out <- hmm_output[translator$proper_name]
names(bigger_hmm_out) = translator$seq

bigger_hmm_out <- lapply(1:length(bigger_hmm_out), function(x){
  tmp <- bigger_hmm_out[[x]]
  tmp$V1 <- names(bigger_hmm_out)[x]
  return(tmp)
})

bigger_hmm_out <- rbindlist(bigger_hmm_out)

bigger_hmm_out$sample <- translator$sample[match(bigger_hmm_out$V1, translator$seq)]
bigger_hmm_out$seq_in_samp <- translator$samp_ID[match(bigger_hmm_out$V1, translator$seq)]

bigger_hmm_out_uniq <- unique(bigger_hmm_out$V1)
bigger_hmm_out <- bigger_hmm_out[match(bigger_hmm_out_uniq, bigger_hmm_out$V1),]


#printing



colnames(interpro_output)[1:13] = c("key", "sequence_MD5", "sequence_representative_length", "best_hit_database", "accession", "description", "start", "stop", "e-value_score", "match_status", "date_of_run", "interpro_accession", "interpro_description")

print("IP_done")


deeparg_output <- deeparg_output[,c(1:3, 11, 13:14)]
colnames(deeparg_output)[4] = c("key")





colnames(eggnog_output)[1] = "key"

signalP_output <- signalP_output[,c(1,3,6,9, 10,11)]
colnames(signalP_output)[c(1:4)] = c("key", "annotation", "signal_peptide_score", "is_tat")

hmm_output <- bigger_hmm_out[,c(1,3,4,5,6)]
colnames(hmm_output)[c(1:3)] = c("key", "annotation", "gene_window")


dir.create("interpro_out")
dir.create("deeparg_out")
dir.create("eggnog_out")
dir.create("signalp_out")
dir.create("tmhmm_out")

split_and_print <- function(dataframe, output_dir, type){
  
  samps <- unique(dataframe$sample)
  
  lapply(samps, function(x){
    
    write.table(dataframe[dataframe$sample == x,], paste0(output_dir, "/", x, "_", type, "_", "functional_annotation.tsv"), sep = "\t", row.names = F, quote = F)
    
  })
  
  
}


split_and_print(interpro_output, "interpro_out", "interpro")
split_and_print(deeparg_output, "deeparg_out", "deeparg")
split_and_print(signalP_output, "signalp_out", "signalp")
split_and_print(hmm_output, "tmhmm_out", "tmhmm")
split_and_print(eggnog_output, "eggnog_out", "eggnog")




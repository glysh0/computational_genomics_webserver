options <- commandArgs()

epidata <- options[which(grepl("-epidat", options))+1]
mlst <- options[which(grepl("-mlst", options))+1]
argpath <- options[which(grepl("-deeparg", options))+1]

if(!file.exists(epidata)){
  print("provide epidata file absolute path with -epidat")
}
if(!file.exists(mlst)){
  print("provide mlst file absolute path with -mlst")
}

libloc <- "/home/projects/group-a/bin/miniconda3/lib/R/library"

libloc <- .libPaths()

#Really just needs epidata and MLST tsv and deeparg out dir...

#No idea how to handle clustering.

library(data.table, lib.loc = libloc)
library(gplots, lib.loc = libloc)
library(lubridate, lib.loc = libloc)
library(ggplot2, lib.loc = libloc)


chewbacca_data <- fread(mlst, sep = "\t", header = T)


rows <- chewbacca_data$FILE

rows <- unname(substr(rows, 1, 7))

chewbacca_matrix <- as.matrix(chewbacca_data[,2:ncol(chewbacca_data)])

rownames(chewbacca_matrix) = rows

colnames(chewbacca_matrix) = rep("", ncol(chewbacca_matrix))

#hclust_bacca <- heatmap.2(chewbacca_matrix)

pdf("MLST_Clusters_Per_Sample.pdf")
heatmap(chewbacca_matrix, Colv =  NA)
dev.off()


epigen_data <- fread(epidata, sep = "\t")

epigen_data$SampleDate <-  as.Date(epigen_data$SampleDate)

food1 <- ggplot(epigen_data, aes(x = Location, y = SampleDate))+
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_text(aes(label = Food1), hjust = 1, nudge_x = -0.07) +
  ylab("Date of Sample") +
  ggtitle("Outbreak Timeline", subtitle = "First food") 

food2 <- ggplot(epigen_data, aes(x = Location, y = SampleDate))+
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_text(aes(label = Food2), hjust = 1, nudge_x = -0.07) +
  ylab("Date of Sample") +
  ggtitle("Outbreak Timeline", subtitle = "Second food")

food3 <- ggplot(epigen_data, aes(x = Location, y = SampleDate))+
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_text(aes(label = Food3), hjust = 1, nudge_x = -0.07) +
  ylab("Date of Sample") +
  ggtitle("Outbreak Timeline", subtitle = "Third food")

food4 <- ggplot(epigen_data, aes(x = Location, y = SampleDate))+
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 45)) +
  geom_text(aes(label = Food4), hjust = 1, nudge_x = -0.07) +
  ylab("Date of Sample") +
  ggtitle("Outbreak Timeline", subtitle = "Fourth food")

pdf("Outbreak_epidata_timeline.pdf", width = 16, height = 9)
print(food1)
print(food2)
print(food3)
print(food4)
dev.off()



arg_data <- list.files(path = argpath, full.names = T)

arg_data <- lapply(arg_data, function(x){
  tmp <- fread(x, sep = "\t")
  tmp$source <- substr(x, nchar(x)-67, nchar(x)-61)
  tmp <- tmp[identity >= 40 & alignment.evalue <= 0.0001,]
  return(tmp)
})

arg_data <- rbindlist(arg_data)

arg_class_by_sample <- arg_data[, list(list(sort(unique(predicted_ARG.class)))), by = source]

classes <- unique(arg_class_by_sample$V1)

arg_class_by_sample$class <- match(arg_class_by_sample$V1, classes)

epigen_data$arg_class <- arg_class_by_sample$class[match(epigen_data$ID, arg_class_by_sample$source)]
epigen_data$args <- arg_class_by_sample$V1[match(epigen_data$ID, arg_class_by_sample$source)]

p <- ggplot(epigen_data, aes(x = Location, y = SampleDate, color = factor(arg_class)))+
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 45)) +
  ylab("Date of Sample") +
  ggtitle("Outbreak Timeline", subtitle = "ARGs")+
  guides(color = guide_legend(title = "Identical ARG\nprofile group"))

pdf("Epidata_and_ARG_profiles.pdf", width = 16, height= 9)
print(p)
dev.off()


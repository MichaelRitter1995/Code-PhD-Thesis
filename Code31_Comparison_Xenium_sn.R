#Run Comparison_Xenium_single_cell.R first.

###Analyse Tumour distribution

fin_Tum <- data.frame()
Diff_fin <- data.frame()

for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "Tumour" & subset$Type == "sn"),2])
  sn_Tum <-   sum(subset[which(subset$Annotation == "Tumour" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "Tumour" & subset$Type == "Xenium"),2])
  Xenium_Tum <-   sum(subset[which(subset$Annotation == "Tumour" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("Tumour", "Healthy"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_Tum/(sn_Tum+sn_sum), sn_sum/(sn_Tum+sn_sum), Xenium_Tum/(Xenium_Tum+Xenium_sum), Xenium_sum/(Xenium_Tum+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_Tum <- rbind(fin_Tum, subset)
  Diff <- data.frame(sample, Xenium_Tum/(Xenium_Tum+Xenium_sum)-sn_Tum/(sn_Tum+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin <- rbind(Diff_fin, Diff)
}

p01 <- ggplot(fin_Tum, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("grey", "#6B3F98"))+ 
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave(plot = p01, filename = "Tumour_barplot.pdf")

Xenium_ss <- subset(fin_Tum,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_Tum,  Type == "sn", Count,
                drop = TRUE)
# Plot paired data

fin_Tum1 <- fin_Tum[fin_Tum$Annotation=="Tumour",]

p02 <-ggplot(fin_Tum1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot=p02, filename = "Tumour_change.pdf")


wilcox.test(fin_Tum[which(fin_Tum$Type=="sn" & fin_Tum$Annotation=="Tumour"),2], fin_Tum[which(fin_Tum$Type=="Xenium" & fin_Tum$Annotation=="Tumour"),2], paired = TRUE, alternative = "two.sided")


ggplot(data=Diff_fin, aes(x=Sample, y=Difference))+geom_bar(stat='identity')+theme_classic()



###Analyse Neuron distribution

fin_Neur <- data.frame()
Diff_fin_Neur <- data.frame()
for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "Neurons" & subset$Type == "sn"),2])
  sn_Neur <-   sum(subset[which(subset$Annotation == "Neurons" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "Neurons" & subset$Type == "Xenium"),2])
  Xenium_Neur <-   sum(subset[which(subset$Annotation == "Neurons" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("Neurons", "Other"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_Neur/(sn_Neur+sn_sum), sn_sum/(sn_Neur+sn_sum), Xenium_Neur/(Xenium_Neur+Xenium_sum), Xenium_sum/(Xenium_Neur+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_Neur <- rbind(fin_Neur, subset)
  Diff <- data.frame(sample, Xenium_Neur/(Xenium_Neur+Xenium_sum)-sn_Neur/(sn_Neur+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin_Neur <- rbind(Diff_fin_Neur, Diff)
}

p1 <-ggplot(fin_Neur, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("#1FAD4B", "grey"))+
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave(plot = p1, filename = "Neurons_barplot.pdf")

Xenium_ss <- subset(fin_Neur,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_Neur,  Type == "sn", Count,
                drop = TRUE)

fin_Neur1 <- fin_Neur[fin_Neur$Annotation=="Neurons",]

p2 <-ggplot(fin_Neur1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot=p2, filename = "Neurons_change.pdf")

wilcox.test(fin_Neur[which(fin_Neur$Type=="sn" & fin_Neur$Annotation=="Neurons"),2], fin_Neur[which(fin_Neur$Type=="Xenium" & fin_Neur$Annotation=="Neurons"),2], paired = TRUE, alternative = "two.sided")



###Analyse OPC distribution

fin_OPC <- data.frame()
Diff_fin_OPC <- data.frame()
for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "OPC" & subset$Type == "sn"),2])
  sn_OPC <-   sum(subset[which(subset$Annotation == "OPC" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "OPC" & subset$Type == "Xenium"),2])
  Xenium_OPC <-   sum(subset[which(subset$Annotation == "OPC" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("OPC", "Other"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_OPC/(sn_OPC+sn_sum), sn_sum/(sn_OPC+sn_sum), Xenium_OPC/(Xenium_OPC+Xenium_sum), Xenium_sum/(Xenium_OPC+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_OPC <- rbind(fin_OPC, subset)
  Diff <- data.frame(sample, Xenium_OPC/(Xenium_OPC+Xenium_sum)-sn_OPC/(sn_OPC+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin_OPC <- rbind(Diff_fin_OPC, Diff)
}

p3 <- ggplot(fin_OPC, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("#F57F20", "grey"))+
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave(plot = p3, filename = "OPC_barplot.pdf")

Xenium_ss <- subset(fin_OPC,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_OPC,  Type == "sn", Count,
                drop = TRUE)

fin_OPC1 <- fin_OPC[fin_OPC$Annotation=="OPC",]

p4 <-ggplot(fin_OPC1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot=p4, filename = "OPC_change.pdf")

wilcox.test(fin_OPC[which(fin_OPC$Type=="sn" & fin_OPC$Annotation=="OPC"),2], fin_OPC[which(fin_OPC$Type=="Xenium" & fin_OPC$Annotation=="OPC"),2], paired = TRUE, alternative = "two.sided")


###Analyse Neuron distribution

fin_Olig <- data.frame()
Diff_fin_Olig <- data.frame()
for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "Oligodendrocytes" & subset$Type == "sn"),2])
  sn_Olig <-   sum(subset[which(subset$Annotation == "Oligodendrocytes" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "Oligodendrocytes" & subset$Type == "Xenium"),2])
  Xenium_Olig <-   sum(subset[which(subset$Annotation == "Oligodendrocytes" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("Oligodendrocytes", "Other"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_Olig/(sn_Olig+sn_sum), sn_sum/(sn_Olig+sn_sum), Xenium_Olig/(Xenium_Olig+Xenium_sum), Xenium_sum/(Xenium_Olig+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_Olig <- rbind(fin_Olig, subset)
  Diff <- data.frame(sample, Xenium_Olig/(Xenium_Olig+Xenium_sum)-sn_Olig/(sn_Olig+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin_Olig <- rbind(Diff_fin_Olig, Diff)
}

p5 <- ggplot(fin_Olig, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("#FDBF6F", "grey"))+
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave(plot=p5, filename = "Oligos_barplot.pdf")

Xenium_ss <- subset(fin_Olig,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_Olig,  Type == "sn", Count,
                drop = TRUE)

fin_Olig1 <- fin_Olig[fin_Olig$Annotation=="Oligodendrocytes",]

p6 <- ggplot(fin_Olig1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot = p6, filename = "Oligo_change.pdf")

wilcox.test(fin_Olig[which(fin_Olig$Type=="sn" & fin_Olig$Annotation=="Oligodendrocytes"),2], fin_Olig[which(fin_Olig$Type=="Xenium" & fin_Olig$Annotation=="Oligodendrocytes"),2], paired = TRUE, alternative = "two.sided")


###Analyse Astrocytes distribution

fin_Astro <- data.frame()
Diff_fin_Astro <- data.frame()
for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "Astrocytes" & subset$Type == "sn"),2])
  sn_Astro <-   sum(subset[which(subset$Annotation == "Astrocytes" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "Astrocytes" & subset$Type == "Xenium"),2])
  Xenium_Astro <-   sum(subset[which(subset$Annotation == "Astrocytes" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("Astrocytes", "Other"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_Astro/(sn_Astro+sn_sum), sn_sum/(sn_Astro+sn_sum), Xenium_Astro/(Xenium_Astro+Xenium_sum), Xenium_sum/(Xenium_Astro+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_Astro <- rbind(fin_Astro, subset)
  Diff <- data.frame(sample, Xenium_Astro/(Xenium_Astro+Xenium_sum)-sn_Astro/(sn_Astro+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin_Astro <- rbind(Diff_fin_Astro, Diff)
}

p7 <- ggplot(fin_Astro, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("#A6CEE2", "grey"))+ 
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave(plot = p7, filename = "Astro_barplot.pdf")

Xenium_ss <- subset(fin_Astro,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_Astro,  Type == "sn", Count,
                drop = TRUE)

fin_Astro1 <- fin_Astro[fin_Astro$Annotation=="Astrocytes",]

p8 <- ggplot(fin_Astro1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot = p8, filename = "Astro_change.pdf")

wilcox.test(fin_Astro[which(fin_Astro$Type=="sn" & fin_Astro$Annotation=="Astrocytes"),2], fin_Astro[which(fin_Astro$Type=="Xenium" & fin_Astro$Annotation=="Astrocytes"),2], paired = TRUE, alternative = "two.sided")


###Analyse Endotehlia distribution

fin_Endo <- data.frame()
Diff_fin_Endo <- data.frame()
for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "Endothelia" & subset$Type == "sn"),2])
  sn_Endo <-   sum(subset[which(subset$Annotation == "Endothelia" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "Endothelia" & subset$Type == "Xenium"),2])
  Xenium_Endo <-   sum(subset[which(subset$Annotation == "Endothelia" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("Endothelia", "Other"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_Endo/(sn_Endo+sn_sum), sn_sum/(sn_Endo+sn_sum), Xenium_Endo/(Xenium_Endo+Xenium_sum), Xenium_sum/(Xenium_Endo+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_Endo <- rbind(fin_Endo, subset)
  Diff <- data.frame(sample, Xenium_Endo/(Xenium_Endo+Xenium_sum)-sn_Endo/(sn_Endo+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin_Endo <- rbind(Diff_fin_Endo, Diff)
}

p9 <- ggplot(fin_Endo, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("#2179B4", "grey"))+
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave(plot = p9, filename = "Endo_barplot.pdf")

Xenium_ss <- subset(fin_Endo,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_Endo,  Type == "sn", Count,
                drop = TRUE)

fin_Endo1 <- fin_Endo[fin_Endo$Annotation=="Endothelia",]

p10 <- ggplot(fin_Endo1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot = p10, filename = "Endo_change.pdf")

wilcox.test(fin_Endo[which(fin_Endo$Type=="sn" & fin_Endo$Annotation=="Endothelia"),2], fin_Endo[which(fin_Endo$Type=="Xenium" & fin_Endo$Annotation=="Endothelia"),2], paired = TRUE, alternative = "two.sided")


###Analyse Lymphon distribution

fin_Lymph <- data.frame()
Diff_fin_Lymph <- data.frame()
for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "Lymphocytes" & subset$Type == "sn"),2])
  sn_Lymph <-   sum(subset[which(subset$Annotation == "Lymphocytes" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "Lymphocytes" & subset$Type == "Xenium"),2])
  Xenium_Lymph <-   sum(subset[which(subset$Annotation == "Lymphocytes" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("Lymphocytes", "Other"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_Lymph/(sn_Lymph+sn_sum), sn_sum/(sn_Lymph+sn_sum), Xenium_Lymph/(Xenium_Lymph+Xenium_sum), Xenium_sum/(Xenium_Lymph+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_Lymph <- rbind(fin_Lymph, subset)
  Diff <- data.frame(sample, Xenium_Lymph/(Xenium_Lymph+Xenium_sum)-sn_Lymph/(sn_Lymph+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin_Lymph <- rbind(Diff_fin_Lymph, Diff)
}

p11 <- ggplot(fin_Lymph, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("#EA8EAC", "grey"))+
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave(plot = p11, filename = "Lympho_barplot.pdf")

Xenium_ss <- subset(fin_Lymph,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_Lymph,  Type == "sn", Count,
                drop = TRUE)

fin_Lymph1 <- fin_Lymph[fin_Lymph$Annotation=="Lymphocytes",]

p12 <- ggplot(fin_Lymph1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot = p12, filename = "Lymph_change.pdf")

wilcox.test(fin_Lymph[which(fin_Lymph$Type=="sn" & fin_Lymph$Annotation=="Lymphocytes"),2], fin_Lymph[which(fin_Lymph$Type=="Xenium" & fin_Lymph$Annotation=="Lymphocytes"),2], paired = TRUE, alternative = "two.sided")


###Analyse Microglia distribution

fin_M <- data.frame()
Diff_fin_M <- data.frame()
for(i in 1:length(unique(fin$Alias))){
  sample <- unique(fin$Alias)[i]
  subset <- fin1[which(fin1$Alias == sample),]
  sn_sum <- sum(subset[which(subset$Annotation != "Macrophages/Microglia" & subset$Type == "sn"),2])
  sn_M <-   sum(subset[which(subset$Annotation == "Macrophages/Microglia" & subset$Type == "sn"),2])
  Xenium_sum <- sum(subset[which(subset$Annotation != "Macrophages/Microglia" & subset$Type == "Xenium"),2])
  Xenium_M <-   sum(subset[which(subset$Annotation == "Macrophages/Microglia" & subset$Type == "Xenium"),2])
  subset <- subset[1:4,]
  subset$Annotation <- rep(c("Macrophages/Microglia", "Other"), 2)
  subset$Type <- c("sn", "sn", "Xenium", "Xenium")
  subset$Count<- c(sn_M/(sn_M+sn_sum), sn_sum/(sn_M+sn_sum), Xenium_M/(Xenium_M+Xenium_sum), Xenium_sum/(Xenium_M+Xenium_sum))
  subset$Alias1 <- subset$Alias1 <- paste0(subset$Alias,"_", subset$Type)
  fin_M <- rbind(fin_M, subset)
  Diff <- data.frame(sample, Xenium_M/(Xenium_M+Xenium_sum)-sn_M/(sn_M+sn_sum))
  colnames(Diff) <- c("Sample", "Difference")
  Diff_fin_M <- rbind(Diff_fin_M, Diff)
}

p13 <- ggplot(fin_M, aes(fill=Annotation, y=Count, x=Alias1)) +scale_fill_manual(values=c("#E21F26", "grey"))+ 
  geom_bar(position="fill", stat="identity")+
  guides(x =  guide_axis(angle = 90))+
  theme_classic()
ggsave
ggsave(plot = p13, filename = "Microglia_barplot.pdf")

Xenium_ss <- subset(fin_M,  Type == "Xenium", Count,
                    drop = TRUE)
# subset weight data after treatment
sn_ss <- subset(fin_M,  Type == "sn", Count,
                drop = TRUE)

fin_M1 <- fin_M[fin_M$Annotation=="Macrophages/Microglia",]

p14 <- ggplot(fin_M1, aes(x = Type, y = Count)) + 
  geom_boxplot(aes(fill = Type), alpha = .2)+
  geom_line(aes(group = Alias)) + 
  geom_point(size = 2)+
  theme_classic()
ggsave(plot = p14, filename = "Microglia_change.pdf")

wilcox.test(fin_M[which(fin_M$Type=="sn" & fin_M$Annotation=="Macrophages/Microglia"),2], fin_M[which(fin_M$Type=="Xenium" & fin_M$Annotation=="Macrophages/Microglia"),2], paired = TRUE, alternative = "two.sided")




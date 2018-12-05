################################################################################
#stool_analysis.R
################################################################################
#Imports sequence data and creates:
#Alpha and beta diversity graphs, top 5 relative abundant ASVs over time.
#Finds chamber's coverage of stool community via raw counts and rarified species 
#By Lucas Flett

################################################################################
#Import Libraries
################################################################################
library('ggplot2') #used for graphics 
library('ggsignif') #p values on box plot
library('tidyr') #used to tidy data
library('dplyr') #used to manipulate data
library('vegan') #used for statistics
library('readr') #open seq .txt file from contaminate controls

################################################################################
#Import data
################################################################################
dat <- read.csv("SPE_SV2.csv", header = T, row.names = 1)
control1 <- read_file('control2.txt') #Full sequence from contaminated control
#control1 <- readChar('control1.txt', file.info('control1.txt')$size)
taxa = read.csv("taxa_Lucasrun272_gg2013.csv", header = T)

################################################################################
#Filter and clean data
################################################################################
control1 <- toupper(control1)
sample <- c("JCS159", "JCS160", "JCS161", "JCS162", "JCS163", "JCS164", 
            "JCS165", "JCS166","JCS167", "JCS168", "JCS169", "JCS170", "JCS187")
dat <- dat[sample,]
row.names(dat) = c("S1C1", "S1C2", "S1C3", "S1C6", "S1C7", "S1C8",
                   "S2C1", "S2C2", "S2C3", "S2C6", "S2C7", "S2C8", "Stool")
dat_t <- as.data.frame(t(dat), drop = F)
taxa$X <- as.character(taxa$X)

#Getting relative abundances from counts
rel_abun <- function(x)
{
  total <- sum(x)
  rel <- x/sum(x)
  return(rel)
}
dat_rel <- apply(dat, 1, rel_abun)
dat_rel_t <- as.data.frame(dat_rel)
dat_rel <- as.data.frame(t(dat_rel_t), drop = F)


################################################################################
#Data analysis 
################################################################################

#Beta diversity box plot between chambers and stool
bray_mat <- vegdist(dat_rel, method = "bray")
mat <- as.matrix(bray_mat)

mat_4.5 <- mat[1:6,1:6]
bray_4.5 <- mat_4.5[lower.tri(mat_4.5, diag = F)]

mat_28.5 <- mat[7:12,7:12]
bray_28.5 <- mat_28.5[lower.tri(mat_28.5, diag = F)]

bray_4.5stool <- 1:15
bray_4.5stool <- NA*bray_4.5stool
bray_4.5stool[1:6] <- mat[13,1:6]

bray_28.5stool <- 1:15
bray_28.5stool <- NA*bray_28.5stool
bray_28.5stool[1:6] <- mat[13,7:12] 

bray_bind <- as.data.frame(cbind(bray_4.5, bray_28.5, bray_4.5stool, 
                                 bray_28.5stool))
bray_gather <- gather(bray_bind, sample, bray)
bray_gather$sample <- factor(bray_gather$sample,
                            levels = c('bray_4.5', 'bray_4.5stool', 'bray_28.5', 
                                        'bray_28.5stool'), ordered = T)

box_plot <- ggplot(bray_gather, aes(x = sample, y = bray)) + geom_boxplot() + 
  labs(y = "Bray-Curtis", x = "") + scale_x_discrete(labels = 
  c("Chambers at 4.5hrs", "Chambers at 4.5hrs vs stool", "Chambers at 28.5hrs", 
    "Chambers at 28.5hrs vs stool")) + theme(axis.text.x = 
                                               element_text(face = "bold"))
#add this for sig diff.
#geom_signif(comparisons = list(c("bray_4.5stool", "bray_28.5stool")),
#                                                   map_signif_level = T) 


#Alpha diversity within chambers comparing 4.5hrs to 28.5hrs.
D <- diversity(dat)
t <- c(rep(4.5,6), rep(28.5,6))
chamber <- c("C1", "C2", "C3", "C6", "C7", "C8")
t_D <- cbind(t, chamber, D[1:12])
colnames(t_D)[3] <- "shannon"
t_D <- as.data.frame(t_D)
t_D$t <- as.numeric(t_D$t)
t_D$t[1:6] <- 4.5
t_D$t[7:12] <- 28.5
t_D$chamber <- as.character(t_D$chamber)
t_D$shannon <- as.numeric(as.character(t_D$shannon))
rownames(t_D) <- c()
shan_plot <- ggplot(t_D, aes(x = t, y = shannon, colour = chamber)) + 
  geom_point() + geom_line() + ylim(0.325,2.75) + labs(y = "Shannon Index", 
                                                       x = "Time (hours)") + 
  geom_hline(yintercept = 2.6165) + annotate("text", x = 10, y = 2.7, 
                                             label = "stool = 2.62")
  


#Top 5 most abundant microbes:
goo <- as.data.frame(apply(dat[1:12,], 2, sum))
colnames(goo) <- "count"
goo <- as.data.frame(goo[order(-goo$count), , drop = F])
goo <- goo[1:10, , drop = F]

E.Shi_2 <- dat_rel$s_Escherichia.Shigella_gg_Escherichia_2
Entero_6 <- dat_rel$s_Enterococcus_gg_Enterococcus_6
Pepto_29 <- dat_rel$s_Peptostreptococcaceae_gg_SMB53_29
Kleb_33 <- dat_rel$s_Klebsiella_gg_Klebsiella_33
Entero_18 <- dat_rel$s_Escherichia.Shigella_gg_Enterobacteriaceae_18
Bifido_3 <- dat_rel$s_Bifidobacterium_gg_Bifidobacterium_3 #10th most abundant

################################################################################
#Plotting most abundant strains over time
################################################################################

#s_Escherichia.Shigella_gg_Escherichia_
t_E.Shi <- cbind(t, chamber, E.Shi_2[1:12])
colnames(t_E.Shi)[3] <- "Relative"
t_E.Shi <- as.data.frame(t_E.Shi)
t_E.Shi$t <- as.numeric(t_E.Shi$t)
t_E.Shi$t[1:6] <- 4.5
t_E.Shi$t[7:12] <- 28.5
t_E.Shi$chamber <- as.character(t_E.Shi$chamber)
t_E.Shi$Relative <- as.numeric(as.character(t_E.Shi$Relative))
rownames(t_E.Shi) <- c()
E.Shi_plot <- ggplot(t_E.Shi, aes(x = t, y = Relative, colour = chamber)) + 
  geom_point() + geom_line() + geom_hline(yintercept = 0.0636) + 
  annotate("text", x = 11, y = 0.10, label = "stool = 6.36e-2") + 
  labs(y = "Relative Abundance", x = "Time (hours)",title = 
         "s_Escherichia.Shigella_gg_Escherichia_")

#s_Enterococcus_gg_Enterococcus_6
t_Entero_6 <- cbind(t, chamber, Entero_6[1:12])
colnames(t_Entero_6)[3] <- "Relative"
t_Entero_6 <- as.data.frame(t_Entero_6)
t_Entero_6$t <- as.numeric(t_Entero_6$t)
t_Entero_6$t[1:6] <- 4.5
t_Entero_6$t[7:12] <- 28.5
t_Entero_6$chamber <- as.character(t_Entero_6$chamber)
t_Entero_6$Relative <- as.numeric(as.character(t_Entero_6$Relative))
rownames(t_Entero_6) <- c()
Entero_6_plot <- ggplot(t_Entero_6, aes(x = t, y = Relative, colour =chamber)) + 
  geom_point() + geom_line() + geom_hline(yintercept = 0.000815) + 
  annotate("text", x = 12, y = 0.025, label = "stool = 8.15e-4") + 
  labs(y = "Relative Abundance", x = "Time (hours)",title = 
         "s_Enterococcus_gg_Enterococcus_6")

#s_Peptostreptococcaceae_gg_SMB53_29
t_Pepto_29 <- cbind(t, chamber, Pepto_29[1:12])
colnames(t_Pepto_29)[3] <- "Relative"
t_Pepto_29 <- as.data.frame(t_Pepto_29)
t_Pepto_29$t <- as.numeric(t_Pepto_29$t)
t_Pepto_29$t[1:6] <- 4.5
t_Pepto_29$t[7:12] <- 28.5
t_Pepto_29$chamber <- as.character(t_Pepto_29$chamber)
t_Pepto_29$Relative <- as.numeric(as.character(t_Pepto_29$Relative))
rownames(t_Pepto_29) <- c()
Pepto_29_plot <- ggplot(t_Pepto_29, aes(x = t, y = Relative, colour =chamber)) + 
  geom_point() + geom_line() + geom_hline(yintercept = 0.0125) + 
  annotate("text", x = 12, y = 0.02, label = "stool = 1.25e-2") + 
  labs(y = "Relative Abundance", x = "Time (hours)",title = 
         "s_Peptostreptococcaceae_gg_SMB53_29")

#s_Klebsiella_gg_Klebsiella_33
t_Kleb_33 <- cbind(t, chamber, Kleb_33[1:12])
colnames(t_Kleb_33)[3] <- "Relative"
t_Kleb_33 <- as.data.frame(t_Kleb_33)
t_Kleb_33$t <- as.numeric(t_Kleb_33$t)
t_Kleb_33$t[1:6] <- 4.5
t_Kleb_33$t[7:12] <- 28.5
t_Kleb_33$chamber <- as.character(t_Kleb_33$chamber)
t_Kleb_33$Relative <- as.numeric(as.character(t_Kleb_33$Relative))
rownames(t_Kleb_33) <- c()
Kleb_33_plot <- ggplot(t_Kleb_33, aes(x = t, y = Relative, colour = chamber)) + 
  geom_point() + geom_line() + geom_hline(yintercept = 0.00248) + 
  annotate("text", x = 12, y = 0.004, label = "stool = 2.48e-3") + 
  labs(y = "Relative Abundance", x = "Time (hours)",title = 
         "s_Klebsiella_gg_Klebsiella_33")

#s_Escherichia.Shigella_gg_Enterobacteriaceae_18
t_Entero_18 <- cbind(t, chamber, Entero_18[1:12])
colnames(t_Entero_18)[3] <- "Relative"
t_Entero_18 <- as.data.frame(t_Entero_18)
t_Entero_18$t <- as.numeric(t_Entero_18$t)
t_Entero_18$t[1:6] <- 4.5
t_Entero_18$t[7:12] <- 28.5
t_Entero_18$chamber <- as.character(t_Entero_18$chamber)
t_Entero_18$Relative <- as.numeric(as.character(t_Entero_18$Relative))
rownames(t_Entero_18) <- c()
Entero_18_plot <- ggplot(t_Entero_18, aes(x = t, y = Relative,colour=chamber)) + 
  geom_point() + geom_line() + geom_hline(yintercept = 0.000735) + 
  annotate("text", x = 12, y = 0.002, label = "stool = 7.35e-4") + 
  labs(y = "Relative Abundance", x = "Time (hours)",title = 
         "s_Escherichia.Shigella_gg_Enterobacteriaceae_18")

#s_Bifidobacterium_gg_Bifidobacterium_3
t_Bifido_3 <- cbind(t, chamber, Bifido_3[1:12])
colnames(t_Bifido_3)[3] <- "Relative"
t_Bifido_3 <- as.data.frame(t_Bifido_3)
t_Bifido_3$t <- as.numeric(t_Bifido_3$t)
t_Bifido_3$t[1:6] <- 4.5
t_Bifido_3$t[7:12] <- 28.5
t_Bifido_3$chamber <- as.character(t_Bifido_3$chamber)
t_Bifido_3$Relative <- as.numeric(as.character(t_Bifido_3$Relative))
rownames(t_Bifido_3) <- c()
Bifido_3_plot <- ggplot(t_Bifido_3, aes(x = t, y = Relative, colour =chamber)) + 
  geom_line() + geom_hline(yintercept = 0.001) + 
  annotate("text", x = 10, y = 0.0015, label = "stool = 0.01%") + 
  labs(y = "Relative Abundance", x = "Time (hours)",title = 
         "s_Bifidobacterium_gg_Bifidobacterium_3")


#Facet all species plots
abun1 <- cbind(t_E.Shi, rep("Escherichia_2",12))
colnames(abun1) <- c("Time", "Chamber", "Relative", "Species")

abun2 <- cbind(t_Entero_6, rep("Enterococcus_6",12))
colnames(abun2) <- c("Time", "Chamber", "Relative", "Species")

abun3 <- cbind(t_Pepto_29, rep("SMB53_29",12))
colnames(abun3) <- c("Time", "Chamber", "Relative", "Species")

abun4 <- cbind(t_Kleb_33, rep("Klebsiella_33",12))
colnames(abun4) <- c("Time", "Chamber", "Relative", "Species")

abun5 <- cbind(t_Entero_18, rep("Enterobacteriaceae_18",12))
colnames(abun5) <- c("Time", "Chamber", "Relative", "Species")

abun <- rbind(abun1, abun2, abun3, abun4, abun5)
abun$Species <- factor(abun$Species, levels =c("Escherichia_2","Enterococcus_6",
                            "SMB53_29","Klebsiella_33","Enterobacteriaceae_18"))

boo <- ggplot(abun, aes(x = Time, y = Relative, colour = Chamber)) + 
  geom_point() + geom_line() + facet_wrap(~Species, scales = "free") + 
  labs(x = "Time (hours)", y = "Relative Abundances")



#How many raw counts are in each sample?
stool_num <- length(which(dat_t$Stool > 0))
S1C1_num <- length(which(dat_t$S1C1 > 0))
S1C2_num <- length(which(dat_t$S1C2 > 0))
S1C3_num <- length(which(dat_t$S1C3 > 0))
S1C6_num <- length(which(dat_t$S1C6 > 0))
S1C7_num <- length(which(dat_t$S1C7 > 0))
S1C8_num <- length(which(dat_t$S1C8 > 0))
S2C1_num <- length(which(dat_t$S2C1 > 0))
S2C2_num <- length(which(dat_t$S2C2 > 0))
S2C3_num <- length(which(dat_t$S2C3 > 0))
S2C6_num <- length(which(dat_t$S2C6 > 0))
S2C7_num <- length(which(dat_t$S2C7 > 0))
S2C8_num <- length(which(dat_t$S2C8 > 0))

S1 <- rbind(S1C1_num, S1C2_num, S1C3_num, S1C6_num, S1C7_num, S1C8_num)
S2 <- rbind(S2C1_num, S2C2_num, S2C3_num, S2C6_num, S2C7_num, S2C8_num)

#How many ASVs are recovered in each sample?
both <- 0*(1:12)
for (i in 1:12)
{
  both[i] <- length(which(dat_t$Stool > 0 & dat_t[,i] > 0))
}
strain_count <- as.data.frame(cbind(S1,both[1:6],S2,both[7:12]))
colnames(strain_count) <- c("hr4.5", "hr4.5_stool", "hr28.5", "hr28.5_stool")
rownames(strain_count) <- c()
percent <- as.data.frame(cbind(strain_count$hr4.5_stool/62, 
                               strain_count$hr28.5_stool/62))
avg_4.5 <- mean(strain_count$hr4.5)
avg_28.5 <- mean(strain_count$hr28.5)
avg_4.5percent <- mean(percent[,1])
avg_28.5percent <- mean(percent[,2])

#The above species nunbers aren't very accurate, should rarify instead:
N <- rowSums(dat) 
S.rar <-rarefy(dat, min(N))
strain_count_rar <- cbind(S.rar[1:6],S.rar[1:6]/46.98, S.rar[7:12],
                          S.rar[7:12]/46.98)
avg_4.5.r <- mean(strain_count_rar[,1])
avg_28.5.r <- mean(strain_count_rar[,3])
avg_4.5percent.r <- mean(strain_count_rar[,2])
avg_28.5percent.r <- mean(strain_count_rar[,4])
avg_rar <- cbind(avg_4.5.r,avg_28.5.r,avg_4.5percent.r,avg_28.5percent.r)

#Find contaminated control sequences in my other chambers
needle <- 0
for (i in 1:length(taxa$X)){
  needle[i] <- grepl(taxa$X[i], control1, fixed = T)
}


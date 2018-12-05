################################################################################
#ODvsTime.R
################################################################################
#Imports log files from turbidostat and creates:
#OD and dilution vs time graphs 
#Quantifies the reproducibility between chambers
#By Lucas Flett

################################################################################
#Import Libraries needed
################################################################################
library('ggplot2') #used for graphics 
library('tidyr') #used to tidy data
library('dplyr') #used to manipulate data
library('zoo') #used for rolling averages, smoothing curves

################################################################################
#Import data
################################################################################
dilution = read.csv("dilution.csv", header = F)
od = read.csv("od.csv", header = F)
t = read.csv("time.csv", header = F)

################################################################################
#Cleaning data
################################################################################
names(t)[names(t) == 'V1'] <- 'time'
od <- rename(od,chamber_1=V1,chamber_2=V2,chamber_3=V3,chamber_4=V4,
             chamber_5=V5,chamber_6=V6,chamber_7=V7,chamber_8=V8)
dilution <- rename(dilution,chamber_1=V1,chamber_2=V2,chamber_3=V3,chamber_4=V4,
                   chamber_5=V5,chamber_6=V6,chamber_7=V7,chamber_8=V8)

t_od <- cbind(t,od)
t_dilution <- cbind(t,dilution)

#Rolling averages
#Include this code if you want smooth curves
#
#t_od <- rollmean(t_od,100)
#t_od <- data.frame(t_od)
#t_dilution <- rollmean(t_dilution,100)
#t_dilution <- data.frame(t_dilution)

#Transform data into acceptable format for ggplot
t_od_gather <- gather(t_od, chamber, od, -time)
t_dilution_gather <- gather(t_dilution, chamber, dilution, -time)


################################################################################
#Ploting all chambers 
################################################################################
all_ods <- ggplot(filter(t_od_gather, time < 70), aes(x = time, y = od, 
              colour = chamber)) + geom_line() + labs(x = "time (hours)")
all_dilutions <- ggplot(t_dilution_gather, aes(x = time, y = dilution, 
                    colour = chamber)) + geom_line() + labs(x = "time (hours)")


################################################################################
#Quantifying Reproducibility between chambers: 
#Pairwise Euclidean distance between OD growth curves
################################################################################
repo <- t_od[,-c(1,5,6)] #for bacteroides -c(1,5,6), for E.coli -c(1,8,9)
colnames(repo) <- c("chamber_1","chamber_2","chamber_3","chamber_4","chamber_5",
                    "chamber_6")
#Compare up to 0.6 OD as threshold OD's where different between experiments
repo <- filter(repo, chamber_1 < 0.6 & chamber_2 < 0.6 & chamber_3 < 0.6, 
               chamber_4 < 0.6 & chamber_5 < 0.6 & chamber_6 < 0.6) 

y12 <- sqrt(sum((repo$chamber_1 - repo$chamber_2)^2))
y13 <- sqrt(sum((repo$chamber_1 - repo$chamber_3)^2))
y14 <- sqrt(sum((repo$chamber_1 - repo$chamber_4)^2))
y15 <- sqrt(sum((repo$chamber_1 - repo$chamber_5)^2))
y16 <- sqrt(sum((repo$chamber_1 - repo$chamber_6)^2))
y23 <- sqrt(sum((repo$chamber_2 - repo$chamber_3)^2))
y24 <- sqrt(sum((repo$chamber_2 - repo$chamber_4)^2))
y25 <- sqrt(sum((repo$chamber_2 - repo$chamber_5)^2))
y26 <- sqrt(sum((repo$chamber_2 - repo$chamber_6)^2))
y34 <- sqrt(sum((repo$chamber_3 - repo$chamber_4)^2))
y35 <- sqrt(sum((repo$chamber_3 - repo$chamber_5)^2))
y36 <- sqrt(sum((repo$chamber_3 - repo$chamber_6)^2))
y45 <- sqrt(sum((repo$chamber_4 - repo$chamber_5)^2))
y46 <- sqrt(sum((repo$chamber_4 - repo$chamber_6)^2))
y56 <- sqrt(sum((repo$chamber_5 - repo$chamber_6)^2))

repo_num <- sum(y12,y13,y14,y15,y16,y23,y24,y25,y26,y34,y35,y36,y45,y46,y56)
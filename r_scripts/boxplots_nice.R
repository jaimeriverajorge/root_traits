# This R-script creates nice box plots of different traits compared by functional
# groups: specific root length (SRL) and root tissue density (RTD)

#CLEAN UP
rm(list=ls(all=TRUE))

# read in root traits file
file_path <- "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\data\\root_traits_summarized_noN.csv"
roots.df <- read.csv(file_path)

# since we are comparing by functional group, we want to get rid of the shrubs
roots.df <- subset(roots.df, func_group!='SS' & func_group != 'SSN')

#vector of SRL trait optima
srl_optima <- c(3631.481, 8307.738, 1527.34)

# create boxplot
SRL<- ggplot(roots.df, aes(x = func_group, y = SRL, fill=func_group)) +
  theme_linedraw() +
  geom_boxplot(notch=F,
               lwd=.5,
               stat="boxplot",
               outlier.shape=NA)+
  geom_text(aes(x=1, y=3631.481), label="*", size=5)+
  geom_text(aes(x=2, y=8307.738), label="*", size=5)+
  geom_text(aes(x=3, y=1527.34), label="*", size=5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=.7),
        axis.title.x = element_text(margin = margin(t = 10, b=5), size=16),
        axis.title.y = element_text(margin = margin(l = 5, r=5), size=16),
        axis.text.x= element_text(margin = margin(t = 10), size=12),
        axis.text.y=element_text(margin = margin(r = 10), size=12),
        legend.title = element_text(colour="black", size=16),
        legend.text = element_text(size=14))+
  labs(y = expression(atop("Specific Root Length", paste(~(cm/g)))),
       x = 'Functional group')+
  scale_x_discrete (limits = c("F", "G", "N"),
                    labels = c("Forbs", "Graminoids", "N-fixers"))+
  scale_fill_manual(name= "functional_groups",
                    values=c("#E69F00", "#009E73", "#56B4E9"),
                    limits = c("F", "G", "N"),
                    labels = c("F", "G", "N"))+
  theme(legend.position="none")
SRL
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\boxplot_SRL_func_group.png"
png(file=file_name,height=3.5,width=4,units="in",res=1200)
plot(SRL)
dev.off()

#vector of RTD trait optima
rtd_optima <- c(0.3027799, 0.3454367, 0.8135943)
rtd_optima[1]
# create boxplot
RTD<- ggplot(roots.df, aes(x = func_group, y = RTD, fill=func_group)) +
  theme_linedraw() +
  geom_boxplot(notch=F,
               lwd=.5,
               stat="boxplot",
               outlier.shape=NA)+
  geom_text(aes(x=1, y=rtd_optima[1]), label="*", size=5)+
  geom_text(aes(x=2, y=rtd_optima[2]), label="*", size=5)+
  geom_text(aes(x=3, y=rtd_optima[3]), label="*", size=5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black", size=.7),
        axis.title.x = element_text(margin = margin(t = 10, b=5), size=16),
        axis.title.y = element_text(margin = margin(l = 5, r=5), size=16),
        axis.text.x= element_text(margin = margin(t = 10), size=12),
        axis.text.y=element_text(margin = margin(r = 10), size=12),
        legend.title = element_text(colour="black", size=16),
        legend.text = element_text(size=14))+
  labs(y = expression(atop("Root Tissue Density", paste(~(g/cm^{3})))),
       x = 'Functional group')+
  scale_x_discrete (limits = c("F", "G", "N"),
                    labels = c("Forbs", "Graminoids", "N-fixers"))+
  scale_fill_manual(name= "functional_groups",
                    values=c("#E69F00", "#009E73", "#56B4E9"),
                    limits = c("F", "G", "N"),
                    labels = c("F", "G", "N"))+
  theme(legend.position="none")
RTD
file_name = "C:\\Users\\jjaime-rivera\\Dropbox\\PC\\Documents\\figures\\boxplot_RTD_func_group.png"
png(file=file_name,height=3.8,width=3.8,units="in",res=1200)
plot(RTD)
dev.off()









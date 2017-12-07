# Figures for ISEEC 2017 conference


# trinucleotide ---- 


load("output_data/trinuc_data_LUAD.RData")
library(ggplot2)
p <- ggplot(data=trinuc.mutation_data, aes(Downstream, Upstream)) +
  geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent")
p <- p + facet_grid(.~section_labels, labeller = label_parsed) 
p <- p +  geom_text(aes(label = round(proportion, 4)*100),size=3)
# p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.ticks = element_blank(),
                            strip.text=element_text(size=15),
                            axis.title.x = element_text(size=15),
                            axis.title.y = element_text(size=15),
                            axis.text.x = element_text(size=12),
                            axis.text.y=element_text(size=12),plot.title = element_text(hjust = 0.5)) 
# ggtitle(paste("Trinucleotide profile for ",tumor.name,sep=""))
p
ggsave(paste("figures/",tumor.name,"_trinuc_heatmap_noname_ISEEC.pdf",sep=""),height = 2.5,width = 10)

# nuc mutations ----

p <- ggplot(normalized.mut.matrix.df.m[normalized.mut.matrix.df.m$Pos %in% positions.to.plot,],aes(Pos, value, fill=variable)) 
p <- p + geom_bar(stat="identity",position="dodge",color="black") 
# p <- p + theme_bw() 
p <- p + theme(panel.background = element_blank())
p <- p + theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey"))
p <- p + theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p <- p + scale_y_continuous(labels=fancy_scientific)
p <- p + geom_text(aes(label=round(value*1e6,2)), vjust=-0.2, position=position_dodge(width=0.9),size=5)
p <- p + geom_text(aes(label=variable,y=-1e-7),position=position_dodge(width = 0.9),size=5)
p <- p + xlab("Nucleotide position")
p <- p + ylab("Mutation rate")
p <- p + scale_fill_discrete(name="Mutation")
p <- p + theme(axis.text.x= element_text(size=30))
p <- p + theme(axis.text.y= element_text(size=25))
p <- p + theme(axis.title.y = element_text(size = 25 , angle = 90))
p <- p + theme(axis.title.x = element_text(size = 25, angle = 00))
p <- p + theme(legend.position = c(0.92, .85),legend.text=element_text(size=25))
# p <- p + scale_x_discrete(labels=paste(Pos))
p
ggsave(p, file="figures/nuc_mut_rates_NCIupdated.png",units = "in",height=7,width = 10)


# amino acid rate ----

just.12 <- KRAS.mutation.output$amino_acid_mutation_rates[,12]
ordered.12 <- sort(just.12,decreasing = T)[1:8] #Order it high --> low
ordered.12 <- ordered.12[c("Phe","Ser","Tyr","Arg","Trp","STOP","Cys","Gly")] #Change the ordering for the figure
message("Total mutation rate for oncogenic KRAS G12C mutations: ");print(sum(ordered.12[c("Phe","Tyr","Ser","Arg","Trp")])) #Total mutation rate 
message("Total mutation rate for all KRAS G12C mutations: ");print(sum(ordered.12)) #Total mutation rate 


KRAS.mutation.output$all_mutations[which(KRAS.mutation.output$all_mutations$AA_Pos==61 & KRAS.mutation.output$all_mutations$freq>1),]
message("Total rate of resistant mutations within KRAS") 
sum(ordered.12[1:5],KRAS.mutation.output$all_mutations$mu[which(KRAS.mutation.output$all_mutations$AA_Pos==61 &KRAS.mutation.output$all_mutations$freq>1)])

message("Sum of position 61:")
sum(KRAS.mutation.output$all_mutations$mu[which(KRAS.mutation.output$all_mutations$AA_Pos==61 &KRAS.mutation.output$all_mutations$freq>1)])

# ordered.12

ordered.12.df <- as.data.frame(matrix(NA,nrow=length(ordered.12),ncol=2))
colnames(ordered.12.df) <- c("Mutation","Rate")
ordered.12.df$Mutation <- names(ordered.12)
ordered.12.df$Rate <- as.numeric(ordered.12)
ordered.12.df
ordered.12.df$Mutation <- factor(ordered.12.df$Mutation, levels = c("Phe","Ser","Tyr","Arg","Trp","STOP","Cys","Gly"))#ordered.12.df$Mutation) #this makes Mutation an ordered factor, which ggplot uses. (http://stackoverflow.com/questions/20041136/avoid-ggplot-sorting-the-x-axis-while-plotting-geom-bar) 
ordered.12.df$alpha_vals <- c(rep(1,5),rep(1,3))
ordered.12.df$color_vals <- c(rep("black",5),rep("white",3))
# scale_fill_manual(values=cbPalette)


fig3 <- ggplot(ordered.12.df,aes(x=Mutation, y=Rate, fill=Mutation))
fig3 <- fig3 + geom_bar(stat="identity",position="dodge",color="black")#,aes(alpha=factor(alpha_vals)))
fig3 <- fig3 + theme(panel.background = element_blank())
fig3 <- fig3 + scale_fill_manual(values=ordered.12.df$color_vals)
fig3 <- fig3 + theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey"))
fig3 <- fig3 + theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
fig3 <- fig3 + scale_y_continuous(labels=fancy_scientific)
fig3 <- fig3 + geom_text(aes(label=round(Rate*1e6,2)), vjust=-0.2, position=position_dodge(width=0.9),size=5)
fig3 <- fig3 + xlab("Amino acid mutation")
fig3 <- fig3 + ylab("Mutation rate")
fig3 <- fig3 + theme(axis.text.x= element_text(size=30))
fig3 <- fig3 + theme(axis.text.y= element_text(size=25))
fig3 <- fig3 + theme(axis.title.y = element_text(size = 25, angle = 90))
fig3 <- fig3 + theme(axis.title.x = element_text(size = 25, angle = 00))
fig3 <- fig3 + theme(legend.position = 'none')
fig3


ggsave(fig3, file="figures/AA_mutation_rates_NCIupdateddata.png",units = "in",height=7,width = 10)


###
#Heatmaps to show silent mutation rate and prevalence ---- 
###

# nucleotide level ---- 

load("output_data/KRAS_NO_mutation_output.RData")

head(KRAS.no.mut.output$norm_mut[,1:5])




# amino acid level ---- 

load("output_data/KRAS_mutation_output.RData")
load("output_data/KRAS_NO_mutation_output.RData")
#mutation rate 
rownames(KRAS.mutation.output$amino_acid_mutation_rates)


KRAS.nuc.mut.df <- as.data.frame(matrix(data = NA,nrow=4,ncol=ncol(KRAS.no.mut.output$norm_mut)+1))

KRAS.nuc.mut.df[,1] <- rownames(KRAS.no.mut.output$norm_mut)
KRAS.nuc.mut.df[,2:ncol(KRAS.nuc.mut.df)] <- KRAS.no.mut.output$norm_mut
colnames(KRAS.nuc.mut.df)[1] <- "Mutation"
colnames(KRAS.nuc.mut.df)[2:ncol(KRAS.nuc.mut.df)] <- colnames(KRAS.no.mut.output$norm_mut)

library(reshape2)
source("R/fancy_scientific_code.R")


KRAS.nuc.mut.df.melt <- melt(KRAS.nuc.mut.df)
# order.we.want <- unique(KRAS.nuc.mut.df.melt$Mutation)[order(unique(KRAS.nuc.mut.df.melt$Mutation))]
# order.we.want <- order.we.want[c(1:16,18,19,20,21,17)]
# KRAS.nuc.mut.df.melt$Mutation <- factor(KRAS.nuc.mut.df.melt$Mutation, levels = rev(order.we.want))


KRAS.mut.tiles.nuc <- ggplot(data = KRAS.nuc.mut.df.melt[1:(4*18),], aes(x=variable, y= Mutation)) 
KRAS.mut.tiles.nuc <- KRAS.mut.tiles.nuc + geom_tile(aes(fill=value)) +  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
KRAS.mut.tiles.nuc <- KRAS.mut.tiles.nuc +  scale_fill_gradient(low = "gray100", high = "red",labels=fancy_scientific,name="Mutation rate\n") + labs(x="Nucleotide position") + theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold")) #+ ggtitle("Silent mutation rate")
KRAS.mut.tiles.nuc <- KRAS.mut.tiles.nuc + theme(axis.text.x= element_text(size=30))
KRAS.mut.tiles.nuc <- KRAS.mut.tiles.nuc + theme(axis.text.y= element_text(size=25))
KRAS.mut.tiles.nuc <- KRAS.mut.tiles.nuc + theme(axis.title.y = element_text(size = 25, angle = 90))
KRAS.mut.tiles.nuc <- KRAS.mut.tiles.nuc + theme(axis.title.x = element_text(size = 25, angle = 00)) +  theme(legend.text=element_text(size=20),legend.title = element_text(size=20))
KRAS.mut.tiles.nuc

ggsave(filename = "figures/_ISEEC2017_NUC_mutation_heatmap.png",plot = KRAS.mut.tiles.nuc)


# KRAS.aa.mut.df <- as.data.frame(KRAS.mutation.output$amino_acid_mutation_rates)
# rownames(KRAS.aa.mut.df)

KRAS.aa.mut.df <- as.data.frame(matrix(data=NA,nrow=21,ncol=ncol(KRAS.no.mut.output$amino_acid_mutation_rates)+1))
KRAS.aa.mut.df[,1] <- rownames(KRAS.no.mut.output$amino_acid_mutation_rates)
KRAS.aa.mut.df[,2:191] <- KRAS.no.mut.output$amino_acid_mutation_rates
colnames(KRAS.aa.mut.df)[1] <- "Mutation"
colnames(KRAS.aa.mut.df)[2:191] <- colnames(KRAS.no.mut.output$amino_acid_mutation_rates)


library(reshape2)
source("R/fancy_scientific_code.R")

KRAS.aa.mut.df.melt <- melt(KRAS.aa.mut.df)
order.we.want <- unique(KRAS.aa.mut.df.melt$Mutation)[order(unique(KRAS.aa.mut.df.melt$Mutation))]
order.we.want <- order.we.want[c(1:16,18,19,20,21,17)]
KRAS.aa.mut.df.melt$Mutation <- factor(KRAS.aa.mut.df.melt$Mutation, levels = rev(order.we.want))


KRAS.mut.tiles <- ggplot(data = KRAS.aa.mut.df.melt[1:(21*20),], aes(x=variable, y= Mutation)) 
KRAS.mut.tiles <- KRAS.mut.tiles + geom_tile(aes(fill=value)) +  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
KRAS.mut.tiles <- KRAS.mut.tiles +  scale_fill_gradient(low = "gray100", high = "red",labels=fancy_scientific,name="Mutation rate\n") + labs(x="Amino acid position") + theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))# + ggtitle("Silent mutation rate")
KRAS.mut.tiles <- KRAS.mut.tiles + theme(axis.text.x= element_text(size=22))
KRAS.mut.tiles <- KRAS.mut.tiles + theme(axis.text.y= element_text(size=17))
KRAS.mut.tiles <- KRAS.mut.tiles + theme(axis.title.y = element_text(size = 25, angle = 90))
KRAS.mut.tiles <- KRAS.mut.tiles + theme(axis.title.x = element_text(size = 25, angle = 00)) +  theme(legend.text=element_text(size=20),legend.title = element_text(size=20))
KRAS.mut.tiles

KRAS.mut.tiles

ggsave(filename = "figures/_ISEEC2017_mutation_heatmap.png",plot = KRAS.mut.tiles, width = 16,height = 9)




#prevalence 
head(KRAS.no.mut.output$amino_acid_tally)


KRAS.aa.tally.df <- as.data.frame(matrix(data=NA,nrow=21,ncol=ncol(KRAS.no.mut.output$amino_acid_tally)+1))
KRAS.aa.tally.df[,1] <- rownames(KRAS.no.mut.output$amino_acid_tally)
KRAS.aa.tally.df[,2:191] <- KRAS.no.mut.output$amino_acid_tally
colnames(KRAS.aa.tally.df)[1] <- "Mutation"
colnames(KRAS.aa.tally.df)[2:191] <- colnames(KRAS.no.mut.output$amino_acid_tally)


library(reshape2)

KRAS.aa.tally.df.melt <- melt(KRAS.aa.tally.df)
order.we.want <- unique(KRAS.aa.tally.df.melt$Mutation)[order(unique(KRAS.aa.tally.df.melt$Mutation))]
order.we.want <- order.we.want[c(1:16,18,19,20,21,17)]
KRAS.aa.tally.df.melt$Mutation <- factor(KRAS.aa.tally.df.melt$Mutation, levels = rev(order.we.want))


KRAS.tally.tiles <- ggplot(data = KRAS.aa.tally.df.melt[1:(21*20),], aes(x=variable, y= Mutation)) 
KRAS.tally.tiles <- KRAS.tally.tiles + geom_tile(aes(fill=value)) + geom_text(aes(x=variable, y= Mutation,label=value))
KRAS.tally.tiles <- KRAS.tally.tiles + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  scale_fill_gradient(low = "gray100", high = "red",name="Mutation number") + labs(x="Amino acid position") + theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold")) #+ ggtitle("Observed mutation prevalence")
KRAS.tally.tiles <- KRAS.tally.tiles + theme(axis.text.x= element_text(size=22))
KRAS.tally.tiles <- KRAS.tally.tiles + theme(axis.text.y= element_text(size=17))
KRAS.tally.tiles <- KRAS.tally.tiles + theme(axis.title.y = element_text(size = 25, angle = 90))
KRAS.tally.tiles <- KRAS.tally.tiles + theme(axis.title.x = element_text(size = 25, angle = 00)) +  theme(legend.text=element_text(size=20),legend.title = element_text(size=20))

KRAS.tally.tiles

ggsave(filename = "figures/_ISEEC2017_tally_heatmap.png",plot = KRAS.tally.tiles, width = 16,height = 9)

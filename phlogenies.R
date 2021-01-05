#libraries
library(ape)
library(seqinr)
library(phangorn)
library(tidyverse)
library(ggtext)

#read in aligned data and grab some sequence info- could probably use an API function call
#but tbh idk how and this works just fine
#setwd("~/Desktop/genomics/final_project")
covid<-read.aa("aligned.txt", format="fasta")
#update labels to just be the accession number for easier graphing later on
old=as.matrix(names(covid))
new<-read.delim("states.txt", header=FALSE)
names(new)<-c("accession","state","date")
covid=updateLabel(covid, old, as.character(new$accession))

#subset based off period 1 Jan1-July1 vs period 2 July1-Dec1
timeline<-new%>%
  na.omit()%>%
  mutate(date= zoo::as.yearmon(date, "%Y-%m"))%>%
  mutate(date=format(date, "%m"))%>%
  mutate(period=ifelse(as.numeric(date)<6,"pre June","postJune"))%>%
  mutate(region=case_when(
    state %in% c("CA","NM","WA") ~ "West",
    state %in% c("FL","GA","TX","SC") ~ "South",
    state %in% c("NY","MD","VA","NC") ~ "East",
    state %in% c("MN","WI","IL","IN") ~ "Midwest",
    TRUE~ "Other"
  ))%>%
  mutate(state=as.factor(state))%>%
  na.omit()

timeline%>%ggplot(aes(date,fill=region))+
  geom_bar()+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.5),
        text=element_text(size=7))+
  labs(
    title = "**COVID spike protein**  
    <span style='font-size:11pt'> Sequence counts by region
    <span style='color:#F8766D;'>East</span>, 
    <span style='color:#7CAE00;'>Midwest</span>, 
    <span style='color:#00BFC4;'>South</span>, and
    <span style='color:#C77CFF;'>West</span> split into time intervals 
    <span style='color:#B79F00;'>pre June</span>(n=75) and
    <span style='color:#F564E3;'>post June</span>(n=75)</span>",
    x = "Month in 2020", y = "sequence count"
  ) +
  geom_vline(xintercept = 6.5,color="darkgrey")+
  theme_minimal() +
  geom_text(x=6.3, y=34.5, label="pre June", col="#B79F00",size=4, angle=90)+
  geom_text(x=6.7, y=35, label="post June", col="#F564E3",size=4,angle=90)+
  theme(
    plot.title = element_markdown(lineheight = 1.1),
    legend.text = element_markdown(size = 11))
  #text(x=)
 #facet_wrap(~regionWest,ncol = 2)

#change alignement to phyDat object for use in the phanghon package.
covid_phyDat <- phyDat(covid, type = "AA", levels = NULL)

#run some model testing to see which distance matrix is best for our data.
#We can use mt list to pick three and compare trees later on.
mt <- modelTest(covid_phyDat, model="all")
dna_dist <- dist.ml(covid, model="mtmam")

covid_UPGMA <- upgma(dna_dist)
covid_NJ  <- NJ(dna_dist)
covid_fastme<-fastme.bal(dna_dist)
par(mfrow=c(1,3))
plot(covid_UPGMA,main="UPGMA",show.tip=FALSE)
title("UPGMA")
tiplabels(pch=15, cex=.5, 
          col = as.factor(timeline$period[match(covid_fastme$tip.label, timeline$accession)]))
#legend("right", legend=c("Midwest","East","South","West"),fill = c("1","2","3","4"))
legend("topleft", legend=c("Pre June", "Post June"),fill = c("1","2"))
plot(covid_NJ,main = "Neighbor Joining",show.tip=FALSE)
plot(covid_fastme,,main = "FastME", show.tip=FALSE)

covid_optim <- optim.parsimony(covid_NJ, covid_phyDat)
covid_pratchet <- pratchet(covid_phyDat)

fit <- pml(covid_NJ, covid10)
print(fit)
fitMTMAM <- optim.pml(fit, model = "mtmam", rearrangement = "stochastic")
logLik(fitMTMAM)
bs <- bootstrap.pml(fitMTMAM, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitMTMAM$tree), bs, p = 50, type="p")
write.tree(bs, file="bootstrap_MTMAM.tre")


#plotting states
library(usmap)
monthly<-timeline%>%group_by(state,date)%>%count()
periodly<-timeline%>%group_by(state,period)%>%count()

plot_usmap(data=monthly,values="n", 
           include = c("CA","NM","WA","FL","GA","TX",
                       "SC","NY","MD","VA","NC","MN",
                       "WI","IL","IN"),
           color="purple") +
  scale_fill_continuous(
    low = "white", high = "purple", name = "seq count", 
    label = scales::comma
  ) + theme(legend.position = "right")+
  labs(title="Sequence counts per month in 2020",
       subtitle = )+
  facet_wrap(~date)
  

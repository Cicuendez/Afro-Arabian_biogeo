# PLOT CUMULATIVE NUMBER OF EVENTS (EMPIRICAL AND SIMULATED)
setwd()


sim_event_all_cons <- readRDS("objects/sim_event_all_cons.rds")
event_all_cons <- readRDS("objects/event_all_cons.rds")

nsim <- 1000
myr <- 60

# Import empirical
event_list_cons <- readRDS("objects/event_list_cons.rds")
event_df <- data.frame(node=0, genus=NA, height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)

for (i in 1:length(event_list_cons)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      event_df <- rbind(event_df, event_list_cons[[i]][j,])
    }
  }
}

event_df <- event_df[-1,]
head(event_df)

#sort(event_df$min)
event_df_sorted <- event_df[rev(order(event_df$max)),]
xx <- event_df_sorted

# Import simulations
sim_event_cons <- readRDS("objects/sim_event_cons.rds")
sim_event_df <- vector("list", length(nsim))
for (i in 1:nsim){
  sim_event_df[[i]] <- data.frame(node=0, genus=NA, height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)
}

head(sim_event_df)

for (i in 1:length(sim_event_cons)){
  for (j in 1:nsim){
    if (nrow(sim_event_cons[[i]][[j]]) > 0){
      for (r in 1:nrow(sim_event_cons[[i]][[j]])){
        sim_event_df[[j]] <- rbind(sim_event_df[[j]], sim_event_cons[[i]][[j]][r,])
      }
    }
  }
}

for (i in 1:length(sim_event_df)){
  sim_event_df[[i]] <- sim_event_df[[i]][-1,]
}

# sort each dataframe
event_df_sorted <- event_df[rev(order(event_df$max)),]

sim_event_df_sorted <- vector("list", nsim)
for (i in 1:nsim){
  sim_event_df_sorted[[i]] <- sim_event_df[[i]][rev(order(sim_event_df[[i]]$max)),]
}

head(sim_event_df_sorted)
sim_xx <- sim_event_df_sorted


##

time_vec <- seq(from=60, to=1, by=-1)

rmat <- matrix(NA, ncol = 2, nrow=length(time_vec))
colnames(rmat) <- c("Ma", "Ncum")
rmat[, "Ma"] <- time_vec

sim_rmat <- vector("list", nsim)
for (i in 1:nsim){
  sim_rmat[[i]] <- matrix(NA, ncol=2, nrow=length(time_vec))
  colnames(sim_rmat[[i]]) <- c("Ma", "Ncum")
  sim_rmat[[i]][, "Ma"] <- time_vec
}


nrow(xx[xx$max >= time_vec[20], ])
for (ii in 1:length(time_vec)){
  
  tmp <- xx[xx$max >= time_vec[ii], ]
  rmat[ii, "Ncum"] <- nrow(tmp)
}

nrow(sim_xx[[1]][sim_xx[[1]]$max >= time_vec[20], ])
for (s in 1:nsim){
  for (ii in 1:length(time_vec)){
    tmp <- sim_xx[[s]][sim_xx[[s]]$max >= time_vec[ii], ]
    sim_rmat[[s]][ii, "Ncum"] <- nrow(tmp)
  }
}

head(sim_rmat)

saveRDS(sim_rmat, "objects/sim_rmat.rds")
saveRDS(rmat, "objects/rmat.rds")

## PLOT CUMULATIVE NUMBER OF EVENTS (OBSERVED AND SIM) ----
sim_rmat <- readRDS("objects/sim_rmat.rds")
rmat <- readRDS("objects/rmat.rds")

##### 95% CI CALCULATION #####
# quantiles and mean of each My.
prob_qup <- 0.975
prob_qlow <- 0.025
nsim <- 1000
myr <- 60

##### 95% CI ALL EVENTS #####
# Let's do a list where each element is a vector with the
# number of transitions per simulation for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist <- vector("list", myr)
for (i in 1:myr){
  dlist[[i]] <- vector("numeric", nsim)
}
#dlist[[1]][2]

for (i in 1:nsim){
  for (j in 1:myr){
    dlist[[j]][i] <- as.data.frame(sim_rmat[[i]])$Ncum[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df <- data.frame(Ma=c(myr:1), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df$mean[i] <- mean(dlist[[i]])
  Q_df$qlow[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df$qup[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line <- spline(Q_df$Ma-1, Q_df$mean)
qup_line <- spline(Q_df$Ma-1, Q_df$qup)
qlow_line <- spline(Q_df$Ma-1, Q_df$qlow)


# Set colors for observed and simulated lines ----
color_all <- "#EE6A50"
color_sim <- "gray93"


pdf("plots/cumulative_events_OK.pdf", paper="a4", height=20, width=10)
plot(1,type='n',xlim=c(myr,0),ylim=c(0,80),xlab='Ma', ylab='Cumulative N', main="Cumulative biogeographic events")
for (s in 1:nsim){
  lines(spline(sim_rmat[[s]][, "Ma"]-1, sim_rmat[[s]][, "Ncum"]), type="l", col=color_sim, lwd=0.25, pch=16, cex=0.7)
}
lines(spline(rmat[, "Ma"]-1, rmat[, "Ncum"]), type="l", col=color_all, lwd=3)
lines(mean_line, type="l", col="black", lwd=1, pch=16, cex=1)
lines(qlow_line, type="l", col="black", lwd=2, pch=16, cex=1)
lines(qup_line, type="l", col="black", lwd=2, pch=16, cex=1)

axis(side=4)
dev.off()


# Transition nodes ----
event_list_cons <- readRDS("objects/event_list_cons.rds")
event_df <- data.frame(node=0, genus=NA, height=0, min=0, max=0, Af2Ar=0, Ar2Af=0, Vic=0, ExtAf=0, ExtAr=0)

for (i in 1:length(event_list_cons)){
  if (nrow(event_list_cons[[i]]) > 0){
    for (j in 1:nrow(event_list_cons[[i]])){
      event_df <- rbind(event_df, event_list_cons[[i]][j,])
    }
  }
}

event_df <- event_df[-1,]
head(event_df)

myr <- 60

# Set color
color_all <- "#EE6A50"

# Plot all events in same color ----
pdf("plots/event_nodes.pdf", paper="a4r")
plot(1,type='n',xlim=c(myr,0),ylim=c(0,nrow(event_df)+5),xlab='Ma', ylab= "genus", main="Transition nodes", yaxt="n")
axis(2, at=c(nrow(event_df):1), labels=event_df$genus, col.axis="black", las=1, cex.axis=0.5)
j <- nrow(event_df)
for (i in 1:nrow(event_df)){
  segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_all, lwd=3)
  j <- j-1
}
dev.off()


# Define colors for different types of event ----
color_Af2Ar <- "#A8CB66"
color_Ar2Af <- "#2E8B57"
color_Vic <- "#FFC125"
color_ExtAf <- "#27408B"
color_ExtAr <- "#5CACEE"
event_cols <- c(color_Af2Ar, color_Ar2Af, color_Vic, color_ExtAf, color_ExtAr)
names(event_cols) <- c("Af2Ar", "Ar2Af", "Vic", "ExtAf", "ExtAr")

plot(c(1:5), c(1:5), col=event_cols, pch=16, cex=5)
?segments

# Plot with different colors for different types of event ----
pdf("plots/event_nodes_colors.pdf", paper="a4", height=20, width=10)
plot(1,type='n',xlim=c(myr,0),ylim=c(0,nrow(event_df)+5),xlab='Ma', ylab= "genus", main="Transition nodes", yaxt="n")
axis(2, at=c(nrow(event_df):1), labels=event_df$genus, col.axis="black", las=1, cex.axis=0.5)
j <- nrow(event_df)
for (i in 1:nrow(event_df)){
  if (event_df[i,]$Af2Ar==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Af2Ar, lwd=3)
  }
  if (event_df[i,]$Ar2Af==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Ar2Af, lwd=3)
  }
  if (event_df[i,]$Vic==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Vic, lwd=3)
  }
  if (event_df[i,]$ExtAf==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_ExtAf, lwd=3)
  }
  if (event_df[i,]$ExtAr==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_ExtAr, lwd=3)
  }
  j <- j-1
}

segments(x0=60, x1=55, y0=c(70,75,80,85,90)-10, y1=c(70, 75, 80, 85, 90)-10, col=event_cols, lwd=4)
text(x=55, y=c(70, 75, 80, 85, 90)-10, labels = names(event_cols), cex=0.7, pos = 4)


dev.off()

saveRDS(event_df, "objects/event_df.rds")
event_df <- readRDS("objects/event_df.rds")
colSums(event_df[,6:10])

# Plot cumulative in the background and the nodes on top (base R) ----
pdf("plots/cumulative_nodes.pdf", paper="a4", height=20, width=10)

# Plot cumulative number
par(mar=c(5,5,5,5))
plot(1,type='s',xlim=c(myr,0),ylim=c(0,nrow(event_df)+5),xlab='Ma', ylab= "genus", 
     main="Biogeographic events through time", yaxt="n", 
     frame = FALSE)
axis(side=4)
mtext("Cumulative events", side = 4, padj = 5)


for (s in 1:nsim){
  lines(spline(sim_rmat[[s]][, "Ma"]-1, sim_rmat[[s]][, "Ncum"]), type="l", col=color_sim, lwd=0.25, pch=16, cex=0.7)
}
lines(spline(rmat[, "Ma"]-1, rmat[, "Ncum"]), type="l", col=color_all, lwd=3)
#lines(mean_line, type="l", col="black", lwd=1, pch=16, cex=1)
lines(qlow_line, type="l", col="black", lwd=0.5, pch=16, cex=1)
lines(qup_line, type="l", col="black", lwd=0.5, pch=16, cex=1)


# Plot event nodes with different colors
axis(2, at=c(nrow(event_df):1), labels=event_df$genus, col.axis="black", las=1, cex.axis=0.5)
j <- nrow(event_df)
for (i in 1:nrow(event_df)){
  if (event_df[i,]$Af2Ar==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Af2Ar, lwd=3)
  }
  if (event_df[i,]$Ar2Af==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Ar2Af, lwd=3)
  }
  if (event_df[i,]$Vic==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_Vic, lwd=3)
  }
  if (event_df[i,]$ExtAf==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_ExtAf, lwd=3)
  }
  if (event_df[i,]$ExtAr==1){
    segments(x0=event_df[i,]$min, y0=j, x1=event_df[i,]$max, y1=j, col=color_ExtAr, lwd=3)
  }
  j <- j-1
}

segments(x0=60, x1=55, y0=c(70,75,80,85,90)-10, y1=c(70, 75, 80, 85, 90)-10, col=event_cols, lwd=4)
text(x=55, y=c(70, 75, 80, 85, 90)-10, labels = names(event_cols), cex=0.7, pos = 4)


dev.off()







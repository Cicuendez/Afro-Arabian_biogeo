



libs <- c("tidyverse", "deeptime", "here", "cowplot")
lapply(libs, require, character.only = TRUE)
# https://github.com/willgearty/deeptime

setwd()

# Import event objects ----
sim_event_all_cons <- readRDS("objects/sim_event_all_cons.rds")
event_all_cons <- readRDS("objects/event_all_cons.rds")
Af2Ar_sim_all_cons <- readRDS("objects/Af2Ar_sim_all_cons.rds")
Af2Ar_all_cons <- readRDS("objects/Af2Ar_all_cons.rds")
Ar2Af_sim_all_cons <- readRDS("objects/Ar2Af_sim_all_cons.rds")
Ar2Af_all_cons <- readRDS("objects/Ar2Af_all_cons.rds")
vicariance_sim_all_cons <- readRDS("objects/vicariance_sim_all_cons.rds")
vicariance_all_cons <- readRDS("objects/vicariance_all_cons.rds")
extirp_from_Af_sim_all_cons <- readRDS("objects/extirp_from_Af_sim_all_cons.rds")
extirp_from_Af_all_cons <- readRDS("objects/extirp_from_Af_all_cons.rds")
extirp_from_Ar_sim_all_cons <- readRDS("objects/extirp_from_Ar_sim_all_cons.rds")
extirp_from_Ar_all_cons <- readRDS("objects/extirp_from_Ar_all_cons.rds")


# Theme ----
# Set a customized theme for all the plots
theme_htc <- function(){
  theme_bw() +
    theme(panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          #        axis.text.x = element_blank(),
          #        axis.text.y = element_text(size = 5),
          #        axis.title.x = element_blank(),
          #        axis.title.y = element_blank(),
          plot.title = element_text(size = 13, vjust = 1, hjust = 0.5)
    )
}

# Geologic timescale ----
# Set the geologic period information 
data(periods)
data(epochs)

periods_htc <- periods
periods_htc$name[1] <- "Q"

epochs_htc <- epochs
epochs_htc$abbr[epochs_htc$abbr == "Plicn"] <- "Pli"
epochs_htc$abbr[epochs_htc$abbr == "Pls"] <- "Ple"
epochs_htc$abbr[epochs_htc$abbr == "Mc"] <- "M"
epochs_htc$abbr[epochs_htc$abbr == "Palcn"] <- "P"

# Quantiles information ----
prob_qup <- 0.975
prob_qlow <- 0.025
nsim <- 1000
myr <- 60


# Parameters (color, size, transparency) ----
# Colors
color_sim <- "gray93"
col_mean <- 'black'
col_q <- 'black'
color_Af2Ar <- "#A8CB66"
color_Ar2Af <- "#2E8B57"
color_Vic <- "#FFC125"
color_ExtAf <- "#27408B"
color_ExtAr <- "#5CACEE"
color_all <- "#EE6A50"

# Line width (size)
lwd_q <- 0.5
lwd_mean <- 0.3
lwd_events <- 2
lwd_sim <- 0.3

# Transparency (alpha)
alpha_sim <- 0.9

######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# ALL BIOGEOGRAPHIC EVENTS ----

# Create the lines to plot ----
# spline consensus
line_all <- data.frame(spline(event_all_cons$Ma-1, event_all_cons$N))
colnames(line_all) <- c("time", "events")

# spline simulations
lines_sim_all <- vector("list", nsim)
for (i in 1:length(sim_event_all_cons)){
  lines_sim_all[[i]] <- data.frame(spline(sim_event_all_cons[[i]]$Ma-1,
                                          sim_event_all_cons[[i]]$N))
  colnames(lines_sim_all[[i]]) <- c("time", "events")
}
#head(lines_sim_all)

# Quantiles and mean ----
# 95% CI ALL EVENTS
# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist_all <- vector("list", myr)
for (i in 1:myr){
  dlist_all[[i]] <- vector("numeric", nsim)
}

for (i in 1:nsim){
  for (j in 1:myr){
    dlist_all[[j]][i] <- sim_event_all_cons[[i]]$N[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df_all <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df_all$mean[i] <- mean(dlist_all[[i]])
  Q_df_all$qlow[i] <- quantile(dlist_all[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df_all$qup[i] <- quantile(dlist_all[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line_all <- data.frame(spline(Q_df_all$Ma-1, Q_df_all$mean))
qup_line_all <- data.frame(spline(Q_df_all$Ma-1, Q_df_all$qup))
qlow_line_all <- data.frame(spline(Q_df_all$Ma-1, Q_df_all$qlow))
colnames(mean_line_all) <- colnames(qup_line_all) <- colnames(qlow_line_all) <- 
  c("time", "events")

# Plot! ----
(plot_all <- ggplot() +
  
  # simulations
  geom_line(bind_rows(lines_sim_all, .id="sim"), 
            mapping = aes(x=time, y=events, group=sim), 
            color=color_sim, alpha=alpha_sim, size=lwd_sim) +
  #observed events
  geom_line(data = line_all, aes(x = time, y = events), color = color_all,
             size = lwd_events) + 
  
  # Mean and quantiles
  geom_line(mean_line_all, mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
  geom_line(qup_line_all, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
  geom_line(qlow_line_all, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
  
  xlim(60,0) +
  labs(x = "Time before present (Ma)", y = "Number of events") +
  ggtitle("All biogeographic events") +
  
  # Insert geologic scale
  coord_geo(xlim = c(60, 0), ylim = c(0,40), pos = as.list(rep("bottom", 2)),
            dat = list(epochs_htc, periods_htc),
            height = list(unit(1, "lines"), unit(1, "line")),
            rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
            skip = c('Holocene'), 
            lab = TRUE) +
  
  # Set the theme 
  theme_htc() +
  ggsave("plots/plot_all.pdf")
)

######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# DISPERSAL FROM AFRICA TO ARABIA ----

# Create the lines to plot ----
# spline consensus
line_Af2Ar <- data.frame(spline(Af2Ar_all_cons$Ma-1, Af2Ar_all_cons$N))
colnames(line_Af2Ar) <- c("time", "events")

# spline simulations
lines_sim_Af2Ar <- vector("list", nsim)
for (i in 1:length(Af2Ar_sim_all_cons)){
  lines_sim_Af2Ar[[i]] <- data.frame(spline(Af2Ar_sim_all_cons[[i]]$Ma-1,
                                          Af2Ar_sim_all_cons[[i]]$N))
  colnames(lines_sim_Af2Ar[[i]]) <- c("time", "events")
}
#head(lines_sim_Af2Ar)

# Quantiles and mean ----
# 95% CI ALL EVENTS
# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist_Af2Ar <- vector("list", myr)
for (i in 1:myr){
  dlist_Af2Ar[[i]] <- vector("numeric", nsim)
}

for (i in 1:nsim){
  for (j in 1:myr){
    dlist_Af2Ar[[j]][i] <- Af2Ar_sim_all_cons[[i]]$N[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df_Af2Ar <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df_Af2Ar$mean[i] <- mean(dlist_Af2Ar[[i]])
  Q_df_Af2Ar$qlow[i] <- quantile(dlist_Af2Ar[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df_Af2Ar$qup[i] <- quantile(dlist_Af2Ar[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line_Af2Ar <- data.frame(spline(Q_df_Af2Ar$Ma-1, Q_df_Af2Ar$mean))
qup_line_Af2Ar <- data.frame(spline(Q_df_Af2Ar$Ma-1, Q_df_Af2Ar$qup))
qlow_line_Af2Ar <- data.frame(spline(Q_df_Af2Ar$Ma-1, Q_df_Af2Ar$qlow))
colnames(mean_line_Af2Ar) <- colnames(qup_line_Af2Ar) <- colnames(qlow_line_Af2Ar) <- 
  c("time", "events")

# Plot! ----
(plot_Af2Ar <- ggplot() +
   
   # simulations
   geom_line(bind_rows(lines_sim_Af2Ar, .id="sim"), 
             mapping = aes(x=time, y=events, group=sim), 
             color=color_sim, alpha=alpha_sim, size=lwd_sim) +
   
   #observed events
   geom_line(data = line_Af2Ar, aes(x = time, y = events), color = color_Af2Ar,
             size = lwd_events) +
   
   # Mean and quantiles
   geom_line(mean_line_Af2Ar, mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
   geom_line(qup_line_Af2Ar, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   geom_line(qlow_line_Af2Ar, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   
   xlim(60,0) +
   labs(x = "Time before present (Ma)", y = "Number of events") +
   ggtitle("Dispersal from Africa to Arabia") +
   
   # Insert geologic scale
   coord_geo(xlim = c(60, 0), ylim = c(0,27), pos = as.list(rep("bottom", 2)),
             dat = list(epochs_htc, periods_htc),
             height = list(unit(1, "lines"), unit(1, "line")),
             rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
             skip = c('Holocene'), 
             lab = TRUE) +
   
   # Set the theme 
   theme_htc() +
   
   # Save it
   ggsave("plots/plot_Af2Ar.pdf")
)

######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# DISPERSAL FROM ARABIA TO AFRICA ----

# Create the lines to plot ----
# spline consensus
line_Ar2Af <- data.frame(spline(Ar2Af_all_cons$Ma-1, Ar2Af_all_cons$N))
colnames(line_Ar2Af) <- c("time", "events")

# spline simulations
lines_sim_Ar2Af <- vector("list", nsim)
for (i in 1:length(Ar2Af_sim_all_cons)){
  lines_sim_Ar2Af[[i]] <- data.frame(spline(Ar2Af_sim_all_cons[[i]]$Ma-1,
                                            Ar2Af_sim_all_cons[[i]]$N))
  colnames(lines_sim_Ar2Af[[i]]) <- c("time", "events")
}
#head(lines_sim_Ar2Af)

# Quantiles and mean ----
# 95% CI ALL EVENTS
# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist_Ar2Af <- vector("list", myr)
for (i in 1:myr){
  dlist_Ar2Af[[i]] <- vector("numeric", nsim)
}

for (i in 1:nsim){
  for (j in 1:myr){
    dlist_Ar2Af[[j]][i] <- Ar2Af_sim_all_cons[[i]]$N[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df_Ar2Af <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df_Ar2Af$mean[i] <- mean(dlist_Ar2Af[[i]])
  Q_df_Ar2Af$qlow[i] <- quantile(dlist_Ar2Af[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df_Ar2Af$qup[i] <- quantile(dlist_Ar2Af[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line_Ar2Af <- data.frame(spline(Q_df_Ar2Af$Ma-1, Q_df_Ar2Af$mean))
qup_line_Ar2Af <- data.frame(spline(Q_df_Ar2Af$Ma-1, Q_df_Ar2Af$qup))
qlow_line_Ar2Af <- data.frame(spline(Q_df_Ar2Af$Ma-1, Q_df_Ar2Af$qlow))
colnames(mean_line_Ar2Af) <- colnames(qup_line_Ar2Af) <- colnames(qlow_line_Ar2Af) <- 
  c("time", "events")

# Plot! ----
(plot_Ar2Af <- ggplot() +
   
   # simulations
   geom_line(bind_rows(lines_sim_Ar2Af, .id="sim"), 
             mapping = aes(x=time, y=events, group=sim), 
             color=color_sim, alpha=alpha_sim, size=lwd_sim) +
   
   # observed events
   geom_line(data = line_Ar2Af, aes(x = time, y = events), color = color_Ar2Af,
             size = lwd_events) +
   
   # Mean and quantiles
   geom_line(mean_line_Ar2Af, mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
   geom_line(qup_line_Ar2Af, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   geom_line(qlow_line_Ar2Af, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   
   xlim(60,0) +
   labs(x = "Time before present (Ma)", y = "Number of events") +
   ggtitle("Dispersal from Arabia to Africa") +
   
   # Insert geologic scale
   coord_geo(xlim = c(60, 0), ylim = c(0,15), pos = as.list(rep("bottom", 2)),
             dat = list(epochs_htc, periods_htc),
             height = list(unit(1, "lines"), unit(1, "line")),
             rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
             skip = c('Holocene'), 
             lab = TRUE) +
   
   # Set the theme 
   theme_htc() +
   
   # Save it
   ggsave("plots/plot_Ar2Af.pdf")
)

######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# VICARIANCE ----

# Create the lines to plot ----
# spline consensus
line_vicariance <- data.frame(spline(vicariance_all_cons$Ma-1, vicariance_all_cons$N))
colnames(line_vicariance) <- c("time", "events")

# spline simulations
lines_sim_vicariance <- vector("list", nsim)
for (i in 1:length(vicariance_sim_all_cons)){
  lines_sim_vicariance[[i]] <- data.frame(spline(vicariance_sim_all_cons[[i]]$Ma-1,
                                                 vicariance_sim_all_cons[[i]]$N))
  colnames(lines_sim_vicariance[[i]]) <- c("time", "events")
}
#head(lines_sim_vicariance)

# Quantiles and mean ----
# 95% CI ALL EVENTS
# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist_vicariance <- vector("list", myr)
for (i in 1:myr){
  dlist_vicariance[[i]] <- vector("numeric", nsim)
}

for (i in 1:nsim){
  for (j in 1:myr){
    dlist_vicariance[[j]][i] <- vicariance_sim_all_cons[[i]]$N[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df_vicariance <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df_vicariance$mean[i] <- mean(dlist_vicariance[[i]])
  Q_df_vicariance$qlow[i] <- quantile(dlist_vicariance[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df_vicariance$qup[i] <- quantile(dlist_vicariance[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line_vicariance <- data.frame(spline(Q_df_vicariance$Ma-1, Q_df_vicariance$mean))
qup_line_vicariance <- data.frame(spline(Q_df_vicariance$Ma-1, Q_df_vicariance$qup))
qlow_line_vicariance <- data.frame(spline(Q_df_vicariance$Ma-1, Q_df_vicariance$qlow))
colnames(mean_line_vicariance) <- colnames(qup_line_vicariance) <- colnames(qlow_line_vicariance) <- 
  c("time", "events")

# Plot! ----
(plot_vicariance <- ggplot() +
   
   # simulations
   geom_line(bind_rows(lines_sim_vicariance, .id="sim"), 
             mapping = aes(x=time, y=events, group=sim), 
             color=color_sim, alpha=alpha_sim, size=lwd_sim) +
   
   # observed events
   geom_line(data = line_vicariance, aes(x = time, y = events), color = color_Vic,
             size = lwd_events) +
   
   # Mean and quantiles
   geom_line(mean_line_vicariance, mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
   geom_line(qup_line_vicariance, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   geom_line(qlow_line_vicariance, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   
   xlim(60,0) +
   labs(x = "Time before present (Ma)", y = "Number of events") +
   ggtitle("Vicariance") +
   
   # Insert geologic scale
   coord_geo(xlim = c(60, 0), ylim = c(0,12), pos = as.list(rep("bottom", 2)),
             dat = list(epochs_htc, periods_htc),
             height = list(unit(1, "lines"), unit(1, "line")),
             rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
             skip = c('Holocene'), 
             lab = TRUE) +
   
   # Set the theme 
   theme_htc() +
   
   # Save it
   ggsave("plots/plot_vicariance.pdf")
)

######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# EXTIRPATION FROM AFRICA ----

# Create the lines to plot ----
# spline consensus
line_extAf <- data.frame(spline(extirp_from_Af_all_cons$Ma-1, extirp_from_Af_all_cons$N))
colnames(line_extAf) <- c("time", "events")

# spline simulations
lines_sim_extAf <- vector("list", nsim)
for (i in 1:length(extirp_from_Af_sim_all_cons)){
  lines_sim_extAf[[i]] <- data.frame(spline(extirp_from_Af_sim_all_cons[[i]]$Ma-1,
                                            extirp_from_Af_sim_all_cons[[i]]$N))
  colnames(lines_sim_extAf[[i]]) <- c("time", "events")
}
#head(lines_sim_vicariance)

# Quantiles and mean ----
# 95% CI ALL EVENTS
# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist_extAf <- vector("list", myr)
for (i in 1:myr){
  dlist_extAf[[i]] <- vector("numeric", nsim)
}

for (i in 1:nsim){
  for (j in 1:myr){
    dlist_extAf[[j]][i] <- extirp_from_Af_sim_all_cons[[i]]$N[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df_extAf <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df_extAf$mean[i] <- mean(dlist_extAf[[i]])
  Q_df_extAf$qlow[i] <- quantile(dlist_extAf[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df_extAf$qup[i] <- quantile(dlist_extAf[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line_extAf <- data.frame(spline(Q_df_extAf$Ma-1, Q_df_extAf$mean))
qup_line_extAf <- data.frame(spline(Q_df_extAf$Ma-1, Q_df_extAf$qup))
qlow_line_extAf <- data.frame(spline(Q_df_extAf$Ma-1, Q_df_extAf$qlow))
colnames(mean_line_extAf) <- colnames(qup_line_extAf) <- colnames(qlow_line_extAf) <- 
  c("time", "events")

# Plot! ----
(plot_extAf <- ggplot() +
   
   # simulations
   geom_line(bind_rows(lines_sim_extAf, .id="sim"), 
             mapping = aes(x=time, y=events, group=sim), 
             color=color_sim, alpha=alpha_sim, size=lwd_sim) +
   
   # observed events
   geom_line(data = line_extAf, aes(x = time, y = events), color = color_ExtAf,
             size = lwd_events) +
   
   # Mean and quantiles
   geom_line(mean_line_extAf, mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
   geom_line(qup_line_extAf, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   geom_line(qlow_line_extAf, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   
   xlim(60,0) +
   labs(x = "Time before present (Ma)", y = "Number of events") +
   ggtitle("Extirpation from Africa") +
   
   # Insert geologic scale
   coord_geo(xlim = c(60, 0), ylim = c(0,12), pos = as.list(rep("bottom", 2)),
             dat = list(epochs_htc, periods_htc),
             height = list(unit(1, "lines"), unit(1, "line")),
             rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
             skip = c('Holocene'), 
             lab = TRUE) +
   
   # Set the theme 
   theme_htc() +
   
   # Save it
   ggsave("plots/plot_extAf.pdf")
)

######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# EXTIRPATION FROM ARABIA ----

# Create the lines to plot ----
# spline consensus
line_extAr <- data.frame(spline(extirp_from_Ar_all_cons$Ma-1, extirp_from_Ar_all_cons$N))
colnames(line_extAr) <- c("time", "events")

# spline simulations
lines_sim_extAr <- vector("list", nsim)
for (i in 1:length(extirp_from_Ar_sim_all_cons)){
  lines_sim_extAr[[i]] <- data.frame(spline(extirp_from_Ar_sim_all_cons[[i]]$Ma-1,
                                            extirp_from_Ar_sim_all_cons[[i]]$N))
  colnames(lines_sim_extAr[[i]]) <- c("time", "events")
}
#head(lines_sim_vicariance)

# Quantiles and mean ----
# 95% CI ALL EVENTS
# Let's do a list where each element is a vector with the
# number of transitions per simulations for each Ma.
# 60 vectors (60 Ma), and 1000 numbers in each vector.
dlist_extAr <- vector("list", myr)
for (i in 1:myr){
  dlist_extAr[[i]] <- vector("numeric", nsim)
}

for (i in 1:nsim){
  for (j in 1:myr){
    dlist_extAr[[j]][i] <- extirp_from_Ar_sim_all_cons[[i]]$N[j]
  }
}

# Now we can create a dataframe with the quantiles and 
# the mean per Ma, which we will calculate with the different
# elements of the list of distributions per Ma (dlist).
Q_df_extAr <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df_extAr$mean[i] <- mean(dlist_extAr[[i]])
  Q_df_extAr$qlow[i] <- quantile(dlist_extAr[[i]], probs=c(prob_qlow, prob_qup))[1]
  Q_df_extAr$qup[i] <- quantile(dlist_extAr[[i]], probs=c(prob_qlow, prob_qup))[2]
}

mean_line_extAr <- data.frame(spline(Q_df_extAr$Ma-1, Q_df_extAr$mean))
qup_line_extAr <- data.frame(spline(Q_df_extAr$Ma-1, Q_df_extAr$qup))
qlow_line_extAr <- data.frame(spline(Q_df_extAr$Ma-1, Q_df_extAr$qlow))
colnames(mean_line_extAr) <- colnames(qup_line_extAr) <- colnames(qlow_line_extAr) <- 
  c("time", "events")

# Plot! ----
(plot_extAr <- ggplot() +
   
   # simulations
   geom_line(bind_rows(lines_sim_extAr, .id="sim"), 
             mapping = aes(x=time, y=events, group=sim), 
             color=color_sim, alpha=alpha_sim, size=lwd_sim) +
   
   # observed events
   geom_line(data = line_extAr, aes(x = time, y = events), color = color_ExtAr,
             size = lwd_events) +
   
   # Mean and quantiles
   geom_line(mean_line_extAr, mapping=aes(x=time, y=events), color=col_mean, size=lwd_mean) +
   geom_line(qup_line_extAr, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   geom_line(qlow_line_extAr, mapping=aes(x=time, y=events), color=col_q, size=lwd_q) +
   
   xlim(60,0) +
   labs(x = "Time before present (Ma)", y = "Number of events") +
   ggtitle("Extirpation from Arabia") +
   
   # Insert geologic scale
   coord_geo(xlim = c(60, 0), ylim = c(0,12), pos = as.list(rep("bottom", 2)),
             dat = list(epochs_htc, periods_htc),
             height = list(unit(1, "lines"), unit(1, "line")),
             rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
             skip = c('Holocene'), 
             lab = TRUE) +
   
   # Set the theme 
   theme_htc() +
   
   # Save it
   ggsave("plots/plot_extAr.pdf")
)



######### :::::::::::::::::::::::::::::::::::::::::::::::::::::: #############
# PLOT GRID ----

event_plotlist <- list(plot_all, plot_Ar2Af, plot_extAf, 
                       plot_vicariance, plot_Af2Ar, plot_extAr)
event_plotgrid <- plot_grid(plotlist = event_plotlist,
                            ncol = 3, 
                            labels = 'auto')
ggsave("plots/all_events_plotgrid.pdf", 
       event_plotgrid, height = 7, width = 10)
ggsave("plots/all_events_plotgrid.png", 
       event_plotgrid, height = 7, width = 10)







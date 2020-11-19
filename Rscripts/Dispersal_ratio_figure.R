
# Dispersal ratio figure ----
# Plot the percentage of dispersal in both directions.

setwd()

# load libraries ----
libs <- c("tidyverse", "deeptime", "here", "cowplot")
lapply(libs, require, character.only = TRUE)


# import data ----
Af2Ar_sim_all_cons <- readRDS("objects/Af2Ar_sim_all_cons.rds")
Af2Ar_all_cons <- readRDS("objects/Af2Ar_all_cons.rds")
Ar2Af_sim_all_cons <- readRDS("objects/Ar2Af_sim_all_cons.rds")
Ar2Af_all_cons <- readRDS("objects/Ar2Af_all_cons.rds")

time_interval <- 1:60

# set colors ----
color_Af2Ar <- "#A8CB66"
color_Ar2Af <- "#2E8B57"
color_sim <- "grey80"
color_ratio <- "purple"
disp_colors <- c(Af2Ar="#A8CB66", Ar2Af="#2E8B57")

# customize consensus dataframe
Ar2Af_all_cons$direction <- "Ar2Af"
Af2Ar_all_cons$direction <- "Af2Ar"

both_disp <- rbind(Ar2Af_all_cons, Af2Ar_all_cons)

# Compute percentages with dplyr
data <- both_disp  %>%
  group_by(Ma, direction) %>%
  summarise(n = sum(N)) %>%
  mutate(percentage = n / sum(n))

ggplot() + 
  geom_area(data=data, mapping=aes(x=-Ma-1, y=percentage, fill=direction), alpha=0.6 , size=1, colour="transparent") +
  scale_fill_manual(values=disp_colors) +
  geom_line(mapping=aes(x=-time_interval-1, y=0.5), color="black") +
  geom_line(mapping=aes(x=-time_interval-1, y=data$percentage[data$direction=="Ar2Af"]), color="purple") +
  theme_classic()


# SIMULATIONS ----

# customize consensus dataframe
for (i in 1:length(Ar2Af_sim_all_cons)){
  Ar2Af_sim_all_cons[[i]]$direction <- "Ar2Af"
  Af2Ar_sim_all_cons[[i]]$direction <- "Af2Ar"
}

both_disp_sim <- vector("list", 1000)
for (i in 1:length(both_disp_sim)){
  both_disp_sim[[i]] <- rbind(Ar2Af_sim_all_cons[[i]], Af2Ar_sim_all_cons[[i]])
}
both_disp_sim[[1]]

# Compute percentages with dplyr
data_sim <- vector("list", 1000)

for (i in 1:length(data_sim)){
  data_sim[[i]] <- both_disp_sim[[i]] %>%
    group_by(Ma, direction) %>%
    summarise(n=sum(N)) %>%
    mutate(percentage=n/sum(n))
}

head(data_sim)

data_sim_Ar2Af <- vector("list", 1000)
for (i in 1:length(data_sim_Ar2Af)){
  data_sim_Ar2Af[[i]] <- data_sim[[i]] %>%
    filter(direction=="Ar2Af")
}
head(data_sim_Ar2Af)



ggplot() +
  
  scale_fill_manual(values=disp_colors) +
  geom_line(bind_rows(data_sim_Ar2Af, .id="sim"), 
            mapping=aes(x=-Ma-1, y=percentage, group=sim), color=color_sim, alpha=0.9, size=0.5) +
  geom_line(mapping=aes(x=-time_interval-1, y=data$percentage[data$direction=="Ar2Af"]), 
            color=color_ratio, size=1) +
  geom_line(mapping=aes(x=-time_interval-1, y=0.5), color="black") +
  geom_area(data=data, mapping=aes(x=-Ma-1, y=percentage, fill=direction),
            alpha=0.5, size=1, color="transparent") +
  
  
  theme_classic()


# 95%

prob_qup <- 0.975
prob_qlow <- 0.025
nsim <- 1000
myr <- 60

data_sim_Ar2Af[[1]]

# Create a list of 60 vectors (one per Ma), and each vector with 1000 elements (1000 sims)
dlist <- vector("list", myr)
for (i in 1:myr){
  dlist[[i]] <- vector("numeric", nsim)
}


# In each element, put the percentage of dispersals Ar2Af
for (i in 1:nsim){
  for (j in 1:myr){
    dlist[[j]][i] <- data_sim_Ar2Af[[i]]$percentage[j]
  }
}

Q_df <- data.frame(Ma=c(1:myr), qlow=0, qup=0, mean=0)
for (i in 1:myr){
  Q_df$mean[i] <- mean(dlist[[i]])
  Q_df$qlow[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup), na.rm = T)[1]
  Q_df$qup[i] <- quantile(dlist[[i]], probs=c(prob_qlow, prob_qup), na.rm = T)[2]
}

mean_line <- spline(Q_df$Ma-1, Q_df$mean)
qup_line <- spline(Q_df$Ma-1, Q_df$qup)
qlow_line <- spline(Q_df$Ma-1, Q_df$qlow)


percent_Ar2Af <- data %>%
  filter(direction=="Ar2Af")

line_percent_Ar2Af <- spline(percent_Ar2Af$Ma-1, percent_Ar2Af$percentage)

line_list_sim_Ar2Af <- vector("list", nsim)
for (i in 1:nsim){
  line_list_sim_Ar2Af[[i]] <- data.frame(spline(data_sim_Ar2Af[[i]]$Ma-1, data_sim_Ar2Af[[i]]$percentage))
}


(dispersal_ratio_plot <- ggplot() +
  geom_line(bind_rows(line_list_sim_Ar2Af, .id="sim"), 
            mapping=aes(x=-(x-1), y=y, group=sim), color=color_sim, alpha=0.9, size=0.5) +
  geom_line(mapping=aes(x=-(time_interval-1), y=0.5), color="black") +
  geom_line(data=data.frame(qlow_line), aes(x=-x, y=y), color="black") +
  geom_line(data=data.frame(qup_line), aes(x=-x, y=y), color="black") +
  geom_area(data=data, mapping=aes(x=-(Ma-1), y=percentage, fill=direction),
            alpha=0.5, size=1, color="transparent") +
  scale_fill_manual(values=disp_colors) +
  geom_line(data=data.frame(line_percent_Ar2Af), aes(x=-x, y=y), color=color_Ar2Af, size=2) +
  theme_classic()
)

pdf("plots/dispersal_ratio_plot2.pdf", paper="a4r", height=10, width=20)
dispersal_ratio_plot
dev.off()

# Plot without grey lines in the background
(dispersal_ratio_plot_no_sims <- ggplot() +
    geom_line(mapping=aes(x=-(time_interval-1), y=0.5), color="black") +
    geom_line(data=data.frame(qlow_line), aes(x=-x, y=y), color="black") +
    geom_line(data=data.frame(qup_line), aes(x=-x, y=y), color="black") +
    geom_area(data=data, mapping=aes(x=-(Ma-1), y=percentage, fill=direction),
              alpha=0.5, size=1, color="transparent") +
    scale_fill_manual(values=disp_colors) +
    geom_line(data=data.frame(line_percent_Ar2Af), aes(x=-x, y=y), color=color_Ar2Af, size=2) +
    
    
    theme_classic()
)


pdf("plots/dispersal_ratio_plot_no_sims.pdf", paper="a4r", height=10, width=20)
dispersal_ratio_plot_no_sims
dev.off()



qup_qlow_df <- data.frame(qupx = qup_line$x, qupy = qup_line$y,
                          qlowx = qlow_line$x, qlowy = qlow_line$y)


(dispersal_ratio_plot_no_sims2 <- ggplot() +
    geom_ribbon(data=qup_qlow_df, aes(ymin=qlowy, ymax=qupy, x=-qupx), fill="gray80") +
    geom_area(data=data, mapping=aes(x=-(Ma-1), y=percentage, fill=direction),
              alpha=0.5, size=1, color="transparent") +
    scale_fill_manual(values=disp_colors) +
    geom_line(data=data.frame(line_percent_Ar2Af), aes(x=-x, y=y), color=color_Ar2Af, size=2) +
    geom_line(mapping=aes(x=-(time_interval-1), y=0.5), color="black", alpha = 0.5) +
    xlab("Ma") +
    ylab("Percentage of dispersal from Arabia to Africa") +
    ggtitle("Dispersal asymmetry") +
    theme_classic()
)

pdf("plots/dispersal_ratio_plot_no_sims2.pdf", paper="a4r", height=10, width=20)
dispersal_ratio_plot_no_sims2
dev.off()


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



# Gray shade 
dispersal_ratio_plot_no_sims3 <- ggplot() +
    geom_ribbon(data=qup_qlow_df, aes(ymin=qlowy, ymax=qupy, x=qupx), fill="gray80") +
    geom_line(data=data.frame(qlow_line), aes(x=x, y=y), color="gray50") +
    geom_line(data=data.frame(qup_line), aes(x=x, y=y), color="gray50") +
    geom_area(data=data, mapping=aes(x=(Ma-1), y=percentage, fill=direction),
              alpha=0.5, size=1, color="transparent") +
    scale_fill_manual(values=disp_colors) +
    geom_line(data=data.frame(line_percent_Ar2Af), aes(x=x, y=y), color=color_Ar2Af, size=2) +
    geom_line(mapping=aes(x=(time_interval-1), y=0.5), color="black", alpha = 0.5) +
    xlab("Ma") +
    ylab("Percentage of dispersal from Arabia to Africa") +
    ggtitle("Dispersal asymmetry") +
    xlim(40, 0) +
    
    # Insert geologic scale
    coord_geo(xlim = c(40, 0), ylim = c(0,1), pos = as.list(rep("bottom", 2)),
              dat = list(epochs_htc, periods_htc),
              height = list(unit(1, "lines"), unit(1, "line")),
              rot = list(0, 0), size = list(2, 3), abbrv = list(TRUE, FALSE), 
              skip = c('Holocene'), 
              lab = TRUE) +
    theme_htc() +
    ggsave("plots/dispersal_ratio_plot_no_sims3_2.pdf", paper="a4r", height=10, width=20)


pdf("plots/dispersal_ratio_plot_no_sims3_2.pdf", paper="a4r", height=10, width=20)
dispersal_ratio_plot_no_sims3
dev.off()


library(sf)
library(geosphere)
library(tidyverse)
data=read.csv("C:/Users/qu_la/Desktop/hkqm/af_bio.csv")

target=c("HiC_scaffold_5.29849774","HiC_scaffold_5.29977431","HiC_scaffold_5.44065996","HiC_scaffold_19.55471519")
data=data[,c("long","lat","pop","site_number",target)]

#set up the transect which follows the elevation and climate gradient
point_A <- c(70.818, 48.744) 
point_B <- c(118.842, 22.35) 
samples = data[,c("long","lat","site_number")]
samples_sf <- sf::st_as_sf(samples, coords = c("long", "lat"), crs = 4326)
transect_line <- st_sfc(st_linestring(rbind(point_A, point_B)), crs = 4326)

#project sampling site onto transect line
projected_points <-sf:: st_nearest_points(samples_sf, transect_line) %>%
  st_cast("POINT") %>% .[seq(2, length(.), by = 2)] 

distances_km <- geosphere::distGeo(matrix(rep(point_A, length(projected_points)), 
                                          ncol = 2, byrow = TRUE),
                                   st_coordinates(projected_points)) / 1000
samples$transect_distance_km <- distances_km
samples=arrange(samples,transect_distance_km)
samples_info=dplyr::left_join(samples,data,"site_number")

samples_info$label <- samples_info$transect_distance_km
samples_info$label <- as.character(round(samples_info$label, 0))
samples_info$HiC_scaffold_5.29977431=1-samples_info$HiC_scaffold_5.29977431
samples_info_long <- tidyr::pivot_longer(samples_info,
                                         cols=target,
                          # cols = colnames(samples_info)[9:708],
                          names_to = "Variable",
                          values_to = "Value")

###plot cline
samples_info_long$label=as.numeric(samples_info_long$label)
ggplot()+
stat_smooth(
    data=filter(samples_info_long,pop!="western"),
    aes(y=Value,x=label,
        group = Variable, color= Variable  
       ),
    method = "glm",
    method.args = list(family = "binomial"),
    formula = y ~ x,
    se = F,
    size = 0.8) +
  geom_vline(xintercept = c(3924.423,4382.152), linetype = "dashed", color = "grey90", size = 1)+
  theme_bw()+
  theme(aspect.ratio = 1.1,text = element_text(size = 25))+
  scale_x_continuous(
    breaks = as.numeric(filter(samples_info_long,pop!="western")$label)%>%unique(),
    labels = c(1:11)
  )  +
  theme(
    panel.grid = element_blank()) +
  ylab("Allele frequency")+xlab("Site")


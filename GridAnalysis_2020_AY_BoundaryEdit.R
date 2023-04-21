#####     0. Load Libraries         ###########################################################
require(rgdal);require(sf)
require(data.table);require(dplyr)
require(ggplot2); require(reshape2)
require(tidyverse); require(gapminder)
setwd("C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/ImpactedAreaCalculations")
## Make a list of regex search terms for taxa
rens <- c('B.*Coral','Bolt','Bryo','L.*Gorg','SeaP','S.*Gorg','Sponge')

##### Link to database ########################################################################

# The input file geodatabase
fgdb <- "C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/ImpactedAreaCalculations/SAI_2020.gdb"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
fc_list <- fc_list[grep("Biomass", fc_list, fixed=T)]
print(fc_list)
fc_list <- fc_list[order(fc_list)]


#####     3. Get KDE data     #################################################################

kde <- st_read(fgdb, 
               layer = 'KDE_Grid')

#####     1. Threshold Fishing effort     #####################################################

### Load cutoffs 
#cutoffs <- as.data.table(openxlsx::read.xlsx('NEREIDA2020_Results.xlsx',sheet='Cutoffs4GridAY'))
cutoffs <- read.csv("C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/ImpactedAreaCalculations/impacted_area_thresholds.csv")
head(cutoffs)
###level might be counter-intuitive - it means the proporion of biomass that is exposed to risk, not how much is protected 
##hence '0' has the lowest value of trawl intensity

#names(cutoffs)[3:4] <- c('AY','AYRO')
#cutoffs <- cutoffs[1:7,]
# Taxonlist
taxlist       <- unique(cutoffs$Tax)
thresholdlist <- unique(cutoffs$level)
mlist         <- "AYRO_byLevel"

# Load Fishig effort Grid
eff <- st_read(fgdb, 
               layer = "Fishingeffort_10_19")

# Create a new column for each taxon - cutoff method combination
for (i in taxlist) {
  for (v in thresholdlist){
    
    flab <- paste0(i,"_",v)

    #eff$threshold = v

    eff[eff$TrawlKmYr <  cutoffs$Value[cutoffs$Tax == i & cutoffs$level == v], flab] <- 'Impacted (below cut-off)'
    eff[eff$TrawlKmYr >= cutoffs$Value[cutoffs$Tax == i & cutoffs$level == v], flab] <- 'Impacted (above cut-off)'
    eff[eff$TrawlKmYr == 0, flab] <- 'At Risk'
    eff[eff$Closure == 1 & eff$Footprint == 1, flab] <- 'Protected - CIF'
    eff[eff$Closure == 1 & eff$Footprint == 0, flab] <- 'Protected - COF'
    eff[eff$Closure == 0 & eff$Footprint == 0, flab] <- 'Protected - OFF'
 }
}

#counts of different cells per level

#blc <- eff[,c(1:10,12:29)]
#bol <- eff[,c(1:9,30:48)]
#bry <- eff[,c(1:9,49:67)]
#lrggor <- eff[,c(1:9,68:86)]
#spn    <- eff[,c(1:9,87:105)]
#smggor <- eff[,c(1:9,106:124)]
#spg    <- eff[,c(1:9,125:143)]

# Check categories
summary(as.factor(eff$BlackCoral_0.3))
summary(as.factor(eff$BlackCoral_0.5))
summary(as.factor(eff$BlackCoral_0.7))
summary(as.factor(eff$BlackCoral_0.9)) #as level increases, count of cells below threshold should increase

#summarise data

area_counts <- data.frame()
for(a in seq(10,ncol(eff),1)){
  name <- names(eff[,a])[1]
  tmp <- as.data.frame(eff[,which(colnames(eff) == name)])
  area_counts <- rbind(area_counts, summary(as.factor(tmp[,1])))
}
colnames(area_counts) <- c("At Risk", "Impacted (below cut-off)", "Impacted (above cut-off)", "Protected - CIF", "Protected - COF", "Protected - OFF")

area_counts <- cbind(data.frame(Tax = cutoffs$Tax,
                            level = cutoffs$level),
                 area_counts)

area_counts$Protected <- area_counts$`Protected - CIF`+area_counts$`Protected - COF`+area_counts$`Protected - OFF`

area_counts <- melt(area_counts, id.vars = c("Tax", "level"), value.name = "count")
area_counts <- area_counts[area_counts$variable != "Protected - CIF" &
                           area_counts$variable != "Protected - COF" &
                           area_counts$variable != "Protected - OFF",]

area_counts <- area_counts %>%
  group_by(Tax, level, variable) %>%
   summarise(n = sum(count)) %>%
     mutate(percentage = n / sum(n))

area_counts$variable <- factor(area_counts$variable, levels = c("Protected", "At Risk", "Impacted (below cut-off)", "Impacted (above cut-off)"))

imp_area_plot <- 
    ggplot(area_counts, aes(x = level*100, y = percentage*100, fill = variable)) + 
           geom_area() + 
           scale_fill_manual(values = c("blue3", "yellow", "darkred", "Orange"),
                             breaks = c("Protected", "At Risk", "Impacted (above cut-off)", "Impacted (below cut-off)"),
                             labels = c("Protected", "At Risk", "Impacted (above cut-off)", "Impacted (below cut-off)")) +
           scale_x_continuous(breaks = seq(0,100,20)) +
           scale_y_continuous(breaks = seq(0,100,20)) +
           facet_wrap(~Tax, nrow = 2) +
           theme_bw() + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
           labs(x = "Extent of habitat exposed to SAI (%)", y = "Area within threat category (%)", fill = "Threat Category")
imp_area_plot
ggsave(imp_area_plot, filename = "impacted_area_by_threshold_global.pdf", height = 12, width = 18, units = "cm", device = "pdf", dpi = 150)
### Write to shapefile

st_write(eff,"effort_thresholds_globalKDE.shp",delete_dsn = TRUE)
#rm(eff)

### shapefiles for each VME KDE area separately
# Load Fishig effort Grid
eff2 <- st_read(fgdb, layer = "Fishingeffort_10_19")

# Create a new column for each taxon - cutoff method combination
for (i in 1:length(taxlist)) {
  
  # Determine colum name
  ft <- grep(rens[i],names(kde),value = TRUE)
  # Select the KDE polygons and convert to centroid
  kdet <- st_centroid(filter(kde, get(ft) == 1))
  
      for (v in thresholdlist){
    
            flab <- paste0(taxlist[i],"_",v,"_",mlist)
            print(flab)
            
            efft <- eff2[kdet,]
    
            efft[efft$TrawlKmYr <  cutoffs$Value[cutoffs$Tax == taxlist[i] & cutoffs$level == v], flab] <- 'Impacted (below cut-off)'
            efft[efft$TrawlKmYr >= cutoffs$Value[cutoffs$Tax == taxlist[i] & cutoffs$level == v], flab] <- 'Impacted (above cut-off)'
            efft[efft$TrawlKmYr == 0 & efft$Closure == 0 & efft$Footprint == 1, flab] <- 'At Risk'
            efft[efft$Closure == 1 | efft$Footprint == 0,flab] <- 'Protected or outside FF'

           # st_write(obj = efft,
            #         dsn = paste0("P:/C8408_NAFO_R_and_D_Support/Working_Area/ImpactedAreaCalculations/areas_by_threshold/ImpactedArea_",flab,".shp"),
             #        delete_dsn = TRUE, quiet = T)
            
            st_write(obj = efft,
                     dsn = paste0("C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/ImpactedAreaCalculations/areas_by_threshold/ImpactedArea_",flab,".shp"),
                     delete_dsn = TRUE, quiet = T)
  }
}  

#get protected/impacted/at risk area sums
shp_list_long <- list.files(path = "C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/ImpactedAreaCalculations/areas_by_threshold/", pattern = ".shp", full.names = T)
shp_list_short <- list.files(path = "C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/ImpactedAreaCalculations/areas_by_threshold/", pattern = ".shp", full.names = F)

area_sums.df <- data.frame()

for(s in 1:length(shp_list_long)){
  tmp.shp <- readOGR(shp_list_long[s], verbose = F)
  name <- colsplit(shp_list_short[s], pattern = "_", names = c("ImpArea", "Tax", "level", "AYRO", "byLevel"))
  counts <- as.data.frame(table(tmp.shp@data[,9]))
  area_sums.df <- rbind(area_sums.df, 
                        data.frame(Tax = rep(name$Tax, times = 4),
                                   level = rep(name$level, times = 4),
                                   cat = counts[,1],
                                   sum = counts[,2])) 
  print(s)
}

area_sums.df2 <- area_sums.df %>%
  group_by(Tax, level, cat) %>%
   summarise(n = sum(sum)) %>%
     mutate(percentage = n / sum(n))

area_sums.df2$cat <- factor(area_sums.df2$cat, levels = c("Protected or outside FF", "At Risk", "Impacted (below cut-off)", "Impacted (above cut-off)"))

imp_area_plot_byKDE_perc <- 
    ggplot(area_sums.df2, aes(x = level*100, y = percentage*100, fill = cat)) + 
           geom_area() + 
           scale_fill_manual(values = c("blue3", "darkred", "yellow",  "Orange"),
                             breaks = c("Protected or outside FF", "Impacted (above cut-off)", "At Risk",  "Impacted (below cut-off)"),
                             labels = c("Protected", "Impacted (within FF, above effort cut-off)",
                                        "At Risk (within FF, no fishing)",  "Impacted (within FF, below effort cut-off)")) +
           scale_x_continuous(breaks = seq(0,100,20)) +
           scale_y_continuous(breaks = seq(0,100,20)) +
           facet_wrap(~Tax, nrow = 2) +
           theme_bw() + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
           labs(x = "Cumulative biomass cut-off (%)", y = "Area within threat category (%)", fill = "Threat Category")
imp_area_plot_byKDE_perc
ggsave(imp_area_plot_byKDE_perc, filename = "impacted_area_by_threshold_by_kde.pdf", height = 12, width = 18, units = "cm", device = "pdf", dpi = 150)

imp_area_plot_byKDE_abs <- 
    ggplot(area_sums.df2, aes(x = level*100, y = n, fill = cat)) + 
           geom_area() + 
           scale_fill_manual(values = c("blue3", "darkred", "yellow",  "Orange"),
                             breaks = c("Protected or outside FF", "Impacted (above cut-off)", "At Risk",  "Impacted (below cut-off)"),
                             labels = c("Protected", "Impacted (within FF, above effort cut-off)",
                                        "At Risk (within FF, no fishing)",  "Impacted (within FF, below effort cut-off)")) +
           scale_x_continuous(breaks = seq(0,100,20)) +
           facet_wrap(~Tax, nrow = 2, scales = "free_y") +
           theme_bw() + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
           labs(x = "Cumulative biomass cut-off (%)", y = "Area within threat category (* 1km^2 cells)", fill = "Threat Category")
imp_area_plot_byKDE_abs
ggsave(imp_area_plot_byKDE_abs, filename = "impacted_area_by_threshold_by_kde_aboslute.pdf", height = 12, width = 18, units = "cm", device = "pdf", dpi = 150)

write.csv(area_sums.df2, "NAFO_ImpactedAreas-by-threshold.csv", row.names = F)

### plot fishing intensity by taxa and cutoff
kde_cols <- c("black", "dodgerblue1", "darkkhaki","darkolivegreen3", "hotpink3", "darkorange2", "navajowhite3")

sai_refpoints <- data.frame(Tax       = unique(cutoffs$Tax),
                            refpoint  = c(0.105, 0.533, 6.615, 0.094, 1.273, 0.455, 0.068), # fishing intensity equivalent
                            SAI_level = c(75, 80, 80, 75, 80, 65, 80)) # % loss cutoff

cutoffs2 <- plyr::join(cutoffs, sai_refpoints, by = "Tax")

cutoff_curves <- ggplot(cutoffs2, aes(x = level*100, y = Value)) + 
                        geom_segment(aes(y = refpoint, yend = refpoint, x = -Inf, xend = SAI_level), linewidth = 0.7, colour = "darkgrey") +
                        geom_segment(aes(x = SAI_level, xend = SAI_level, y = 5*10^-3, yend = refpoint ), linewidth = 0.7, colour = "darkgrey") +  
                        geom_line(aes(colour = Tax), show.legend = F, linewidth = 1.4) + 
                        scale_colour_manual(values = kde_cols) +
                        scale_x_continuous(breaks=seq(0,100,20)) +
                        scale_y_log10(labels = scales::number_format(accuracy = 0.01)) +
                        facet_wrap(~Tax, scales = "free_y", nrow = 2) +
                        theme_bw() + 
                        labs(x = "Cumulative Biomass (%)", y = bquote('Log. Fishing Intensity (km/'*km^2*'/yr)'))
  
cutoff_curves
ggsave(cutoff_curves, filename = "cutoff_curves.pdf", height = 10, width = 20, units = "cm", device = "pdf", dpi = 150)

cutoff_curves_inv <- ggplot(cutoffs, aes(x = 100-(level*100), y = Value)) + 
  geom_line(aes(colour = Tax), show.legend = F, size = 1.4) + 
  facet_wrap(~Tax, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = kde_cols) +
  #scale_x_continuous() +
  scale_x_reverse(breaks=seq(0,100,20)) +
  scale_y_log10(labels = scales::number_format(accuracy = 0.01)) +
  theme_bw() + labs(x = "Biomass remaining (%)", y = bquote('Log. Fishing Intensity (km/'*km^2*'/yr)'))

cutoff_curves_inv
ggsave(cutoff_curves_inv, filename = "cutoff_curves_inverted.pdf", height = 10, width = 20, units = "cm", device = "pdf", dpi = 150)

save.image(paste0("NAFO_impacted_areas", Sys.Date(), ".RData"))

##### plot map by status, taxon and protection level #####

conts <- readOGR("C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/shp/NAFO_GEBCO_Contour.shp")
nafodivs <- readOGR("C:/Users/jb26/OneDrive - CEFAS/New laptop/Projects/NAFO NEREIDA/shp/NAFO_Divisions.shp")

conts.f <- fortify(conts)
divs.f <- fortify(nafodivs)
#in loop
for (t in 1:length(taxlist)){
    for (u in 1:length(thresholdlist)){

  shp <- readOGR(paste0("areas_by_threshold/ImpactedArea_", taxlist[t], "_", thresholdlist[u], "_AYRO_byLevel.shp"))
  shp <- spTransform(shp, proj4string(conts))
  colnames(shp@data)[9] <- "status"
  shp.f <- fortify(shp)
  shp@data$id <- rownames(shp@data)
  shp.f <- plyr::join(shp.f, shp@data, by = "id")
  plot <- ggplot() + geom_path(data = conts.f, aes(x = long, y = lat, group = group)) +
                     #geom_path(data = divs.f, aes(x = long, y = lat, group = group)) +
                     geom_polygon(data = shp.f, aes(x = long, y = lat, fill = status, group = group)) +
                     scale_fill_manual(values = c("blue3", "darkred", "yellow",  "Orange"),
                                       breaks = c("Protected or outside FF", "Impacted (above cut-off)", "At Risk",  "Impacted (below cut-off)"),
                                       labels = c("Protected", "Impacted (within FF, above cut-off)",
                                                  "At Risk (within FF, no fishing)",  "Impacted (within FF, below cut-off)")) +
                     labs(fill = "Protection\nStatus", x = "Longitude", y = "Latitude", title = paste0(taxlist[t],": ", thresholdlist[u]))+ 
                     theme_bw() + theme(legend.position = "bottom") + guides(fill=guide_legend(nrow=2, byrow=TRUE))
  plot
  
  dev.new()
  ggsave(plot, filename = paste0("threshold_maps/", taxlist[t], "_", thresholdlist[u], "_AYRO_byLevel.pdf"), device = "pdf", width = 6, height = 6, units = "in")
  dev.off()
  
  dev.new() #flemish cap cut out
  ggsave(plot + coord_cartesian(xlim = c(-48,-42.5), ylim = c(45,49.5)), 
         filename = paste0("threshold_maps/flemcap/", taxlist[t], "_", thresholdlist[u], "_AYRO_byLevel_FlemishCap.pdf"), 
         device = "pdf", width = 6, height = 6, units = "in")
  dev.off()
  
  dev.new() #tail of banks cut out
  ggsave(plot + coord_cartesian(xlim = c(-51,-48), ylim = c(42,45)), 
         filename = paste0("threshold_maps/banktail/", taxlist[t], "_", thresholdlist[u], "_AYRO_byLevel_BankTail.pdf"),
         device = "pdf", width = 6, height = 6, units = "in")
  dev.off()  
  
  print(paste0(taxlist[t], ": ", thresholdlist[u]))
  
}}#out loop

#####     2. Get Taxon Biomasses        #####################################################

# New field ames for the output
fldnames <- paste(taxlist,'Kg',sep ="")

### Make base output feature
bmf <- eff[,1]

for (i in 1:length(fc_list)){
  
  bmss <- st_read(fgdb, 
                 layer = fc_list[i])
  bmsp <- st_centroid(bmss)[,1]
  names(bmsp)[1] <- fldnames[i]
  bmf <- st_join(bmf, bmsp,join=st_intersects)
  
  
}

# Check data
summary(bmf)


st_write(bmf,"SAI/AllYearsAverageBoundaryEdit/BiomassAllVME_010321.shp",delete_dsn = TRUE)



#####   3. Calculate Areas and Biomasses     #################################################

### Combine all required data into one table for calculating statistics

## First select which one of the biomass curve methods to use
nl <- names(eff)
fsel <- c('OBJECTID', nl[grep(mlist[1],nl)])

# Join the Impact, KDE and biomass tables together
adata <-left_join(st_drop_geometry(eff[,fsel]), st_drop_geometry(kde[,c(1,3:9,11)]),by='OBJECTID')
adata <- left_join(adata,st_drop_geometry(bmf))

# Check for Null values
adata[!complete.cases(adata),]
# None of the rows with some biomass data in the incomplete rows are inside their respective KDE
# Removing incomplete rows
adata <- adata[complete.cases(adata),]


##

eff_selected_cutoffs <- eff[,c(1:9, 25, 47, 68, 88, 110, 128, 152)]

cell_status_biomass.df <- dplyr::full_join(adata, eff_selected_cutoffs, by = "OBJECTID")

biomass_per_status <- data.frame()

for(a in 1:length(taxlist)){
      tmp.biom <- aggregate(cell_status_biomass.df[,a+9]~
                            cell_status_biomass.df[,a+24], FUN = "sum")
      names(tmp.biom) <- c("SAI_status", "Biomass_kg")
      tmp.biom$Biomass_perc <- (tmp.biom$Biomass_kg/sum(tmp.biom$Biomass_kg))*100
      tmp.biom <- bind_rows(tmp.biom,
                            data.frame(SAI_status = "protected_all",
                                       Biomass_kg = sum(tmp.biom$Biomass_kg[4:6]),
                                       Biomass_perc = sum(tmp.biom$Biomass_perc[4:6])))
      tmp.biom$Tax <- taxlist[a]
      
      biomass_per_status <- bind_rows(biomass_per_status, tmp.biom)
}
write.csv(biomass_per_status, file = "biomass_per_cell_status.csv", row.names = F)
## Object to copy values into
output <- NULL

## Cycle through taxa and calculate area and biomass by impact and KDE membership

for (i in 1:length(rens)) {
  
  selnames <- grep(rens[i],names(adata),value = TRUE)
  selnames <- c(selnames, grep(paste(rens[i],'.*KDE',sep = ''),names(adata),value = TRUE))
  selnames <- c(selnames, grep(paste(rens[i],'.*Kg',sep = ''),names(adata),value = TRUE))
  
  dt <- as.data.table(adata %>%
                        group_by(get(selnames[2]),get(selnames[1])) %>%
                        summarise(Area = sum(Shape_Area)/1000000,
                                  Biomass = sum(get(selnames[3]))))
  dt[,Tax := taxlist[i]]
  dt[,Data := 'AYRO']
  
  output <- rbind(output,dt)
  
  }
  
  
names(output)[1:2] <- c('KDE','Biomass')

## Add closure type
output[Name %in% oldcs, ClosureType := 'Existing']
output[Name %in% newcs, ClosureType := 'New']
output[Name =='Open', ClosureType := 'Open']

### Make output into tables
require(dplyr)

# SAI Area
TblA <- dcast(output[KDE==1,],Impact~Tax,value.var = 'Area')
cols <- names(TblA)[-1]
TblA[,(cols) := round(.SD,0), .SDcols=cols]
TblA
TblA[is.na(TblA)] <- 0
TblA
out_cols = paste("PCT", cols, sep = ".")
TblA[,c(out_cols) := lapply(.SD, function(x) round(x/sum(x)*100,1)),.SDcols=cols]
TblA
setcolorder(TblA, c(1,2,9,3,10,4,11,5,12,6,13,7,14,8,15))

# SAI Biomass
TblB <- dcast(output[KDE==1,],Impact~Tax,value.var = 'Biomass')
cols <- names(TblB)[-1]
TblB[,(cols) := round(.SD,0), .SDcols=cols]
TblB
TblB[is.na(TblB)] <- 0
TblB
out_cols = paste("PCT", cols, sep = ".")
TblB[,c(out_cols) := lapply(.SD, function(x) round(x/sum(x)*100,1)),.SDcols=cols]
TblB
setcolorder(TblB, c(1,2,9,3,10,4,11,5,12,6,13,7,14,8,15))

# Management options Area
TblTemp <- dcast(output[KDE==1,],Tax~Impact,value.var = 'Area')
TblTemp
TblTemp[is.na(TblTemp)] <- 0
TblTemp

cols <- names(TblTemp)[-1]
TblTemp[,VMEtot:=rowSums(.SD),.SDcols=cols]
TblTemp

TblMOA <- TblTemp[,c(1,7)]
TblMOA

TblMOA[,Closed := TblTemp[,rowSums(.SD), .SDcols=4:5]]
TblMOA[,Conditional := TblTemp[,6]]
TblMOA[,ProtectOvrl := TblTemp[,rowSums(.SD), .SDcols=4:6]]
TblMOA[,Unprotected := TblTemp[,rowSums(.SD), .SDcols=2:3]]

cols <- names(TblMOA)[-1]
TblMOA[,(cols) := round(.SD,2), .SDcols=cols]

# Management options Biomass

TblTemp <- dcast(output[KDE==1,],Tax~Impact,value.var = 'Biomass')
TblTemp
TblTemp[is.na(TblTemp)] <- 0
TblTemp

cols <- names(TblTemp)[-1]
TblTemp[,VMEtot:=rowSums(.SD),.SDcols=cols]
TblTemp

TblMOB <- TblTemp[,c(1,7)]
TblMOB

TblMOB[,Closed := TblTemp[,rowSums(.SD), .SDcols=4:5]]
TblMOB[,Conditional := TblTemp[,6]]
TblMOB[,ProtectOvrl := TblTemp[,rowSums(.SD), .SDcols=4:6]]
TblMOB[,Unprotected := TblTemp[,rowSums(.SD), .SDcols=2:3]]

cols <- names(TblMOB)[-1]
TblMOB[,(cols) := round(.SD,2), .SDcols=cols]

### Write tables into Excel

require(openxlsx)

wb = createWorkbook()

addWorksheet(wb, "GridAnalysisOutput")
addWorksheet(wb, "TblAreaSAI")
addWorksheet(wb, "TblBiomassSAI")
addWorksheet(wb, "TblAreaMO")
addWorksheet(wb, "TblBiomassMO")


writeData(wb,output, sheet="GridAnalysisOutput")
writeData(wb,TblA, sheet="TblAreaSAI")
writeData(wb,TblB, sheet="TblBiomassSAI")
writeData(wb,TblMOA, sheet="TblAreaMO")
writeData(wb,TblMOB, sheet="TblBiomassMO")

saveWorkbook(wb, 'SAI/AllYearsAverageBoundaryEdit/GridAnalysis_AYRO_050321.xlsx')

### Write table into Excel
openxlsx::write.xlsx(output,file = 'SAI/AllYearsAverageBoundaryEdit/GridAnalysis_AYRO_010321.xlsx')

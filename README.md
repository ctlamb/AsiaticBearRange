Asiatic Bear Range Simulations
================
Clayton Lamb
03 January, 2022

## LOAD PACKAGES AND DATA

``` r
##load packages
library(raster)
library(mapview)
library(here)
library(lwgeom)
library(units)
library(sf)
library(fasterize)
library(stars)
library(hrbrthemes)
library(kableExtra)
library(rnaturalearth)
library(rnaturalearthhires)
library(tidyverse)

##define fn
st_erase=function(x, y) st_difference(x, st_union(st_combine(y)))

##load data
bear <- list()
iucn <- list()
points <- list()
hab <- list()
check <- list()
hr <- list()
ranges <- list()
grids <- list()

##PANDAS
bear[[1]] <- "GiantPanda"

iucn[[1]] <- st_read(here::here("data", "Pandas Final", "IUCN Map", "BearRange_Zone47.shp"))%>%
  st_transform(crs=32648)%>%
  group_by()%>%
  summarize()

points[[1]] <-st_read(here::here("data","Pandas Final", "Points", "3rd_census_merge.shp"))%>%
  mutate(grp="Census3rd")%>%
  select(grp)%>%
  rbind(
    st_read(here::here("data","Pandas Final", "Points", "3rd_National_Panda_Survey_Points_Sichuan_Province.shp"))%>%
      mutate(grp="Census3rd_Sichuan")%>%
      select(grp)
  )%>%
  cbind(st_coordinates(.))%>%
  filter(X>0)%>%
  st_transform(crs=32648)

hab[[1]] <- st_read(here::here("data","Pandas Final", "Suitable Habitat","Panda_habitat.shp"))%>%
  st_transform(crs=32648)


check[[1]] <- st_read(here::here("data", "Pandas Final", "Modelled Ranges", "GP", "Buffer_Export.shp"))%>%
  st_transform(crs=32648)

hr[[1]] <- 10*1E6

###ABB

bear[[2]] <- "Asian Black Bear"


iucn[[2]] <- st_read(here::here("data", "ABB Final", "IUCN Map", "BearRange_Zone47.shp"))%>%
  st_transform(crs=32648)%>%
  st_simplify(dTolerance=100)%>%
  group_by()%>%
  summarize()

points[[2]] <-st_read(here::here("data","ABB Final", "Points", "ABB_New_Merged.shp"))%>%
  mutate(grp="New")%>%
  select(grp)%>%
  st_zm()%>%
  rbind(
    st_read(here::here("data","ABB Final", "Points", "ABB_Old_Merged.shp"))%>%
      mutate(grp="Old")%>%
      select(grp)
  )%>%
  cbind(st_coordinates(.))%>%
  filter(X>0)%>%
  st_transform(crs=32648)


hab[[2]] <- raster(here::here("data","ABB Final", "Suitable Habitat","Suitable_Habitiat_Boundary.tif"))
values(hab[[2]])[values(hab[[2]])==0]<-NA
hab[[2]] <- st_as_stars(hab[[2]])
hab[[2]]=st_as_sf(hab[[2]][1], as_points = FALSE, merge = TRUE)%>%
  st_transform(crs=32648)


check[[2]] <- st_read(here::here("data", "ABB Final", "Modelled Ranges", "ABB", "Buffer_Merge_Dissolve.shp"))%>%
  st_transform(crs=32648)

hr[[2]] <- 50*1E6


##PLOT
# mapview(iucn[[1]])+
#   mapview(hab[[1]])+
#   mapview(points[[1]])
# 
# mapview(iucn[[2]])+
#   mapview(hab[[2]])+
#   mapview(points[[2]])
# 
# 
# ###ABB Sichuan only
sich <- rnaturalearth::ne_states("china")%>%
  st_as_sf()%>%
  filter(name%in%"Sichuan")%>%
  st_transform(crs=32648)

# 
# mapview(sich)



bear[[3]] <- "Asian Black Bear-Sichuan Only"


iucn[[3]] <- iucn[[2]]%>%
  st_intersection(sich)
#mapview(sich, alpha=0.3)+mapview(iucn[[3]])

points[[3]] <- points[[2]]%>%
  st_intersection(sich)
#mapview(sich, alpha=0.3)+mapview(points[[3]])

hab[[3]] <- hab[[2]]%>%st_make_valid()%>%
  st_intersection(sich)
#mapview(sich, alpha=0.3)+mapview(hab[[3]])


check[[3]] <- check[[2]]%>%
  st_intersection(sich)

hr[[3]] <- 50*1E6
```

## LOOP TO ASSESS THE MULTIPLICATION FACTOR

``` r
df <- data.frame()
for(b in 1:3){
factor.x <- c(5,10,30,50,75,150)
for(i in 1:length(factor.x)){
  ##make grid
  # We use our species-specific unit of measure (Table 1), to establish a grid within which we evaluate each cell for historic range (range at last assessment), suitable habitat, and number of credible observations.  
  # Cells outside of the historical range and lacking suitable habitat will not be considered without 2 or more credible observations. 
  # The remaining cells either contain suitable habitat or are within the historical range and may, or may not, contain credible observations.  
  factor <- factor.x[i]
  grid.width <- (sqrt(hr[[b]]*factor))%>%round(0)
  aoi.radius <- sqrt((hr[[b]]*factor)/3.14)
  
  grid <- st_make_grid(iucn[[b]]%>%st_buffer(300000),
                         cellsize=grid.width)%>%
    st_sf()
  
  grid$pt_count <- lengths(st_intersects(grid, points[[b]]))
  grid$hist <- lengths(st_intersects(grid, iucn[[b]]))
  grid$hab <- lengths(st_intersects(grid, hab[[b]]))
  
  grid <- grid%>%
    mutate(keep=case_when(pt_count>1~"Y",
                          hist>0~"Y",
                          hab>0~"Y",
                          TRUE~"N"))%>%
    filter(keep%in%"Y")
  

  # mapview(iucn[[b]])+mapview(grid,alpha.regions=0)
  
  ###update range using Bill's rules from the MS
  ##Add to range rules
  #R1 All suitable habitat within cells with 2 or more observations is considered occupied
  grid.r1 <- grid%>%filter(pt_count>=2)%>%st_intersection(hab[[b]]%>%st_make_valid())
  # mapview(p.grid.r1) + mapview(p.grid%>%filter(hab>=2))
  
  #R2 Internal cells that lack observations during the required period will be considered part of range, as long as they overlap with suitable habitat layer and have an observation within one adjoining cell.  (this is is due to the 'area of influence extending into this cell) 
  aoi <- points[[b]]%>%st_buffer(aoi.radius)
  #mapview(p.aoi)
  
  grid$touches <- lengths(st_is_within_distance(grid%>%st_centroid(),dist=grid.width*1.05))
  # mapview(grid,alpha.regions=0)+mapview(grid%>%filter(touches>=5))
  
  grid.r2 <-  grid%>%filter(touches>=5)
  grid.r2$hasaoi <- lengths(st_intersects(grid.r2, aoi))
  grid.r2 <-  grid.r2%>%filter(hasaoi>0)
  grid.r2$hab <- lengths(st_intersects(grid.r2, hab[[b]]))
  grid.r2 <-  grid.r2%>%filter(hab>0)
  
  grid.r2.a <- points[[b]]%>%st_intersection(grid.r2)
  grid.r2.a <- grid.r2.a%>%st_buffer(aoi.radius)
  
  grid.r2.b <- grid.r2%>%st_intersection(hab[[b]]%>%st_make_valid())
  
  # mapview(grid,alpha.regions=0) + mapview(grid.r2.a)
  # mapview(grid,alpha.regions=0) + mapview(grid.r2.b)
  
  
  #R3 All edge cells with observations will have their area of influence extended beyond the cell boundary regardless of suitable habitat status. 
  grid.r3 <-  grid%>%filter(touches<5)
  # mapview(grid,alpha.regions=0) + mapview(grid.r3)
  
  grid.r3 <- st_intersection(grid.r3, points[[b]])
  # mapview(grid,alpha.regions=0) + mapview(grid.r3)
  
  grid.r3 <- grid.r3%>%st_buffer(aoi.radius)
  # mapview(grid,alpha.regions=0) + mapview(grid.r3)
  
  grid.r3 <- grid.r3%>%st_erase(grid%>%group_by()%>%summarize()%>%st_make_valid())
  # mapview(grid,alpha.regions=0) + mapview(grid.r3)
  
  #All edge cells without observations rely on adjoining cells for their status (as with internal cells). 
  grid.r4 <-  grid%>%filter(touches<5)
  grid.r4$touches <- lengths(st_intersects(grid.r4, aoi))
  grid.r4 <- grid.r4%>%filter(touches>0)
  
  grid.r4.a <-  points[[b]]%>%st_intersection(grid.r4)
  grid.r4.a  <- grid.r4.a %>%st_buffer(aoi.radius)
  
  grid.r4.b <- grid.r4%>%st_intersection(hab[[b]]%>%st_make_valid())
  
  # mapview(grid,alpha.regions=0) + mapview(grid.r4.a)
  # mapview(grid,alpha.regions=0) + mapview(grid.r4.b)
  
  ##clean up
  grid.r2.a <- grid.r2.a%>%group_by()%>%summarize()
  grid.r3 <- grid.r3%>%group_by()%>%summarize()
  grid.r4.a <- grid.r4.a%>%group_by()%>%summarize()
  
  #All disjunct patches of suitable habitat lacking observations will be considered outside the current range.
  #All internal areas with intermediate suitable habitat and lacking observations with be subject to expert opinion.  Secondary inquiries and field investigations should be conducted at all suitable habitat along range edge and within potential disjunct populations prior to final decisions. 
  
  ###merge together
  range <- rbind(grid.r1%>%select(pt_count),
                 grid.r2.a%>%mutate(pt_count=1)%>%select(pt_count), 
                 grid.r2.b%>%mutate(pt_count=1)%>%select(pt_count), 
                 grid.r3%>%mutate(pt_count=1)%>%select(pt_count),
                 grid.r4.a%>%mutate(pt_count=1)%>%select(pt_count),
                 grid.r4.b%>%mutate(pt_count=1)%>%select(pt_count))%>%
    group_by()%>%
    summarise()
  
  #mapview(range)
  ranges[[((b-1)*length(factor.x))+i]]<-range
  grids[[((b-1)*length(factor.x))+i]]<-grid
  
  # mapview(range) + mapview(check[[b]])
  # mapview(range) + mapview(p.iucn)
  # mapview(range)+mapview(p.points)
  
  ###summary stats
  iucn.area <- st_area(iucn[[b]])%>%set_units(km^2)
  range.area <- st_area(range)%>%set_units(km^2)
  new.area <- st_erase(range,iucn[[b]])%>%st_area()%>%set_units(km^2)
  lost.area <- st_erase(iucn[[b]],range)%>%st_area()%>%set_units(km^2)
  overlap.area.p <- ((st_intersection(iucn[[b]],range)%>%st_area()%>%set_units(km^2))/iucn.area)*100
  extra.area.p <- ((st_erase(range,iucn[[b]])%>%st_area()%>%set_units(km^2))/iucn.area)*100
  
  df<- rbind(df,
             data.frame(sp=bear[[b]],
                        factor=factor.x[i],
                        iucn.area=iucn.area%>%as.numeric()%>%round(0),
                        range.area=range.area%>%as.numeric()%>%round(0),
                        new.area=new.area%>%as.numeric()%>%round(0),
                        lost.area=lost.area%>%as.numeric()%>%round(0),
                        overlap.area.p=overlap.area.p%>%as.numeric()%>%round(0),
                        extra.area.p=extra.area.p%>%as.numeric()%>%round(0)))
  #print(factor.x[i])
}
}

kable(df)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

sp

</th>

<th style="text-align:right;">

factor

</th>

<th style="text-align:right;">

iucn.area

</th>

<th style="text-align:right;">

range.area

</th>

<th style="text-align:right;">

new.area

</th>

<th style="text-align:right;">

lost.area

</th>

<th style="text-align:right;">

overlap.area.p

</th>

<th style="text-align:right;">

extra.area.p

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

27457

</td>

<td style="text-align:right;">

5179

</td>

<td style="text-align:right;">

12327

</td>

<td style="text-align:right;">

64

</td>

<td style="text-align:right;">

15

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

33688

</td>

<td style="text-align:right;">

8684

</td>

<td style="text-align:right;">

9601

</td>

<td style="text-align:right;">

72

</td>

<td style="text-align:right;">

25

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

47381

</td>

<td style="text-align:right;">

18445

</td>

<td style="text-align:right;">

5669

</td>

<td style="text-align:right;">

84

</td>

<td style="text-align:right;">

53

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

56402

</td>

<td style="text-align:right;">

25901

</td>

<td style="text-align:right;">

4105

</td>

<td style="text-align:right;">

88

</td>

<td style="text-align:right;">

75

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

75

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

65098

</td>

<td style="text-align:right;">

33619

</td>

<td style="text-align:right;">

3126

</td>

<td style="text-align:right;">

91

</td>

<td style="text-align:right;">

97

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

83642

</td>

<td style="text-align:right;">

50687

</td>

<td style="text-align:right;">

1650

</td>

<td style="text-align:right;">

95

</td>

<td style="text-align:right;">

146

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

534574

</td>

<td style="text-align:right;">

68288

</td>

<td style="text-align:right;">

4386374

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

755366

</td>

<td style="text-align:right;">

108377

</td>

<td style="text-align:right;">

4205669

</td>

<td style="text-align:right;">

13

</td>

<td style="text-align:right;">

2

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

1281991

</td>

<td style="text-align:right;">

240128

</td>

<td style="text-align:right;">

3810796

</td>

<td style="text-align:right;">

21

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

1650709

</td>

<td style="text-align:right;">

345780

</td>

<td style="text-align:right;">

3547729

</td>

<td style="text-align:right;">

27

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

75

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

1997194

</td>

<td style="text-align:right;">

457401

</td>

<td style="text-align:right;">

3312865

</td>

<td style="text-align:right;">

32

</td>

<td style="text-align:right;">

9

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

2716659

</td>

<td style="text-align:right;">

718970

</td>

<td style="text-align:right;">

2854970

</td>

<td style="text-align:right;">

41

</td>

<td style="text-align:right;">

15

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear-Sichuan Only

</td>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

327584

</td>

<td style="text-align:right;">

186473

</td>

<td style="text-align:right;">

14545

</td>

<td style="text-align:right;">

155657

</td>

<td style="text-align:right;">

52

</td>

<td style="text-align:right;">

4

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear-Sichuan Only

</td>

<td style="text-align:right;">

10

</td>

<td style="text-align:right;">

327584

</td>

<td style="text-align:right;">

232053

</td>

<td style="text-align:right;">

24552

</td>

<td style="text-align:right;">

120083

</td>

<td style="text-align:right;">

63

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear-Sichuan Only

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:right;">

327584

</td>

<td style="text-align:right;">

319108

</td>

<td style="text-align:right;">

56260

</td>

<td style="text-align:right;">

64736

</td>

<td style="text-align:right;">

80

</td>

<td style="text-align:right;">

17

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear-Sichuan Only

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

327584

</td>

<td style="text-align:right;">

365009

</td>

<td style="text-align:right;">

81956

</td>

<td style="text-align:right;">

44532

</td>

<td style="text-align:right;">

86

</td>

<td style="text-align:right;">

25

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear-Sichuan Only

</td>

<td style="text-align:right;">

75

</td>

<td style="text-align:right;">

327584

</td>

<td style="text-align:right;">

410034

</td>

<td style="text-align:right;">

111554

</td>

<td style="text-align:right;">

29104

</td>

<td style="text-align:right;">

91

</td>

<td style="text-align:right;">

34

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear-Sichuan Only

</td>

<td style="text-align:right;">

150

</td>

<td style="text-align:right;">

327584

</td>

<td style="text-align:right;">

481466

</td>

<td style="text-align:right;">

171388

</td>

<td style="text-align:right;">

17507

</td>

<td style="text-align:right;">

95

</td>

<td style="text-align:right;">

52

</td>

</tr>

</tbody>

</table>

## PLOT RESULTS OF MULTIPLICATION FACTOR

``` r
df%>%
  pivot_longer(c(-factor,-sp))%>%
  filter(name%in%c("overlap.area.p","extra.area.p"))%>%
  mutate(sp=as.character(sp))%>%
  mutate(sp=case_when(sp%in%"Asian Black Bear-Sichuan Only"~"Asian Black Bear-Sichuan",
                      TRUE~sp))%>%
  mutate(sp=fct_relevel(sp,"GiantPanda", "Asian Black Bear", "Asian Black Bear-Sichuan"))%>%
  ggplot(aes(x=factor,y=value, color=name))+
  geom_line()+
  geom_point()+
  theme_bw()+
  labs(x="Multiplication Factor (x)",
       y="Percent of IUCN Range",
       title="How Sensitive is the Range to Mulitplication Factor?",
       subtitle="Range Estimates Compared to IUCN")+
  facet_wrap(vars(sp))+
  geom_hline(yintercept = 100, linetype="dashed")+
  theme_ipsum()
```

![](README_files/figure-gfm/PLOT%20RESULTS%20OF%20MULTIPLICATION%20FACTOR-1.png)<!-- -->

``` r
ggsave(here::here("plots", "mult_factor.png"), height=4, width=10)
```

## LOOP TO ASSESS ABILITY TO DETECT RANGE CHANGES

``` r
set.seed(2020)
df.change <- data.frame()

for(sim.count in 1:50){
# factor.x <- c(5,10,30,75)
 factor.x <- c(50)
p.change.y <- c(1,5,10,20,50)
#ranges.sub <- ranges[c(1:3,5,7:9,11,13:15,17)]
ranges.sub <- ranges[c(4,10,16)]
range.lost <-list()
count<-0

    for(b in 1:3){
      
      for(i in 1:length(factor.x)){
        range.orig <- ranges.sub[[((b-1)*length(factor.x))+i]]
        
        for(y in 1:length(p.change.y)){
          count <- count+1
          p.change <- 1-(p.change.y[y]/100)
          
        ##make grid
        # We use our species-specific unit of measure (Table 1), to establish a grid within which we evaluate each cell for historic range (range at last assessment), suitable habitat, and number of credible observations.  
        # Cells outside of the historical range and lacking suitable habitat will not be considered without 2 or more credible observations. 
        # The remaining cells either contain suitable habitat or are within the historical range and may, or may not, contain credible observations.  
          if(b==1){
          iucn.r.large <- raster(range.orig, res = 15000)
          r.large.change <- fasterize(range.orig, iucn.r.large)
          #plot(r.large)
          }
          
          
          if(b%in%c(2,3)){
            iucn.r.large <- raster(range.orig, res = 30000)
            r.large.change <- fasterize(range.orig, iucn.r.large)
            #plot(r.large.change)
          }
          

          values(r.large.change)[!is.na(values(r.large.change))] <- rbinom(length(values(r.large.change)[!is.na(values(r.large.change))]), 1, p.change)
          values(r.large.change)[values(r.large.change)==1] <- NA
          values(r.large.change)[values(r.large.change)==0] <- 1
          #plot(r.large.change)
          
          
          
          r.large.change <- st_as_stars(r.large.change)
          r.large.change =st_as_sf(r.large.change, as_points = FALSE, merge = TRUE)%>%
            st_transform(crs=32648)%>%
            st_make_valid()
          
         #mapview(range.orig)+mapview(r.large.change)
          
          loss.poly <- st_erase(range.orig, r.large.change)
          
          loss.per <- (((st_area(range.orig)-st_area(loss.poly))/st_area(range.orig))*100)%>%round(0)%>%as.numeric()
          
          points.change <- st_erase(points[[b]],r.large.change)%>%
            st_make_valid()
          
          
          #mapview(range.orig) + mapview(points.change)
          
          loss.per.points <- (((nrow(points[[b]])-nrow(points.change))/nrow(points[[b]]))*100)%>%round(0)
        
          
          
        factor <- factor.x[i]
        grid.width <- (sqrt(hr[[b]]*factor))%>%round(0)
        aoi.radius <- sqrt((hr[[b]]*factor)/3.14)
        
        grid <- st_make_grid(iucn[[b]]%>%st_buffer(300000),
                             cellsize=grid.width)%>%
                st_sf()
        
        grid$pt_count <- lengths(st_intersects(grid, points.change))
        grid$hist <- lengths(st_intersects(grid, iucn[[b]]))
        grid$hab <- lengths(st_intersects(grid, hab[[b]]))
        
        grid <- grid%>%
          mutate(keep=case_when(pt_count>1~"Y",
                                hist>0~"Y",
                                hab>0~"Y",
                                TRUE~"N"))%>%
          filter(keep%in%"Y")
        
        
        # mapview(iucn[[b]])+mapview(grid,alpha.regions=0)
        
        ###update range using Bill's rules from the MS
        ##Add to range rules
        #R1 All suitable habitat within cells with 2 or more observations is considered occupied
        grid.r1 <- grid%>%filter(pt_count>=2)%>%st_intersection(hab[[b]]%>%st_make_valid())
        # mapview(p.grid.r1) + mapview(p.grid%>%filter(hab>=2))
        
        #R2 Internal cells that lack observations during the required period will be considered part of range, as long as they overlap with suitable habitat layer and have an observation within one adjoining cell.  (this is is due to the 'area of influence extending into this cell) 
        aoi <- points.change%>%st_buffer(aoi.radius)
        #mapview(p.aoi)
        
        grid$touches <- lengths(st_is_within_distance(grid%>%st_centroid(),dist=grid.width*1.05))
        # mapview(grid,alpha.regions=0)+mapview(grid%>%filter(touches>=5))
        
        grid.r2 <-  grid%>%filter(touches>=5)
        grid.r2$hasaoi <- lengths(st_intersects(grid.r2, aoi))
        grid.r2 <-  grid.r2%>%filter(hasaoi>0)
        grid.r2$hab <- lengths(st_intersects(grid.r2, hab[[b]]))
        grid.r2 <-  grid.r2%>%filter(hab>0)
        
        grid.r2.a <- points.change%>%st_intersection(grid.r2)
        grid.r2.a <- grid.r2.a%>%st_buffer(aoi.radius)
        
        grid.r2.b <- grid.r2%>%st_intersection(hab[[b]]%>%st_make_valid())
        
        # mapview(grid,alpha.regions=0) + mapview(grid.r2.a)
        # mapview(grid,alpha.regions=0) + mapview(grid.r2.b)
        
        
        #R3 All edge cells with observations will have their area of influence extended beyond the cell boundary regardless of suitable habitat status. 
        grid.r3 <-  grid%>%filter(touches<5)
        # mapview(grid,alpha.regions=0) + mapview(grid.r3)
        
        grid.r3 <- st_intersection(grid.r3, points.change)
        # mapview(grid,alpha.regions=0) + mapview(grid.r3)
        
        grid.r3 <- grid.r3%>%st_buffer(aoi.radius)
        # mapview(grid,alpha.regions=0) + mapview(grid.r3)
        
        grid.r3 <- grid.r3%>%st_erase(grid%>%group_by()%>%summarize()%>%st_make_valid())
        # mapview(grid,alpha.regions=0) + mapview(grid.r3)
        
        #All edge cells without observations rely on adjoining cells for their status (as with internal cells). 
        grid.r4 <-  grid%>%filter(touches<5)
        grid.r4$touches <- lengths(st_intersects(grid.r4, aoi))
        grid.r4 <- grid.r4%>%filter(touches>0)
        
        grid.r4.a <-  points.change%>%st_intersection(grid.r4)
        grid.r4.a  <- grid.r4.a %>%st_buffer(aoi.radius)
        
        grid.r4.b <- grid.r4%>%st_intersection(hab[[b]]%>%st_make_valid())
        
        # mapview(grid,alpha.regions=0) + mapview(grid.r4.a)
        # mapview(grid,alpha.regions=0) + mapview(grid.r4.b)
        
        ##clean up
        grid.r2.a <- grid.r2.a%>%group_by()%>%summarize()
        grid.r3 <- grid.r3%>%group_by()%>%summarize()
        grid.r4.a <- grid.r4.a%>%group_by()%>%summarize()
        
        #All disjunct patches of suitable habitat lacking observations will be considered outside the current range.
        #All internal areas with intermediate suitable habitat and lacking observations with be subject to expert opinion.  Secondary inquiries and field investigations should be conducted at all suitable habitat along range edge and within potential disjunct populations prior to final decisions. 
        
        ###merge together
        range <- rbind(grid.r1%>%select(pt_count),
                       grid.r2.a%>%mutate(pt_count=1)%>%select(pt_count), 
                       grid.r2.b%>%mutate(pt_count=1)%>%select(pt_count), 
                       grid.r3%>%mutate(pt_count=1)%>%select(pt_count),
                       grid.r4.a%>%mutate(pt_count=1)%>%select(pt_count),
                       grid.r4.b%>%mutate(pt_count=1)%>%select(pt_count))%>%
          group_by()%>%
          summarise()
        
        #mapview(range)

        
        # mapview(range) + mapview(check[[b]])
        # mapview(range) + mapview(p.iucn)
        # mapview(range)+mapview(p.points)
        
    
        ###summary stats
        range.orig.area <- st_area(range.orig)%>%set_units(km^2)
        range.area <- st_area(range)%>%set_units(km^2)
        area.lost <- ((range.orig.area-range.area)/range.orig.area)*100
        
        df.change <- rbind(df.change,
                           data.frame(sim=sim.count,
                                      sp=bear[[b]],
                                      factor=factor.x[i],
                                      change.loop=p.change.y[y],
                                      loss.percent.range=loss.per%>%abs(),
                                      loss.percent.points=loss.per.points%>%abs(),
                                      area.lost=area.lost%>%as.numeric()%>%round(0)))
        
        range.lost[[count]] <- range
        
        print(paste(bear[[b]], ", factor=",factor.x[i],"x", ", decline=",p.change.y[y],"%", sep=""))
      }
    }
    }
#print(sim.count)
}
```

    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"
    ## [1] "GiantPanda, factor=50x, decline=1%"
    ## [1] "GiantPanda, factor=50x, decline=5%"
    ## [1] "GiantPanda, factor=50x, decline=10%"
    ## [1] "GiantPanda, factor=50x, decline=20%"
    ## [1] "GiantPanda, factor=50x, decline=50%"
    ## [1] "Asian Black Bear, factor=50x, decline=1%"
    ## [1] "Asian Black Bear, factor=50x, decline=5%"
    ## [1] "Asian Black Bear, factor=50x, decline=10%"
    ## [1] "Asian Black Bear, factor=50x, decline=20%"
    ## [1] "Asian Black Bear, factor=50x, decline=50%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=1%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=5%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=10%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=20%"
    ## [1] "Asian Black Bear-Sichuan Only, factor=50x, decline=50%"

``` r
df.change$factor <- as.character(df.change$factor)%>%
  fct_reorder(as.numeric(.))

df.change$sp <- as.character(df.change$sp)

df.change <- df.change%>%
  mutate(sp=case_when(sp%in%"Asian Black Bear-Sichuan Only"~as.character("Asian Black Bear-Sichuan"),
                      TRUE~as.character(sp)))
df.change$sp <- fct_relevel(df.change$sp%>%as.factor(), "GiantPanda", "Asian Black Bear", "Asian Black Bear-Sichuan")
```

## PLOT ABILITY TO DETECT RANGE CHANGES

``` r
df.change%>%
  select(sp, sim,loss.percent.range,Points=loss.percent.points,Area=area.lost)%>%
  pivot_longer(-c(sp,sim,loss.percent.range))%>%
  arrange(sp,name,loss.percent.range,sim)%>%
  mutate(Approach=forcats::fct_rev(factor(name)))%>%
  ggplot(aes(x=loss.percent.range,y=value, color=Approach, group=sim))+
  geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha=0.5)+
  geom_path(alpha=0.15)+
  theme_ipsum()+
  labs(x="Actual area Lost (%)",
       y="Predicted Loss (%)",
       title="Detecting Simulated Range Losses")+
  coord_fixed(ratio = 1)+
  facet_wrap(vars(sp))+
  xlim(0,55)+
  ylim(0,55)+
  guides(colour = guide_legend(override.aes = list(alpha=1)))
```

![](README_files/figure-gfm/PLOT%20ABILITY%20TO%20DETECT%20RANGE%20CHANGES-1.png)<!-- -->

``` r
ggsave(here::here("plots", "range_loss_sim50.png"), height=4, width=10)

df.change%>%
  select(sp, sim,loss.percent.range,Points=loss.percent.points,Area=area.lost)%>%
  group_by(sp,loss.percent.range)%>%
  summarise(Points=median(Points),
            Area=median(Area))%>%
  write_csv(here::here("tables", "range_loss_sim50_summary.csv"))
```

## LOOP TO ASSESS HOW MUCH SAMPLING IS NEEDED FOR ABB

``` r
####LOOP TO ASSESS THE MULTIPLICATION FACTOR
df.samp.dens2 <- data.frame()
ranges.sd <- list()
for(b in 1:2){
  samp.dens <- c(100,500,1000,3000,5000,10000)
  for(i in 1:length(samp.dens)){
    
    points.i <- st_sample(iucn[[b]], size=samp.dens[i], type="random")%>%
                        st_as_sf()%>%
                        mutate(grp="rnd")%>%
                        cbind(st_coordinates(.))%>%
                        rename(geometry=x)

    
    
    
    
    ##make grid
    # We use our species-specific unit of measure (Table 1), to establish a grid within which we evaluate each cell for historic range (range at last assessment), suitable habitat, and number of credible observations.  
    # Cells outside of the historical range and lacking suitable habitat will not be considered without 2 or more credible observations. 
    # The remaining cells either contain suitable habitat or are within the historical range and may, or may not, contain credible observations.  
    factor <- 50
    grid.width <- (sqrt(hr[[b]]*factor))%>%round(0)
    aoi.radius <- sqrt((hr[[b]]*factor)/3.14)
    
    grid <- st_make_grid(iucn[[b]]%>%st_buffer(300000),
                         cellsize=grid.width)%>%
      st_sf()
    
    grid$pt_count <- lengths(st_intersects(grid, points.i))
    grid$hist <- lengths(st_intersects(grid, iucn[[b]]))
    grid$hab <- lengths(st_intersects(grid, hab[[b]]))
    
    grid <- grid%>%
      mutate(keep=case_when(pt_count>1~"Y",
                            hist>0~"Y",
                            hab>0~"Y",
                            TRUE~"N"))%>%
      filter(keep%in%"Y")
    
    
    # mapview(iucn[[b]])+mapview(grid,alpha.regions=0)
    
    ###update range using Bill's rules from the MS
    ##Add to range rules
    #R1 All suitable habitat within cells with 2 or more observations is considered occupied
    grid.r1 <- grid%>%filter(pt_count>=2)%>%st_intersection(hab[[b]]%>%st_make_valid())
    # mapview(p.grid.r1) + mapview(p.grid%>%filter(hab>=2))
    
    #R2 Internal cells that lack observations during the required period will be considered part of range, as long as they overlap with suitable habitat layer and have an observation within one adjoining cell.  (this is is due to the 'area of influence extending into this cell) 
    aoi <- points.i%>%st_buffer(aoi.radius)
    #mapview(p.aoi)
    
    grid$touches <- lengths(st_is_within_distance(grid%>%st_centroid(),dist=grid.width*1.05))
    # mapview(grid,alpha.regions=0)+mapview(grid%>%filter(touches>=5))
    
    grid.r2 <-  grid%>%filter(touches>=5)
    grid.r2$hasaoi <- lengths(st_intersects(grid.r2, aoi))
    grid.r2 <-  grid.r2%>%filter(hasaoi>0)
    grid.r2$hab <- lengths(st_intersects(grid.r2, hab[[b]]))
    grid.r2 <-  grid.r2%>%filter(hab>0)
    
    grid.r2.a <- points.i%>%st_intersection(grid.r2)
    grid.r2.a <- grid.r2.a%>%st_buffer(aoi.radius)
    
    grid.r2.b <- grid.r2%>%st_intersection(hab[[b]]%>%st_make_valid())
    
    # mapview(grid,alpha.regions=0) + mapview(grid.r2.a)
    # mapview(grid,alpha.regions=0) + mapview(grid.r2.b)
    
    
    #R3 All edge cells with observations will have their area of influence extended beyond the cell boundary regardless of suitable habitat status. 
    grid.r3 <-  grid%>%filter(touches<5)
    # mapview(grid,alpha.regions=0) + mapview(grid.r3)
    
    grid.r3 <- st_intersection(grid.r3, points.i)
    # mapview(grid,alpha.regions=0) + mapview(grid.r3)
    
    grid.r3 <- grid.r3%>%st_buffer(aoi.radius)
    # mapview(grid,alpha.regions=0) + mapview(grid.r3)
    
    grid.r3 <- grid.r3%>%st_erase(grid%>%group_by()%>%summarize()%>%st_make_valid())
    # mapview(grid,alpha.regions=0) + mapview(grid.r3)
    
    #All edge cells without observations rely on adjoining cells for their status (as with internal cells). 
    grid.r4 <-  grid%>%filter(touches<5)
    grid.r4$touches <- lengths(st_intersects(grid.r4, aoi))
    grid.r4 <- grid.r4%>%filter(touches>0)
    
    grid.r4.a <-  points.i%>%st_intersection(grid.r4)
    grid.r4.a  <- grid.r4.a %>%st_buffer(aoi.radius)
    
    grid.r4.b <- grid.r4%>%st_intersection(hab[[b]]%>%st_make_valid())
    
    # mapview(grid,alpha.regions=0) + mapview(grid.r4.a)
    # mapview(grid,alpha.regions=0) + mapview(grid.r4.b)
    
    ##clean up
    grid.r2.a <- grid.r2.a%>%group_by()%>%summarize()
    grid.r3 <- grid.r3%>%group_by()%>%summarize()
    grid.r4.a <- grid.r4.a%>%group_by()%>%summarize()
    
    #All disjunct patches of suitable habitat lacking observations will be considered outside the current range.
    #All internal areas with intermediate suitable habitat and lacking observations with be subject to expert opinion.  Secondary inquiries and field investigations should be conducted at all suitable habitat along range edge and within potential disjunct populations prior to final decisions. 
    
    ###merge together
    range <- rbind(grid.r1%>%select(pt_count),
                   grid.r2.a%>%mutate(pt_count=1)%>%select(pt_count), 
                   grid.r2.b%>%mutate(pt_count=1)%>%select(pt_count), 
                   grid.r3%>%mutate(pt_count=1)%>%select(pt_count),
                   grid.r4.a%>%mutate(pt_count=1)%>%select(pt_count),
                   grid.r4.b%>%mutate(pt_count=1)%>%select(pt_count))%>%
      group_by()%>%
      summarise()
    
    #mapview(range)
    ranges.sd[[i]]<-range
    
    
    # mapview(range) + mapview(check[[b]])
    # mapview(range) + mapview(iucn[[b]])
    # mapview(range)+mapview(p.points)
    
    ###summary stats
    iucn.area <- st_area(iucn[[b]])%>%set_units(km^2)
    range.area <- st_area(range)%>%set_units(km^2)
    new.area <- st_erase(range,iucn[[b]])%>%st_area()%>%set_units(km^2)
    lost.area <- st_erase(iucn[[b]],range)%>%st_area()%>%set_units(km^2)
    overlap.area.p <- ((st_intersection(iucn[[b]],range)%>%st_area()%>%set_units(km^2))/iucn.area)*100
    extra.area.p <- ((st_erase(range,iucn[[b]])%>%st_area()%>%set_units(km^2))/iucn.area)*100
    
    if(length(lost.area)==0){
      lost.area<-0
    }
    df.samp.dens2<- rbind(df.samp.dens2,
                         data.frame(sp=bear[[b]],
                                    samp.dens=samp.dens[i],
                                    factor=50,
                                    iucn.area=iucn.area%>%as.numeric()%>%round(0),
                                    range.area=range.area%>%as.numeric()%>%round(0),
                                    new.area=new.area%>%as.numeric()%>%round(0),
                                    lost.area=lost.area%>%as.numeric()%>%round(0),
                                    overlap.area.p=overlap.area.p%>%as.numeric()%>%round(0),
                                    extra.area.p=extra.area.p%>%as.numeric()%>%round(0)))
    print(samp.dens[i])
  }
}
```

    ## [1] 100
    ## [1] 500
    ## [1] 1000
    ## [1] 3000
    ## [1] 5000
    ## [1] 10000
    ## [1] 100
    ## [1] 500
    ## [1] 1000
    ## [1] 3000
    ## [1] 5000
    ## [1] 10000

``` r
kable(df.samp.dens2)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

sp

</th>

<th style="text-align:right;">

samp.dens

</th>

<th style="text-align:right;">

factor

</th>

<th style="text-align:right;">

iucn.area

</th>

<th style="text-align:right;">

range.area

</th>

<th style="text-align:right;">

new.area

</th>

<th style="text-align:right;">

lost.area

</th>

<th style="text-align:right;">

overlap.area.p

</th>

<th style="text-align:right;">

extra.area.p

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

41062

</td>

<td style="text-align:right;">

15436

</td>

<td style="text-align:right;">

8979

</td>

<td style="text-align:right;">

74

</td>

<td style="text-align:right;">

45

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

500

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

70054

</td>

<td style="text-align:right;">

36186

</td>

<td style="text-align:right;">

738

</td>

<td style="text-align:right;">

98

</td>

<td style="text-align:right;">

105

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

1000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

77038

</td>

<td style="text-align:right;">

42749

</td>

<td style="text-align:right;">

316

</td>

<td style="text-align:right;">

99

</td>

<td style="text-align:right;">

124

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

3000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

86973

</td>

<td style="text-align:right;">

52435

</td>

<td style="text-align:right;">

67

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

152

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

5000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

88554

</td>

<td style="text-align:right;">

53979

</td>

<td style="text-align:right;">

31

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

156

</td>

</tr>

<tr>

<td style="text-align:left;">

GiantPanda

</td>

<td style="text-align:right;">

10000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

34605

</td>

<td style="text-align:right;">

90812

</td>

<td style="text-align:right;">

56210

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

162

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

100

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

646016

</td>

<td style="text-align:right;">

61488

</td>

<td style="text-align:right;">

4268131

</td>

<td style="text-align:right;">

12

</td>

<td style="text-align:right;">

1

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

500

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

2260107

</td>

<td style="text-align:right;">

224729

</td>

<td style="text-align:right;">

2817280

</td>

<td style="text-align:right;">

42

</td>

<td style="text-align:right;">

5

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

1000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

3321054

</td>

<td style="text-align:right;">

361482

</td>

<td style="text-align:right;">

1893088

</td>

<td style="text-align:right;">

61

</td>

<td style="text-align:right;">

7

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

3000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

4621852

</td>

<td style="text-align:right;">

666442

</td>

<td style="text-align:right;">

897249

</td>

<td style="text-align:right;">

82

</td>

<td style="text-align:right;">

14

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

5000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

5106281

</td>

<td style="text-align:right;">

864159

</td>

<td style="text-align:right;">

610536

</td>

<td style="text-align:right;">

87

</td>

<td style="text-align:right;">

18

</td>

</tr>

<tr>

<td style="text-align:left;">

Asian Black Bear

</td>

<td style="text-align:right;">

10000

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

4852659

</td>

<td style="text-align:right;">

5647000

</td>

<td style="text-align:right;">

1163340

</td>

<td style="text-align:right;">

368999

</td>

<td style="text-align:right;">

92

</td>

<td style="text-align:right;">

24

</td>

</tr>

</tbody>

</table>

## PLOT MUCH SAMPLING IS NEEDED FOR ABB

``` r
df.samp.dens2%>%
  mutate(loss=(1-overlap.area.p)+100)%>%
  pivot_longer(c(-samp.dens,-sp, -factor))%>%
  filter(name%in%c("loss","extra.area.p"))%>%
  
  ggplot(aes(x=samp.dens,y=value, color=name))+
  geom_line()+
  geom_point()+
  theme_bw()+
  labs(x="Randomly distributed presences (n)",
       y="Percent",
       title="How Much Sampling is Needed?",
       subtitle="Range Estimates Compared to IUCN")+
  facet_wrap(vars(sp))+
  theme_ipsum()
```

![](README_files/figure-gfm/PLOT%20HOW%20MUCH%20SAMPLING%20IS%20NEEDED%20FOR%20ABB-1.png)<!-- -->

``` r
df.samp.dens2%>%
  pivot_longer(c(-samp.dens, -factor,-sp))%>%
  filter(name%in%c("overlap.area.p","extra.area.p"))%>%
  mutate(sp=as.factor(sp))%>%
  ggplot(aes(x=samp.dens,y=value, color=name))+
  geom_line()+
  geom_point()+
  theme_bw()+
  labs(x="Randomly distributed presences (n)",
       y="Percent of IUCN Range",
       title="How Much Sampling is Needed?",
       subtitle="Range Estimates Compared to IUCN")+
  facet_wrap(vars(sp))+
  geom_hline(yintercept = 100, linetype="dashed")+
  theme_ipsum()
```

![](README_files/figure-gfm/PLOT%20HOW%20MUCH%20SAMPLING%20IS%20NEEDED%20FOR%20ABB-2.png)<!-- -->

``` r
ggsave(here::here("plots", "sampling_points.png"), height=4, width=8)
```

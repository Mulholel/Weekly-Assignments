---
title: "R Assignment, November 5th"
author: "Emma Mulholland"
date: "November 4, 2018"
output: 
  html_document: 
    keep_md: yes
---



**Question 2**
Import RNA-seq file: 

```r
library(readr)
rna_counts<-read_csv("eXpress_dm_counts.csv")
```

```
## Warning: Missing column names filled in: 'X1' [1]
```

```
## Parsed with column specification:
## cols(
##   .default = col_integer(),
##   X1 = col_character()
## )
```

```
## See spec(...) for full column specifications.
```

Creation of mean expression function. I changed any 0 values to NA and then removed the NA values from the mean calculation. This means that the mean expression values will be slightly higher (fewer observations for sum to be divided by), but it makes the log transformations possible. I included the `na.rm = TRUE` condition even for the non-log transformed option to be consistent across each option. 

```r
mean_express<-function(x, y, l2 = FALSE){
  x[x == 0]<- NA
  if(l2 == TRUE){
    mean(log2(x[[y]]), na.rm = TRUE)
  } else {
    mean(x[[y]], na.rm = TRUE)
  }
}
```
Sample use of function on columns from `rna_counts`: 

```r
mean_express(rna_counts, 6)
```

```
## [1] 1458.417
```

```r
mean_express(rna_counts, 6, l2 = TRUE)
```

```
## [1] 8.392313
```

```r
mean_express(rna_counts, 12)
```

```
## [1] 1469.37
```

```r
mean_express(rna_counts, 12, l2 = TRUE)
```

```
## [1] 8.385078
```

```r
mean_express(rna_counts, 9)
```

```
## [1] 2323.996
```

```r
mean_express(rna_counts, 9, l2 = TRUE)
```

```
## [1] 9.385257
```

**Question 3** 
Function for calculating the mean of each column and then placing it into a vector with the associated column names:

```r
mean_express_v<-function(x, l2 = FALSE){
  x[x == 0]<- NA
  vect_mean_express<-numeric()
  if(l2 == TRUE){
    for(i in 2:length(x)){
     z<-mean(log2(x[[i]]), na.rm = TRUE)
     vect_mean_express<-append(vect_mean_express, z)
    } 
    } else {
        for(i in 2:length(x)){
        z<-mean(x[[i]], na.rm = TRUE)
        vect_mean_express<-append(vect_mean_express, z)
        }
    }
  names(vect_mean_express)<-colnames(x[2:length(x)])
  return(vect_mean_express)
}
```
Test with `rna_counts`

```r
mean_express_v(rna_counts)
```

```
##  F101_lg_female_hdhorn F101_lg_female_thxhorn   F101_lg_female_wings 
##               1995.265               2000.166               1611.530 
##  F105_lg_female_hdhorn F105_lg_female_thxhorn   F105_lg_female_wings 
##               2117.814               1458.417               1901.251 
##  F131_lg_female_hdhorn F131_lg_female_thxhorn   F131_lg_female_wings 
##               2131.489               2323.996               2304.294 
##   F135_sm_female_wings  F135_sm_female_hdhorn F135_sm_female_thxhorn 
##               1755.773               1469.370               1794.356 
##  F136_sm_female_hdhorn F136_sm_female_thxhorn   F136_sm_female_wings 
##               2085.804               1792.106               2004.460 
##  F196_sm_female_hdhorn F196_sm_female_thxhorn   F196_sm_female_wings 
##               1360.403               1044.642               3084.915 
##  F197_sm_female_hdhorn F197_sm_female_thxhorn   F197_sm_female_wings 
##               2647.018               2065.087               2118.694 
##  F218_lg_female_hdhorn F218_lg_female_thxhorn   F218_lg_female_wings 
##               2342.413               1968.106               2105.309 
## M120_sm_male_genitalia    M120_sm_male_hdhorn   M120_sm_male_thxhorn 
##               1846.710               2136.398               2119.085 
##     M120_sm_male_wings M125_lg_male_genitalia    M125_lg_male_hdhorn 
##               2565.656               2117.618               2402.461 
##     M125_lg_male_wings M160_lg_male_genitalia    M160_lg_male_hdhorn 
##               2607.358               1745.491               2146.665 
##   M160_lg_male_thxhorn     M160_lg_male_wings M171_sm_male_genitalia 
##               2112.693               2219.074               2056.718 
##    M171_sm_male_hdhorn   M171_sm_male_thxhorn     M171_sm_male_wings 
##               1624.177               1649.561               1851.155 
## M172_sm_male_genitalia    M172_sm_male_hdhorn   M172_sm_male_thxhorn 
##               2209.739               1736.940               1361.131 
##     M172_sm_male_wings M180_lg_male_genitalia    M180_lg_male_hdhorn 
##               2628.183               1940.822               2704.497 
##   M180_lg_male_thxhorn     M180_lg_male_wings M200_sm_male_genitalia 
##               2031.149               3241.669               2429.810 
##    M200_sm_male_hdhorn   M200_sm_male_thxhorn     M200_sm_male_wings 
##               2055.577               2847.835               2237.569 
## M257_lg_male_genitalia    M257_lg_male_hdhorn   M257_lg_male_thxhorn 
##               2190.283               2395.309               2775.140 
##     M257_lg_male_wings 
##               1348.493
```

```r
mean_express_v(rna_counts, l2 = TRUE)
```

```
##  F101_lg_female_hdhorn F101_lg_female_thxhorn   F101_lg_female_wings 
##               9.033518               9.045379               8.314916 
##  F105_lg_female_hdhorn F105_lg_female_thxhorn   F105_lg_female_wings 
##               9.100038               8.392313               8.581792 
##  F131_lg_female_hdhorn F131_lg_female_thxhorn   F131_lg_female_wings 
##               9.009665               9.385257               8.869657 
##   F135_sm_female_wings  F135_sm_female_hdhorn F135_sm_female_thxhorn 
##               8.460318               8.385078               8.943893 
##  F136_sm_female_hdhorn F136_sm_female_thxhorn   F136_sm_female_wings 
##               8.820358               8.853017               8.856807 
##  F196_sm_female_hdhorn F196_sm_female_thxhorn   F196_sm_female_wings 
##               8.564862               8.371856               9.674086 
##  F197_sm_female_hdhorn F197_sm_female_thxhorn   F197_sm_female_wings 
##               9.601041               9.008145               8.730967 
##  F218_lg_female_hdhorn F218_lg_female_thxhorn   F218_lg_female_wings 
##               9.371433               9.182880               8.745546 
## M120_sm_male_genitalia    M120_sm_male_hdhorn   M120_sm_male_thxhorn 
##               8.945391               8.807874               9.315651 
##     M120_sm_male_wings M125_lg_male_genitalia    M125_lg_male_hdhorn 
##               9.203030               8.859041               8.740225 
##     M125_lg_male_wings M160_lg_male_genitalia    M160_lg_male_hdhorn 
##               8.981995               8.925549               8.924142 
##   M160_lg_male_thxhorn     M160_lg_male_wings M171_sm_male_genitalia 
##               9.107550               8.916383               9.148976 
##    M171_sm_male_hdhorn   M171_sm_male_thxhorn     M171_sm_male_wings 
##               8.599081               8.698897               8.526702 
## M172_sm_male_genitalia    M172_sm_male_hdhorn   M172_sm_male_thxhorn 
##               9.200235               8.546968               8.351981 
##     M172_sm_male_wings M180_lg_male_genitalia    M180_lg_male_hdhorn 
##               9.117877               8.950337               8.967102 
##   M180_lg_male_thxhorn     M180_lg_male_wings M200_sm_male_genitalia 
##               8.690568               9.261206               9.378247 
##    M200_sm_male_hdhorn   M200_sm_male_thxhorn     M200_sm_male_wings 
##               8.567874               9.361805               8.784611 
## M257_lg_male_genitalia    M257_lg_male_hdhorn   M257_lg_male_thxhorn 
##               9.213571               8.962720               9.397812 
##     M257_lg_male_wings 
##               8.329061
```
Comparing the outcomes for column 6, 9, 12 shows the same results as those seen in Q2 (in vector from mean_express_v, columns 6,9,12 correspond to elements 5, 8, and 11 because the first column is ignored).
Using the `sort` function, it appears that several of the highest mean expression values are seen in tissue samples from either thorax or head horns. It also appears that the individuals F197, M180, M257, M200 all have high mean expression values in multiple samples. 

**Question 4**
I wrote the following function to accomplish the same task using the apply family of functions ( the "_vapp" in the function name is for "vapply"): 

```r
mean_express_vapp<- function(x, l2 = FALSE){
  if(l2 == TRUE){
  mean_ex_vect<-vapply(x, function(x) {mean(x, na.rm = TRUE)}, FUN.VALUE = 1)
  return(log2(mean_ex_vect))
   } else {
  mean_ex_vect<-vapply(x, function(x) {mean(x, na.rm = TRUE)}, FUN.VALUE = 1)
  return(mean_ex_vect)
  }
}
```
Comparing the speed of each function using `system.time`

```r
system.time(mean_express_v(rna_counts))
```

```
##    user  system elapsed 
##    0.02    0.00    0.01
```

```r
system.time(mean_express_vapp(rna_counts))
```

```
## Warning in mean.default(x, na.rm = TRUE): argument is not numeric or
## logical: returning NA
```

```
##    user  system elapsed 
##    0.02    0.00    0.02
```
The function that uses vapply runs faster (elapsed=0 sec) compared to the function using a for loop. 

**Question 5**
You could use the `colNames` function to calculate and output all the column means with the appropriate names: 

```r
colMeans(rna_counts[2:56], na.rm = TRUE)
```

```
##  F101_lg_female_hdhorn F101_lg_female_thxhorn   F101_lg_female_wings 
##               1978.847               1983.250               1583.904 
##  F105_lg_female_hdhorn F105_lg_female_thxhorn   F105_lg_female_wings 
##               2105.712               1433.749               1869.962 
##  F131_lg_female_hdhorn F131_lg_female_thxhorn   F131_lg_female_wings 
##               2117.847               2307.529               2272.692 
##   F135_sm_female_wings  F135_sm_female_hdhorn F135_sm_female_thxhorn 
##               1728.483               1452.913               1776.309 
##  F136_sm_female_hdhorn F136_sm_female_thxhorn   F136_sm_female_wings 
##               2065.780               1777.769               1988.882 
##  F196_sm_female_hdhorn F196_sm_female_thxhorn   F196_sm_female_wings 
##               1348.898               1025.301               3067.287 
##  F197_sm_female_hdhorn F197_sm_female_thxhorn   F197_sm_female_wings 
##               2639.152               2047.151               2081.889 
##  F218_lg_female_hdhorn F218_lg_female_thxhorn   F218_lg_female_wings 
##               2329.563               1950.561               2074.992 
## M120_sm_male_genitalia    M120_sm_male_hdhorn   M120_sm_male_thxhorn 
##               1832.780               2105.145               2101.163 
##     M120_sm_male_wings M125_lg_male_genitalia    M125_lg_male_hdhorn 
##               2536.920               2088.092               2372.259 
##     M125_lg_male_wings M160_lg_male_genitalia    M160_lg_male_hdhorn 
##               2559.085               1727.538               2111.337 
##   M160_lg_male_thxhorn     M160_lg_male_wings M171_sm_male_genitalia 
##               2087.583               2184.076               2035.093 
##    M171_sm_male_hdhorn   M171_sm_male_thxhorn     M171_sm_male_wings 
##               1598.190               1621.659               1825.344 
## M172_sm_male_genitalia    M172_sm_male_hdhorn   M172_sm_male_thxhorn 
##               2196.101               1713.119               1344.019 
##     M172_sm_male_wings M180_lg_male_genitalia    M180_lg_male_hdhorn 
##               2602.351               1922.634               2670.498 
##   M180_lg_male_thxhorn     M180_lg_male_wings M200_sm_male_genitalia 
##               2003.293               3216.476               2412.038 
##    M200_sm_male_hdhorn   M200_sm_male_thxhorn     M200_sm_male_wings 
##               2032.085               2820.495               2203.813 
## M257_lg_male_genitalia    M257_lg_male_hdhorn   M257_lg_male_thxhorn 
##               2170.258               2361.912               2749.767 
##     M257_lg_male_wings 
##               1325.684
```

```r
#for log2 values: 
log2(colMeans(rna_counts[2:56], na.rm = TRUE))
```

```
##  F101_lg_female_hdhorn F101_lg_female_thxhorn   F101_lg_female_wings 
##               10.95044               10.95365               10.62927 
##  F105_lg_female_hdhorn F105_lg_female_thxhorn   F105_lg_female_wings 
##               11.04009               10.48558               10.86879 
##  F131_lg_female_hdhorn F131_lg_female_thxhorn   F131_lg_female_wings 
##               11.04838               11.17213               11.15019 
##   F135_sm_female_wings  F135_sm_female_hdhorn F135_sm_female_thxhorn 
##               10.75529               10.50473               10.79467 
##  F136_sm_female_hdhorn F136_sm_female_thxhorn   F136_sm_female_wings 
##               11.01247               10.79585               10.95774 
##  F196_sm_female_hdhorn F196_sm_female_thxhorn   F196_sm_female_wings 
##               10.39757               10.00183               11.58275 
##  F197_sm_female_hdhorn F197_sm_female_thxhorn   F197_sm_female_wings 
##               11.36586               10.99940               11.02368 
##  F218_lg_female_hdhorn F218_lg_female_thxhorn   F218_lg_female_wings 
##               11.18584               10.92967               11.01889 
## M120_sm_male_genitalia    M120_sm_male_hdhorn   M120_sm_male_thxhorn 
##               10.83982               11.03970               11.03697 
##     M120_sm_male_wings M125_lg_male_genitalia    M125_lg_male_hdhorn 
##               11.30886               11.02797               11.21205 
##     M125_lg_male_wings M160_lg_male_genitalia    M160_lg_male_hdhorn 
##               11.32141               10.75450               11.04394 
##   M160_lg_male_thxhorn     M160_lg_male_wings M171_sm_male_genitalia 
##               11.02762               11.09281               10.99088 
##    M171_sm_male_hdhorn   M171_sm_male_thxhorn     M171_sm_male_wings 
##               10.64222               10.66326               10.83395 
## M172_sm_male_genitalia    M172_sm_male_hdhorn   M172_sm_male_thxhorn 
##               11.10073               10.74241               10.39234 
##     M172_sm_male_wings M180_lg_male_genitalia    M180_lg_male_hdhorn 
##               11.34560               10.90887               11.38289 
##   M180_lg_male_thxhorn     M180_lg_male_wings M200_sm_male_genitalia 
##               10.96816               11.65127               11.23604 
##    M200_sm_male_hdhorn   M200_sm_male_thxhorn     M200_sm_male_wings 
##               10.98874               11.46173               11.10579 
## M257_lg_male_genitalia    M257_lg_male_hdhorn   M257_lg_male_thxhorn 
##               11.08365               11.20574               11.42509 
##     M257_lg_male_wings 
##               10.37252
```

**Question 6**
You could use the `rowNames` function 

```r
gene_means<- rowMeans(rna_counts[1:nrow(rna_counts), 2:56])
names(gene_means)<-rna_counts$X1
head(gene_means)
```

```
## FBpp0087248 FBpp0293785 FBpp0080383 FBpp0077879 FBpp0311746 FBpp0289081 
##    23.45455  3446.90909    79.54545   139.21818   145.09091  1485.90909
```

```r
#for log2 transformation: 
gene_means_l2<- log(rowMeans(rna_counts[1:nrow(rna_counts), 2:56]))
names(gene_means_l2)<-rna_counts$X1
head(gene_means_l2)
```

```
## FBpp0087248 FBpp0293785 FBpp0080383 FBpp0077879 FBpp0311746 FBpp0289081 
##    3.155064    8.145233    4.376329    4.936042    4.977361    7.303782
```


**Question 7 **
Subsetting data into male head horn samples (lg and sm treatments): 

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
male_samples<- rna_counts %>% select(starts_with("M"))
male_samples_hdhorns<-male_samples %>% select(contains("hdhorn"))
male_samples_hdhorns_sm<- male_samples_hdhorns %>% select(contains("sm"))  
male_samples_hdhorns_lg<- male_samples_hdhorns %>% select(contains("lg"))
```

Calculating mean expression values for each gene in sm and lg treatments, and the mean differences (abs values): 

```r
rMeans_male_samples_hdhorns<-rowMeans(male_samples_hdhorns, na.rm = TRUE)
rMeans_male_samples_hdhorns_lg<-rowMeans(male_samples_hdhorns_lg, na.rm = TRUE)
rMeans_male_samples_hdhorns_sm<-rowMeans(male_samples_hdhorns_sm, na.rm = TRUE)
mean_diff_lg_sm_hdhords<-abs(rMeans_male_samples_hdhorns_lg - rMeans_male_samples_hdhorns_sm)
```


**Question 8**

```r
#Plot of male sample hdhorns means against means difference
plot(rMeans_male_samples_hdhorns, mean_diff_lg_sm_hdhords, type = "p")
```

![](November5_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
#Plot of log2 transformed male sample hdhorns means and log2 transformed means difference
plot(log2(rMeans_male_samples_hdhorns), log2(mean_diff_lg_sm_hdhords), type = "p")
```

![](November5_files/figure-html/unnamed-chunk-12-2.png)<!-- -->


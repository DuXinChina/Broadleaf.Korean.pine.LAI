---
title: "Handbook of Broadleaf.Korean.pine.LAI"
author: ""
date: ''
geometry: "left=0.5cm,right=0.5cm,top=0.5cm,bottom=0.5cm"
output:
  html_document: default
  word_document: default
editor_options:
  chunk_output_type: inline
  fontsize: 9pt
---

## 软件包简介

叶面积指数（LAI）指单位水平地面面积上总叶面积的一半，为植物冠层结构的最常用参数之一。我团队提出了一种基于分层泰森多边形估算复层林叶面积指数的方法，相较于传统基于公式的叶面积指数估计方法，该法可通过泰森多边形估算林分中各林木的树冠下叶面积指数，并在此基础上较为准确的计算林分中任意位点处的若干与叶面积指数相关的冠层结构参数（局部叶面积指数、局部叶面积指数标注差、各类树种对局部叶面积指数的贡献率，叶面积指数垂直结构特征等）。

我团队将基于分层泰森多边形估算复层林叶面积指数的方法用于阔叶红松林叶面积指数的估算，开发了R包Broadleaf.Korean.pine.LAI。

## 依赖关系

本包依赖于 ape；sp ；sf；raster；grDevices；rgeos；gstat ；tcltk ；ggplot2；plot3D ；deldir ；grDevices；stats 运行前请先安装

## 软件包安装方式

library(devtools);install_github("DuXinChina/Broadleaf.Korean.pine.LAI/Broadleaf.Korean.pine.LAI")

## 1 示例数据展示

示例数据为虚拟数据，数据为50m×50m样地，样地中包含色木、紫椴、冷杉、红松四种树种的空间坐标、胸径、树高信息。依据树高，将阔叶红松林分为不同林层。其中，下木层树高小于10m；亚林层树高为10m---16m；主林层树高为16m---30m；超冠层树高大于30m。虚拟数据中，超冠层包含两株红松，为随机分布；主林层中包含15株红松、10株冷杉、10株紫椴、5株色木，均为随机分布；亚林层中包含50株色木，为以（10,40）和（40,10）为中心的各25株的聚集分布；林下层包含四类树种各20株植株，各树种在小尺度为聚集分布，在大尺度为随机分布。超冠层中，林木胸径大于60cm；主林层中，林木胸径为30cm---60 cm；亚林层中林木胸径为10 cm-30cm；下木层胸径小于10 cm。

数据中x与y列分别为林木横纵坐标、DBH为胸径、H为树高、Species列为树种种名。

Species列可包含的常见树种种名为HS 红松、YS 云杉、LS 冷杉、ZD 紫椴、KD 糠椴、MGL 蒙古栎、SQL 水曲柳、HTQ 胡桃楸、HBL 黄檗、SMQ 色木槭、QKQ 青楷槭。HKQ 花楷槭、JSQ 假色槭、NJQ 拧筋槭、BNQ 白牛槭、CY 春榆、BH 白桦、HH 怀槐、FH 枫桦、LYY 裂叶榆、QT 其他树种。

```{r}
data = Broadleaf.Korean.pine.LAI::b
data
```

## 2 示例

### 2.1 LAI.single(a,b,r)

2.1.1功能介绍： 基于传统方法，计算阔叶红松林内以任意位点为中心，一定半径的圆形区域内的叶面积指数，并给出不同树种在叶面积指数中的占比，及针叶树种与阔叶树种叶面积指数的占比。

2.1.2参数意义： a---需要计算叶面积指数的位点坐标 b---样地中林木的坐标、树种、胸径、树高信息 r---以位点中心圆的半径

2.1.3例：

计算以（25,25）为中心5m为半径的圆内的叶面积

```{r}
a=matrix(c(25,25),1,2)
colnames(a)=c("x","y")
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b[,-4]
head(b)
```

```{r}
Broadleaf.Korean.pine.LAI::LAI.single(a=a,b=b,r=5)
```

总叶面积指数为12.23283、树种色木槭叶面积为12.23283、阔叶树种叶面积12.23283、阔叶占比100%

### 2.2 LAI.mult(a,b,r)

2.2.1功能介绍： 计算阔叶红松林内的多个位点为中心，一定半径的圆形区域内的叶面积指数，针叶树种与阔叶树种叶面积指数以及针阔叶面积分别的占比。

2.2.2参数意义：

a---需要计算叶面积指数的多个位点坐标

b---样地中林木的坐标、树种、胸径、树高信息

r---以位点中心圆的半径

2.2.3例：

计算多个样点的叶面积指数，针叶叶面积指数，阔叶叶面积指数以及针阔叶树种叶面积指数的占比

```{r}
a=matrix(runif(10,5,45),5,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b[,-4]
head(b)
```

```{r}
Broadleaf.Korean.pine.LAI::LAI.mult(a=a,b=b,r=5)
```

### 2.3 Plot.LAI.Krig (minx,maxx,miny,maxy,b,seq,r)

2.3.1功能介绍： 绘制阔叶红松林叶面积指数的克里格插值图。计算阔叶红松林叶面积指数半变异函数

2.3.2参数意义：

minx---绘制克里格插值图的最小横坐标

maxx---绘制克里格插值图的最大横坐标

miny---绘制克里格插值图的最小纵坐标

maxy---绘制克里格插值图的最大纵坐标

b---样地中的林木坐标、树种、胸径

seq---绘制克里格图的空间分辨率

r---以位点中心圆的半径

2.3.3例：

\####

```{r }
b=Broadleaf.Korean.pine.LAI::b[,-4]
head(b)

```

```{r }
head(b)
p=Broadleaf.Korean.pine.LAI::Plot.LAI.Krig(5,45,5,45,b=b,seq=20,r=3)
p1=p+geom_vline(xintercept = c(5,45),linetype=2)+geom_hline(yintercept =c(5,45),linetype=2)
p1=p+scale_x_continuous(expand= c(0, 5))+scale_y_continuous(expand= c(0, 5))
p2=p1+geom_point(data=b,aes(x=x,y=y),size=b$DBH/8,col="grey4",alpha=0.2)+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
p2
```

### **2.4 Standing_Forest_LAI(minx,maxx,miny,maxy,b)**

2.4.1功能介绍：

计算阔叶红松林全林分的叶面积指数，并给出不同树种叶面积指数的量与占比，以及针、阔叶树种叶面积指数的量与占比。

2.4.2 参数意义：

minx---计算叶面积指数样地横坐标最小边界

maxx---计算叶面积指数样地横坐标最大边界

miny---计算叶面积指数样地纵坐标最小边界

maxy---计算叶面积指数样地纵坐标最大边界

b---林木分布信息

2.4.3例：

```{r}
b=Broadleaf.Korean.pine.LAI::b[,-4]
head(b)

```

```{r}
Broadleaf.Korean.pine.LAI::Standing_Forest_LAI(minx=5,maxx=45,miny=5,maxy=45,b)

```

林分叶面积指数11.896。红松叶面积指数3.948；冷杉叶面积指数1.242；紫椴叶面积指数1.353；色木叶面积指数5.352。针叶树叶面积指5.191；阔叶树叶面积指数6.705。针叶树种叶面积指数占比43.636%；阔叶树种叶面积指数占比56.364%。

### **2.5 LSD_LAI(a,b,r,Lr)**

2.5.1功能介绍：

在以一定半径圆的面积为分辨率计算阔叶红松林叶面指数后。计算以某一位点为中心，一定半径圆内叶面积指数的局部标准差。通常情况下，林缘处及大树小树的冠层交错处的叶面积指数具有较大的局部标准差。

2.5.2参数意义：

a---需要计算叶面积指数局部标准差的位点坐标

b---林木分布信息及胸径、树高

r---计算林木内叶面积指数时所用圆的半径

Lr---计算位点叶面积指数局部标准差时所用圆的半径

2.5.3例：

```{r}
a=data.frame(28,25)
colnames(a)=c("x","y")
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b[,-4]
head(b)

```

```{r}
result=Broadleaf.Korean.pine.LAI::LSD_LAI(a=a,b=b,r=3,Lr=1.5)
result
```

### **2.6 LSD_LAI_mult(a,b,r,Lr)**

2.6.1功能介绍：

计算多个位点叶面积指数的局部标准差

2.6.2参数意义：

a---需要计算叶面积指数局部标准差的位点坐标

b---林木分布信息及胸径、树高

r---计算林木内叶面积指数时所用圆的半径

Lr---计算位点叶面积指数局部标准差时所用圆的半径

2.6.3例：

```{r}
a=matrix(runif(10,5,45),5,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b[,-4]
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::LSD_LAI_mult(a=a,b=b,r=3,Lr=1.5)
formattable::formattable(result)
```

### **2.7 Plot.LSD_LAI.Krig(minx,maxx,miny,maxy,b,seq,r,Lr)**

2.7.1功能介绍：

绘制阔叶红松林叶面积指数局部标准差的克里格插值图。计算阔叶红松林叶面积指数局部标准差的半变异函数。

2.7.2参数意义：

minx---绘制克里格插值图的最小横坐标

maxx---绘制克里格插值图的最大横坐标

miny---绘制克里格插值图的最小纵坐标

maxy---绘制克里格插值图的最大纵坐标

b---样地中的林木坐标、树种、胸径

seq---绘制克里格图的空间分辨率

r---计算林分叶面积指数时以位点中心圆的半径

Lr---计算林分叶面积指数局部标准差时取样圆的半径

2.7.3例：

```{r}
p=Broadleaf.Korean.pine.LAI::Plot.LSD_LAI.Krig(minx=5,maxx=45,miny=5,maxy=45,b=b,seq=20,r=3,Lr=1)
p1=p+geom_vline(xintercept = c(5,45),linetype=2)+geom_hline(yintercept =c(5,45),linetype=2)
p1=p+scale_x_continuous(expand= c(0, 5))+scale_y_continuous(expand= c(0, 5))
p2=p1+geom_point(data=b,aes(x=x,y=y),size=b$DBH/8,col="grey4",alpha=0.2)+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
p2
```

### **2.8 Voronoi.LAI(minx,maxx,miny,maxy,boundary,b,r)**

2.8.1功能介绍：

基于泰森多边形法计算林地内（单一林层）各株树木下的叶面积指数。

2.8.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

2.8.3例：

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H<10)
b=b[,-4]
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary=5,b=b,r=2)
formattable::formattable(result)
```

### **2.9 Plot.Voronoi.LAI(minx,maxx,miny,maxy,boundary,b,r)**

2.9.1功能介绍：

基于泰森多边形法计算单层林或复层林内单一林层的林地内各位置的叶面积指数，并绘图。

2.9.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

2.9.3例：

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H>20 & b$H<30)
b=b[,-4]
head(b)
```

```{r}
p=Broadleaf.Korean.pine.LAI::Plot.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary=5,b=b,r=3.5)
p1=p+geom_point(data=b,aes(x=x,y=y),size=b$DBH/8,col="grey4",alpha=0.2)
p1
```

### 2.10 Single.point.Voronoi.LAI(minx,maxx,miny,maxy,boundary,a,b,r)

2.10.1功能介绍：

基于泰森多边形法计算单层林或复层林内单一林层的林地内任意单一位点的叶面积指数。

2.10.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要计算叶面积指数的位点的坐标

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

2.10.3例：

```{r}
a=matrix(runif(2,5,45),1,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H>20 & b$H<30)
b=b[,-4]
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::Single.point.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,r=3.5)
formattable::formattable(result)
```

### 2.11 mult.point.Voronoi.LAI(minx,maxx,miny,maxy,boundary,a,b,r)

2.11.1功能介绍：

基于泰森多边形法计算单层林或复层林内单一林层的林地内任意多个位点的叶面积指数。

2.11.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要计算叶面积指数的位点坐标

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

2.11.3例：

```{r}
set.seed(100)
a=matrix(runif(40,5,45),20,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
head(a)
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H>20 & b$H<30)
b=b[,-4]
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::mult.point.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,r=3.5)
formattable::formattable(result)
```

### 2.12 Local.single.point.Voronoi.LAI(minx,maxx,miny,maxy,boundary,a,b,r,Lr)

2.12.1功能介绍：

基于泰森多边形法计算单层林或复层林内单一林层的林地内任意单一位点的与叶面积有关的局部冠层结构参数。

2.12.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要模拟局部冠层结构参数的位置坐标

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

Lr---提取局部冠层结构参数时的圆的半径

2.12.3例：

```{r}
set.seed(100)
a=matrix(runif(2,5,45),1,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H>20 & b$H<30)
b=b[,-4]
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::Local.single.point.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,r=3.5,Lr=1.5)
formattable::formattable(result)
```

### 2.13 Local.mult.point.Voronoi.LAI(minx,maxx,miny,maxy,boundary,a,b,r,Lr)

2.13.1功能介绍：

基于泰森多边形法计算单层林或复层林内单一林层的林地内任意多个位点的与叶面积有关的局部冠层结构参数。

2.13.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要模拟局部冠层结构参数的位置坐标

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

Lr---提取局部冠层结构参数时的圆的半径

2.13.3例：

```{r}
set.seed(100)
a=matrix(runif(40,5,45),20,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
head(a)
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H>20 & b$H<30)
b=b[,-4]
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::Local.mult.point.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,r=3.5,Lr=1.5)
formattable::formattable(result)
```

### 2.14 Semivariogram.Voronoi.LAI.Single(minx,maxx,miny,maxy,boundary,b,r,seq)

2.14.1功能介绍：

基于泰森多边形法模拟单层林或复层林内单一林层的叶面积指数的半变异函数

2.14.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

seq---模拟半变异函数时的空间分辨率

2.14.3例：

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H>20 & b$H<30)
b=b[,-4]
head(b)
```

```{r}
Broadleaf.Korean.pine.LAI::Semivariogram.Voronoi.LAI.Single(minx=0,maxx=50,miny=0,maxy=50,boundary=5,b=b,r=3.5,seq=20)
```

### 2.15 Plot.Voronoi.LAI.Sum(minx,maxx,miny,maxy,boundary,b,strata,r)

2.15.1功能介绍：

基于分层泰森多边形法绘制复层林林分叶面积指数图

2.15.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

b---样地中的林木坐标、树种、胸径

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

2.15.3例：

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r}
p=Broadleaf.Korean.pine.LAI::Plot.Voronoi.LAI.Sum(minx=0,maxx=50,miny=0,maxy=50,boundary=5,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4))
p1=p+geom_vline(xintercept = c(5,45),linetype=2)+geom_hline(yintercept =c(5,45),linetype=2)
p1=p+scale_x_continuous(expand= c(0, 5))+scale_y_continuous(expand= c(0, 5))
p2=p1+geom_point(data=b,aes(x=x,y=y),size=b$DBH/8,col="grey4",alpha=0.2)+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
p2
```

### 2.16 Single.point.Voronoi.LAI.sum(minx,maxx,miny,maxy,boundary,a,b,strata,r)

2.16.1功能介绍：

基于分层泰森多边形法计算复层林内单一任意位点的叶面积指数

2.16.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要计算叶面积指数的位点的坐标

b---样地中的林木坐标、树种、胸径

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

2.16.3例：

```{r}
set.seed(100)
a=matrix(runif(2,5,45),1,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::Single.point.Voronoi.LAI.sum(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4))
result
```

### 2.17 mult.point.Voronoi.LAI.sum(minx,maxx,miny,maxy,boundary,a,b,strata,r)

2.17.1功能介绍：

基于分层泰森多边形法计算复层林内多个任意位点的叶面积指数

2.17.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要计算叶面积指数的位点的坐标

b---样地中的林木坐标、树种、胸径

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

2.17.3 例：

```{r}
set.seed(100)
a=matrix(runif(40,5,45),20,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::mult.point.Voronoi.LAI.sum(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4))
head(result)
```

### 2.18 Local.single.point.Voronoi.LAI.sum(minx,maxx,miny,maxy,boundary,a,b,strata,r,Lr)

2.18.1功能介绍：

基于分层泰森多边形法计算复层林内单一任意位点的与叶面积指数相关的局部冠层结构特征

2.18.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要模拟局部冠层结构参数的位点的坐标

b---样地中的林木坐标、树种、胸径

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

Lr---模拟局部冠层结构参数的圆的半径

2.18.3 例：

```{r}
set.seed(100)
a=matrix(runif(2,5,45),1,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::Local.single.point.Voronoi.LAI.sum(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4),Lr=1.5)
formattable::formattable(result)
```

### 2.19 Local.mult.point.Voronoi.LAI.sum(minx,maxx,miny,maxy,boundary,a,b,strata,r,Lr)

2.19.1功能介绍：

基于分层泰森多边形法计算复层林内多个任意位点的与叶面积指数相关的局部冠层结构特征

2.19.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---样地边界缓冲区宽度

a---需要模拟局部冠层结构参数的位点的坐标

b---样地中的林木坐标、树种、胸径

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

Lr---模拟局部冠层结构参数的圆的半径

2.19.3 例：

```{r}
set.seed(100)
a=matrix(runif(40,5,45),20,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r}
result=Broadleaf.Korean.pine.LAI::Local.mult.point.Voronoi.LAI.sum(minx=0,maxx=50,miny=0,maxy=50,boundary=5,a=a,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4),Lr=1.5)
formattable::formattable(result)
```

### 2.20 Semivariogram.Voronoi.LAI(minx,maxx,miny,maxy,boundary,b,strata,r,seq )

2.20.1功能介绍：

基于分层泰森多边形法模拟复层林的叶面积指数的半变异函数

2.20.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---缓冲区宽度

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

seq---模拟半变异函数时的空间分辨率

2.20.3例：

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r}
Broadleaf.Korean.pine.LAI::Semivariogram.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary =5,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4),seq = 20)
```
### 2.21 plot.Local.single.point.Voronoi.LAI(minx,maxx,miny,maxy,a,b,r,Lr)

2.21.1功能介绍：

基于泰森多边形法计算单层林或复层林内单一林层的林地内任意单一位点的与叶面积有关的局部冠层结构参数,并绘图验证。

2.21.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

a---需要模拟局部冠层结构参数的位置坐标

b---样地中的林木坐标、树种、胸径

r---将在林分中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为添加样点构成泰森多边形外接圆的半径

Lr---提取局部冠层结构参数时的圆的半径

2.21.3例：

```{r}
set.seed(100)
a=matrix(runif(2,5,45),1,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
b=subset(b,b$H>20 & b$H<30)
b=b[,-4]
head(b)
```

```{r,message=FALSE, warning=FALSE}
Broadleaf.Korean.pine.LAI::plot.Local.single.point.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,a=a,b=b,r=3.5,Lr=1.5)

```

### 2.22 plot.Local.single.point.Voronoi.LAI.sum(minx,maxx,miny,maxy,a,b,strata,r,Lr)

2.22.1功能介绍：

基于分层泰森多边形法计算复层林内单一任意位点的与叶面积指数相关的局部冠层结构特征，并绘图

2.22.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

a---需要模拟局部冠层结构参数的位点的坐标

b---样地中的林木坐标、树种、胸径

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

Lr---模拟局部冠层结构参数的圆的半径

2.22.3 例：

```{r}
set.seed(100)
a=matrix(runif(2,5,45),1,2)
colnames(a)=c("x","y")
a=as.data.frame(a)
a
```

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r,message=FALSE, warning=FALSE}
gc()
Broadleaf.Korean.pine.LAI::plot.Local.single.point.Voronoi.LAI.sum(minx=0,maxx=50,miny=0,maxy=50,a=a,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4),Lr=1.5)

```

### 2.23 Voronoi.LAI.sum.ISAA(minx, maxx, miny, maxy, boundary, b, strata, r, indis, lag)

2.23.1功能介绍：

基于分层泰森多边形法进行复层林内叶面积指数的增量空间自相关分析

2.23.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---缓冲区宽度

b---样地中的林木坐标、树种、胸径

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

indis---初始距离

lag---滞后距增量

2.23.3 例：

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r,message=FALSE, warning=FALSE}
gc()
Broadleaf.Korean.pine.LAI::Voronoi.LAI.sum.ISAA(minx=0, maxx=50, miny=0, maxy=50, boundary=5, b=b, strata=c(10,16,30),r=c(2,2.5,3.5,4), indis=1, lag=1)

```


### 2.24 Voronoi.pointcloud(minx, maxx, miny, maxy, boundary, b, seq, strata, r,
    S, theta, phi)

2.24.1功能介绍：

基于分层泰森多边形反演叶面积的三维空间分布，并以点云绘图

2.24.2参数意义：

minx---样地横坐标最小范围

maxx---样地横坐标最大范围

miny---样地纵坐标最小范围

maxy---样地纵坐标最大范围

boundary---缓冲区宽度

b---样地中的林木坐标、树种、胸径

seq---绘制点云时的空间分辨率

strata---不同林层间的分割高度

r---将在林分各林层中没有林木分布的位置添加样点，以避免以林木为中心的泰森多边形过大。r为不同林层添加样点构成泰森多边形外接圆的半径

S---林分内平均冠形率，冠长与冠幅的比例

theta---绘图水平翻转角度

phi---绘图垂直翻转角度

2.24.3 例：

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```


```{r,message=FALSE, warning=FALSE}

Broadleaf.Korean.pine.LAI::Voronoi.pointcloud(minx=0, maxx=50, miny=0, maxy=50, boundary=5, b=b, seq=100, strata=c(10,16,30),r=c(2,2.5,3.5,4), S=1.5, theta=120, phi=20)
Broadleaf.Korean.pine.LAI::Voronoi.pointcloud(minx=0, maxx=50, miny=0, maxy=50, boundary=5, b=b, seq=100, strata=c(10,16,30),r=c(2,2.5,3.5,4), S=1.5, theta=0, phi=90)
```

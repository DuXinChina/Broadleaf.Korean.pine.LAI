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

## ��������

Ҷ���ָ����LAI��ָ��λˮƽ�����������Ҷ�����һ�룬Ϊֲ��ڲ�ṹ����ò���֮һ�����Ŷ������һ�ֻ��ڷֲ�̩ɭ����ι��㸴����Ҷ���ָ���ķ���������ڴ�ͳ���ڹ�ʽ��Ҷ���ָ�����Ʒ������÷���ͨ��̩ɭ����ι����ַ��и���ľ��������Ҷ���ָ�������ڴ˻����Ͻ�Ϊ׼ȷ�ļ����ַ�������λ�㴦��������Ҷ���ָ����صĹڲ�ṹ�������ֲ�Ҷ���ָ�����ֲ�Ҷ���ָ����ע��������ֶԾֲ�Ҷ���ָ���Ĺ����ʣ�Ҷ���ָ����ֱ�ṹ�����ȣ���

���Ŷӽ����ڷֲ�̩ɭ����ι��㸴����Ҷ���ָ���ķ���������Ҷ������Ҷ���ָ���Ĺ��㣬������R��Broadleaf.Korean.pine.LAI��

## ������ϵ

���������� ape��sp ��sf��raster��grDevices��rgeos��gstat ��tcltk ��ggplot2��plot3D ��deldir ��grDevices��stats ����ǰ���Ȱ�װ

## �������װ��ʽ

library(devtools);install_github("DuXinChina/Broadleaf.Korean.pine.LAI/Broadleaf.Korean.pine.LAI")

## 1 ʾ������չʾ

ʾ������Ϊ�������ݣ�����Ϊ50m��50m���أ������а���ɫľ����鲡���ɼ�������������ֵĿռ����ꡢ�ؾ���������Ϣ���������ߣ�����Ҷ�����ַ�Ϊ��ͬ�ֲ㡣���У���ľ������С��10m�����ֲ�����Ϊ10m---16m�����ֲ�����Ϊ16m---30m�����ڲ����ߴ���30m�����������У����ڲ����������ɣ�Ϊ����ֲ������ֲ��а���15����ɡ�10����ɼ��10����鲡�5��ɫľ����Ϊ����ֲ������ֲ��а���50��ɫľ��Ϊ�ԣ�10,40���ͣ�40,10��Ϊ���ĵĸ�25��ľۼ��ֲ������²�����������ָ�20��ֲ�꣬��������С�߶�Ϊ�ۼ��ֲ����ڴ�߶�Ϊ����ֲ������ڲ��У���ľ�ؾ�����60cm�����ֲ��У���ľ�ؾ�Ϊ30cm---60 cm�����ֲ�����ľ�ؾ�Ϊ10 cm-30cm����ľ���ؾ�С��10 cm��

������x��y�зֱ�Ϊ��ľ�������ꡢDBHΪ�ؾ���HΪ���ߡ�Species��Ϊ����������

Species�пɰ����ĳ�����������ΪHS ���ɡ�YS ��ɼ��LS ��ɼ��ZD ��鲡�KD ��鲡�MGL �ɹ��ݡ�SQL ˮ������HTQ ����鱡�HBL ���ޡ�SMQ ɫľ�ʡ�QKQ �࿬�ʡ�HKQ �����ʡ�JSQ ��ɫ�ʡ�NJQ š���ʡ�BNQ ��ţ�ʡ�CY ���ܡ�BH ���롢HH ������FH ���롢LYY ��Ҷ�ܡ�QT �������֡�

```{r}
data = Broadleaf.Korean.pine.LAI::b
data
```

## 2 ʾ��

### 2.1 LAI.single(a,b,r)

2.1.1���ܽ��ܣ� ���ڴ�ͳ������������Ҷ��������������λ��Ϊ���ģ�һ���뾶��Բ�������ڵ�Ҷ���ָ������������ͬ������Ҷ���ָ���е�ռ�ȣ�����Ҷ��������Ҷ����Ҷ���ָ����ռ�ȡ�

2.1.2�������壺 a---��Ҫ����Ҷ���ָ����λ������ b---��������ľ�����ꡢ���֡��ؾ���������Ϣ r---��λ������Բ�İ뾶

2.1.3����

�����ԣ�25,25��Ϊ����5mΪ�뾶��Բ�ڵ�Ҷ���

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

��Ҷ���ָ��Ϊ12.23283������ɫľ��Ҷ���Ϊ12.23283����Ҷ����Ҷ���12.23283����Ҷռ��100%

### 2.2 LAI.mult(a,b,r)

2.2.1���ܽ��ܣ� ������Ҷ�������ڵĶ��λ��Ϊ���ģ�һ���뾶��Բ�������ڵ�Ҷ���ָ������Ҷ��������Ҷ����Ҷ���ָ���Լ�����Ҷ����ֱ��ռ�ȡ�

2.2.2�������壺

a---��Ҫ����Ҷ���ָ���Ķ��λ������

b---��������ľ�����ꡢ���֡��ؾ���������Ϣ

r---��λ������Բ�İ뾶

2.2.3����

�����������Ҷ���ָ������ҶҶ���ָ������ҶҶ���ָ���Լ�����Ҷ����Ҷ���ָ����ռ��

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

2.3.1���ܽ��ܣ� ������Ҷ������Ҷ���ָ���Ŀ�����ֵͼ��������Ҷ������Ҷ���ָ������캯��

2.3.2�������壺

minx---���ƿ�����ֵͼ����С������

maxx---���ƿ�����ֵͼ����������

miny---���ƿ�����ֵͼ����С������

maxy---���ƿ�����ֵͼ�����������

b---�����е���ľ���ꡢ���֡��ؾ�

seq---���ƿ����ͼ�Ŀռ�ֱ���

r---��λ������Բ�İ뾶

2.3.3����

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

2.4.1���ܽ��ܣ�

������Ҷ������ȫ�ֵַ�Ҷ���ָ������������ͬ����Ҷ���ָ��������ռ�ȣ��Լ��롢��Ҷ����Ҷ���ָ��������ռ�ȡ�

2.4.2 �������壺

minx---����Ҷ���ָ�����غ�������С�߽�

maxx---����Ҷ���ָ�����غ��������߽�

miny---����Ҷ���ָ��������������С�߽�

maxy---����Ҷ���ָ���������������߽�

b---��ľ�ֲ���Ϣ

2.4.3����

```{r}
b=Broadleaf.Korean.pine.LAI::b[,-4]
head(b)

```

```{r}
Broadleaf.Korean.pine.LAI::Standing_Forest_LAI(minx=5,maxx=45,miny=5,maxy=45,b)

```

�ַ�Ҷ���ָ��11.896������Ҷ���ָ��3.948����ɼҶ���ָ��1.242�����Ҷ���ָ��1.353��ɫľҶ���ָ��5.352����Ҷ��Ҷ���ָ5.191����Ҷ��Ҷ���ָ��6.705����Ҷ����Ҷ���ָ��ռ��43.636%����Ҷ����Ҷ���ָ��ռ��56.364%��

### **2.5 LSD_LAI(a,b,r,Lr)**

2.5.1���ܽ��ܣ�

����һ���뾶Բ�����Ϊ�ֱ��ʼ�����Ҷ������Ҷ��ָ���󡣼�����ĳһλ��Ϊ���ģ�һ���뾶Բ��Ҷ���ָ���ľֲ���׼�ͨ������£���Ե��������С���Ĺڲ㽻����Ҷ���ָ�����нϴ�ľֲ���׼�

2.5.2�������壺

a---��Ҫ����Ҷ���ָ���ֲ���׼���λ������

b---��ľ�ֲ���Ϣ���ؾ�������

r---������ľ��Ҷ���ָ��ʱ����Բ�İ뾶

Lr---����λ��Ҷ���ָ���ֲ���׼��ʱ����Բ�İ뾶

2.5.3����

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

2.6.1���ܽ��ܣ�

������λ��Ҷ���ָ���ľֲ���׼��

2.6.2�������壺

a---��Ҫ����Ҷ���ָ���ֲ���׼���λ������

b---��ľ�ֲ���Ϣ���ؾ�������

r---������ľ��Ҷ���ָ��ʱ����Բ�İ뾶

Lr---����λ��Ҷ���ָ���ֲ���׼��ʱ����Բ�İ뾶

2.6.3����

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

2.7.1���ܽ��ܣ�

������Ҷ������Ҷ���ָ���ֲ���׼��Ŀ�����ֵͼ��������Ҷ������Ҷ���ָ���ֲ���׼��İ���캯����

2.7.2�������壺

minx---���ƿ�����ֵͼ����С������

maxx---���ƿ�����ֵͼ����������

miny---���ƿ�����ֵͼ����С������

maxy---���ƿ�����ֵͼ�����������

b---�����е���ľ���ꡢ���֡��ؾ�

seq---���ƿ����ͼ�Ŀռ�ֱ���

r---�����ַ�Ҷ���ָ��ʱ��λ������Բ�İ뾶

Lr---�����ַ�Ҷ���ָ���ֲ���׼��ʱȡ��Բ�İ뾶

2.7.3����

```{r}
p=Broadleaf.Korean.pine.LAI::Plot.LSD_LAI.Krig(minx=5,maxx=45,miny=5,maxy=45,b=b,seq=20,r=3,Lr=1)
p1=p+geom_vline(xintercept = c(5,45),linetype=2)+geom_hline(yintercept =c(5,45),linetype=2)
p1=p+scale_x_continuous(expand= c(0, 5))+scale_y_continuous(expand= c(0, 5))
p2=p1+geom_point(data=b,aes(x=x,y=y),size=b$DBH/8,col="grey4",alpha=0.2)+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
p2
```

### **2.8 Voronoi.LAI(minx,maxx,miny,maxy,boundary,b,r)**

2.8.1���ܽ��ܣ�

����̩ɭ����η������ֵ��ڣ���һ�ֲ㣩������ľ�µ�Ҷ���ָ����

2.8.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

2.8.3����

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

2.9.1���ܽ��ܣ�

����̩ɭ����η����㵥���ֻ򸴲����ڵ�һ�ֲ���ֵ��ڸ�λ�õ�Ҷ���ָ��������ͼ��

2.9.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

2.9.3����

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

2.10.1���ܽ��ܣ�

����̩ɭ����η����㵥���ֻ򸴲����ڵ�һ�ֲ���ֵ������ⵥһλ���Ҷ���ָ����

2.10.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫ����Ҷ���ָ����λ�������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

2.10.3����

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

2.11.1���ܽ��ܣ�

����̩ɭ����η����㵥���ֻ򸴲����ڵ�һ�ֲ���ֵ���������λ���Ҷ���ָ����

2.11.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫ����Ҷ���ָ����λ������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

2.11.3����

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

2.12.1���ܽ��ܣ�

����̩ɭ����η����㵥���ֻ򸴲����ڵ�һ�ֲ���ֵ������ⵥһλ�����Ҷ����йصľֲ��ڲ�ṹ������

2.12.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫģ��ֲ��ڲ�ṹ������λ������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

Lr---��ȡ�ֲ��ڲ�ṹ����ʱ��Բ�İ뾶

2.12.3����

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

2.13.1���ܽ��ܣ�

����̩ɭ����η����㵥���ֻ򸴲����ڵ�һ�ֲ���ֵ���������λ�����Ҷ����йصľֲ��ڲ�ṹ������

2.13.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫģ��ֲ��ڲ�ṹ������λ������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

Lr---��ȡ�ֲ��ڲ�ṹ����ʱ��Բ�İ뾶

2.13.3����

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

2.14.1���ܽ��ܣ�

����̩ɭ����η�ģ�ⵥ���ֻ򸴲����ڵ�һ�ֲ��Ҷ���ָ���İ���캯��

2.14.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

seq---ģ�����캯��ʱ�Ŀռ�ֱ���

2.14.3����

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

2.15.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η����Ƹ������ַ�Ҷ���ָ��ͼ

2.15.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

b---�����е���ľ���ꡢ���֡��ؾ�

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

2.15.3����

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

2.16.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η����㸴�����ڵ�һ����λ���Ҷ���ָ��

2.16.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫ����Ҷ���ָ����λ�������

b---�����е���ľ���ꡢ���֡��ؾ�

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

2.16.3����

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

2.17.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η����㸴�����ڶ������λ���Ҷ���ָ��

2.17.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫ����Ҷ���ָ����λ�������

b---�����е���ľ���ꡢ���֡��ؾ�

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

2.17.3 ����

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

2.18.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η����㸴�����ڵ�һ����λ�����Ҷ���ָ����صľֲ��ڲ�ṹ����

2.18.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫģ��ֲ��ڲ�ṹ������λ�������

b---�����е���ľ���ꡢ���֡��ؾ�

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

Lr---ģ��ֲ��ڲ�ṹ������Բ�İ뾶

2.18.3 ����

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

2.19.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η����㸴�����ڶ������λ�����Ҷ���ָ����صľֲ��ڲ�ṹ����

2.19.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���ر߽绺�������

a---��Ҫģ��ֲ��ڲ�ṹ������λ�������

b---�����е���ľ���ꡢ���֡��ؾ�

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

Lr---ģ��ֲ��ڲ�ṹ������Բ�İ뾶

2.19.3 ����

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

2.20.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η�ģ�⸴���ֵ�Ҷ���ָ���İ���캯��

2.20.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

seq---ģ�����캯��ʱ�Ŀռ�ֱ���

2.20.3����

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```

```{r}
Broadleaf.Korean.pine.LAI::Semivariogram.Voronoi.LAI(minx=0,maxx=50,miny=0,maxy=50,boundary =5,b=b,strata=c(10,16,30),r=c(2,2.5,3.5,4),seq = 20)
```
### 2.21 plot.Local.single.point.Voronoi.LAI(minx,maxx,miny,maxy,a,b,r,Lr)

2.21.1���ܽ��ܣ�

����̩ɭ����η����㵥���ֻ򸴲����ڵ�һ�ֲ���ֵ������ⵥһλ�����Ҷ����йصľֲ��ڲ�ṹ����,����ͼ��֤��

2.21.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

a---��Ҫģ��ֲ��ڲ�ṹ������λ������

b---�����е���ľ���ꡢ���֡��ؾ�

r---�����ַ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ������㹹��̩ɭ��������Բ�İ뾶

Lr---��ȡ�ֲ��ڲ�ṹ����ʱ��Բ�İ뾶

2.21.3����

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

2.22.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η����㸴�����ڵ�һ����λ�����Ҷ���ָ����صľֲ��ڲ�ṹ����������ͼ

2.22.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

a---��Ҫģ��ֲ��ڲ�ṹ������λ�������

b---�����е���ľ���ꡢ���֡��ؾ�

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

Lr---ģ��ֲ��ڲ�ṹ������Բ�İ뾶

2.22.3 ����

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

2.23.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η����и�������Ҷ���ָ���������ռ�����ط���

2.23.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���������

b---�����е���ľ���ꡢ���֡��ؾ�

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

indis---��ʼ����

lag---�ͺ������

2.23.3 ����

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

2.24.1���ܽ��ܣ�

���ڷֲ�̩ɭ����η���Ҷ�������ά�ռ�ֲ������Ե��ƻ�ͼ

2.24.2�������壺

minx---���غ�������С��Χ

maxx---���غ��������Χ

miny---������������С��Χ

maxy---�������������Χ

boundary---���������

b---�����е���ľ���ꡢ���֡��ؾ�

seq---���Ƶ���ʱ�Ŀռ�ֱ���

strata---��ͬ�ֲ��ķָ�߶�

r---�����ַָ��ֲ���û����ľ�ֲ���λ��������㣬�Ա�������ľΪ���ĵ�̩ɭ����ι���rΪ��ͬ�ֲ�������㹹��̩ɭ��������Բ�İ뾶

S---�ַ���ƽ�������ʣ��ڳ���ڷ��ı���

theta---��ͼˮƽ��ת�Ƕ�

phi---��ͼ��ֱ��ת�Ƕ�

2.24.3 ����

```{r}
b=Broadleaf.Korean.pine.LAI::b
head(b)
```


```{r,message=FALSE, warning=FALSE}

Broadleaf.Korean.pine.LAI::Voronoi.pointcloud(minx=0, maxx=50, miny=0, maxy=50, boundary=5, b=b, seq=100, strata=c(10,16,30),r=c(2,2.5,3.5,4), S=1.5, theta=120, phi=20)
Broadleaf.Korean.pine.LAI::Voronoi.pointcloud(minx=0, maxx=50, miny=0, maxy=50, boundary=5, b=b, seq=100, strata=c(10,16,30),r=c(2,2.5,3.5,4), S=1.5, theta=0, phi=90)
```

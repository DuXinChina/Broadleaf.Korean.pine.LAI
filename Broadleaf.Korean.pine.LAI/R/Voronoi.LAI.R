
Voronoi.LAI=function(minx,maxx,miny,maxy,boundary,b,r)
{
  library(deldir)####deldir命令计算泰森多边形用
  dis=r/2*sqrt(3)
  ####���ɾ��ȷֲ��� 
  xpoint=seq(minx,maxx,2*dis)
  ypoint=seq(miny,maxy,1.5*r)
  ly=length(ypoint)
  lx=length(xpoint)
  xpoint=rep(xpoint,each=ly)
  ypoint=rep(ypoint,lx)
  uy=unique(ypoint)
  point=cbind(xpoint,ypoint,1:length(xpoint))
  for (i in 1:(0.5*length(uy)))
  {
    point[which(point[,2]==uy[2*i]),1]=point[which(point[,2]==uy[2*i]),1]+dis
  }
  colnames(point)=c("x","y","ID")
  
  ####计算以单一林木为中心，r为半径范围外的全部样�?
  pointoutside.single=function(a)
  {
    ####删除林木a周围半径r以内的样�?
    pointinside=subset(point,point[,1]<a[,1]+r&point[,1]>a[,1]-r&point[,2]<a[,2]+r&point[,2]>a[,2]-r)
    pointinside=as.data.frame(pointinside)
    if(nrow(pointinside)>0)
    {
      c=pointinside[,1:2]
      for (i in 1:nrow(pointinside))
      {
        c[i,]=(pointinside[i,1:2]-a[1,1:2])^2
        d=(c[,1]+c[,2])^(1/2)
      }
      d=cbind(pointinside,d)
      d=subset(d,d>0)
      colnames(d) = c("x","y","ID","Distance")
      d
      pointinside=subset(d,d[,4]<r)
      pointinside=pointinside[,1:3]
      ###计算圈内点的补集�?
      pointoutside_ID=setdiff(point[,3],pointinside[,3])
      pointoutside=point[pointoutside_ID,1:3]
    }
    if(nrow(pointinside)==0)
    {
      pointoutside=point
    }
    pointoutside
  }
  ###以循环语句计算各林木为中心，r为半径范围外的全部样点的交集
  b_for=pointoutside.single(b[1,1:2])[,3]
  for(i in 2:nrow(b))
  {
    bn=pointoutside.single(b[i,1:2])[,3]
    b_for=intersect(b_for,bn)
  }
  point=point[b_for,1:2]
  ####将均匀分布的数据点导为林木数据b的相同形�?
  point=cbind(point,0)
  point=as.data.frame(point)
  point$Species=c("LX")
  colnames(point)=c("x","y","DBH","Species")
  ####将随机分布点与原始数据b结合
  bnew=rbind(point,b)
  
  ####删除样地中重复样�?
  bnew=bnew[deldir(bnew[,1],bnew[,2])$ind.orig,]
  ####计算林木的胸高断面积与泰森多边形的面�?
  BasalA_m2=pi*(bnew[,3]/(2*100))^2
  deldir_area=deldir(bnew[,1],bnew[,2])$summary$dir.area
  ###将泰森多边形面积由平方米转换为公�?
  deldir_area=deldir_area/10000
  ####求个林木胸高断面积与泰森多边形面积的比�?
  BasalA_m2hm=BasalA_m2/deldir_area
  ####将数据bnew合并形成新的bnew
  bnew=cbind(bnew,deldir_area,BasalA_m2,BasalA_m2hm)
  
  ####按照物种样地内的林木分类
  HS=subset(bnew,bnew$Species=="HS")
  LS=subset(bnew,bnew$Species=="LS")
  ZD=subset(bnew,bnew$Species=="ZD")
  SM=subset(bnew,bnew$Species=="SM")
  FH=subset(bnew,bnew$Species=="FH")
  LYY=subset(bnew,bnew$Species=="LYY")
  QT=subset(bnew,bnew$Species=="QT")
  LX=subset(bnew,bnew$Species=="LX")
  ####导入刘志理的叶面积指数方�?
  ####导入刘志理文章中各类树种的叶面积指数方程，HS红松、LS冷杉、ZD紫椴、SM色木、FH枫桦、LYY裂叶榆、QT其他；云杉叶面积使用冷杉方程，其他阔叶树种使用QT方程
  HS$LAI=0.3431*HS$BasalA_m2hm^0.7972
  LS$LAI=0.1995*LS$BasalA_m2hm^0.9539
  ZD$LAI=0.2584*ZD$BasalA_m2hm^0.6361
  SM$LAI=0.4575*SM$BasalA_m2hm^0.5524
  FH$LAI=0.3369*FH$BasalA_m2hm^0.541
  LYY$LAI=0.2743*LYY$BasalA_m2hm^0.6814
  QT$LAI=0.3004*QT$BasalA_m2hm^0.6298
  LX$LAI=0
  ####将各个树种的计算结果合并
  bnew=rbind(HS,LS,ZD,SM,FH,LYY,QT,LX)
  ####删除缓冲区内的样�?
  bnew=subset(bnew,bnew[,1]>(minx+boundary)&bnew[,1]<(maxx-boundary)&bnew[,2]>(miny+boundary)&bnew[,2]<(maxy-boundary))
  bnew
}


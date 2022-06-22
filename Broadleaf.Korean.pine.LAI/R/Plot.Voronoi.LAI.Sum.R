
Plot.Voronoi.LAI.Sum=function(minx,maxx,miny,maxy,boundary,b,strata,r)
{
  
  library(ggplot2)
  library(tcltk)
  ####加载计算泰森多边形叶面积指数的function
  Voronoi.LAI.mult=function(minx,maxx,miny,maxy,boundary,b,r)
  {
    library(deldir)####deldir命令计算泰森多边形用
    ###在样地内均匀布点，以生成正六边形的蜂窝状泰森多边形
    ###基于正六变形外接圆半径r去推算样地中布点行数
    ###基于外接圆半径r去推算正六边形的边心距
    dis=r/2*sqrt(3)
    ####生成均匀分布点
    xpoint=seq(minx,maxx,2*dis)
    ypoint=seq(miny,maxy,1.5*r)
    ly=length(ypoint)
    lx=length(xpoint)
    xpoint=rep(xpoint,each=ly)
    ypoint=rep(ypoint,lx)
    for (i in 1:(0.5*length(xpoint)))
    {
      xpoint[2*i]=xpoint[2*i]+dis
    }
    point=cbind(xpoint,ypoint,1:length(xpoint))
    colnames(point)=c("x","y","ID")
    
    ####计算以单一林木为中心，r为半径范围外的全部样点
    pointoutside.single=function(a)
    {
      ####删除林木a周围半径r以内的样点
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
        ###计算圈内点的补集合
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
    ####将均匀分布的数据点导为林木数据b的相同形式
    point=cbind(point,0)
    point=as.data.frame(point)
    point$Species=c("LX")
    colnames(point)=c("x","y","DBH","Species")
    ####将随机分布点与原始数据b结合
    bnew=rbind(point,b)
    
    ####删除样地中重复样点
    bnew=bnew[deldir(bnew[,1],bnew[,2])$ind.orig,]
    ####计算林木的胸高断面积与泰森多边形的面积
    BasalA_m2=pi*(bnew[,3]/(2*100))^2
    deldir_area=deldir(bnew[,1],bnew[,2])$summary$dir.area
    ###将泰森多边形面积由平方米转换为公顷
    deldir_area=deldir_area/10000
    ####求个林木胸高断面积与泰森多边形面积的比值
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
    ####导入刘志理的叶面积指数方程
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
    ####删除缓冲区内的样点
    bnew=subset(bnew,bnew[,1]>(minx+boundary)&bnew[,1]<(maxx-boundary)&bnew[,2]>(miny+boundary)&bnew[,2]<(maxy-boundary))
    bnew
  }
  
  
  strata=c(0,strata,Inf)
  Forest_strata=list()
  for (i in 1:(length(strata)-1))
  {
    Forest_strata[[i]]=subset(b,b$H>strata[i]&b$H<strata[i+1])
    Forest_strata[[i]]=Voronoi.LAI.mult(minx,maxx,miny,maxy,0,Forest_strata[[i]][,c(1:3,5)],r[i])
  }
  
  ####生成样点以备林分尺度地统计分析
  pointx=seq(minx+boundary,maxx-boundary,length.out=200)
  pointy=seq(miny+boundary,maxy-boundary,length.out=200)
  point=expand.grid(pointx, pointy)
  
  
  
  ##提取第一个林层
  stra_single=Forest_strata[[1]]
  ####对各样点所属泰森多边形进行分类,计算第一个林层各点的叶面积指数
  Lbnew1=subset(stra_single,stra_single[,1]>=(point[1,1]-1.5*r[1])&stra_single[,1]<=(point[1,1]+1.5*r[1])&stra_single[,2]>=(point[1,2]-1.5*r[1])&stra_single[,2]<=(point[1,2]+1.5*r[1]))
  Lbnew1$d=(point[1,1]-Lbnew1[,1])^2+(point[1,2]-Lbnew1[,2])^2
  Lbnew1=Lbnew1[which.min(Lbnew1$d),]
  #Class_point=as.data.frame(Class_point)
  pb=tkProgressBar("林层1","已完成 %", 0, 100)
  for(i in 2:nrow(point))
  {
    Lbnew=subset(stra_single,stra_single[,1]>=(point[i,1]-1.5*r[1])&stra_single[,1]<=(point[i,1]+1.5*r[1])&stra_single[,2]>=(point[i,2]-1.5*r[1])&stra_single[,2]<=(point[i,2]+1.5*r[1]))
    Lbnew$d=(point[i,1]-Lbnew[,1])^2+(point[i,2]-Lbnew[,2])^2
    Lbnew=Lbnew[which.min(Lbnew$d),]
    Lbnew1=rbind(Lbnew1,Lbnew)
    info=sprintf("已完成 %d%%", round(i*100/nrow(point)))  ## 设置进度条的完成度
    setTkProgressBar(pb, i*100/nrow(point), "林层1" ,info)  ## 设置进度条
  }
  Lbnew1=Lbnew1$LAI
  close(pb)
  
  
  ###提取其余林层的叶面积指数
  for (j in 2:length(Forest_strata))
  {
    stra_single=Forest_strata[[j]]
    Lbnewn=subset(stra_single,stra_single[,1]>=(point[1,1]-1.5*r[j])&stra_single[,1]<=(point[1,1]+1.5*r[j])&stra_single[,2]>=(point[1,2]-1.5*r[j])&stra_single[,2]<=(point[1,2]+1.5*r[j]))
    Lbnewn$d=(point[1,1]-Lbnewn[,1])^2+(point[1,2]-Lbnewn[,2])^2
    Lbnewn=Lbnewn[which.min(Lbnewn$d),]
    pb=tkProgressBar(paste("林层",j),"已完成 %", 0, 100)
    
    for(i in 2:nrow(point))
    {
      Lbnew=subset(stra_single,stra_single[,1]>=(point[i,1]-1.5*r[j])&stra_single[,1]<=(point[i,1]+1.5*r[j])&stra_single[,2]>=(point[i,2]-1.5*r[j])&stra_single[,2]<=(point[i,2]+1.5*r[j]))
      Lbnew$d=(point[i,1]-Lbnew[,1])^2+(point[i,2]-Lbnew[,2])^2
      Lbnew=Lbnew[which.min(Lbnew$d),]
      Lbnewn=rbind(Lbnewn,Lbnew)
      info=sprintf("已完成 %d%%", round(i*100/nrow(point)))  ## 设置进度条的完成度
      setTkProgressBar(pb, i*100/nrow(point), paste("林层",j),info)  ## 设置进度条
    }
    Lbnewn=Lbnewn$LAI
    Lbnew1=Lbnew1+Lbnewn
    close(pb)
  }
  point=cbind(point,Lbnew1)
  colnames(point)=c("x","y","LAI")
  point=as.data.frame(point)
  
  ####画图
  ggplot() + geom_raster(data=point, aes(x=x,y=y,fill=LAI))+theme_bw()+
    scale_fill_gradientn(colours =c("white","yellow2","green","green3","green4"))+
    scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  
}




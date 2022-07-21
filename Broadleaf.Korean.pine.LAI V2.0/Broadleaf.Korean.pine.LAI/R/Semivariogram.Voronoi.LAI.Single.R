
Semivariogram.Voronoi.LAI.Single=function(minx,maxx,miny,maxy,boundary,b,r,seq)
{
  library(sp)
  library(gstat)
  library(deldir)
    Voronoi.LAI.mult=function(minx,maxx,miny,maxy,boundary,b,r)
  {
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
    uy=unique(ypoint)
    point=cbind(xpoint,ypoint,1:length(xpoint))
    for (i in 1:(0.5*length(uy)))
    {
      point[which(point[,2]==uy[2*i]),1]=point[which(point[,2]==uy[2*i]),1]+dis
      
    }
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
        pointinside=pointinside[which.min(pointinside[,4]),]
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
    ####计算泰森多边形的面积
    deldir_area=deldir(bnew[,1],bnew[,2])$summary$dir.area
    ####将数据bnew合并形成新的bnew
    bnew=cbind(bnew,deldir_area)
    
    ####按照物种样地内的林木分类
    HS=subset(bnew,bnew$Species=="HS") 
    #红松
    ZD=subset(bnew,bnew$Species=="ZD")
    #紫椴
    KD=subset(bnew,bnew$Species=="KD")
    #糠椴
    MGL=subset(bnew,bnew$Species=="MGL")
    #蒙古栎
    SQL=subset(bnew,bnew$Species=="SQL")
    #水曲柳
    HTQ=subset(bnew,bnew$Species=="HTQ")
    #胡桃楸
    HBL=subset(bnew,bnew$Species=="HBL")
    #黄菠萝
    SMQ=subset(bnew,bnew$Species=="SMQ")
    #色木槭
    QKQ=subset(bnew,bnew$Species=="QKQ")
    #青楷槭
    JSQ=subset(bnew,bnew$Species=="JSQ")
    #假色槭
    NJQ=subset(bnew,bnew$Species=="NJQ")
    #拧筋槭
    BNQ=subset(bnew,bnew$Species=="BNQ")
    #白牛槭
    HKQ=subset(HKQ,HKQ$Species=="HKQ")
    #花楷槭
    CY=subset(bnew,bnew$Species=="CY")
    #春榆
    BH=subset(bnew,bnew$Species=="BH")
    #白桦
    HH=subset(bnew,bnew$Species=="HH")
    #坏槐
    LS=subset(bnew,bnew$Species=="LS")
    #冷杉
    YS=subset(bnew,bnew$Species=="YS")
    #云杉
    FH=subset(bnew,bnew$Species=="FH")
    #枫桦
    LYY=subset(bnew,bnew$Species=="LYY")
    #裂叶榆
    QT=subset(bnew,bnew$Species=="QT")
    #其他
    LX=subset(bnew,bnew$Species=="LX")
    #林隙
    
    ####导入叶面积指数方程
    HS$LAI= (0.04321* HS$DBH^1.831)/(108.786/1000*HS$deldir_area)
    ZD$LAI=(0.02540* ZD$DBH^1.632)/(33.583/1000*ZD$deldir_area)
    KD$LAI=(0.01034* KD$DBH^1.797)/(27.086 /1000*KD$deldir_area)
    MGL$LAI=( 0.05095* MGL$DBH^1.639)/(45.413 /1000*MGL$deldir_area)
    SQL$LAI=( 0.10465* SQL$DBH^1.417)/(57.241 /1000*SQL$deldir_area)
    HTQ$LAI=( 0.01144* HTQ$DBH^1.834)/(54.843  /1000*HTQ$deldir_area)
    HBL$LAI=( 0.01546* HBL$DBH^1.758)/(38.184 /1000*HBL$deldir_area)
    SMQ$LAI=( 0.02938*SMQ$DBH^1.681)/(30.300 /1000*SMQ$deldir_area)
    QKQ$LAI=(0.07652*QKQ$DBH^1.511)/(34.425 /1000*QKQ$deldir_area)
    HKQ$LAI=0.0081*HKQ$DBH^2.3418*37.899/HKQ$deldir_area
    JSQ$LAI=( 0.00933 *JSQ$DBH^1.939)/(23.117 /1000*JSQ$deldir_area)
    NJQ$LAI=( 0.01117 *NJQ$DBH^1.873)/(30.084 /1000*NJQ$deldir_area)
    BNQ$LAI=( 0.04882 *BNQ$DBH^1.597)/(35.058 /1000*BNQ$deldir_area)
    CY$LAI=( 0.04997 *CY$DBH^1.530)/(61.778 /1000*CY$deldir_area)
    BH$LAI=( 0.02491 *BH$DBH^1.913)/(47.379 /1000*BH$deldir_area)
    HH$LAI=( 0.04255 *HH$DBH^1.635)/(38.858 /1000*HH$deldir_area)
    LS$LAI=exp(-4.1730+2.0713*log(LS$DBH))*(7.096)/LS$deldir_area
    YS$LAI=exp(-3.5764+1.9801*log(YS$DBH))*(4.984)/YS$deldir_area
    FH$LAI=0.0081*FH$DBH^2.3418*19.744/FH$deldir_area
    LYY$LAI=0.0081*LYY$DBH^2.3418*30.016/LYY$deldir_area
    QT$LAI=0.0081*QT$DBH^2.3418*26.63/QT$deldir_area
    LX$LAI=0
    
    ####将各个树种的计算结果合并
    bnew=rbind(HS,LS,YS,ZD,KD,MGL,SQL,HTQ,HBL,SMQ,QKQ,HKQ,JSQ,NJQ,BNQ,CY,BH,HH,FH,LYY,QT,LX)
    ####删除缓冲区内的样点
    bnew=subset(bnew,bnew[,1]>(minx+boundary)&bnew[,1]<(maxx-boundary)&bnew[,2]>(miny+boundary)&bnew[,2]<(maxy-boundary))
    bnew
  }
  
  
  
  Voronoi_LAI=Voronoi.LAI.mult(minx,maxx,miny,maxy,0,b,r)
  pointx=seq(minx+boundary,maxx-boundary,length.out=seq+1)
  pointy=seq(miny+boundary,maxy-boundary,length.out=seq+1)
  point=expand.grid(pointx, pointy)
  
  
  minxb=minx+boundary
  maxxb=maxx-boundary
  minyb=miny+boundary
  maxyb=maxy-boundary
  
  xgrid.cen=seq(minxb+0.5*(maxxb-minxb)-0.5/seq*(maxxb-minxb),minxb+0.5*(maxxb-minxb)+0.5/seq*(maxxb-minxb), length.out = 3)
  ygrid.cen=seq(minyb+0.5*(maxyb-minyb)-0.5/seq*(maxyb-minyb),minyb+0.5*(maxyb-minyb)+0.5/seq*(maxyb-minyb), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)
  
  xgrid.left.top=seq(minxb+0.25*(maxxb-minxb)-0.5/seq*(maxxb-minxb),minxb+0.25*(maxxb-minxb)+0.5/seq*(maxxb-minxb), length.out = 3)
  ygrid.left.top=seq(minyb+0.75*(maxyb-minyb)-0.5/seq*(maxyb-minyb),minyb+0.75*(maxyb-minyb)+0.5/seq*(maxyb-minyb), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)
  
  xgrid.left.bottom=seq(minxb+0.25*(maxxb-minxb)-0.5/seq*(maxxb-minxb),minxb+0.25*(maxxb-minxb)+0.5/seq*(maxxb-minxb), length.out = 3)
  ygrid.left.bottom=seq(minyb+0.25*(maxyb-minyb)-0.5/seq*(maxyb-minyb),minyb+0.25*(maxyb-minyb)+0.5/seq*(maxyb-minyb), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)
  
  xgrid.right.top=seq(minxb+0.75*(maxxb-minxb)-0.5/seq*(maxxb-minxb),minxb+0.75*(maxxb-minxb)+0.5/seq*(maxxb-minxb), length.out = 3)
  ygrid.right.top=seq(minyb+0.75*(maxyb-minyb)-0.5/seq*(maxyb-minyb),minyb+0.75*(maxyb-minyb)+0.5/seq*(maxyb-minyb), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)
  
  xgrid.right.bottom=seq(minxb+0.75*(maxxb-minxb)-0.5/seq*(maxxb-minxb),minxb+0.75*(maxxb-minxb)+0.5/seq*(maxxb-minxb), length.out = 3)
  ygrid.right.bottom=seq(minyb+0.25*(maxyb-minyb)-0.5/seq*(maxyb-minyb),minyb+0.25*(maxyb-minyb)+0.5/seq*(maxyb-minyb), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  
  point=rbind(point,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  stra_single=Voronoi_LAI
  Lbnew1=subset(stra_single,stra_single[,1]>=(point[1,1]-1.5*r[1])&stra_single[,1]<=(point[1,1]+1.5*r[1])&stra_single[,2]>=(point[1,2]-1.5*r[1])&stra_single[,2]<=(point[1,2]+1.5*r[1]))
  Lbnew1$d=(point[1,1]-Lbnew1[,1])^2+(point[1,2]-Lbnew1[,2])^2
  Lbnew1=Lbnew1[which.min(Lbnew1$d),]
  for(i in 2:nrow(point))
  {
    Lbnew=subset(stra_single,stra_single[,1]>=(point[i,1]-1.5*r[1])&stra_single[,1]<=(point[i,1]+1.5*r[1])&stra_single[,2]>=(point[i,2]-1.5*r[1])&stra_single[,2]<=(point[i,2]+1.5*r[1]))
    Lbnew$d=(point[i,1]-Lbnew[,1])^2+(point[i,2]-Lbnew[,2])^2
    Lbnew=Lbnew[which.min(Lbnew$d),]
    Lbnew1=rbind(Lbnew1,Lbnew)
  }
  Lbnew1=Lbnew1$LAI
  
  

  point=cbind(point,Lbnew1)
  colnames(point)=c("x","y","LAI")
  point=as.data.frame(point)
  
  data1=point
  coordinates(point) <- c("x","y")
  spplot(point,"LAI")
  vgm1 <- variogram(LAI~1, point)
  plot(vgm1, plot.numbers = TRUE)
  m <- fit.variogram(vgm1,vgm("Sph"))
  print(m)
  p1=plot(vgm1, model=m)
  print(p1)
  sd=subset(vgm1$dist,vgm1$dist<m$range[2])
  spre=m$psill[1]+(m$psill[2])*((3*sd)/(2*m$range[2])-sd^3/(2*(m$range[2])^3))
  bd=subset(vgm1$dist,vgm1$dist>m$range[2])
  bpre=rep((m$psill[1]+m$psill[2]),length(bd))
  pre=rbind(as.matrix(spre),as.matrix(bpre))
  Coefficient_of_Determination=1- sum((pre-vgm1$gamma)^2)/sum((vgm1$gamma-mean(vgm1$gamma))^2)
  print(paste("Coefficient_of_Determination=",Coefficient_of_Determination))
}

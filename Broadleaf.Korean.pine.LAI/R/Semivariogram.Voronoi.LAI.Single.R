
Semivariogram.Voronoi.LAI.Single=function(minx,maxx,miny,maxy,boundary,b,r,seq)
{
  library(sp)
  library(gstat)
  Voronoi.LAI.mult=function(minx,maxx,miny,maxy,boundary,b,r)
  {
    library(deldir)
    dis=r/2*sqrt(3)
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
    pointoutside.single=function(a)
    {
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
        pointoutside_ID=setdiff(point[,3],pointinside[,3])
        pointoutside=point[pointoutside_ID,1:3]
      }
      if(nrow(pointinside)==0)
      {
        pointoutside=point
      }
      pointoutside
    }
    b_for=pointoutside.single(b[1,1:2])[,3]
    for(i in 2:nrow(b))
    {
      bn=pointoutside.single(b[i,1:2])[,3]
      b_for=intersect(b_for,bn)
    }
    point=point[b_for,1:2]
    point=cbind(point,0)
    point=as.data.frame(point)
    point$Species=c("LX")
    colnames(point)=c("x","y","DBH","Species")
    bnew=rbind(point,b)

    bnew=bnew[deldir(bnew[,1],bnew[,2])$ind.orig,]
    BasalA_m2=pi*(bnew[,3]/(2*100))^2
    deldir_area=deldir(bnew[,1],bnew[,2])$summary$dir.area
    deldir_area=deldir_area/10000
    BasalA_m2hm=BasalA_m2/deldir_area
    bnew=cbind(bnew,deldir_area,BasalA_m2,BasalA_m2hm)
    

    HS=subset(bnew,bnew$Species=="HS")
    LS=subset(bnew,bnew$Species=="LS")
    ZD=subset(bnew,bnew$Species=="ZD")
    SM=subset(bnew,bnew$Species=="SM")
    FH=subset(bnew,bnew$Species=="FH")
    LYY=subset(bnew,bnew$Species=="LYY")
    QT=subset(bnew,bnew$Species=="QT")
    LX=subset(bnew,bnew$Species=="LX")
    HS$LAI=0.3431*HS$BasalA_m2hm^0.7972
    LS$LAI=0.1995*LS$BasalA_m2hm^0.9539
    ZD$LAI=0.2584*ZD$BasalA_m2hm^0.6361
    SM$LAI=0.4575*SM$BasalA_m2hm^0.5524
    FH$LAI=0.3369*FH$BasalA_m2hm^0.541
    LYY$LAI=0.2743*LYY$BasalA_m2hm^0.6814
    QT$LAI=0.3004*QT$BasalA_m2hm^0.6298
    LX$LAI=0
    bnew=rbind(HS,LS,ZD,SM,FH,LYY,QT,LX)
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


Plot.Voronoi.LAI=function(minx,maxx,miny,maxy,boundary,b,r)
{
  library(ggplot2)
  library(grDevices)
  library(deldir)
  ####基于泰森多边形计算样地内不同斑块的叶面积指数
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
  bnew=Voronoi.LAI.mult(minx,maxx,miny,maxy,0,b,r)
  bnew$group=1:nrow(bnew)
  
  
  a=deldir(bnew[,1],bnew[,2])
  a=a$dirsgs
  a1=cbind(a[,1:2],a[,5])
  a2=cbind(a[,3:4],a[,5])
  colnames(a1)=c("x","y","group")
  colnames(a2)=c("x","y","group")
  a12=rbind(a1,a2)
  b1=cbind(a[,1:2],a[,6])
  b2=cbind(a[,3:4],a[,6])
  colnames(b1)=c("x","y","group")
  colnames(b2)=c("x","y","group")
  b12=rbind(b1,b2)
  ab12=rbind(a12,b12)
  #ab12=ab12[order(ab12[,1]),]
  ab12=ab12[order(ab12[,3]),]
  
  ab12_1=subset(ab12,ab12[,3]==1)
  order=chull(ab12_1[,1],ab12_1[,2])
  ab12_1=ab12_1[order,]  
  for(i in 2:max(ab12[,3]))
  {
    ab12_n=subset(ab12,ab12[,3]==i)
    order=chull(ab12_n[,1],ab12_n[,2])
    ab12_n=ab12_n[order,]
    ab12_1=rbind(ab12_1,ab12_n)
  }
  
  for(i in 1:nrow(ab12_1))
  {
    ab12_1[i,4]=bnew[which(bnew$group==ab12_1[i,3]),8]
  }
  colnames(ab12_1)=c("x","y","group","LAI")
  ab12_1[which(ab12_1[,4]==0),4]=0.001
  expandx=(maxx-minx)*0.095
  expandy=(maxy-miny)*0.095
  aex=1.053*expandx
  aey=1.053*expandy
  x1=c(minx-aex,minx-aex,minx+boundary,minx+boundary)
  y1=c(miny-aey,maxy+aey,maxy+aey,miny-aey)
  xy1=cbind(x1,y1,1)
  xy1=as.data.frame(xy1)
  colnames(xy1)=c("x","y","group")
  x2=c(minx-aex,minx-aex,maxx+aex,maxx+aex)
  y2=c(miny-aey,miny+boundary,miny+boundary,miny-aey)
  xy2=cbind(x2,y2,2)
  xy2=as.data.frame(xy2)
  colnames(xy2)=c("x","y","group")
  x3=c(maxx-boundary,maxx-boundary,maxx+aex,maxx+aex)
  y3=c(miny-aey,maxy+aex,maxy+aex,miny-aey)
  xy3=cbind(x3,y3,3)
  xy3=as.data.frame(xy3)
  colnames(xy3)=c("x","y","group")
  x4=c(minx-aex,minx-aex,maxx+aex,maxx+aex)
  y4=c(maxy-boundary,maxy+aey,maxy+aey,maxy-boundary)
  xy4=cbind(x4,y4,4)
  xy4=as.data.frame(xy4)
  colnames(xy4)=c("x","y","group")
  xy=rbind(xy1,xy2,xy3,xy4)
  ggplot()+
    geom_polygon(data=ab12_1,aes(x=x,y=y,group=group,fill=LAI),colour="black")+
    scale_fill_gradientn(colours =c("white","yellow2","green","green3","green4"))+
    geom_polygon(data=xy,aes(x=x,y=y,group=group),colour="white",fill="white")+
    #xlim(0,50)+ylim(0,50)+
    scale_x_continuous(expand= c(0, -expandx))+scale_y_continuous(expand= c(0, -expandy))+
    geom_vline(xintercept = c(minx+boundary,maxx-boundary),linetype=2,size=1)+geom_hline(yintercept =c(miny+boundary,maxy-boundary),linetype=2,size=1)+
    theme_bw()
  
}



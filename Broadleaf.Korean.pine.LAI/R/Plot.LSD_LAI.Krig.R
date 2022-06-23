
###绘制林分叶面积指数局域标准差的克里格插值图
Plot.LSD_LAI.Krig=function(minx,maxx,miny,maxy,b,seq,r,Lr)
{
  library(sp)
  library(gstat)
  library(tcltk)  
  library(ggplot2)
  
  ###计算多个位点的叶面积的局部标准差
  LSD_LAI_mult=function(a,b,r,Lr)
  {
    library(sp)
    library(gstat)
    library(tcltk)  
    library(ggplot2)
    ####计算单一位点叶面积的局部标准差
    LSD_LAI=function(a,b,r,Lr)
    {
      ###以a点为中心划定半径为Lr的正方形范围
      aminx=a[,1]-Lr
      amaxx=a[,1]+Lr
      aminy=a[,2]-Lr
      amaxy=a[,2]+Lr
      
      ###在方形范围内均匀生成坐标点
      ax=seq(aminx,amaxx,length.out=5)
      ay=seq(aminy,amaxy,length.out=5)
      ax=rep(ax,5)
      ay=rep(ay,each=5)
      acenterpoint=cbind(ax,ay)
      acenterpoint=as.data.frame(acenterpoint)
      ###筛选正方形点阵中与中心a点距离小于Lr的样点
      Neighbourhood.single=function(a,acenterpoint,Lr)
      {
        ###计算正方形范围内到各样点到中心点a的距离
        
        c=acenterpoint
        for (i in 1:nrow(acenterpoint))
        {
          c[i,]=(acenterpoint[i,1:2]-a[1,1:2])^2
          d=(c[,1]+c[,2])^(1/2)
        }
        d=cbind(acenterpoint,d)
        d=subset(d,d>0)
        colnames(d) = c("x","y","Distance")####x,y为样点坐标，Distance为样点与中心点a间的距离
        d
        
        ###绘制半径为r的尺度圆，以找出样地内位于以a为圆心Lr为半径尺度圆内的样点
        Neighbourhood.single=subset(d,d$Distance<=Lr)
        Neighbourhood.single
      }
      
      ###生成局域位点阵
      Local_point=rbind(a,Neighbourhood.single(a,acenterpoint,Lr)[,1:2])
      
      ###计算点阵内各位点的叶面积指数
      LAI.mult=function(a,b,r)
      {
        
        ####计算以随机选择的位点为中心，一定半径尺度圆内叶面积指数
        LAI.single=function(a,b,r)
        {
          ###找出以a点为中心，r为半径的尺度圆内的林木
          ###数据a为某一样点，第一列x，第二列y；数据b为林分中各林木的分布坐标，胸径与种名，第一列x，第二列y，第三列DBH，第四列Species；r为尺度圆半径
          Neighbourhood.single=function(a,b,r)
          {
            ###计算样地内林木到a点的距离
            Neighbourhood.single1=function(a,b)###计算样地内林木到a点的距离
            {
              c=b[,1:2]
              for (i in 1:nrow(b))
              {
                c[i,]=(b[i,1:2]-a[1,1:2])^2
                d=(c[,1]+c[,2])^(1/2)
              }
              d=cbind(b,d)
              d=subset(d,d>0)
              colnames(d) = c("x","y","DBH","Species","Distance")####x,y为林木坐标Species为树种,DBH为胸径，Distance为林木与中心点a1间的距离
              d
            }
            
            ###绘制半径为r的尺度圆，以找出样地内位于以a为圆心r为半径尺度圆内的林木
            Neighbourhood.single1=Neighbourhood.single1(a,b)
            Neighbourhood.single=subset(Neighbourhood.single1,Neighbourhood.single1$Distance<r)
            Neighbourhood.single
          }
          
          ####筛选尺度圆内的植株
          Neighbourhood.single=Neighbourhood.single(a,b,r)
          
          ####计算树种胸高断面积与尺度圆面积的比值
          ####计算尺度圆面积
          scale_circleA=pi*r^2
          ####将尺度圆面积由平方米换算为公顷
          scale_circleA=scale_circleA/10000
          ####求尺度圆内算各树种的胸高断面积
          BasalA=pi*(Neighbourhood.single[,3]/(2*100))^2
          ####求比值
          BasalA=BasalA/scale_circleA
          e=cbind(Neighbourhood.single,BasalA)
          colnames(e) = c("x","y","DBH","Species","Distance","Basal_Area")
          
          ####按照物种将尺度圆内的林木分类
          HS=subset(e,e$Species=="HS")
          LS=subset(e,e$Species=="LS")
          ZD=subset(e,e$Species=="ZD")
          SM=subset(e,e$Species=="SM")
          FH=subset(e,e$Species=="FH")
          LYY=subset(e,e$Species=="LYY")
          QT=subset(e,e$Species=="QT")
          
          ####求尺度圆内各个树种的胸高断面积
          HSBA=sum(HS[,6])
          LSBA=sum(LS[,6])
          ZDBA=sum(ZD[,6])
          SMBA=sum(SM[,6])
          FHBA=sum(FH[,6])
          LYYBA=sum(LYY[,6])
          QTBA=sum(QT[,6])
          
          ####导入刘志理文章中各类树种的叶面积指数方程，HS红松、LS冷杉、ZD紫椴、SM色木、FH枫桦、LYY裂叶榆、QT其他；云杉叶面积使用冷杉方程，其他阔叶树种使用QT方程
          HS_LAI=0.3431*HSBA^0.7972
          LS_LAI=0.1995*LSBA^0.9539
          ZD_LAI=0.2584*ZDBA^0.6361
          SM_LAI=0.4575*SMBA^0.5524
          FH_LAI=0.3369*FHBA^0.541
          LYY_LAI=0.2743*LYYBA^0.6814
          QT_LAI=0.3004*QTBA^0.6298
          
          ###计算尺度圆内的总叶面积
          LAI=sum(HS_LAI,LS_LAI,ZD_LAI,SM_LAI,FH_LAI,LYY_LAI,QT_LAI)
          
          LAI
        }
        
        
        d=matrix(NA,nrow(a),3)
        
        for(j in 1:nrow(a))
        {
          d[j,]=cbind(as.matrix(a[j,]),as.matrix(LAI.single(a[j,],b,r)))
        }
        colnames(d)=c("x","y","LAI")
        rownames(d)=1:nrow(a)
        d=as.data.frame(d)
        d
      }
      Lsd=cbind(a,sd(LAI.mult(Local_point,b,r)[,3]))
      colnames(Lsd)=c("x","y","LSD_LAI")
      Lsd
    }
    
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("Progress","Percent complete %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(LSD_LAI(a[j,],b,r,Lr)$LSD_LAI))
      info=sprintf("Percent complete %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("Progress (%s)", info),info)  ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)  
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","LSD_LAI")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  
  ####生成插值点
  xgrid=seq(minx,maxx, length.out = seq+1)
  ygrid=seq(miny,maxy, length.out = seq+1)
  basexy=expand.grid(xgrid, ygrid)
  xgrid.cen=seq(minx+0.5*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.5*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.cen=seq(miny+0.5*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.5*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.cen=expand.grid(xgrid.cen, ygrid.cen)
  
  xgrid.left.top=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.top=expand.grid(xgrid.left.top, ygrid.left.top)
  
  xgrid.left.bottom=seq(minx+0.25*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.25*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.left.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.left.bottom=expand.grid(xgrid.left.bottom, ygrid.left.bottom)
  
  xgrid.right.top=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.top=seq(miny+0.75*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.75*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.top=expand.grid(xgrid.right.top, ygrid.right.top)
  
  xgrid.right.bottom=seq(minx+0.75*(maxx-minx)-0.5/seq*(maxx-minx),minx+0.75*(maxx-minx)+0.5/seq*(maxx-minx), length.out = 3)
  ygrid.right.bottom=seq(miny+0.25*(maxy-miny)-0.5/seq*(maxy-miny),miny+0.25*(maxy-miny)+0.5/seq*(maxy-miny), length.out = 3)
  basexy.right.bottom=expand.grid(xgrid.right.bottom, ygrid.right.bottom)
  #basexy1=basexy
  xgrid=seq(minx,maxx, length.out = 200)
  ygrid=seq(miny,maxy, length.out = 200)
  basexy1=expand.grid(xgrid, ygrid)
  Basexy1=basexy1
  basexy=rbind(basexy,basexy.cen,basexy.left.top,basexy.left.bottom,basexy.right.top,basexy.right.bottom)
  colnames(basexy1) <- c("x", "y")
  coordinates(basexy1) <- ~x+y
  gridded(basexy1) <- TRUE
  colnames(basexy) <- c("x", "y")
  basexy= dplyr::distinct(basexy)####删除重复样点
  
  ######计算插值点的LSD_LAI
  data=LSD_LAI_mult(basexy,b,r,Lr)
  data1=data
  coordinates(data) <- c("x","y")#定义坐标 
  spplot(data,"LSD_LAI")
  vgm1 <- variogram(LSD_LAI~1, data)
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
  krige_res <- krige(LSD_LAI~1, data, basexy1, model = m)
  
  ### 查看克里格插值的结果
  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours =c("red2","orange2","yellow2","lightgreen","green","green3"))
  p2=p2+labs(title = "LSD_LAI")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}

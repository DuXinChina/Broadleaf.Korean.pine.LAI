
####绘制林分的叶面积指数克里格插值图
Plot.LAI.Krig=function(minx,maxx,miny,maxy,b,seq,r)
{
  library(sp)
  library(gstat)
  library(tcltk)  
  library(ggplot2)
  ####计算以随机选择的位点为中心，一定半径尺度圆内叶面积指数
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
      colnames(d) = c("x","y","DBH","Species","Distance")
      ####x,y为林木坐标Species为树种,DBH为胸径，Distance为林木与中心点a1间的距离
      d
    }
    
    ###绘制半径为r的尺度圆，以找出样地内位于以a为圆心r为半径尺度圆内的林木
    Neighbourhood.single1=Neighbourhood.single1(a,b)
    Neighbourhood.single=subset(Neighbourhood.single1,Neighbourhood.single1$Distance<r)
    Neighbourhood.single
  }
  
  ####筛选尺度圆内的植株
  Neighbourhood.single=Neighbourhood.single(a,b,r)
  
    ####计算尺度圆面积
    scale_circleA=pi*r^2
     
   # colnames(e) = c("x","y","DBH","Species","Distance","Basal_Area")
    
    ####按照物种将尺度圆内的林木分类
    HS=subset(Neighbourhood.single,Neighbourhood.single$Species=="HS") 
#红松
    ZD=subset(Neighbourhood.single,Neighbourhood.single$Species=="ZD")
#紫椴
    KD=subset(Neighbourhood.single,Neighbourhood.single$Species=="KD")
#糠椴
    MGL=subset(Neighbourhood.single,Neighbourhood.single$Species=="MGL")
#蒙古栎
    SQL=subset(Neighbourhood.single,Neighbourhood.single$Species=="SQL")
#水曲柳
    HTQ=subset(Neighbourhood.single,Neighbourhood.single$Species=="HTQ")
#胡桃楸
    HBL=subset(Neighbourhood.single,Neighbourhood.single$Species=="HBL")
#黄菠萝
    SMQ=subset(Neighbourhood.single,Neighbourhood.single$Species=="SMQ")
#色木槭
    QKQ=subset(Neighbourhood.single,Neighbourhood.single$Species=="QKQ")
#青楷槭
    JSQ=subset(Neighbourhood.single,Neighbourhood.single$Species=="JSQ")
#假色槭
    NJQ=subset(Neighbourhood.single,Neighbourhood.single$Species=="NJQ")
#拧筋槭
    BNQ=subset(Neighbourhood.single,Neighbourhood.single$Species=="BNQ")
#白牛槭
    HKQ=subset(Neighbourhood.single,Neighbourhood.single$Species=="HKQ")
#花楷槭
    CY=subset(Neighbourhood.single,Neighbourhood.single$Species=="CY")
#春榆
    BH=subset(Neighbourhood.single,Neighbourhood.single$Species=="BH")
#白桦
    HH=subset(Neighbourhood.single,Neighbourhood.single$Species=="HH")
#坏槐
    LS=subset(Neighbourhood.single,Neighbourhood.single$Species=="LS")
#冷杉
    YS=subset(Neighbourhood.single,Neighbourhood.single$Species=="YS")
#云杉
    FH=subset(Neighbourhood.single,Neighbourhood.single$Species=="FH")
#枫桦
    LYY=subset(Neighbourhood.single,Neighbourhood.single$Species=="LYY")
#裂叶榆
    QT=subset(Neighbourhood.single,Neighbourhood.single$Species=="QT")
    

    
    ####导入叶面积指数方程
    HS_LAI= sum((0.04321* HS$DBH^1.831)/(108.786/1000*scale_circleA))
    ZD_LAI=sum((0.02540* ZD$DBH^1.632)/(33.583/1000*scale_circleA))
    KD_LAI=sum((0.01034* KD$DBH^1.797)/(27.086 /1000*scale_circleA))
    MGL_LAI=sum(( 0.05095* MGL$DBH^1.639)/(45.413 /1000*scale_circleA))
    SQL_LAI=sum(( 0.10465* SQL$DBH^1.417)/(57.241 /1000*scale_circleA))
    HTQ_LAI=sum(( 0.01144* HTQ$DBH^1.834)/(54.843  /1000*scale_circleA))
    HBL_LAI=sum(( 0.01546* HBL$DBH^1.758)/(38.184 /1000*scale_circleA))
    SMQ_LAI=sum(( 0.02938*SMQ$DBH^1.681)/(30.300 /1000*scale_circleA))
    QKQ_LAI=sum((0.07652*QKQ$DBH^1.511)/(34.425 /1000*scale_circleA))
    HKQ_LAI=sum(0.0081*HKQ$DBH^2.3418*37.899/scale_circleA)
    JSQ_LAI=sum(( 0.00933 *JSQ$DBH^1.939)/(23.117 /1000*scale_circleA))
    NJQ_LAI=sum(( 0.01117 *NJQ$DBH^1.873)/(30.084 /1000*scale_circleA))
    BNQ_LAI=sum(( 0.04882 *BNQ$DBH^1.597)/(35.058 /1000*scale_circleA))
    CY_LAI=sum(( 0.04997 *CY$DBH^1.530)/(61.778 /1000*scale_circleA))
    BH_LAI=sum(( 0.02491 *BH$DBH^1.913)/(47.379 /1000*scale_circleA))
    HH_LAI=sum(( 0.04255 *HH$DBH^1.635)/(38.858 /1000*scale_circleA))
    LS_LAI=sum(exp(-4.1730+2.0713*log(LS$DBH))*(7.096)/scale_circleA)
    YS_LAI=sum(exp(-3.5764+1.9801*log(YS$DBH))*(4.984)/scale_circleA)
    FH_LAI=sum(0.0081*FH$DBH^2.3418*19.744/scale_circleA)
    LYY_LAI=sum(0.0081*LYY$DBH^2.3418*30.016/scale_circleA)
    QT_LAI=sum(0.0081*QT$DBH^2.3418*26.63/scale_circleA)

    
    ###计算尺度圆内的总叶面积
    LAI=sum(HS_LAI,LS_LAI,YS_LAI,ZD_LAI,KD_LAI,MGL_LAI,SQL_LAI,HTQ_LAI,HBL_LAI,SMQ_LAI,QKQ_LAI,HKQ_LAI,JSQ_LAI,NJQ_LAI,BNQ_LAI,CY_LAI,BH_LAI,HH_LAI,FH_LAI,LYY_LAI,QT_LAI)
    
    ###计算尺度圆内针叶与阔叶树的叶面积
    Needles_LAI=sum(HS_LAI,LS_LAI,YS_LAI)
    Broadleaf_LAI=sum(ZD_LAI,KD_LAI,MGL_LAI,SQL_LAI,HTQ_LAI,HBL_LAI,SMQ_LAI,QKQ_LAI,JSQ_LAI,NJQ_LAI,BNQ_LAI,CY_LAI,BH_LAI,HH_LAI,FH_LAI,LYY_LAI,QT_LAI)
    
    ###计算尺度圆内针叶树叶面积与阔叶树叶面积所占的比例
    N_L_percent=Needles_LAI/LAI*100
    B_L_percent=Broadleaf_LAI/LAI*100
    
    ####整理计算结果
    Species_LAI=c(HS_LAI,LS_LAI,YS_LAI,ZD_LAI,KD_LAI,MGL_LAI,SQL_LAI,HTQ_LAI,HBL_LAI,SMQ_LAI,QKQ_LAI,JSQ_LAI,NJQ_LAI,BNQ_LAI,CY_LAI,BH_LAI,HH_LAI,FH_LAI,LYY_LAI,QT_LAI)
    Species_LAI=matrix(Species_LAI,1,)
    colnames(Species_LAI)=c("HS_LAI","LS_LAI","YS_LAI","ZD_LAI","KD_LAI","MGL_LAI","SQL_LAI","HTQ_LAI","HBL_LAI","SMQ_LAI","QKQ_LAI","JSQ_LAI","NJQ_LAI","BNQ_LAI","CY_LAI","BH_LAI","HH_LAI","FH_LAI","LYY_LAI","QT_LAI")
  N_B_LAI=c(Needles_LAI,Broadleaf_LAI,N_L_percent,B_L_percent)
  N_B_LAI=matrix(N_B_LAI,1,)
  colnames(N_B_LAI)=c("Needles_LAI","Broadleaf_LAI","N_L_percent","B_L_percent")
  list(LAI=LAI,Species_LAI=Species_LAI,N_B_LAI=N_B_LAI)
}

  
  ####以循环语句计算多个点的叶面积指数，并给出计算进度条
  LAI.mult=function(a,b,r)
  {
    d=matrix(NA,nrow(a),3)
    pb=tkProgressBar("Progress","Percent complete  %", 0, 100)
    star_time=Sys.time()
     ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(LAI.single(a[j,],b,r)$LAI))
      info=sprintf("Percent complete %d%%", round(j*100/nrow(a)))
       ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("Progress (%s)", info),info)
       ## 设置进度条
    }
    end_time=Sys.time()
     ## 记录程序结束时间
    close(pb)  
    run_time=end_time - star_time
    ## 计算程序运行时间
    colnames(d)=c("x","y","LAI")
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
  basexy= dplyr::distinct(basexy)
  ####删除重复样点
  
  ######计算插值点的LAI
  data=LAI.mult(basexy,b,r)
  data1=data
  coordinates(data) <- c("x","y")
  #定义坐标 
  spplot(data,"LAI")
  vgm1 <- variogram(LAI~1, data)
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
  krige_res <- krige(LAI~1, data, basexy1, model = m)
  
  ### 查看克里格插值的结果
  z=krige_res$var1.pred
  Basexyz=cbind(Basexy1,z)
  colnames(Basexyz)=c("x","y","Vaule")
  p2=ggplot() + geom_raster(data=Basexyz, aes(x=x,y=y,fill=Vaule))+theme_bw()+ scale_fill_gradientn(colours =c("white","green1","green2","green3","green4"))
  p2=p2+labs(title = "LAI")+scale_x_continuous(expand= c(0, 0))+scale_y_continuous(expand= c(0, 0))
  print(p2)
}
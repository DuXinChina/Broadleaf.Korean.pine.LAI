
LSD_LAI=function(a,b,r,Lr)
{
  aminx=a[,1]-Lr
  amaxx=a[,1]+Lr
  aminy=a[,2]-Lr
  amaxy=a[,2]+Lr
  ax=seq(aminx,amaxx,length.out=5)
  ay=seq(aminy,amaxy,length.out=5)
  ax=rep(ax,5)
  ay=rep(ay,each=5)
  acenterpoint=cbind(ax,ay)
  acenterpoint=as.data.frame(acenterpoint)
  Neighbourhood.single=function(a,acenterpoint,Lr)
  {
    c=acenterpoint
    for (i in 1:nrow(acenterpoint))
    {
      c[i,]=(acenterpoint[i,1:2]-a[1,1:2])^2
      d=(c[,1]+c[,2])^(1/2)
    }
    d=cbind(acenterpoint,d)
    d=subset(d,d>0)
    colnames(d) = c("x","y","Distance")
    ####x,y为样点坐标，Distance为样点与中心点a间的距离
    d
    Neighbourhood.single=subset(d,d$Distance<=Lr)
    Neighbourhood.single
  }
  
  Local_point=rbind(a,Neighbourhood.single(a,acenterpoint,Lr)[,1:2])
  
  LAI.mult=function(a,b,r)
  {
    library(sp)
    library(gstat)
    library(tcltk)  
    library(ggplot2)
   
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
      
LAI
    }  
    
    
    
    
    d=matrix(NA,nrow(a),3)
    e=matrix(NA,nrow(a),6)
    pb=tkProgressBar("Progress","Percent complete %", 0, 100)
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(LAI.single(a[j,],b,r)))
      info=sprintf("Percent complete %d%%", round(j*100/nrow(a)))
      ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("Progress (%s)", info),info)
      ## 设置进度条
    }
    end_time=Sys.time()  ## 记录程序结束时间
    close(pb)  
    run_time=end_time - star_time  ## 计算程序运行时间
    colnames(d)=c("x","y","LAI")
    rownames(d)=1:nrow(a)
    d=as.data.frame(d)
    d
  }
  Lsd=cbind(a,sd(LAI.mult(Local_point,b,r)[,3]))
  colnames(Lsd)=c("x","y","LSD_LAI")
  Lsd
}
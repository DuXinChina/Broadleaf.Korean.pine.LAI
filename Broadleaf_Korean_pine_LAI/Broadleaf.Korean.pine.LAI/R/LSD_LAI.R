
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
    colnames(d) = c("x","y","Distance")####x,y为样点坐标，Distance为样点与中心点a间的距离
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
      Neighbourhood.single=function(a,b,r)
      {
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
        Neighbourhood.single1=Neighbourhood.single1(a,b)
        Neighbourhood.single=subset(Neighbourhood.single1,Neighbourhood.single1$Distance<r)
        Neighbourhood.single
      }
   
      Neighbourhood.single=Neighbourhood.single(a,b,r)
    
      scale_circleA=pi*r^2
      scale_circleA=scale_circleA/10000
      BasalA=pi*(Neighbourhood.single[,3]/(2*100))^2
      BasalA=BasalA/scale_circleA
      e=cbind(Neighbourhood.single,BasalA)
      colnames(e) = c("x","y","DBH","Species","Distance","Basal_Area")
      HS=subset(e,e$Species=="HS")
      LS=subset(e,e$Species=="LS")
      ZD=subset(e,e$Species=="ZD")
      SM=subset(e,e$Species=="SM")
      FH=subset(e,e$Species=="FH")
      LYY=subset(e,e$Species=="LYY")
      QT=subset(e,e$Species=="QT")
      
      HSBA=sum(HS[,6])
      LSBA=sum(LS[,6])
      ZDBA=sum(ZD[,6])
      SMBA=sum(SM[,6])
      FHBA=sum(FH[,6])
      LYYBA=sum(LYY[,6])
      QTBA=sum(QT[,6])
      HS_LAI=0.3431*HSBA^0.7972
      LS_LAI=0.1995*LSBA^0.9539
      ZD_LAI=0.2584*ZDBA^0.6361
      SM_LAI=0.4575*SMBA^0.5524
      FH_LAI=0.3369*FHBA^0.541
      LYY_LAI=0.2743*LYYBA^0.6814
      QT_LAI=0.3004*QTBA^0.6298
      
      LAI=sum(HS_LAI,LS_LAI,ZD_LAI,SM_LAI,FH_LAI,LYY_LAI,QT_LAI)
      
      LAI
    }
    d=matrix(NA,nrow(a),3)
    e=matrix(NA,nrow(a),6)
    pb=tkProgressBar("进度","已完成 %", 0, 100) 
    star_time=Sys.time() ## 记录程序开始时间
    for(j in 1:nrow(a))
    {
      d[j,]=cbind(as.matrix(a[j,]),as.matrix(LAI.single(a[j,],b,r)))
      info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
      setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
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



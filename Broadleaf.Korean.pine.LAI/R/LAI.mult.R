

####以循环语句计算多个点的叶面积指数，并给出计算进度条
LAI.mult=function(a,b,r)
{
  library(sp)
  library(gstat)
  library(tcltk)  
  library(ggplot2)
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
    
    ###计算尺度圆内针叶与阔叶树的叶面积
    Needles_LAI=sum(HS_LAI,LS_LAI)
    Broadleaf_LAI=sum(ZD_LAI,SM_LAI,FH_LAI,LYY_LAI,QT_LAI)
    
    ###计算尺度圆内针叶树叶面积与阔叶树叶面积所占的比例
    N_L_percent=Needles_LAI/LAI*100
    B_L_percent=Broadleaf_LAI/LAI*100
    
    ####整理计算结果
    Species_LAI=c(HS_LAI,LS_LAI,ZD_LAI,SM_LAI,FH_LAI,LYY_LAI,QT_LAI)
    Species_LAI=matrix(Species_LAI,1,)
    colnames(Species_LAI)=c("HS_LAI","LS_LAI","ZD_LAI","SM_LAI","FH_LAI","LYY_LAI","QT_LAI")
    N_B_LAI=c(Needles_LAI,Broadleaf_LAI,N_L_percent,B_L_percent)
    N_B_LAI=matrix(N_B_LAI,1,)
    colnames(N_B_LAI)=c("Needles_LAI","Broadleaf_LAI","N_L_percent","B_L_percent")
    list(LAI=LAI,Species_LAI=Species_LAI,N_B_LAI=N_B_LAI)
  }
  d=matrix(NA,nrow(a),3)
  e=matrix(NA,nrow(a),6)
  pb=tkProgressBar("进度","已完成 %", 0, 100) 
  star_time=Sys.time() ## 记录程序开始时间
  for(j in 1:nrow(a))
  {
    d[j,]=cbind(as.matrix(a[j,]),as.matrix(LAI.single(a[j,],b,r)$LAI))
    e[j,]=cbind(as.matrix(a[j,]),as.matrix(LAI.single(a[j,],b,r)$N_B_LAI))
    info=sprintf("已完成 %d%%", round(j*100/nrow(a)))  ## 设置进度条的完成度
    setTkProgressBar(pb, j*100/nrow(a), sprintf("进度 (%s)", info),info)  ## 设置进度条
  }
  end_time=Sys.time()  ## 记录程序结束时间
  close(pb)  
  run_time=end_time - star_time  ## 计算程序运行时间
  colnames(d)=c("x","y","LAI")
  rownames(d)=1:nrow(a)
  d=as.data.frame(d)
  colnames(e)=c("x","y","Needles_LAI","Broadleaf_LAI","N_L_percent","B_L_percent")
  rownames(e)=1:nrow(a)  
  list(LAI=d,B_N_LAI=e)
  
}


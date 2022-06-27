
####计算林分尺度的叶面积指数
Standing_Forest_LAI=function(minx,maxx,miny,maxy,b)
{
  ###计算样地面积
  xlength=abs(maxx-minx)
  ylength=abs(maxy-miny)
  plotarea=xlength*ylength
  ###将样地面积换算为公顷
  plotarea=plotarea/10000
  
  ###提取选定样地范围内的林木
  b=subset(b,b[,1]>minx&b[,1]<maxx&b[,2]>miny&b[,2]<maxy)
  
  ####求选定的样地范围内各树种的胸高断面积
  BasalA=pi*(b[,3]/(2*100))^2
  ####求比值
  BasalA=BasalA/plotarea
  e=cbind(b,BasalA)
  colnames(e) = c("x","y","DBH","Species","Basal_Area")
  
  ####按照物种将划定的样地范围内内的林木分类
  HS=subset(e,e$Species=="HS")
  LS=subset(e,e$Species=="LS")
  ZD=subset(e,e$Species=="ZD")
  SM=subset(e,e$Species=="SM")
  FH=subset(e,e$Species=="FH")
  LYY=subset(e,e$Species=="LYY")
  QT=subset(e,e$Species=="QT")
  
  ####求尺度圆内各个树种的胸高断面积
  HSBA=sum(HS[,5])
  LSBA=sum(LS[,5])
  ZDBA=sum(ZD[,5])
  SMBA=sum(SM[,5])
  FHBA=sum(FH[,5])
  LYYBA=sum(LYY[,5])
  QTBA=sum(QT[,5])
  
  ####导入刘志理文章中各类树种的叶面积指数方程，HS红松、LS冷杉、ZD紫椴、SM色木、FH枫桦、LYY裂叶榆、QT其他；云杉叶面积使用冷杉方程，其他阔叶树种使用QT方程
  HS_LAI=0.3431*HSBA^0.7972
  LS_LAI=0.1995*LSBA^0.9539
  ZD_LAI=0.2584*ZDBA^0.6361
  SM_LAI=0.4575*SMBA^0.5524
  FH_LAI=0.3369*FHBA^0.541
  LYY_LAI=0.2743*LYYBA^0.6814
  QT_LAI=0.3004*QTBA^0.6298
  
  ###计算划定样地内的总叶面积指数
  LAI=sum(HS_LAI,LS_LAI,ZD_LAI,SM_LAI,FH_LAI,LYY_LAI,QT_LAI)
  
  ###计算划定样地内针叶与阔叶树的叶面积指数
  Needles_LAI=sum(HS_LAI,LS_LAI)
  Broadleaf_LAI=sum(ZD_LAI,SM_LAI,FH_LAI,LYY_LAI,QT_LAI)
  
  ###计算划定样地内针叶树叶面积与阔叶树叶面积所占的比例
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


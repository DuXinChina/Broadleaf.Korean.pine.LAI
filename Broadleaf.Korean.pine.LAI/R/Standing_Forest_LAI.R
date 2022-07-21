
####计算林分尺度的叶面积指数
Standing_Forest_LAI=function(minx,maxx,miny,maxy,b)
{
  ###计算样地面积
  xlength=abs(maxx-minx)
  ylength=abs(maxy-miny)
  plotarea=xlength*ylength
 
  ###提取选定样地范围内的林木
  b=subset(b,b[,1]>minx&b[,1]<maxx&b[,2]>miny&b[,2]<maxy)

Neighbourhood.single=b  
    
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
    HS_LAI= sum((0.04321* HS$DBH^1.831)/(108.786/1000*plotarea))
    ZD_LAI=sum((0.02540* ZD$DBH^1.632)/(33.583/1000*plotarea))
    KD_LAI=sum((0.01034* KD$DBH^1.797)/(27.086 /1000*plotarea))
    MGL_LAI=sum(( 0.05095* MGL$DBH^1.639)/(45.413 /1000*plotarea))
    SQL_LAI=sum(( 0.10465* SQL$DBH^1.417)/(57.241 /1000*plotarea))
    HTQ_LAI=sum(( 0.01144* HTQ$DBH^1.834)/(54.843  /1000*plotarea))
    HBL_LAI=sum(( 0.01546* HBL$DBH^1.758)/(38.184 /1000*plotarea))
    SMQ_LAI=sum(( 0.02938*SMQ$DBH^1.681)/(30.300 /1000*plotarea))
    QKQ_LAI=sum((0.07652*QKQ$DBH^1.511)/(34.425 /1000*plotarea))
    HKQ_LAI=sum(0.0081*HKQ$DBH^2.3418*37.899/plotarea)
    JSQ_LAI=sum(( 0.00933 *JSQ$DBH^1.939)/(23.117 /1000*plotarea))
    NJQ_LAI=sum(( 0.01117 *NJQ$DBH^1.873)/(30.084 /1000*plotarea))
    BNQ_LAI=sum(( 0.04882 *BNQ$DBH^1.597)/(35.058 /1000*plotarea))
    CY_LAI=sum(( 0.04997 *CY$DBH^1.530)/(61.778 /1000*plotarea))
    BH_LAI=sum(( 0.02491 *BH$DBH^1.913)/(47.379 /1000*plotarea))
    HH_LAI=sum(( 0.04255 *HH$DBH^1.635)/(38.858 /1000*plotarea))
    LS_LAI=sum(exp(-4.1730+2.0713*log(LS$DBH))*(7.096)/plotarea)
    YS_LAI=sum(exp(-3.5764+1.9801*log(YS$DBH))*(4.984)/plotarea)
    FH_LAI=sum(0.0081*FH$DBH^2.3418*19.744/plotarea)
    LYY_LAI=sum(0.0081*LYY$DBH^2.3418*30.016/plotarea)
    QT_LAI=sum(0.0081*QT$DBH^2.3418*26.63/plotarea)
    
    ###计算样地内的总叶面积
    LAI=sum(HS_LAI,LS_LAI,YS_LAI,ZD_LAI,KD_LAI,MGL_LAI,SQL_LAI,HTQ_LAI,HBL_LAI,SMQ_LAI,QKQ_LAI,HKQ_LAI,JSQ_LAI,NJQ_LAI,BNQ_LAI,CY_LAI,BH_LAI,HH_LAI,FH_LAI,LYY_LAI,QT_LAI)
    
    ###计算样地内针叶与阔叶树的叶面积
    Needles_LAI=sum(HS_LAI,LS_LAI,YS_LAI)
    Broadleaf_LAI=sum(ZD_LAI,KD_LAI,MGL_LAI,SQL_LAI,HTQ_LAI,HBL_LAI,SMQ_LAI,QKQ_LAI,HKQ_LAI,JSQ_LAI,NJQ_LAI,BNQ_LAI,CY_LAI,BH_LAI,HH_LAI,FH_LAI,LYY_LAI,QT_LAI)
    

  ###计算划定样地内针叶树叶面积与阔叶树叶面积所占的比例
  N_L_percent=Needles_LAI/LAI*100
  B_L_percent=Broadleaf_LAI/LAI*100
  
  ####整理计算结果
  Species_LAI=c(HS_LAI,LS_LAI,YS_LAI,ZD_LAI,KD_LAI,MGL_LAI,SQL_LAI,HTQ_LAI,HBL_LAI,SMQ_LAI,QKQ_LAI,HKQ_LAI,JSQ_LAI,NJQ_LAI,BNQ_LAI,CY_LAI,BH_LAI,HH_LAI,FH_LAI,LYY_LAI,QT_LAI)
  Species_LAI=matrix(Species_LAI,1,)
  colnames(Species_LAI)=c("HS_LAI","LS_LAI","YS_LAI","ZD_LAI","KD_LAI","MGL_LAI","SQL_LAI","HTQ_LAI","HBL_LAI","SMQ_LAI","QKQ_LAI", "HKQ_LAI","JSQ_LAI","NJQ_LAI","BNQ_LAI","CY_LAI","BH_LAI","HH_LAI","FH_LAI","LYY_LAI","QT_LAI")
  N_B_LAI=c(Needles_LAI,Broadleaf_LAI,N_L_percent,B_L_percent)
  N_B_LAI=matrix(N_B_LAI,1,)
  colnames(N_B_LAI)=c("Needles_LAI","Broadleaf_LAI","N_L_percent","B_L_percent")
  list(LAI=LAI,Species_LAI=Species_LAI,N_B_LAI=N_B_LAI)
}


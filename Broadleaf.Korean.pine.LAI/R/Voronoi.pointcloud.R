Voronoi.pointcloud=function(minx,maxx,miny,maxy,boundary,b,seq,strata,r,S,theta,phi)
{
  library(ape)
  LAI.Vaule = function(minx, maxx, miny, maxy, a, b,  r)
  {
    maxr = max(r)
    pointa = a
    library(raster)
    library(rgeos)
    library(sf)
    library(ggplot2)
    library(ggforce)
    library(grDevices)
    library(deldir)
    library(tcltk)
    library(deldir)
    raster.Voronoi.LAI = function(minx, maxx, miny, maxy, 
                                  b, r) {
      b1=b
      b = b[, -4]
      options(warn = -1)
      Voronoi.LAI.mult = function(minx, maxx, miny, maxy, 
                                  b, r) {
        dis = r/2 * sqrt(3)
        xpoint = seq(minx - 2 * r, maxx + 2 * r, 2 * 
                       dis)
        ypoint = seq(miny - 2 * r, maxy + 2 * r, 1.5 * 
                       r)
        ly = length(ypoint)
        lx = length(xpoint)
        xpoint = rep(xpoint, each = ly)
        ypoint = rep(ypoint, lx)
        uy = unique(ypoint)
        point = cbind(xpoint, ypoint, 1:length(xpoint))
        for (i in 1:(0.5 * length(uy))) {
          point[which(point[, 2] == uy[2 * i]), 1] = point[which(point[, 
                                                                       2] == uy[2 * i]), 1] + dis
        }
        colnames(point) = c("x", "y", "ID")
        pointoutside.single = function(a) {
          pointinside = subset(point, point[, 1] < a[, 
                                                     1] + r & point[, 1] > a[, 1] - r & point[, 
                                                                                              2] < a[, 2] + r & point[, 2] > a[, 2] - r)
          pointinside = as.data.frame(pointinside)
          if (nrow(pointinside) > 0) {
            c = pointinside[, 1:2]
            for (i in 1:nrow(pointinside)) {
              c[i, ] = (pointinside[i, 1:2] - a[1, 1:2])^2
              d = (c[, 1] + c[, 2])^(1/2)
            }
            d = cbind(pointinside, d)
            d = subset(d, d > 0)
            colnames(d) = c("x", "y", "ID", 
                            "Distance")
            d
            pointinside = subset(d, d[, 4] < r)
            pointinside = pointinside[which.min(pointinside[, 
                                                            4]), ]
            pointinside = pointinside[, 1:3]
            pointoutside_ID = setdiff(point[, 3], pointinside[, 
                                                              3])
            pointoutside = point[pointoutside_ID, 1:3]
          }
          if (nrow(pointinside) == 0) {
            pointoutside = point
          }
          pointoutside
        }
        b_for = pointoutside.single(b[1, 1:2])[, 3]
        for (i in 2:nrow(b)) {
          bn = pointoutside.single(b[i, 1:2])[, 3]
          b_for = intersect(b_for, bn)
        }
        point = point[b_for, 1:2]
        point = cbind(point, 0)
        point = cbind(point, 0)
        point = as.data.frame(point)
        point$Species = c("LX")
        
        colnames(point) = c("x", "y", "DBH","H", 
                            "Species")
        bnew = rbind(point, b1)
        bnew = bnew[deldir(bnew[, 1], bnew[, 2])$ind.orig, 
        ]
        deldir_area = deldir(bnew[, 1], bnew[, 2])$summary$dir.area
        bnew = cbind(bnew, deldir_area)
        HS = subset(bnew, bnew$Species == "HS")
        ZD = subset(bnew, bnew$Species == "ZD")
        KD = subset(bnew, bnew$Species == "KD")
        MGL = subset(bnew, bnew$Species == "MGL")
        SQL = subset(bnew, bnew$Species == "SQL")
        HTQ = subset(bnew, bnew$Species == "HTQ")
        HBL = subset(bnew, bnew$Species == "HBL")
        SMQ = subset(bnew, bnew$Species == "SMQ")
        QKQ = subset(bnew, bnew$Species == "QKQ")
        JSQ = subset(bnew, bnew$Species == "JSQ")
        NJQ = subset(bnew, bnew$Species == "NJQ")
        BNQ = subset(bnew, bnew$Species == "BNQ")
        HKQ = subset(bnew, bnew$Species == "HKQ")
        CY = subset(bnew, bnew$Species == "CY")
        BH = subset(bnew, bnew$Species == "BH")
        HH = subset(bnew, bnew$Species == "HH")
        LS = subset(bnew, bnew$Species == "LS")
        YS = subset(bnew, bnew$Species == "YS")
        FH = subset(bnew, bnew$Species == "FH")
        LYY = subset(bnew, bnew$Species == "LYY")
        QT = subset(bnew, bnew$Species == "QT")
        LX = subset(bnew, bnew$Species == "LX")
        HS$LAI = (0.04321 * HS$DBH^1.831)/(108.786/1000 * 
                                             HS$deldir_area)
        ZD$LAI = (0.0254 * ZD$DBH^1.632)/(33.583/1000 * 
                                            ZD$deldir_area)
        KD$LAI = (0.01034 * KD$DBH^1.797)/(27.086/1000 * 
                                             KD$deldir_area)
        MGL$LAI = (0.05095 * MGL$DBH^1.639)/(45.413/1000 * 
                                               MGL$deldir_area)
        SQL$LAI = (0.10465 * SQL$DBH^1.417)/(57.241/1000 * 
                                               SQL$deldir_area)
        HTQ$LAI = (0.01144 * HTQ$DBH^1.834)/(54.843/1000 * 
                                               HTQ$deldir_area)
        HBL$LAI = (0.01546 * HBL$DBH^1.758)/(38.184/1000 * 
                                               HBL$deldir_area)
        SMQ$LAI = (0.02938 * SMQ$DBH^1.681)/(30.3/1000 * 
                                               SMQ$deldir_area)
        QKQ$LAI = (0.07652 * QKQ$DBH^1.511)/(34.425/1000 * 
                                               QKQ$deldir_area)
        HKQ$LAI = 0.0081 * HKQ$DBH^2.3418 * 37.899/HKQ$deldir_area
        JSQ$LAI = (0.00933 * JSQ$DBH^1.939)/(23.117/1000 * 
                                               JSQ$deldir_area)
        NJQ$LAI = (0.01117 * NJQ$DBH^1.873)/(30.084/1000 * 
                                               NJQ$deldir_area)
        BNQ$LAI = (0.04882 * BNQ$DBH^1.597)/(35.058/1000 * 
                                               BNQ$deldir_area)
        CY$LAI = (0.04997 * CY$DBH^1.53)/(61.778/1000 * 
                                            CY$deldir_area)
        BH$LAI = (0.02491 * BH$DBH^1.913)/(47.379/1000 * 
                                             BH$deldir_area)
        HH$LAI = (0.04255 * HH$DBH^1.635)/(38.858/1000 * 
                                             HH$deldir_area)
        LS$LAI = exp(-4.173 + 2.0713 * log(LS$DBH)) * 
          (7.096)/LS$deldir_area
        YS$LAI = exp(-3.5764 + 1.9801 * log(YS$DBH)) * 
          (4.984)/YS$deldir_area
        FH$LAI = 0.0081 * FH$DBH^2.3418 * 19.744/FH$deldir_area
        LYY$LAI = 0.0081 * LYY$DBH^2.3418 * 30.016/LYY$deldir_area
        QT$LAI = 0.0081 * QT$DBH^2.3418 * 26.63/QT$deldir_area
        LX$LAI = 0
        bnew = rbind(HS, LS, YS, ZD, KD, MGL, SQL, HTQ, 
                     HBL, SMQ, QKQ, HKQ, JSQ, NJQ, BNQ, CY, BH, 
                     HH, FH, LYY, QT, LX)
        bnew
      }
      bnew = Voronoi.LAI.mult(minx, maxx, miny, maxy, b, 
                              r)
      bnew$group = 1:nrow(bnew)
      a = deldir(bnew[, 1], bnew[, 2])
      a = a$dirsgs
      a1 = cbind(a[, 1:2], a[, 5])
      a2 = cbind(a[, 3:4], a[, 5])
      colnames(a1) = c("x", "y", "group")
      colnames(a2) = c("x", "y", "group")
      a12 = rbind(a1, a2)
      b1 = cbind(a[, 1:2], a[, 6])
      b2 = cbind(a[, 3:4], a[, 6])
      colnames(b1) = c("x", "y", "group")
      colnames(b2) = c("x", "y", "group")
      b12 = rbind(b1, b2)
      ab12 = rbind(a12, b12)
      ab12 = ab12[order(ab12[, 3]), ]
      ab12_1 = subset(ab12, ab12[, 3] == 1)
      order = chull(ab12_1[, 1], ab12_1[, 2])
      ab12_1 = ab12_1[order, ]
      for (i in 2:max(ab12[, 3])) {
        ab12_n = subset(ab12, ab12[, 3] == i)
        order = chull(ab12_n[, 1], ab12_n[, 2])
        ab12_n = ab12_n[order, ]
        ab12_1 = rbind(ab12_1, ab12_n)
      }
      for (i in 1:nrow(ab12_1)) {
        ab12_1[i, c(4,5)] = bnew[which(bnew$group == ab12_1[i, 
                                                            3]),c(4,7)]
      }
      colnames(ab12_1) = c("x", "y", "group", "H",
                           "LAI")
      ab12_1 = subset(ab12_1, ab12_1$LAI > 0.002)
      n = max(ab12_1$group)
      limitxy = data.frame(x = c(minx - 2 * maxr, maxx + 
                                   2 * maxr), y = c(miny - 2 * maxr, maxy + 2 * 
                                                      maxr), group = c(n + 1, n + 2),H=0, LAI = 0.01)
      ab12_1 = rbind(ab12_1, limitxy)
      point = b[, 1:2]
      listpolygons = list()
      n = n + 2
      for (i in 1:n) {
        polygonsi = subset(ab12_1, ab12_1$group == i)[, 
                                                      1:2]
        polygonsi = Polygon(polygonsi)
        listpolygons[i] = Polygons(list(polygonsi), paste("DI", 
                                                          i))
      }
      p1 = SpatialPolygons(listpolygons, 1:n)
      poldata = createSPComment(p1)
      poldata$H = ab12_1[!duplicated(ab12_1$group), ]$H
      poldata$LAI = ab12_1[!duplicated(ab12_1$group), ]$LAI
      poldata = as(poldata, "SpatialPolygonsDataFrame")
      ra = raster(poldata, res = 0.05)
      shape_r_LAI = rasterize(poldata, ra, "LAI")
      shape_r_LAI = reclassify(shape_r_LAI , cbind(NA, 0), right = F)
      shape_r_H = rasterize(poldata, ra, "H")
      shape_r_H = reclassify(shape_r_H, cbind(NA, 0), right = F)
      result=list(H=shape_r_H,LAI=shape_r_LAI)
      result
    }
    
    shape_r1 = raster.Voronoi.LAI(minx, maxx, miny, maxy, 
                                  b, r)
    
    a1 = a
    coordinates(a1) <- ~x + y
    H = extract(shape_r1$H, a1)
    LAI = extract(shape_r1$LAI, a1)
    result = a
    result$H=H
    result$LAI = LAI
    result
  }
  
  
  
  strata = c(0, strata, Inf)
  Forest_strata = list()
  for (i in 1:(length(strata) - 1)) {
    Forest_strata[[i]] = subset(b, b$H >= strata[i] & 
                                  b$H < strata[i + 1])
  }
  maxx1 = maxx - boundary
  minx1 = minx + boundary
  maxy1 = maxy - boundary
  miny1 = miny + boundary
  seq=seq
  ax=seq(minx1,maxx1,length.out=seq)
  ay=seq(miny1,maxy1,length.out=seq)
  a=expand.grid(ax,ay)
  colnames(a)=c("x","y")
  
  pointcould1=LAI.Vaule(minx, maxx, miny, maxy,a ,Forest_strata[[1]],  r[1])
  
  
  
  for (i in 2:length(Forest_strata)) {
    pointcould=LAI.Vaule(minx, maxx, miny, maxy, a ,Forest_strata[[i]],  r[i])
    pointcould1 = rbind(pointcould1,pointcould)
  }
  
  pointcould1=subset(pointcould1,pointcould1$LAI>0)
  pointcould1$LAI=pointcould1$LAI/max(pointcould1$LAI)*30
  pointcould1$LAI=round(pointcould1$LAI)
  
  point_strata = list()
  for (i in 1:(length(strata) - 1)) {
    point_strata[[i]] = subset(pointcould1, pointcould1$H >= strata[i] & 
                                 pointcould1$H < strata[i + 1])
  }
  ri=r[1]
  p_s1=point_strata[[1]]
  
  P_could=function(p_s1,ri)
  {
    p_s1$H=p_s1$H-S*ri
    p_s1_1=matrix(NA,p_s1[1,4],2)
    p_s1_1[,1]=p_s1[1,1]
    p_s1_1[,2]=p_s1[1,2]
    p_s1_1=as.data.frame(p_s1_1)
    p_s1_1$H=runif(p_s1[1,4],0,S*ri)+p_s1[1,3]
    for(i in 2:nrow(p_s1))
    {
      p_s1_i=matrix(NA,p_s1[i,4],2)
      p_s1_i[,1]=p_s1[i,1]
      p_s1_i[,2]=p_s1[i,2]
      p_s1_i=as.data.frame(p_s1_i)
      p_s1_i$H=runif(p_s1[i,4],0,S*ri)+p_s1[i,3] 
      p_s1_1=rbind(p_s1_1,p_s1_i)
    }
    result=p_s1_1
    result
  }
  
  P_C_1=P_could(p_s1,ri)
  P_C=P_C_1
  for (i in 2:length(r))
  {
    ri=r[i]
    p_si=point_strata[[i]] 
    P_C_i=P_could(p_si,ri)
    P_C=rbind(P_C,P_C_i)
  }
  colnames(P_C)=c("x","y","z")
  
  ri=r[1]
  Tree_T=Forest_strata[[1]]
  Tree_Trunk=function(Tree_T,ri)
  {
    Tree_T$H=Tree_T$H-S*ri
    Tree_T
  }
  Tree_T_1=Tree_Trunk(Tree_T,ri)
  Tree_trunk=Tree_T_1
  for (i in 2:length(r))
  {
    ri=r[i]
    Tree_T=Forest_strata[[i]] 
    Tree_T_i=Tree_Trunk(Tree_T,ri)
    Tree_trunk=rbind(Tree_trunk,Tree_T_i)
  }
  Tree_trunk=subset(Tree_trunk,Tree_trunk$x>minx1 & Tree_trunk$x<maxx1 & Tree_trunk$y>miny1 & Tree_trunk$y<maxy1)
  Tree_trunk=Tree_trunk[,c(1,2,4)]
  colnames(Tree_trunk)=c("x","y","z")
  
  lengthx=maxx1-minx1
  lengthy=maxy1-miny1
  
  library(plot3D)
  scatter3D(x = P_C$x, y = P_C$y, z = P_C$z,zlim=c(0,max(lengthx,lengthy,P_C$z)),
            pch = 21, cex = 0.25,col="lightgreen",bg="green3",ticktype = "detailed",
            
            theta = theta, phi =phi, 
            colkey = FALSE)
  scatter3D(x =Tree_trunk$x, y = Tree_trunk$y, z = Tree_trunk$z,
            pch = 21, cex = 0.25,col="black",bg="black", type = "h",
            theta = theta, phi =phi,
            colkey = FALSE,add=T)
}
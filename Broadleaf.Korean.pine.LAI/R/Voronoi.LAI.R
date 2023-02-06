Voronoi.LAI=function (minx, maxx, miny, maxy, boundary, b, r) 
{
  library(deldir)
  dis = r/2 * sqrt(3)
  xpoint = seq(minx-2*r, maxx+2*r, 2 * dis)
  ypoint = seq(miny-2*r, maxy+2*r, 1.5 * r)
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
    pointinside = subset(point, point[, 1] < a[, 1] + r & 
                           point[, 1] > a[, 1] - r & point[, 2] < a[, 2] + r & 
                           point[, 2] > a[, 2] - r)
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
  point = as.data.frame(point)
  point$Species = c("LX")
  colnames(point) = c("x", "y", "DBH", "Species")
  bnew = rbind(point, b)
  bnew = bnew[deldir(bnew[, 1], bnew[, 2])$ind.orig, ]
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
  HS$LAI = (0.04321 * HS$DBH^1.831)/(108.786/1000 * HS$deldir_area)
  ZD$LAI = (0.0254 * ZD$DBH^1.632)/(33.583/1000 * ZD$deldir_area)
  KD$LAI = (0.01034 * KD$DBH^1.797)/(27.086/1000 * KD$deldir_area)
  MGL$LAI = (0.05095 * MGL$DBH^1.639)/(45.413/1000 * MGL$deldir_area)
  SQL$LAI = (0.10465 * SQL$DBH^1.417)/(57.241/1000 * SQL$deldir_area)
  HTQ$LAI = (0.01144 * HTQ$DBH^1.834)/(54.843/1000 * HTQ$deldir_area)
  HBL$LAI = (0.01546 * HBL$DBH^1.758)/(38.184/1000 * HBL$deldir_area)
  SMQ$LAI = (0.02938 * SMQ$DBH^1.681)/(30.3/1000 * SMQ$deldir_area)
  QKQ$LAI = (0.07652 * QKQ$DBH^1.511)/(34.425/1000 * QKQ$deldir_area)
  HKQ$LAI = 0.0081 * HKQ$DBH^2.3418 * 37.899/HKQ$deldir_area
  JSQ$LAI = (0.00933 * JSQ$DBH^1.939)/(23.117/1000 * JSQ$deldir_area)
  NJQ$LAI = (0.01117 * NJQ$DBH^1.873)/(30.084/1000 * NJQ$deldir_area)
  BNQ$LAI = (0.04882 * BNQ$DBH^1.597)/(35.058/1000 * BNQ$deldir_area)
  CY$LAI = (0.04997 * CY$DBH^1.53)/(61.778/1000 * CY$deldir_area)
  BH$LAI = (0.02491 * BH$DBH^1.913)/(47.379/1000 * BH$deldir_area)
  HH$LAI = (0.04255 * HH$DBH^1.635)/(38.858/1000 * HH$deldir_area)
  LS$LAI = exp(-4.173 + 2.0713 * log(LS$DBH)) * (7.096)/LS$deldir_area
  YS$LAI = exp(-3.5764 + 1.9801 * log(YS$DBH)) * (4.984)/YS$deldir_area
  FH$LAI = 0.0081 * FH$DBH^2.3418 * 19.744/FH$deldir_area
  LYY$LAI = 0.0081 * LYY$DBH^2.3418 * 30.016/LYY$deldir_area
  QT$LAI = 0.0081 * QT$DBH^2.3418 * 26.63/QT$deldir_area
  LX$LAI = 0
  bnew = rbind(HS, LS, YS, ZD, KD, MGL, SQL, HTQ, HBL, SMQ, 
               QKQ, HKQ, JSQ, NJQ, BNQ, CY, BH, HH, FH, LYY, QT, LX)
  bnew=subset(bnew,bnew$LAI>0)
  bnew
}

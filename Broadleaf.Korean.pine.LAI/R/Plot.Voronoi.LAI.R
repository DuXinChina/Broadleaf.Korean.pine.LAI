Plot.Voronoi.LAI=function (minx, maxx, miny, maxy, boundary, b, r) 
{
  library(ggplot2)
  library(grDevices)
  library(deldir)
  Voronoi.LAI.mult = function(minx, maxx, miny, maxy, boundary, 
                              b, r) {
    dis = r/2 * sqrt(3)
    xpoint = seq(minx, maxx, 2 * dis)
    ypoint = seq(miny, maxy, 1.5 * r)
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
      pointinside = subset(point, point[, 1] < a[, 1] + 
                             r & point[, 1] > a[, 1] - r & point[, 2] < a[, 
                                                                          2] + r & point[, 2] > a[, 2] - r)
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
    colnames(point) = c("x", "y", "DBH", 
                        "Species")
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
    bnew = rbind(HS, LS, YS, ZD, KD, MGL, SQL, HTQ, HBL, 
                 SMQ, QKQ, HKQ, JSQ, NJQ, BNQ, CY, BH, HH, FH, LYY, 
                 QT, LX)
    bnew = subset(bnew, bnew[, 1] > (minx + boundary) & bnew[, 
                                                             1] < (maxx - boundary) & bnew[, 2] > (miny + boundary) & 
                    bnew[, 2] < (maxy - boundary))
    bnew
  }
  bnew = Voronoi.LAI.mult(minx, maxx, miny, maxy, 0, b, r)
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
    ab12_1[i, 4] = bnew[which(bnew$group == ab12_1[i, 3]), 
                        6]
  }
  colnames(ab12_1) = c("x", "y", "group", 
                       "LAI")
  ab12_1[which(ab12_1[, 4] == 0), 4] = 0.001
  expandx = (maxx - minx) * 0.095
  expandy = (maxy - miny) * 0.095
  aex = 1.053 * expandx
  aey = 1.053 * expandy
  x1 = c(minx - aex, minx - aex, minx + boundary, minx + boundary)
  y1 = c(miny - aey, maxy + aey, maxy + aey, miny - aey)
  xy1 = cbind(x1, y1, 1)
  xy1 = as.data.frame(xy1)
  colnames(xy1) = c("x", "y", "group")
  x2 = c(minx - aex, minx - aex, maxx + aex, maxx + aex)
  y2 = c(miny - aey, miny + boundary, miny + boundary, miny - 
           aey)
  xy2 = cbind(x2, y2, 2)
  xy2 = as.data.frame(xy2)
  colnames(xy2) = c("x", "y", "group")
  x3 = c(maxx - boundary, maxx - boundary, maxx + aex, maxx + 
           aex)
  y3 = c(miny - aey, maxy + aex, maxy + aex, miny - aey)
  xy3 = cbind(x3, y3, 3)
  xy3 = as.data.frame(xy3)
  colnames(xy3) = c("x", "y", "group")
  x4 = c(minx - aex, minx - aex, maxx + aex, maxx + aex)
  y4 = c(maxy - boundary, maxy + aey, maxy + aey, maxy - boundary)
  xy4 = cbind(x4, y4, 4)
  xy4 = as.data.frame(xy4)
  colnames(xy4) = c("x", "y", "group")
  xy = rbind(xy1, xy2, xy3, xy4)
  ggplot() + geom_polygon(data = ab12_1, aes(x = x, y = y, 
                                             group = group, fill = LAI), colour = "black") + 
    scale_fill_gradientn(colours = c("white", "yellow2", 
                                            "green", "green3", "green4")) + 
                                              geom_polygon(data = xy, aes(x = x, y = y, group = group), 
                                                           colour = "white", fill = "white") + scale_x_continuous(expand = c(0, 
                                                                                                                             -expandx)) + scale_y_continuous(expand = c(0, -expandy)) + 
    geom_vline(xintercept = c(minx + boundary, maxx - boundary), 
               linetype = 2, size = 1) + geom_hline(yintercept = c(miny + 
                                                                     boundary, maxy - boundary), linetype = 2, size = 1) + 
    theme_bw()
}
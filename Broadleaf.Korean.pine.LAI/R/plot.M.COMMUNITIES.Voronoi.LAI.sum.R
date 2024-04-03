
plot.M.COMMUNITIES.Voronoi.LAI.sum=function(minx,maxx,miny,maxy,a,b,strata,r,Min_com_edge)
{
  options(warn = -1)
  library(grDevices)
  M.COMMUNITIES.Voronoi.LAI.sum = function(minx, maxx, miny,maxy, boundary,a,b,strata,r,Min_com_edge)
    {
    library(tcltk)
    library(deldir)
    Voronoi.LAI.mult = function(minx, maxx, miny, maxy, boundary, 
                                b, r) {
      dis = r/2 * sqrt(3)
      xpoint = seq(minx - 2 * r, maxx + 2 * r, 2 * dis)
      ypoint = seq(miny - 2 * r, maxy + 2 * r, 1.5 * r)
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
          colnames(d) = c("x", "y", "ID", "Distance")
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
      ZD$LAI = (0.0254 * ZD$DBH^1.632)/(33.583/1000 * ZD$deldir_area)
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
      CY$LAI = (0.04997 * CY$DBH^1.53)/(61.778/1000 * CY$deldir_area)
      BH$LAI = (0.02491 * BH$DBH^1.913)/(47.379/1000 * 
                                           BH$deldir_area)
      HH$LAI = (0.04255 * HH$DBH^1.635)/(38.858/1000 * 
                                           HH$deldir_area)
      LS$LAI = exp(-4.173 + 2.0713 * log(LS$DBH)) * (7.096)/LS$deldir_area
      YS$LAI = exp(-3.5764 + 1.9801 * log(YS$DBH)) * (4.984)/YS$deldir_area
      FH$LAI = 0.0081 * FH$DBH^2.3418 * 19.744/FH$deldir_area
      LYY$LAI = 0.0081 * LYY$DBH^2.3418 * 30.016/LYY$deldir_area
      QT$LAI = 0.0081 * QT$DBH^2.3418 * 26.63/QT$deldir_area
      LX$LAI = 0
      bnew = rbind(HS, LS, YS, ZD, KD, MGL, SQL, HTQ, HBL, 
                   SMQ, QKQ, HKQ, JSQ, NJQ, BNQ, CY, BH, HH, FH, 
                   LYY, QT, LX)
      bnew
    }
    strata = c(0, strata, Inf)
    Forest_strata = list()
    for (i in 1:(length(strata) - 1)) {
      Forest_strata[[i]] = subset(b, b$H >= strata[i] & 
                                    b$H < strata[i + 1])
      Forest_strata[[i]] = Voronoi.LAI.mult(minx, maxx, 
                                            miny, maxy, 0, Forest_strata[[i]][, c(1:3, 5)], 
                                            r[i])
    }
    ConvexHull_point = function(a, Min_com_edge) {
      if (missing(Min_com_edge)) {
        warning("missing Min_com_edge, automatic set 0.1 ")
        Min_com_edge = 0.1
      }
      if (nrow(a) == 1) {
        warning("just one point")
      }
      if (nrow(a) == 1 & Min_com_edge == 0) {
        stop("Min_com_edge cant be 0 when one point")
      }
      up_point = a
      left_point = a
      right_point = a
      down_point = a
      up_point[, 2] = up_point[, 2] + Min_com_edge
      left_point[, 1] = left_point[, 1] - Min_com_edge
      right_point[, 1] = right_point[, 1] + Min_com_edge
      down_point[, 2] = down_point[, 2] - Min_com_edge
      total_point = rbind(up_point, left_point, right_point, 
                          down_point)
      total_point = unique(total_point)
      total_point
      hull_indices = chull(total_point)
      hull_vertices <<- total_point[hull_indices, ]
      point_x = seq(min(hull_vertices[, 1]), max(hull_vertices[, 
                                                               1]), length.out = 20)
      point_y = seq(min(hull_vertices[, 2]), max(hull_vertices[, 
                                                               2]), length.out = 20)
      point_xy = expand.grid(point_x, point_y)
      colnames(point_xy) = c("x", "y")
      new_point_inside = function(new_point) {
        point = rbind(new_point, hull_vertices)
        pointn = chull(point)
        if (1 %in% pointn) TRUE
        else F
      }
      inside_row = apply(point_xy, 1, new_point_inside)
      point_xy = point_xy[!inside_row, ]
      point_xy
      list(point_xy = point_xy, hull_vertices = hull_vertices)
    }
    point = ConvexHull_point(a, Min_com_edge)$point_xy
    hull_vertices = ConvexHull_point(a, Min_com_edge)$hull_vertices
    stra_single = Forest_strata[[1]]
    Lbnew1 = subset(stra_single, stra_single[, 1] >= (point[1, 
                                                            1] - 1.5 * r[1]) & stra_single[, 1] <= (point[1, 
                                                                                                          1] + 1.5 * r[1]) & stra_single[, 2] >= (point[1, 
                                                                                                                                                        2] - 1.5 * r[1]) & stra_single[, 2] <= (point[1, 
                                                                                                                                                                                                      2] + 1.5 * r[1]))
    Lbnew1$d = (point[1, 1] - Lbnew1[, 1])^2 + (point[1, 
                                                      2] - Lbnew1[, 2])^2
    Lbnew1 = Lbnew1[which.min(Lbnew1$d), ]
    for (i in 2:nrow(point)) {
      Lbnew = subset(stra_single, stra_single[, 1] >= (point[i, 
                                                             1] - 1.5 * r[1]) & stra_single[, 1] <= (point[i, 
                                                                                                           1] + 1.5 * r[1]) & stra_single[, 2] >= (point[i, 
                                                                                                                                                         2] - 1.5 * r[1]) & stra_single[, 2] <= (point[i, 
                                                                                                                                                                                                       2] + 1.5 * r[1]))
      Lbnew$d = (point[i, 1] - Lbnew[, 1])^2 + (point[i, 
                                                      2] - Lbnew[, 2])^2
      Lbnew = Lbnew[which.min(Lbnew$d), ]
      Lbnew1 = rbind(Lbnew1, Lbnew)
    }
    Lbnew1 = Lbnew1[, c(4, 6)]
    for (j in 2:length(Forest_strata)) {
      stra_single = Forest_strata[[j]]
      Lbnewn = subset(stra_single, stra_single[, 1] >= 
                        (point[1, 1] - 1.5 * r[j]) & stra_single[, 1] <= 
                        (point[1, 1] + 1.5 * r[j]) & stra_single[, 2] >= 
                        (point[1, 2] - 1.5 * r[j]) & stra_single[, 2] <= 
                        (point[1, 2] + 1.5 * r[j]))
      Lbnewn$d = (point[1, 1] - Lbnewn[, 1])^2 + (point[1, 
                                                        2] - Lbnewn[, 2])^2
      Lbnewn = Lbnewn[which.min(Lbnewn$d), ]
      for (i in 2:nrow(point)) {
        Lbnew = subset(stra_single, stra_single[, 1] >= 
                         (point[i, 1] - 1.5 * r[j]) & stra_single[, 
                                                                  1] <= (point[i, 1] + 1.5 * r[j]) & stra_single[, 
                                                                                                                 2] >= (point[i, 2] - 1.5 * r[j]) & stra_single[, 
                                                                                                                                                                2] <= (point[i, 2] + 1.5 * r[j]))
        Lbnew$d = (point[i, 1] - Lbnew[, 1])^2 + (point[i, 
                                                        2] - Lbnew[, 2])^2
        Lbnew = Lbnew[which.min(Lbnew$d), ]
        Lbnewn = rbind(Lbnewn, Lbnew)
      }
      Lbnewn = Lbnewn[, c(4, 6)]
      Lbnew1 = cbind(Lbnew1, Lbnewn)
    }
    aresult = Lbnew1[1, ]
    Lbnew1 = Lbnew1[-1, ]
    aresult = aresult[, 1:length(Forest_strata) * 2]
    aresult = sum(aresult)
    multLAI = Lbnew1[, 1:length(Forest_strata) * 2]
    rbLAI = Lbnew1[, 1:2]
    for (i in 2:length(Forest_strata)) {
      rbLAI = rbind(rbLAI, Lbnew1[, c(2 * i - 1, 2 * i)])
    }
    Canopy = rbLAI[which(!rbLAI$LAI == 0), ]
    needle = subset(Canopy, Canopy$Species == "HS" | Canopy$Species == 
                      "YS" | Canopy$Species == "LS")
    N_L_Percent = sum(needle$LAI)/sum(Canopy$LAI)
    broad = Canopy[which(!Canopy$Species == "HS" & !Canopy$Species == 
                           "YS" & !Canopy$Species == "LS"), ]
    B_L_Percent = sum(broad$LAI)/sum(Canopy$LAI)
    LAIsum = sum(multLAI)
    colSumsLAI = colSums(multLAI)
    rowSumsLAI = rowSums(multLAI)
    sdLAI = sd(rowSumsLAI)
    meanLAI = mean(rowSumsLAI)
    minLAI = min(rowSumsLAI)
    maxLAI = max(rowSumsLAI)
    Gap_percent = length(which(rowSumsLAI == 0))/length(rowSumsLAI)
    Canopy_percent = length(which(!rowSumsLAI == 0))/length(rowSumsLAI)
    weightLAI = matrix(NA, length(Forest_strata), 1)
    for (i in 1:length(Forest_strata)) {
      weightLAI[i, ] = colSumsLAI[i] * i
    }
    Strata_cont = sum(weightLAI)/(LAIsum * mean(1:length(Forest_strata)))
    result = data.frame(mean(a$x), mean(a$y), minLAI, 
                        maxLAI, meanLAI, sdLAI, Gap_percent, Canopy_percent, 
                        N_L_Percent, B_L_Percent, Strata_cont)
    colnames(result) = c("x", "y", "Local_min_LAI", 
                         "Local_max_LAI", "Local_mean_LAI", "Local_sd_LAI", 
                         "Gap_percent", "Canopy_percent", "N_L_Percent", "B_L_Percent", 
                         "Strata_cont")
    result
  }
  print(M.COMMUNITIES.Voronoi.LAI.sum(minx, maxx, miny, maxy, boundary, a, b, strata, r, Min_com_edge))
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
  raster.Voronoi.LAI = function(minx, maxx, miny, maxy, b, 
                                r) {
    b = b[, -4]
    options(warn = -1)
    Voronoi.LAI.mult = function(minx, maxx, miny, maxy, b, 
                                r) {
      dis = r/2 * sqrt(3)
      xpoint = seq(minx - 2 * r, maxx + 2 * r, 2 * dis)
      ypoint = seq(miny - 2 * r, maxy + 2 * r, 1.5 * r)
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
          colnames(d) = c("x", "y", "ID", "Distance")
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
      ZD$LAI = (0.0254 * ZD$DBH^1.632)/(33.583/1000 * ZD$deldir_area)
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
      CY$LAI = (0.04997 * CY$DBH^1.53)/(61.778/1000 * CY$deldir_area)
      BH$LAI = (0.02491 * BH$DBH^1.913)/(47.379/1000 * 
                                           BH$deldir_area)
      HH$LAI = (0.04255 * HH$DBH^1.635)/(38.858/1000 * 
                                           HH$deldir_area)
      LS$LAI = exp(-4.173 + 2.0713 * log(LS$DBH)) * (7.096)/LS$deldir_area
      YS$LAI = exp(-3.5764 + 1.9801 * log(YS$DBH)) * (4.984)/YS$deldir_area
      FH$LAI = 0.0081 * FH$DBH^2.3418 * 19.744/FH$deldir_area
      LYY$LAI = 0.0081 * LYY$DBH^2.3418 * 30.016/LYY$deldir_area
      QT$LAI = 0.0081 * QT$DBH^2.3418 * 26.63/QT$deldir_area
      LX$LAI = 0
      bnew = rbind(HS, LS, YS, ZD, KD, MGL, SQL, HTQ, HBL, 
                   SMQ, QKQ, HKQ, JSQ, NJQ, BNQ, CY, BH, HH, FH, 
                   LYY, QT, LX)
      bnew
    }
    bnew = Voronoi.LAI.mult(minx, maxx, miny, maxy, b, r)
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
      ab12_1[i, 4] = bnew[which(bnew$group == ab12_1[i, 
                                                     3]), 6]
    }
    colnames(ab12_1) = c("x", "y", "group", "LAI")
    ab12_1 = subset(ab12_1, ab12_1$LAI > 0.002)
    n = max(ab12_1$group)
    limitxy = data.frame(x = c(minx - 2 * maxr, maxx + 2 * 
                                 maxr), y = c(miny - 2 * maxr, maxy + 2 * maxr), group = c(n + 
                                                                                             1, n + 2), LAI = 0.01)
    ab12_1 = rbind(ab12_1, limitxy)
    point = b[, 1:2]
    listpolygons = list()
    n = n + 2
    for (i in 1:n) {
      polygonsi = subset(ab12_1, ab12_1$group == i)[, 1:2]
      polygonsi = Polygon(polygonsi)
      listpolygons[i] = Polygons(list(polygonsi), paste("DI", 
                                                        i))
    }
    p1 = SpatialPolygons(listpolygons, 1:n)
    poldata = createSPComment(p1)
    poldata$LAI = ab12_1[!duplicated(ab12_1$group), ]$LAI
    poldata = as(poldata, "SpatialPolygonsDataFrame")
    ra = raster(poldata, res = 0.05)
    shape_r = rasterize(poldata, ra, "LAI")
    shape_r = reclassify(shape_r, cbind(NA, 0), right = F)
    shape_r
  }
  strata = c(0, strata, Inf)
  Forest_strata = list()
  for (i in 1:(length(strata) - 1)) {
    Forest_strata[[i]] = subset(b, b$H >= strata[i] & b$H < 
                                  strata[i + 1])
  }
  shape_r1 = raster.Voronoi.LAI(minx, maxx, miny, maxy, Forest_strata[[1]], 
                                r[1])
  pb = tkProgressBar("Finsih", "Percentcomplete%", 0, 100)
  for (i in 2:length(Forest_strata)) {
    shape_r = raster.Voronoi.LAI(minx, maxx, miny, maxy, 
                                 Forest_strata[[i]], r[i])
    shape_r1 = overlay(shape_r1, shape_r, fun = function(x, 
                                                         y) {
      return(x + y)
    })
    info = sprintf("Percentcomplete%d%%", round(i * 100/length(Forest_strata)))
    setTkProgressBar(pb, i * 100/length(Forest_strata), "Finish", 
                     info)
  }
  close(pb)
  shape_r1 = as(shape_r1, "SpatialPixelsDataFrame")
  point = as.data.frame(shape_r1)
  point = data.frame(x = point$x, y = point$y, LAI = point$layer)
  point = subset(point, point$LAI > 0)
  ggplot() + geom_raster(data = point, aes(x = x, y = y, fill = LAI)) + 
    geom_point(data = pointa, aes(x = x, y = y), size = 0.8) + 
    theme_classic() + scale_fill_gradientn(colours = c("white", 
                                                              "yellow2", "green", "green3", "green4")) + scale_x_continuous(expand = c(0, 
                                                              0)) + scale_y_continuous(expand = c(0, 0)) + geom_polygon(data = hull_vertices, 
                                                                                                                        aes(x, y), col = "red4", fill = "white", alpha = 0) 
}

Voronoi.LAI.sum.ISAA=function (minx, maxx, miny, maxy, boundary, b, seq, strata, r, indis, lag) 
{
  library(ape)
  LAI.Vaule = function(minx, maxx, miny, maxy, a, b, strata, 
                       r) {
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
        point = as.data.frame(point)
        point$Species = c("LX")
        colnames(point) = c("x", "y", "DBH", 
                            "Species")
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
        ab12_1[i, 4] = bnew[which(bnew$group == ab12_1[i, 
                                                       3]), 6]
      }
      colnames(ab12_1) = c("x", "y", "group", 
                           "LAI")
      ab12_1 = subset(ab12_1, ab12_1$LAI > 0.002)
      n = max(ab12_1$group)
      limitxy = data.frame(x = c(minx - 2 * maxr, maxx + 
                                   2 * maxr), y = c(miny - 2 * maxr, maxy + 2 * 
                                                      maxr), group = c(n + 1, n + 2), LAI = 0.01)
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
      Forest_strata[[i]] = subset(b, b$H >= strata[i] & 
                                    b$H < strata[i + 1])
    }
    shape_r1 = raster.Voronoi.LAI(minx, maxx, miny, maxy, 
                                  Forest_strata[[1]], r[1])
    for (i in 2:length(Forest_strata)) {
      shape_r = raster.Voronoi.LAI(minx, maxx, miny, maxy, 
                                   Forest_strata[[i]], r[i])
      shape_r1 = overlay(shape_r1, shape_r, fun = function(x, 
                                                           y) {
        return(x + y)
      })
    }
    a1 = a
    coordinates(a1) <- ~x + y
    LAI = extract(shape_r1, a1)
    result = a
    result$LAI = LAI
    result
  }
  maxx1 = maxx - boundary
  minx1 = minx + boundary
  maxy1 = maxy - boundary
  miny1 = miny + boundary
  
  nx = (maxx1 - minx1 - indis)%/%lag
  ny = (maxy1 - miny1 - indis)%/%lag
  n = min(nx, ny)
  I = matrix(NA, n, 1)
  Z = matrix(NA, n, 1)
  distance = matrix(NA, n, 1)
  x = seq(minx1, maxx1, length = seq + 1)
  y = seq(miny1, maxy1, length = seq + 1)
  xy = expand.grid(x, y)
  xy$size = 0
  colnames(xy) = c("x", "y", "size")
  vaule = LAI.Vaule(minx, maxx, miny, maxy, xy, b, strata, r)
  vaule = vaule[complete.cases(vaule), ]
  vaule = vaule[, c(1, 2, 4)]
  colnames(vaule) = c("x", "y", "LAI")
  vaule = as.data.frame(vaule)
  for (i in 1:n) {
    vaule.dists = as.matrix(dist(cbind(vaule$x, vaule$y)))
    vaule.dists = (vaule.dists - indis)%/%(i * lag) * i * 
      lag + 0.5 * (i * lag + indis)
    vaule.dists[vaule.dists <= 0] = 0.5 * (i * lag + indis)
    vaule.dists.inv = 1/vaule.dists
    Moran = Moran.I(vaule$LAI, vaule.dists.inv)
    I[i] = Moran$observed
    Z[i] = (Moran$observed - Moran$expected)/Moran$sd
    distance[i] = lag * (i) + indis
    ISAA = cbind(distance, I, Z)
    colnames(ISAA) = c("lag", "Moran_I", "Z_Score")
  }
  vaule.dists = as.matrix(dist(cbind(vaule$x, vaule$y)))
  vaule.dists[vaule.dists < indis] = 0.5 * indis
  vaule.dists.inv = 1/vaule.dists
  Moran = Moran.I(vaule$LAI, vaule.dists.inv)
  I = Moran$observed
  Z = (Moran$observed - Moran$expected)/Moran$sd
  distance = indis
  ISAA_1 = cbind(distance, I, Z)
  colnames(ISAA_1) = c("lag", "Moran_I", "Z_Score")
  ISAA = rbind(ISAA_1, ISAA)
  yBL = max(ISAA[, 3])/max(ISAA[, 2])
  yBL = 1.1 * yBL
  ISAA = as.data.frame(ISAA)
  ISAAI = ISAA[, c(1, 2)]
  ISAAI$Index = c(" Moran_I")
  colnames(ISAAI) = c("Lag", "Moran_I", "Index")
  ISAAZ = ISAA[, c(1, 3)]
  ISAAZ$Index = c("Z_Score")
  colnames(ISAAZ) = c("Lag", "Moran_I", "Index")
  ISAAZ$Moran_I = ISAAZ$Moran_I/yBL
  ISAAIZ = rbind(ISAAI, ISAAZ)
  Z = ISAA[, 3]
  Z_positive = subset(Z, Z > 0)
  Zbef = matrix(NA, length(ISAA[, 3]) - 2, 1)
  Zbac = matrix(NA, length(ISAA[, 3]) - 2, 1)
  for (i in 1:(length(Z) - 2)) {
    Zbef[i] = Z[i + 1] - Z[i]
    Zbac[i] = Z[i + 2] - Z[i + 1]
  }
  n = length(vaule$LAI)
  ZBB = cbind(Zbef, Zbac)
  ZBB = cbind(ISAA$lag[-c(1, length(Z))], ZBB)
  ZBB = ZBB[which(ZBB[, 2] > 0 & ZBB[, 3] < 0), 1]
  Z_positive = cbind(ISAA$lag[-c(1, length(Z))], Z[-c(1, length(Z))])
  ZBB_positive = intersect(ZBB, Z_positive[which(Z_positive[, 
                                                            2] > 0), 1])
  p = ggplot(ISAAIZ, aes(x = Lag)) + 
    geom_hline(aes(yintercept = 1.96/yBL),linetype = 5, col = "black", size = 1) + 
    geom_hline(aes(yintercept = -1.96/yBL), linetype = 5, col = "black", size = 1) + 
    geom_hline(aes(yintercept = (-1/(n - 1))), linetype = 2, col = "blue", size = 1) + 
    geom_vline(xintercept = ZBB_positive, linetype = 5, col = "red", size = 1) + 
    geom_hline(aes(yintercept = 0), linetype = 1, col = "gray", size = 1) + 
    geom_line(aes(y = Moran_I,group = Index)) + 
    geom_point(aes(y = Moran_I, shape = Index), size = 2)
  p = p + scale_y_continuous(sec.axis = sec_axis(~. * yBL, name = "Z_Score")) + 
    theme_bw()
  p = p + theme(legend.title = element_blank(), legend.position = c(0.8,0.8))
  p
}

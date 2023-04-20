MaXent

library(tidyverse)
library(data.table)
library(cowplot)
library(spocc)
library(raster)
library(dismo)
library(rJava)
library(rgeos)
library(grDevices)
library(mapr)
library(rgdal)
library(boot)
library(ENMeval)
library(ROCR)
library(vcd)
library(sp)
library(maptools)
library(rgbif)
library(geodata)
library(sf)

##################Pruebas#######################

spp <- c("Candelaria concolor")
gb <- data.frame(Especie = spp)
keys <- 
        gb %>%
        pull("Especie") %>%
        name_backbone_checklist()  %>%
        filter(!matchType == "NONE") %>%
        pull(usageKey)

keys1 <- occ_download(
        pred_in("taxonKey", keys),
        pred("country", "CO"),
        pred("hasCoordinate", TRUE),
        pred("occurrenceStatus","PRESENT"),
        format = "SIMPLE_CSV",
        user="javier_biologia",pwd="Biologia20-",email="6Javier7@gmail.com"
)

keys2 <- grep("^",keys1, value = T) #Leptogium "0146333-230224095556074"
keys2 <- "0146333-230224095556074"
keys3 <- occ_download_wait(keys2)
keys4 <- grep("^",keys3, value = T)

status <- keys4[8]


#cambia cada descarga fijate en el nombre del mensaje despues de correr occ_download
keys2 <- "0146333-230224095556074"
d <- 
        occ_download_get(keys2) %>%
        occ_download_import()

d1 <- d[, c(10, 22, 23, 36)]
colnames(d1) <- c("especies", "lat", "lon", "basisOfRecord")
d1 <- subset(d1,!is.na(lat) & !is.na(lon))
d1 <- data.table(d1)
d1 <- d1[!especies == ""]

pru <- with(d1, data.table(especies, cor = paste(lat, lon))) 
pru <- pru %>% group_by(especies) %>% duplicated(by = key(cor))
d2 <- d1[!pru]
dp <- subset(d2, basisOfRecord == "PRESERVED_SPECIMEN")
with(dp, sort(table(especies)))
spp <- with(d1, unique(especies))



setwd("/Users/javiermontanochiriboga/Documents/Javier/Packages/Neutral1")
colombia1 <- getData('GADM' , country = "COL", level = 1)
#colombia1 <- gadm("COL", level = 1, path = tempdir())
#fra <- gadm(country="FRA", level=1, path =tempdir())
setwd("/Users/javiermontanochiriboga/Documents/Javier/Packages/Neutral1")
cli2 <- worldclim_country("Colombia", var = 'bio', res = 0.5, path = getwd())
#cli2 <- worldclim_country("Colombia", var = 'bio', res = 0.5, path = tempdir())
for (i in 1:16){assign(paste0("c", i), cli2[[i]])}

c17 <- elevation_30s("Colombia", path = getwd())
#c17 <- elevation_30s("Colombia", path = tempdir())
#c[[17]] <- c17

#clima2 <- mask(clima2, colombia1)
#climag <- elevation_global(0.5, path = "~/Downloads/")
#clima3 <- elevation_30s("Colombia", path = "~/Downloads/")

datos1 <- d2[, -c(1, 4)]
coordinates(datos1) <- ~lon + lat #espaciales para graficar
src <- CRS("+init=epsg:4326")
crs(datos1) <- src
colombia1 <- spTransform(colombia1, crs(datos1))

c1 <- crop(raster(c1), colombia1)
c2 <- crop(raster(c2), colombia1)
c3 <- crop(raster(c3), colombia1)
c4 <- crop(raster(c4), colombia1)
c5 <- crop(raster(c5), colombia1)
c6 <- crop(raster(c6), colombia1)
c7 <- crop(raster(c7), colombia1)
c8 <- crop(raster(c8), colombia1)
c9 <- crop(raster(c9), colombia1)
c10 <- crop(raster(c10), colombia1)
c11 <- crop(raster(c11), colombia1)
c12 <- crop(raster(c12), colombia1)
c13 <- crop(raster(c13), colombia1)
c14 <- crop(raster(c14), colombia1)
c15 <- crop(raster(c15), colombia1)
c16 <- crop(raster(c16), colombia1)
c17 <- crop(raster(c17), colombia1)


c1 <- mask(c1, colombia1)
c2 <- mask(c2, colombia1)
c3 <- mask(c3, colombia1)
c4 <- mask(c4, colombia1)
c5 <- mask(c5, colombia1)
c6 <- mask(c6, colombia1)
c7 <- mask(c7, colombia1)
c8 <- mask(c8, colombia1)
c9 <- mask(c9, colombia1)
c10 <- mask(c10, colombia1)
c11 <- mask(c11, colombia1)
c12 <- mask(c12, colombia1)
c13 <- mask(c13, colombia1)
c14 <- mask(c14, colombia1)
c15 <- mask(c15, colombia1)
c16 <- mask(c16, colombia1)
c17 <- mask(c17, colombia1)


cb <- list(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17)

clima <- stack(cb)

rm(list = listlist("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "cb", "d", "d1", "cli2"))

condiciones <- extract(clima, datos1)
head(condiciones)

sindatos <- is.na(condiciones[,1])
table(sindatos)

rm(sindatos)
rm(condiciones)

count = 0

dir.create("Graficos")

for (sp in spp[1:9]){
        
        count = count + 1
        ocu1 <- dp[especies == sp]
        #ocu1 <- dp[especies == spp[16]]
        #sp <- spp[16]
        
        if (dim(ocu1)[1] > 10){
                
                #esp[8]
                # "Leptogium phyllocarpum (Pers.) Mont."
                # esp[12]
                # "Ramalina calcarata Krog & Swinscow"
                
                datos1 <- ocu1[, -c(1, 4)]
                ocu1 <- ocu1[,-1]
                
                
                coordinates(datos1) <- ~lon + lat #espaciales para graficar
                src <- CRS("+init=epsg:4326")
                crs(datos1) <- src
                
                mapa <- spTransform(colombia1, crs(datos1))
                sobrelapan <- over(datos1, mapa)
                #countrydata <- subset(datos1, !is.na(sobrelapan$FIPS)) #solo los registros que sobrelapan el poligono
                celdas <- cellFromXY(clima[[1]], datos1)
                datos <- datos1[!duplicated(celdas),]
                
                #Buffer
                buffer <- buffer(datos, width = 100000) #bufer de 10km
                areas <- crop(clima, extent(buffer)) #un rectangulo del area de los puntos
                areas <- mask(areas, buffer) #el area exacta que coge todos los registros
                
                areas1 <- crop(clima, extent(mapa)) #un rectangulo del area de los puntos
                areas1 <- mask(areas1, mapa) #el area exacta que coge todos los registros
                
                #forma2
                k <- if(length(datos1$lat) >= 3){k = 3} else{k = 1}
                
                fold <- kfold(datos1, k = k)
                s <- sample(k, 1)
                ts <- datos1[fold == s, ] # el segundo bloque se usara para testear
                tr <- datos1[fold != s, ] # los otros cuatro bloques se utlizaran en el modelo
                if (length(tr@coords) > length(ts@coords)) {train <- tr
                test <- ts} else {test <- tr
                train <- ts}
                #plot(train, col = "blue")
                #plot(test, col = "red")
                
                
                #background
                set.seed(1)
                u <- runif(1)
                if (u < .6){bg <- sampleRandom(x = areas, size = 10000, na.rm = T, sp = T)} else {bg <- randomPoints(areas, 10000, p = datos)} #hacen lo mismo pero diferente
                plot(areas[[1]])
                plot(bg, add = T)
                plot(datos, add = T, col = "red")
                
                
                #se extraen las condiciones para la prueba, el training y el entorno
                #p <- raster::extract(clima, train[, -3])
                #ptest <- raster::extract(clima, test[, -3])
                p <- raster::extract(clima, train)
                ptest <- raster::extract(clima, test)
                a <- raster::extract(clima, bg)
                
                
                pa <- c(rep(1, nrow(p)), rep(0, nrow(a))) #presente en el test y ausente en el entorno
                pder <- as.data.frame(rbind(p, a))
                #pder[is.na(pder)] <- 0
                
                
                #mod <- dismo::maxent(x = pder, p = pa, factors = "biome", nbg = 10000,, args = c("-J", "-P"))
                mod <- dismo::maxent(x = pder, p = pa, factors = "biome", nbg = 5000,, args = c("-J", "-P"))
                
                response(mod)
                ped1 <- predict(mod, areas1)
                
                # for ggplot, we need the prediction to be a data frame 
                ped1 <- as(ped1, "SpatialPixelsDataFrame")
                dfp <- as.data.frame(ped1)
                xmax <- max(dfp$x)
                xmin <- min(dfp$x)
                ymax <- max(dfp$y)
                ymin <- min(dfp$y)
                
                wrld <- ggplot2::map_data("world")
                
                
                
                
                #evaluaciones
                p1 <- predict(mod, data.frame(ptest))
                a1 <- predict(mod, data.frame(a))
                combinado <- c(p1, a1)
                label <- c(rep(1, length(p1)), rep(0, length(a1)))
                predic <- prediction(combinado, label)
                #perfor <- performance(predic, "sens", "spec")
                perfor <- performance(predic, "tpr","fpr")
                auc <- performance(predic, "auc")
                
                #evaltrain <- dismo::evaluate(p = ptest, a = a, model = mod)
                #plot(evaltrain, "ROC")
                
                #bdf2 <- data.frame("Sensitivity" = perfor@y.values, "Specificity" = perfor@x.values, "Cutoff" = perfor@alpha.values)
                
                #colnames(bdf2) <- c("Sensitivity", "Specificity", "Cutoff")
                
                bdf2 <- data.frame("TPR" = perfor@y.values, "FPR" = perfor@x.values, "Cutoff" = perfor@alpha.values)
                
                colnames(bdf2) <- c("TPR", "TFR", "Cutoff")
                bdf2 <- data.table(bdf2)
                
                
                
                #th3 <- threshold(evaltrain, stat = "sensitivity", sensitivity = 0.95)
                
                AUC <- function(data = list(ptest, a), i) {
                        p1 <- predict(mod, data.frame(ptest[i,]))
                        a1 <- predict(mod, data.frame(a))
                        combinado <- c(p1, a1)
                        label <- c(rep(1, length(p1)), rep(0, length(a1)))
                        predic <- prediction(combinado, label)
                        
                        auc <- performance(predic, "auc")
                        result <- list()
                        result[[1]] <- auc@y.values[[1]]
                        result[[2]] <- sd(perfor@y.values[[1]])
                        return(result[[1]])
                }
                
                #bo <- boot(ptest, AUC, 500)
                bo <- boot(ptest, AUC, 100)
                d <- density(bo$t)
                
                bdf <- data.table(id = 1:length(bo$t), boot = bo$t)
                colnames(bdf) <- c("id", "boot")
                
                breaks <- pretty(range(bo$t), n = nclass.FD(bo$t), min.n = 1)
                bwidth <- breaks[2]-breaks[1]
                
                p2 <- ggplot(bdf, aes(x = boot)) + 
                        geom_histogram(aes(y = ..density..), colour = "grey51", fill = "grey87", binwidth=bwidth) + 
                        geom_density(alpha = .3, linetype = "dashed", color = "coral", linewidth = 1) +
                        geom_vline(aes(xintercept = mean(boot)),
                                   color = "coral", linewidth = 1) +
                        #geom_vline(aes(xintercept = median(boot)),
                        #      color = "orchid", linetype = "dashed", linewidth = 1) +
                        annotate("text", label = paste("Boot Prom. =", with(bo, round(t0, 3))), x = quantile(d$x, .75), y = max(d$y), size = 4, color = "coral") +
                        theme_classic() + labs(x = "AUC", y = "Densidad", title = paste0("Histograma Boot ", sp)) +
                        theme(plot.title = element_text(colour = "coral"), axis.title = element_text(colour = "grey58"))
                
                p2
                ggsave(paste0(paste0(sp, count), " boot hist.png"), path = "/Users/javiermontanochiriboga/Documents/Javier/Packages/Neutral1/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2)
                
                
                p3 <- ggplot(bdf, aes(sample = boot))+ stat_qq_line(colour = "grey58", linewidth = 1.0) + theme_classic() + stat_qq(colour = "coral", size = 1.6) + labs(x = "", y = "", title = paste0("QQ plot boot "), sp) +
                        theme(plot.title = element_text(colour = "coral"), axis.title = element_text(colour = "grey58"))
                
                p3
                ggsave(paste0(paste0(sp, count)," boot qqplot.png"), path = "/Users/javiermontanochiriboga/Documents/Javier/Packages/Neutral1/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2) 
                
                p1 <- ggplot() +
                        geom_polygon(data = wrld, mapping = aes(x = long, y = lat, group = group), fill = "white") + geom_raster(data = dfp, aes(x = x, y = y, fill = layer)) +  scale_fill_gradientn(colors = terrain.colors(10, rev = T))  +  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = F) + # expand = F fixes weird margin
                        scale_size_area() +
                        theme_classic() +
                        borders("state") + labs(title = paste("SDM", sp), x = "Longitud", y = "Latitud", fill = "Pontencial") + theme(legend.box.background = element_blank(), legend.box.margin = margin(5,5,5,5), plot.title = element_text(colour = "coral"), axis.title = element_text(colour = "grey58"), legend.key = element_rect(fill = "white", colour = "white"), legend.title = element_text(colour = "grey58"))
                
                p1
                ggsave(paste0(paste0(sp, count)," potencial.png"), path = "/Users/javiermontanochiriboga/Documents/Javier/Packages/Neutral1/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2)
                
                inter <- boot.ci(bo, type = "basic")
                
                
                p4 <- ggplot() +     
                        geom_line(data = bdf2, aes(x = TFR, y = TPR), linetype = "dashed", color = "coral", linewidth = 0.8) +
                        geom_point(data = bdf2, aes(x = TFR, y = TPR), shape = 16, color = "coral", size = 2) +
                        geom_line(aes(x = c(0, 1), y = c(0, 1)), linetype = "dashed", color = "grey58", linewidth = 1) +
                        annotate("text", label = paste0(paste0(paste0(paste0(paste0("AUC = ", round(auc@y.values[[1]], 3)), "\nStd. error = "), round(with(bo, sd(t)), 4)), "\nIC 95% = (", paste(paste0(round(inter$basic[1,][4], 3), ","), paste(round(inter$basic[1,][5], 3), ")")))), x = 0.775, y = 0.15, size = 4, color = "coral") +
                        labs(x = "Tasa de Falsos Positivos(1-Especificidad)", y = "Tasa de Verdaderos Positivos(Sensitividad)", title = paste("Curva ROC ", sp)) +
                        theme_classic() +
                        theme(plot.title = element_text(colour = "coral"), axis.title = element_text(colour = "grey58"))
                
                p4
                ggsave(paste0(paste0(sp, count)," Curva AUC.png"), path = "/Users/javiermontanochiriboga/Documents/Javier/Packages/Neutral1/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2)
                
                #p <- plot_grid(p1, p4, p2, p3, labels = c("A", "B", "C", "D"), label_size = 12)
                
                #png(paste0(sp, count), height = 2000, width = 3000)
                #plot_grid(p1, p2, p3, labels = c("A", "B", "C"), label_size = 12)
                #dev.off()
                
                #ggsave(paste0(paste0(sp, count),".png"), width = 26, height = 14, dpi = 700, units = "cm", limitsize = FALSE, scale = 2)
                
                #save_plot(paste0(paste0(sp, count),".png"), ncol = 2, nrow = 2,  base_width = 26, base_height = 14, dpi = 700)
                
                #save_plot(paste0(paste0(sp, count),".png"), p, nrow = 2, ncol = 2)}else{
                print(paste0(paste0("La especie ", sp), " tiene menos de 5 registros"))
        }
}

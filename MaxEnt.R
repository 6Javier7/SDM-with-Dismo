#MaXent

#Paquetes

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

#Datos

usuario <- "Biologx" #cuenta de gbif
contras <- "AsiliiPower" #contraseña en gbif
correo <- "Biologx@mail.com" #email en gbif

spp <- c("Candelaria concolor")
gb <- data.frame(Especie = spp)
keys <- 
        gb %>%
        pull("Especie") %>%
        name_backbone_checklist()  %>%
        filter(!matchType == "NONE") %>%
        pull(usageKey) #crea un archivo con los codigos del gbif

keys1 <- occ_download(
        pred_in("taxonKey", keys),
        pred("country", "CO"),
        pred("hasCoordinate", TRUE),
        pred("occurrenceStatus","PRESENT"),
        format = "SIMPLE_CSV",
        user = usuario, pwd = contras, email = correo
) #para descargar los registros

keys2 <- grep("^", keys1, value = T)
keys3 <- occ_download_wait(keys2) #revisar si ya se descargaron
keys4 <- grep("^",keys3, value = T)

status <- keys4[8]


#cambia cada descarga fijate en el nombre del mensaje despues de correr occ_download
d <- 
        occ_download_get(keys2) %>%
        occ_download_import() #archivo con los datos

d1 <- d[, c(10, 22, 23, 36)] #solo dejar el nombre, las coordenadas y el tipo de observacion
colnames(d1) <- c("especies", "lat", "lon", "basisOfRecord")
d1 <- subset(d1,!is.na(lat) & !is.na(lon)) #filtrar datos sin coordenadas en x o en y
d1 <- data.table(d1) #cambia de data.frame a data.table
d1 <- d1[!especies == ""] #quitar los registros sin nobre cinetifico

pru <- with(d1, data.table(especies, cor = paste(lat, lon))) #la idea es crear pru para filtar duplicaods(registros de las mismas coordenads)
pru <- pru %>% group_by(especies) %>% duplicated(by = key(cor))
d2 <- d1[!pru] #filtro de duplicados
dp <- subset(d2, basisOfRecord == "PRESERVED_SPECIMEN") #filtro de especimes de herbario
with(dp, sort(table(especies))) #numero de registros por especie
spp <- with(d1, unique(especies))



colombia1 <- getData('GADM' , country = "COL", level = 1) #shape de Colombia
cli2 <- worldclim_country("Colombia", var = 'bio', res = 0.5, path = getwd()) #datos clima para Colombia
for (i in 1:16){assign(paste0("c", i), cli2[[i]])} #separacion de los raster del clima

c17 <- elevation_30s("Colombia", path = getwd()) #Raster de alturas

datos1 <- d2[, -c(1, 4)] #datos1 solo tiene las coordenadas
coordinates(datos1) <- ~lon + lat #crea una estructura espacial de coordenadas
src <- CRS("+init=epsg:4326") #establece el sistema de referencia
crs(datos1) <- src
colombia1 <- spTransform(colombia1, crs(datos1)) #cambia la proyeccion del shape a la del crs elegido anteriormente

#Se corta cada raster del clima para que tenga la misma extencion que el mapa de Colombia
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

#El mapa cortado anteriormente se refina para que tenga la misma forma del mapa de Colombia
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

#se unen los raster en una lista
cb <- list(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17)

clima <- stack(cb) #se combinan todos los rasters

rm(list = listlist("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "c15", "c16", "c17", "cb", "d", "d1", "cli2"))

condiciones <- extract(clima, datos1) #es una revision para ver si se tienen datos en las coordenadas de los registros
head(condiciones)

sindatos <- is.na(condiciones[,1]) 
table(sindatos) #Analizar el numero de registros sin datos climaticos

rm(sindatos)
rm(condiciones)

count = 0

dir.create("Graficos") #Se crea una carpeta para guardar las graficas que se van a crear

for (sp in spp){ #se repite el mismo proceso por cada especie en spp
        
        count = count + 1
        ocu1 <- dp[especies == sp]
        
        if (dim(ocu1)[1] > 10){ #se filtran las especies con menos de 10 registros
                
                datos1 <- ocu1[, -c(1, 4)]
                ocu1 <- ocu1[,-1]
                coordinates(datos1) <- ~lon + lat #espaciales para graficar
                src <- CRS("+init=epsg:4326")
                crs(datos1) <- src
                
                mapa <- spTransform(colombia1, crs(datos1))
                sobrelapan <- over(datos1, mapa)
                
                #Buffer
                buffer <- buffer(datos, width = 100000) #bufer de 10km
                areas <- crop(clima, extent(buffer)) #un rectangulo del area de los puntos
                areas <- mask(areas, buffer) #el area exacta que coge todos los registros
                
                areas1 <- crop(clima, extent(mapa)) #un rectangulo del area de los puntos
                areas1 <- mask(areas1, mapa) #el area exacta que coge todos los registros
                
                #forma2
                k <- if(length(datos1$lat) >= 3){k = 3} else{k = 1} 
                
                fold <- kfold(datos1, k = k) #se dividen los datos en tres aleatoriamente
                s <- sample(k, 1) #se elige un tercio de los datos de manera aleatoria
                ts <- datos1[fold == s, ] # el segundo bloque se usara para testear
                tr <- datos1[fold != s, ] # los otros cuatro bloques se utlizaran en el modelo
                if (length(tr@coords) > length(ts@coords)) {train <- tr
                test <- ts} else {test <- tr
                train <- ts}
                
                
                #background
                set.seed(1)
                u <- runif(1)
                if (u < .6){bg <- sampleRandom(x = areas, size = 10000, na.rm = T, sp = T)} else {bg <- randomPoints(areas, 10000, p = datos)} #hacen lo mismo pero diferente
                plot(areas[[1]])
                plot(bg, add = T)
                plot(datos, add = T, col = "red")
                
                
                p <- raster::extract(clima, train)
                ptest <- raster::extract(clima, test)
                a <- raster::extract(clima, bg)
                
                
                pa <- c(rep(1, nrow(p)), rep(0, nrow(a))) #presente en el test y ausente en el entorno
                pder <- as.data.frame(rbind(p, a))
                
                mod <- dismo::maxent(x = pder, p = pa, factors = "biome", nbg = 5000, args = c("-J", "-P")) #modelo de distribucion de especies
                
                response(mod) #Ponderacion de la variacion de cada variables
                ped1 <- predict(mod, areas1) #Extrapolacion del modelo en todo Colombia
                
                # Para graficar se cambia la estructura de la prediccion a una data.frame 
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
                
                perfor <- performance(predic, "tpr","fpr") 
                auc <- performance(predic, "auc") #auc del modelo
                
                bdf2 <- data.frame("TPR" = perfor@y.values, "FPR" = perfor@x.values, "Cutoff" = perfor@alpha.values)
                
                colnames(bdf2) <- c("TPR", "TFR", "Cutoff")
                bdf2 <- data.table(bdf2)
                
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
                
                
                bo <- boot(ptest, AUC, 100) #Bootrap del modelo
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
                        annotate("text", label = paste("Boot Prom. =", with(bo, round(t0, 3))), x = quantile(d$x, .75), y = max(d$y), size = 4, color = "coral") +
                        theme_classic() + labs(x = "AUC", y = "Densidad", title = paste0("Histograma Boot ", sp)) +
                        theme(plot.title = element_text(colour = "coral"), axis.title = element_text(colour = "grey58"))
                
                p2
                ggsave(paste0(paste0(sp, count), " boot hist.png"), path = "/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2)
                
                
                p3 <- ggplot(bdf, aes(sample = boot))+ stat_qq_line(colour = "grey58", linewidth = 1.0) + theme_classic() + stat_qq(colour = "coral", size = 1.6) + labs(x = "", y = "", title = paste0("QQ plot boot "), sp) +
                        theme(plot.title = element_text(colour = "coral"), axis.title = element_text(colour = "grey58"))
                
                p3
                ggsave(paste0(paste0(sp, count)," boot qqplot.png"), path = "/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2) 
                
                p1 <- ggplot() +
                        geom_polygon(data = wrld, mapping = aes(x = long, y = lat, group = group), fill = "white") + geom_raster(data = dfp, aes(x = x, y = y, fill = layer)) +  scale_fill_gradientn(colors = terrain.colors(10, rev = T))  +  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = F) + # expand = F fixes weird margin
                        scale_size_area() +
                        theme_classic() +
                        borders("state") + labs(title = paste("SDM", sp), x = "Longitud", y = "Latitud", fill = "Pontencial") + theme(legend.box.background = element_blank(), legend.box.margin = margin(5,5,5,5), plot.title = element_text(colour = "coral"), axis.title = element_text(colour = "grey58"), legend.key = element_rect(fill = "white", colour = "white"), legend.title = element_text(colour = "grey58"))
                
                p1
                ggsave(paste0(paste0(sp, count)," potencial.png"), path = "/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2)
                
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
                ggsave(paste0(paste0(sp, count)," Curva AUC.png"), path = "/Graficos", width = 14, height = 7, dpi = 700, units = "cm", limitsize = FALSE, scale = 2)}else{
                print(paste0(paste0("La especie ", sp), " tiene menos de 5 registros"))
        }
}

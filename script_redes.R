
tr<-read.csv("transectos.csv", sep=",", header=TRUE)
str(tr)


library(bipartite)


sites <- unique(tr$Bosque)


out.site <- data.frame(Site_id = NA, Recorrido=NA, Periodo=NA, Anidamiento=NA, Uniformidad=NA, H2=NA,
                       species.poll=NA, Especies_plantas=NA, 
                       compl_fun_polinizadores=NA, compl_fun_plantas=NA, Robustez=NA)

outsp.site <- data.frame(Site_id = NA, Recorrido=NA, Periodo=NA, Especies=NA,  Grado = NA,
                         Centralidad=NA, Especialización=NA)


outsp.pl.site <- data.frame(Site_id = NA, Recorrido=NA, Periodo=NA, Especies=NA,  Grado = NA,
                            Centralidad=NA, Especialización=NA)



webs <- list()
for(i in 1:length(sites)){
  print(i)
  temp <- subset(tr[,c(1,3:7)], Bosque == sites[i])
  
  recorrido<-unique(tr$Transecto)
  
  
  for(j in 1:length(recorrido)){
    
    temp2<-subset(temp, temp$Transecto==recorrido[j])
    periodo<-unique(temp2$Periodo)
    for(k in 1:length(periodo)){
      temp3<-subset(temp2, temp2$Periodo==periodo[k])
      
      web<-table(temp3$Planta, temp3$Polinizador)
      
      
      plotweb(web, col.interaction = "grey", col.high = "blue", col.low = "yellow", text.rot =90 )
      spntw <- try(specieslevel(web), TRUE)
      ntw <- networklevel(web)
      if(isTRUE(class(spntw)=="try-error")) {next} 
      
      ex <- second.extinct(web, participant="lower", method="random", nrep=100, details=FALSE)
      r<- robustness(ex)
      
      n <- nrow(out.site)
      webs[[n]] <- web  
      n2 <- nrow(outsp.site)
      n3 <- nrow(outsp.pl.site)
      
     
      
      out.site[n + 1,1] <- as.character(sites[i])
      out.site[n + 1,2] <- as.character(recorrido[j])
      out.site[n + 1,3] <- as.character(periodo[k])
      out.site[n + 1,4:10] <- c(ntw[10], ntw[17], ntw[19], ntw[20:21], ntw[42:43])
      out.site[n + 1,11] <-r
      
      
      outsp.site[n2+seq(nrow(spntw$`higher level`)),1] <- as.character(sites[i])
      outsp.site[n2+seq(nrow(spntw$`higher level`)),2] <- as.character(recorrido[j])
      outsp.site[n2+seq(nrow(spntw$`higher level`)),3] <- as.character(periodo[k])
      outsp.site[n2+seq(nrow(spntw$`higher level`)),4] <- try(as.character(rownames(spntw$`higher level`)), TRUE)
      outsp.site[n2+seq(nrow(spntw$`higher level`)),5:7] <- try(c(spntw$`higher level`[2], 
                                                                  spntw$`higher level`[14],
                                                                  spntw$`higher level`[20]), TRUE)
      
      
      
      outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),1] <- as.character(sites[i])
      outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),2] <- as.character(recorrido[j])
      outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),3] <- as.character(periodo[k])
      outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),4] <- try(as.character(rownames(spntw$`lower level`)), TRUE)
      outsp.pl.site[n3+seq(nrow(spntw$`lower level`)),5:7] <- try(c(spntw$`lower level`[2], 
                                                                    spntw$`lower level`[14],
                                                                    spntw$`lower level`[20]), TRUE)
      
    }  }}

str(out.site)
head(out.site)
out.site


head(outsp.site)
str(outsp.site)

head(outsp.pl.site)
str(outsp.pl.site)

write.table(out.site[2:nrow(out.site),], file= "SITE_network_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.site[2:nrow(outsp.site),], file= "SITE_species_level_metrics.csv", row.names= FALSE, sep= ",")
write.table(outsp.pl.site[2:nrow(outsp.pl.site),], file= "SITE_plant_species_level_metrics.csv", row.names= FALSE, sep= ",")


### CREO 3 TABLAS (PARA LA RED, PARA LAS ESPECIES DE
# POLINIZADORES Y PARA LAS ESPECIES DE PLANTAS) Y LAS 
# RENOMBRO COMO SIGUE PARA HACER LOS DIFERENTES ANALISIS.

sitems<-read.csv("SITE_network_level_metrics.csv")
sitems1<-read.csv("SITE_species_level_metrics.csv")
sitems2<-read.csv("SITE_plant_species_level_metrics.csv")

library(Matrix)
library(lme4)
library(carData)
library(effects)
library(performance)
library(see)
library(gridExtra)
library(car)
library(ggplot2)
library(lattice)

### TABLA RED: HAGO FACTOR LO QUE DEDE SER FACTOR
# (SITIO, RECORRIDO Y PERIODO).
sitems$Site_id<-as.factor(sitems$Site_id)
sitems$Recorrido<-as.factor(sitems$Recorrido)
sitems$Periodo<-as.factor(sitems$Periodo)

# HISTOGRAMAS DE MIS VARIABLES.
hist(sitems$species.poll)
hist(sitems$Especies_plantas)
hist(sitems$Anidamiento)
hist(sitems$Uniformidad)
hist(sitems$H2)
hist(sitems$Robustez)
hist(sitems$compl_fun_polinizadores)
hist(sitems$compl_fun_plantas)

# MODELO 1
m1<-lmer(species.poll ~ Recorrido + Periodo + (1|Site_id),  data=sitems)
# COMPRUEBO NORMALIDAD RESIDUOS.
check_model(m1)
# SIGUEN DISTRIBUCION NORMAL.
Anova(m1)
# GRAFICO RESULTADOS
allEffects(m1)
plot(allEffects(m1))

# MODELO 2
m2<-lmer(Especies_plantas ~ Recorrido + Periodo + (1|Site_id), data=sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m2)
# SIGUEN DISTRIBUCION NORMAL
Anova(m2)
# GRAFICO RESULTADOS
allEffects(m2)
plot(allEffects(m2))

# MODELO 3 - ANIDAMIENTO
m3<-lmer(Anidamiento ~ Recorrido + Periodo + (1|Site_id), data = sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m3)
# NO SIGUEN UNA DISTRIBUCION NORMAL

# MODELO 4 - ANIDAMIENTO
m4<-glmer(Anidamiento ~ Recorrido + Periodo + (1|Site_id), family = "poisson", data = sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
qqmath(m4)
# NO SIGUEN UNA DISTRIBUCION NORMAL

# MODELO 5 - ANIDAMIENTO
m5<-glmer.nb(Anidamiento ~ Recorrido + Periodo + (1|Site_id), data = sitems)
Anova(m5)
# RESUMEN MODELO 5 Y GRAFICO RESULTADOS
summary(m5)
allEffects(m5)
plot(allEffects(m5))

# MISMO MODELO QUE EL ANTERIOR PERO CON VARIABLE
# "ESPECIES PLANTAS" COMO VARIABLE EXPLICATIVA.
m6<-glmer.nb(Anidamiento ~ Recorrido + Periodo + Especies_plantas + (1|Site_id), data=sitems)
# RESUMEN MODELO 6
summary(m6)

# MODELO 7 - UNIFORMIDAD
m7<-lmer(Uniformidad ~ Recorrido + Periodo + (1|Site_id), data=sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m7)
# SIGUEN DISTRIBUCION NORMAL
Anova(m7)
# GRAFICO RESULTADOS
allEffects(m7)
plot(allEffects(m7))

# MODELO 8 - ESPECIALIZACION COMP.
m8<-lmer(H2 ~ Recorrido + Periodo + (1|Site_id), data=sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m8)
# SIGUEN DISTRIBUCION NORMAL
Anova(m8)
# GRAFICO RESULTADOS
allEffects(m8)
plot(allEffects(m8))

# MODELO 9 - COMPLEMENTARIEDAD FUNC. POL.
m9<-lmer(compl_fun_polinizadores ~ Recorrido + Periodo + Especies_plantas + (1|Site_id), data=sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m9)
#NO SIGUEN DISTRIBUCION NORMAL

# MODELO 10 - COMPLEMENTARIEDAD FUNC. POL.
m10<-glmer(compl_fun_polinizadores ~ Recorrido + Periodo + Especies_plantas + (1|Site_id), family = "poisson", data=sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
qqmath(m10)
# PARECE QUE SIGUEN UNA NORMALIDAD
# RESUMEN NO DEJA VER NADA
summary(m10)

# MODELO 11 - COMPLEMENTARIEDAD FUNC. POL.
m11<-glmer.nb(compl_fun_polinizadores ~ Recorrido + Periodo + Especies_plantas + (1|Site_id), data=sitems)
# RESIDUOS
qqmath(m11)
Anova(m11)
# RESUMEN MODELO 11 Y GRAFICOS
summary(m11)
allEffects(m11)
plot(allEffects(m11))

# MODELO 12 - ROBUSTEZ A LA PERDIDA DE SPP.
m12<-lmer(Robustez ~ Recorrido + Periodo + (1|Site_id), data = sitems)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m12)
# SIGUEN DISTRIBUCION NORMAL
car::Anova(m12)
# GRAFICOS
allEffects(m12)
plot(allEffects(m12))

### TABLA POLINIZADORES: HAGO FACTOR LO QUE DEDE SER FACTOR
# (SITIO, RECORRIDO Y PERIODO).
sitems1$Site_id<-as.factor(sitems1$Site_id)
sitems1$Recorrido<-as.factor(sitems1$Recorrido)
sitems1$Periodo<-as.factor(sitems1$Periodo)

# HISTOGRAMAS DE MIS VARIABLES
hist(sitems1$Grado)
hist(sitems1$Centralidad)
hist(sitems1$Especialización)

# MODELO 13 - GRADO
m13<-lmer(Grado ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems1)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m13)
# SIGUEN DISTRIBUCION NORMAL
Anova(m13)
# GRAFICOS
allEffects(m13)
plot(allEffects(m13))


sitems1$Especies[sitems1$Especies=="Bombus pascurorum"]<-"Bombus pascuorum"
sitems1$Especies[sitems1$Especies=="Zygaena sp."]<-"Zygaena lonicerae"

## DADO QUE AUN QUEDAN ESPECIES DE POLINIZADORES POR IDENTIFICAR,
## SELECCIONAMOS LAS MAS COMUNES Y CALCULAMOS LAS DIFERENTES 
## METRICAS PARA ESAS.

selected<-c("Apis mellifera", "Bombus pascuorum", 
"Bombus terrestris", "Sphaerophoria scripta",
"Episyrphus balteatus", "Eristalis tenax",
"Bombus pratorum", "Erebia meolans",  "Bombus hortorum")

sitems1_sub<-sitems1[sitems1$Especies %in% selected,]

# MODELO 13b - GRADO SPP. MAS COMUNES
m13b<-lmer(Grado ~ Recorrido + Periodo  + (1|Especies) + (1|Site_id),  data=sitems1_sub)
# COMPRUEBO NORMALIDAD RESIDUOS
check_model(m13b)
qqmath(m13b)
# SIGUEN DISTRIBUCION NORMAL
Anova(m13b)
# GRAFICOS
allEffects(m13b)
plot(allEffects(m13b))
ggplot(sitems1_sub, aes(x=Periodo, y=Grado, fill=Especies)) + 
  geom_boxplot() + theme_classic(base_size = 10)

# MODELO 14 - CENTRALIDAD
m14<-lmer(Centralidad ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems1)
# RESIDUOS
check_model(m14)
# NO SIGUEN NORMALIDAD

# MODELO 15 - CENTRALIDAD
m15<-glmer(Centralidad ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), family = "poisson", data=sitems1)
# RESIDUOS
qqmath(m15)
# NO SIGUEN NORMALIDAD

# MODELO 16 - CENTRALIDAD
m16<-glmer.nb(Centralidad ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems1)
Anova(m16)
# RESUMEN Y GRAFICOS
summary(m16)
allEffects(m16)
plot(allEffects(m16))

# MODELO 16b - CENTRALIDAD SPP. MAS COMUNES
m16b<-glmer.nb(Centralidad ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems1_sub)
# RESIDUOS
qqmath(m16b)
Anova(m16b)
# RESUMEN Y GRAFICOS
summary(m16b)
allEffects(m16b)
plot(allEffects(m16b))
ggplot(sitems1_sub, aes(x=Periodo, y=Centralidad, fill=Especies)) + 
  geom_boxplot()

# MODELO 17 - ESPECIALIZACION
m17<-lmer(Especialización ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems1)
# RESIDUOS
check_model(m17)
# NO SIGUEN NORMALIDAD

# MODELO 18 - ESPECIALIZACION
m18<-glmer(Especialización ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), family = "poisson", data=sitems1)
# RESIDUOS
qqmath(m18)
# PARECE QUE SON NORMALES
# RESUMEN NO ME DEJA VER NADA
summary(m18) 

# MODELO 19 - ESPECIALIZACION
m19<-glmer.nb(Especialización ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems1)
# RESIDUOS
qqmath(m19)

Anova(m19)
# RESUMEN Y GRAFICOS
summary(m19)
allEffects(m19)
plot(allEffects(m19))

# MODELO 19b - ESPECIALIZACION SPP. MAS COMUNES
m19b<-glmer.nb(Especialización ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems1_sub)
# RESIDUOS
qqmath(m19b)

Anova(m19b)
# RESUMEN Y GRAFICOS
summary(m19b)
allEffects(m19b)
plot(allEffects(m19b))
ggplot(sitems1_sub, aes(x=Periodo, y=d, fill=species)) + 
  geom_boxplot()


### TABLA PLANTAS: HAGO FACTOR LO QUE DEDE SER FACTOR
# (SITIO, RECORRIDO Y PERIODO).
sitems2$Site_id<-as.factor(sitems2$Site_id)
sitems2$Recorrido<-as.factor(sitems2$Recorrido)
sitems2$Periodo<-as.factor(sitems2$Periodo)

# HISTOGRAMAS DE MIS VARIABLES
hist(sitems2$Grado)
hist(sitems2$Centralidad)
hist(sitems2$Especialización)

# PLANTAS MAS COMUNES
selected<-c("Bellis perennis", "Potentilla erecta", "Taraxacum sp.", "Thymus sp.")

sitems2_sub<-sitems2[sitems2$Especies %in% selected,]

# MODELO 20 - GRADO
m20<-lmer(Grado ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems2)
# RESIDUOS
check_model(m20)
# NO SON NORMALES

# MODELO 21 - GRADO
m21<-glmer(Grado ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), family = "poisson", data=sitems2)
# RESIDUOS
qqmath(m21)
# SE PUEDE ASUMIR NORMALIDAD
summary(m21)
# RESUMEN NO SE VE

# MODELO 22 - GRADO
m22<-glmer.nb(Grado ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems2)
qqmath(m21)

Anova(m22)
# RESUMEN Y GRAFICOS
summary(m22)
allEffects(m22)
plot(allEffects(m22))

m22b<-glmer.nb(Grado ~ Recorrido + Periodo  + (1|Especies) + (1|Site_id),  data=sitems2_sub)
Anova(m22b)
allEffects(m22b)
plot(allEffects(m22b))
ggplot(sitems2_sub, aes(x=Periodo, y=Grado, fill=Especies)) + 
  geom_boxplot() + theme_classic(base_size = 10)


# MODELO 23 - CENTRALIDAD
m23<-lmer(Centralidad ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems2)
# RESIDUOS
check_model(m23)
# NO SON NORMALES

# MODELO 24 - CENTRALIDAD
m24<-glmer(Centralidad ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), family = "poisson", data=sitems2)
# RESIDUOS
qqmath(m24)
# NO SON NORMALES

# MODELO 25 - CENTRALIDAD
m25<-glmer.nb(Centralidad ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems2)

Anova(m25)
# RESUMEN Y GRAFICOS
summary(m25)
allEffects(m25)
plot(allEffects(m25))

m25b<-glmer.nb(Centralidad ~ Recorrido + Periodo  + (1|Especies) + (1|Site_id),  data=sitems2_sub)
Anova(m25b)

# MODELO 26 - ESPECIALIZACION
m26<-lmer(Especialización ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems2)
# RESIDUOS
check_model(m26)
# NO SON NORMALES

# MODELO 27 - ESPECIALIZACION
m27<-glmer(Especialización ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), family = "poisson", data=sitems2)
# RESIDUOS
qqmath(m27)
# SE PUEDE ASUMIR NORMALIDAD
summary(m27)
# NO SE VE NADA

# MODELO 28 - ESPECIALIZACION
m28<-glmer.nb(Especialización ~ Recorrido + Periodo + (1|Especies) + (1|Site_id), data=sitems2)
# RESIDUOS
qqmath(m28)

Anova(m28)
# RESUMEN Y GRAFICOS
summary(m28)
allEffects(m28)
plot(allEffects(m28))

m28b<-glmer.nb(Especialización ~ Recorrido + Periodo  + (1|Especies) + (1|Site_id),  data=sitems2_sub)
Anova(m28b)

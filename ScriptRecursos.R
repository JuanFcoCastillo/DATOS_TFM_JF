#renombro mi tabla de datos
r<-Rec_Flo2

#convierto mi variable Bosque que son nUmeros en un factor
r$Bosque1<-factor(r$Bosque)

library(dplyr)
#me quedo con las filas para valores medios
r1<-filter(r, filtro == "A")
r1

#Hago 5 nuevas tablas, una para cada bosque
b1<-filter(r, Bosque == 1)
b1

b2<-filter(r, Bosque == 2)
b2

b3<-filter(r, Bosque == 3)
b3

b4<-filter(r, Bosque == 4)
b4

b5<-filter(r, Bosque == 5)
b5



library(ggplot2)

#numero de especies de plantas juntando todos los bosques por fecha
ggplot(r, aes(y=N_especies, x=Fecha, colour=Bosque1, group=Bosque1, shape=Bosque1)) + 
  geom_point() + geom_line() + ylab("N? especies plantas") + xlab("Fecha") + theme_classic(base_size = 15)

#numero medio de inflorescencias por fecha y por bosque
ggplot(r, aes(y=N_inf_medio, x=Fecha, colour=Bosque1, group=Bosque1, shape=Bosque1)) + 
  geom_point() + geom_line() + ylab("N? medio inflorescencias") + xlab("Fecha") + theme_classic(base_size = 15)


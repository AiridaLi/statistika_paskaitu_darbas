#Priskyrimas
x<-c(1, 2, 8, 3)
x
(x<-c(1, 2, 8, 3))
length(x)
assign("y", c(0, -2, 4, 8, 5))
y#kai n?ra priskyrimo
1/y#kai skiriasi ilgiai
length(x)
length(y)
z<- 2*x-3*y
length(z)#sekos sukurimas
1:10
seq(1,10)
seq(-10,-10)
seq(1,10, by=0.5)
seq(-10,10, by=2)
y
(s1<-rep(y, times=3)) #kartoji vektoriu tris kartus
(s2<-rep(y, each=3)) #paimi pirma kartoji kiekviena elementa tris kartus#priskyrimas su pav. skaiciu tiek kiek pavadinimu!
kaina<- c(25,15,3,0.5,47,86,12,4,6,34)
names(kaina)<-c("batai", "kremas", "suris", "druska", "rankine", "tortas", "mesa", "darzoves", "auskarai", "vitaminai")
kainamaistas<-kaina[c("suris","druska","tortas","mesa","darzoves")] #r mato skc
maistas
sum(maistas)attributes(kaina)class(kaina)#logic < integar < double < complex < character(factor) hierarchine struktura, zemesnis tampa aukstesnes h. elementu#logic x
x>7
x
y
x>y 
(x>3)&(x<10) #& ir
(x>3)|(x<10) #| arba
sum((x>3)&(x<10))
as.integer((x>3)&(x<10)) #true 1 false 0u<-c(1,2,NA,3,5)
is.na(u) #ar yra praleistu reiksmiu?
any(is.na(u)) #pasako ar ya bent viena praleista
na.omit(u) #pasalina praleistas reiksmes is vektoriaus0/0 #NaN Not a number
Inf-Infdta<-c(1,2,3,2,3,1,1) #tai kas pazymeta skc ne butinai yra skaicius, gali buti uzkoduota zodziai ir pan
dta
class(dta)cdta<-as.character(dta)
cdta
class(cdta)fdta<-factor(dta) #levelius suranda (pasikartojancias reiksmes), patogu dirbti nes matai skirtingas reiksmes
fdta
class(fdta)ofdta<-factor(dta, levels=1:3, ordered = T) #klase pasikeicia, atsiranda hierarchija
ofdta
class(ofdta)rdta<-factor(dta, labels=c("I", "II", "III"))
rdta
class(rdta)
m1<-matrix(1:6, nrow=2, ncol=3, byrow=T) #byrow skiriasi matricoje skc tvarka
m1
m2<-matrix(1:6, nrow=2, ncol=3, byrow=F)
m2
m3<-array(1:6, c(2,3))
m3mt<-matrix(1:12, nrow=3, ncol=4)
mt
dim(mt)m1[1,2]
m2[2,3]
m1[,2] #visos eilutes antras stulpelis
m1[2,] #anta eilute visi stulpeliai
m4<-array(1:24, c(4,3,2)) #c trecias skaicius rodo kiek lenteliu
m4
m4[,,2] #antra lentele
m4[,2,] #is abieju lenteliu po 2 stuplTitanic
dim(Titanic) #keturmate, 4,2,2,2
ftable(Titanic)
dim(ftable(Titanic))m1
m2
m1+m2 #matricu sudetis
m1-m2 #matricu atimtis
m1*m2 #matricu elementu daugyba
t(m2) #transponuota matrica
m1 %*% m2 #tikroji matricu daugyba!!!!!!
m1 %*% t(m2) solve(m1 %*% t(m2)) #matricos atvirkstine
det(m1 %*% t(m2)) #matricos det
getwd()
setwd("C:/Users/Studentas/Desktop/s1811406")##############Duomenu_importavimas_3########################m1<-matrix(1:6, nrow=2, ncol=3, byrow=T)
m1<-matrix(1:6, nrow=2, ncol=3, byrow=F )
m2<-matrix(1:6, nrow=2, ncol=3, byrow=F)
m2
m5<-rbind(m1,m2) #eilute
m5
m6<-cbind(m1,m2)
m6dat<- data.frame(x=1:3, y=c("a","b","c"), f=factor(c("v","m","v")))
m6dat
class(dat)
str(dat)
summary(dat) #dazniu lenteles y, f; x stulp kvantiliai,medianos
(salis<-c("Prancuzija", "Olandija", "Vokietija", "Anglija", "Ispanija"))
(oro_bendrove<-c("Wizzair", "Lufthansa", "Lufthansa", "Ryanair", "Ryanair"))
(kaina<-c(60,120,80,150,130))
(viesbutis<-c(250,400,200,300,300))
(naktys<-c(3,7,5,4,7))matr<-cbind(salis, oro_bendrove, kaina, viesbutis, naktys) #cbind sukure matrica, character tipas su kabutemis gavosi, nes skc dingo, objektai pasikeite
matrsistema<-data.frame(salis, oro_bendrove, kaina, viesbutis, naktys) #pranasesnis budas nei cbind
sistemamode(sistema) #koks tipas? ats. misrus
class(sistema)
sistema[,2]
class(sistema[,2])
sistema[,4]
class(sistema[,4])sistema[,3]        #kaip rasti duomenis
sistema[,"kaina"] #kaip rasti duomenis
sistema$kaina      #kaip rasti duomenis
attach(sistema) #paima lenteles pavadinimus ir leidzia jais vari.
kaina
oro_bendrove
viesbutis/naktys
kaina+viesbutis
detach(sistema) #atjungti lenteles
#list
srs<-list(c("gerai", "blogai"), c(T, F), c(1,5,2), c(4,6))
srs
srs[[2]]
srs[[3]][[2]]  #reikia dvieju skliausteliu
anketa<-list(vardas=c("Kazys Jonaitis", "Povilas Petraitis"),
pareigos=c("mokytojas", "vairuotojas"),
seima=c("nevedes","vedes"),
vaiku_skaicius=c(2,3),
vaiku_amzius=matrix(c(5,6,NA,1,4,7),
                    ncol=3, byrow=T))
anketamode(anketa)
class(anketa)
attributes(anketa)
dim(anketa)
anketa[[2]]
anketa$pareigos
names(anketa)#############################################################################################################################library(readxl)
nus <- read_excel("Nusikalstamumas.xlsx")
View(nus)nmm<-read.table("Nusikalstamumas.csv", header = T,
                         sep=";")
head(nmm)
summary(nmm)nmm<-read.table("Nusikalstamumas.csv", header = TRUE, sep=";", dec="," )
head(nmm)
summary(nmm)onmm<-read.csv2("Nusikalstamumas.csv", header = TRUE)
head(onmm)
summary(onmm)

install.packages("dplyr")
library(dplyr)data()
head(mtcars) #glava, matomos tik pirmos kelios eilutes
?mtcars
select(mtcars, am:wt)
select(mtcars, am:mpg)
select(mtcars, ends_with("t"))
select(mtcars, contains("m"))mtcars %>% select(am:wt) #antras variantas patogesnisfilter(mtcars, cyl %in% c(4,6)) #ar cilindu stulpelyje yra 4 arba 6, iesko elemento
filter(mtcars, cyl==4 | cyl==6) #same kaip 245 cilindras4 arba cilindas6
filter(mtcars, mpg>=15, mpg<=22) # , & reiskia IR zodeli
filter(mtcars, mpg>=15 & mpg<=22) #tas pats kas pries tai
filter(mtcars, hp>100)mtcars %>% filter(cyl %in% c(4,6))
arrange(mtcars, hp ,mpg) #mpg mazejimo tvarka, hp didejimo tvarka
arrange(mtcars, desc(hp)) #mazejimo hp, nuo didziausio iki maziausiorename(mtcars, svoris=wt) #naujas stulpelis vadinsis svoris, reiksmes imamos is wt
rename(mtcars, svoris=wt, aj=hp)
mtcars2<-mutate(mtcars,
                kml=mpg/3.7854*1.6093)
arrange(mtcars2, desc(kml))grpby <-group_by(mtcars,am)
summarise(grpby, vid= mean(mpg)) #automatas ryja daugiau degaluselect(mtcars, hp) #paima tik hp stulpeli
select(mtcars, -hp) #nepaima tik to stulpeliofilter(mtcars, am==0)
filter(mtcars, am==0 & wt>4)
filter(mtcars, !cyl %in% c(4,6)) #ismesti 4 ir 6, !=ne#== loginis lygu
#!= loginis nelyguarrange(mtcars, wt)
arrange(mtcars, desc(cyl))
mutate(mtcars, wtkg=wt/2.2046)
mtcarsgrp<-group_by(mtcars, am)
summarise(mtcarsgrp, nmb=n(),
          vidsvoris=mean(wt))vidurkiai<-mtcars %>%
        group_by(am) %>%
        summarise_all(list( ~ mean(.)))
vidurkiai #vidurkis vidutine reiksme/ komunizmasmedianos<-mtcars %>%
group_by(am) %>%
        summarise_all(list( ~ median(.)))
medianos #vidurine reiksme#jei mean 2k didesnis uz mediana tai turim isskirciu daug, med=mean, kai yra simetrijasample_n(mtcars, 5) #paimi 5 eilutes, pasirenki atsitiktinaisample_frac(mtcars, 0.1) #nurodo kiek % paimti, atsitiktinaidistinct(mtcars)distinct(mtcars, hp, .keep_all=T)
distinct(mtcars, hp, vs, .keep_all=T) #ismeta pasikartojancias reiksmesmtcars3<-data.frame(rownames(mtcars),
mtcars)
colnames(mtcars3)<-c("Car", colnames(mtcars))
rownames(mtcars3)<-NULL
filter(mtcars3, grepl("er", Car)) #grepl suranda atitikimus 
summarise_at(mtcars, vars(mpg, wt),
             funs(n(),mean, median)) #senas variantast<-mtcars %>%
filter(cyl %in% c(4,8)) %>%
        group_by(am) %>%
        do(head(. , 2)) #do reiskia ivykdyti, taskas reiskia kazkokia manipuliacija
tp <-mtcars %>% select(mpg, cyl, wt, am) %>%
        filter(cyl %in% c(4,8)) %>%
        group_by(am) %>%
        do(arrange(.,desc(wt))) %>% slice(3)
p
#####################################################################################2019-10-01##################################
setwd("C:/Users/Airida/Desktop")
library(ggplot2)
library(dplyr)
dat<-read.csv("nations.csv", header=T)
head(dat)
dt<-select(dat,country, year, gdp_percap, birth_rate, population,region)
dt1<-filter(dt, year=="2013")
head(dt1)
p<-ggplot(data=dt1, mapping=aes(x=gdp_percap, y=birth_rate))+
        geom_point()
p #atvirkstine priklausomybe, BVP rodo salies issivysymo lygi
p<-ggplot(data=dt1, mapping=
                  aes(x=gdp_percap, y=birth_rate,
                      size=population, color=region))
p+geom_point() #salies populiacija tasku dydziu, regionai atskira spalvos
p<-ggplot(data=dt1, mapping=
aes(x=gdp_percap, y=birth_rate))
p+geom_point(color="green", size=2, shape=20, fill="white")+ #21 apskritimas, formos 
        geom_smooth(method="loess")+scale_x_log10()+ #norint pamatyti tendencija suglotnini taskus, suranda bendra tendencija, scale_x_log10() x asyje log pagrindu 10
        theme(axis.text.x=element_text(angle=45, hjust=1, size="12", color="brown"))+ #theme kaip atrodys grafikas
        xlab("LOG_BVP")+ylab("Gimstamumas")+ #asiu pavadinimai
        ggtitle("BVP vs Gimstamumas") #pavadinimas#10e+2 (10^2)
#10e-3 (10^-3)
#e pakelia kvadratu
#asys turi buti matomos, kad nebutu apgauta auditorija? nvm
#################################################################################
ggplot(mtcars, aes(x="", y=mpg)) +
        geom_boxplot(state="boxplot",color="red")
ggplot(dt1, aes(x="", y=gdp_percap)) +
        geom_boxplot(stat="boxplot", outlier.color="blue",
                     outlier.shape = 23, outlier_fill="green") #mediana labai zema, nes turtingu ir issivysciusiu saliu nera daug, jos pakenta i isskirtis, LT patektu ten kur mediana, arba siek tiek auksciau
#isskirciu radimas?
is_outlier<-function(x) {
        return(x<quantile(x, 0.25, na.rm=T) - #NA reiksmes panikina
                       1.5*IQR(x, na.rm=T) |
                       x>quantile(x, 0.75, na.rm=T) +
                       1.5*IQR(x, na.rm=T))
}
dt1 %>%
        group_by(region) %>%
        mutate(outlier = ifelse(is_outlier(gdp_percap), country, NA)) %>% #mutate sukuria nauja stupeli, jei nera isskirties NA, o jei yra tai salies pavadinimas
        ggplot(., aes(x=factor(region), y = gdp_percap)) + 
        geom_boxplot() +
        geom_text(aes(label = outlier), na.rm=T, hjust = -0.3) +
        theme(axis.text.x = element_text(angle = 25))
#####################################didziausia mediana N.America ir issibarstymas itin mazas, 3 kvartilis sutampa su mediana, isskirtys suskaiciuotos pagal regionus
ggplot(mtcars, aes(x=as.factor(cyl))) +
geom_bar(colour = "blue", fill="white") +
        xlab("Cilindr7 skaicius") + ylab("Daznis")
ggplot(mtcars, aes(x=mpg)) +
        geom_histogram(binwidth = 5, colour="blue", fill="white") +
        scale_x_continuous(breaks = seq(10,40,5)) #kai daug skiritngu reiksmiu breziama grupuotu duomenu dazniu lentele
dt1 %>% filter(region == "Europe & Central Asia") %>% #pasiimi dt1 duomenis, isfilturoji duomenis
        mutate(cuts = cut(gdp_percap, breaks = seq(0,155000, by=10000))) %>% #sukuri nauja kintamaji, cut f-ja sugrupuoja duomenis 0-155000 kas 10k
        ggplot(.,aes(x=as.factor(cuts))) + #
        geom_histogram(binwidth = 5, colour="blue", fill="yellow",
                       stat="count") + #statistika count, suskaiciuoja kiek kiekvienam intervale pateko skaiciuku
        theme(axis.text.x = element_text(angle =45)) #pasukam asis
##########################################################################
ggplot(mtcars, aes(x=hp, y=mpg)) +
        geom_point() +
        geom_line()  #kek myliu su vienu galonu degalu, sklaidos diagramu nerekomenduojama jungti eilutemisggplot(mtcars, aes(x=hp, y=mpg)) + #sis variantas daug patogesnis, nes matoma tendecija
geom_point() +
        geom_smooth(methods="loess")ggplot(data=economics, aes(x=date, y=pop)) +
        geom_line(colour="blue", size=2)

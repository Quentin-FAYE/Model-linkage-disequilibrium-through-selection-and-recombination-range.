# List of functions to simulate single locus bi-allelic population genetics model
# Sylvain Glémin 2023
# sylvain.glemin@univ-rennes.fr



# Notation
# Gametes 1:A0B0 2:A0B1 3:A1B0 4:A1B1
# Genotypes noted G11 = A0A0;B0B0 etc...
# Vector of genotypes geno=c(G11,G12,G13,G14,G22,G23,G24,G33,G34,G44)

# VECTEUR de FREQUENCES INITALES ou valeur seectif associe a chacun
# G11 A0A0;B0B0   1
# G12 A0A0;B1B0   1 + hb*sb     # valeur selective sb associe à la frequence à hb ? ???
# G13 A0A1;B0B0   1 + ha*sa   
# G14 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G22 A0A0;B1B1   1 + sb
# G23 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G24 A0A1;B1B1   1 + ha*sa + sb + ehb
# G33 A1A1;B0B0   1 + sa
# G34 A1A1;B0B1   1 + sa + hb*sb + eah
# G44 A1A1;B1B1   1 + sa + sb + eab


# VECTEUR de FREQUENCES INITALES ou valeur seectif associe a chacun
# G11 A0A0;B0B0   
# G12 A0A0;B1B0   
# G13 A0A1;B0B0      
# G14 A0A1;B0B1   
# G22 A0A0;B1B1   
# G23 A0A1;B0B1   
# G24 A0A1;B1B1   
# G33 A1A1;B0B0   
# G34 A1A1;B0B1   
# G44 A1A1;B1B1   

########################### #
# Life cycle functions ######
########################### #



##################### MUTATION  #########################################
#' @title mutation
#' @description Function that return the genotype frequencies after mutation
#' Multiple mutation are not allowed
#' @param geno: a vector of the three genotype frequencies
#' @param mu11: mutation rate from A0 to A1
#' @param mu12: mutation rate from A1 to A0
#' @param mu21: mutation rate from B0 to B1
#' @param mu22: mutation rate from B1 to B0
#' @return the genotype frequencies after mutation
mutation <- function(geno,mu11,mu12,mu21,mu22){
  # Before mutation
  G11 <- geno[1]
  G12 <- geno[2]
  G13 <- geno[3]
  G14 <- geno[4]
  G22 <- geno[5]
  G23 <- geno[6]
  G24 <- geno[7]
  G33 <- geno[8]
  G34 <- geno[9]
  G44 <- geno[10]
  # After mutation
  g11 <- G11*(1-2*mu11-2*mu21) + G12*mu22 + G13*mu12
  g12 <- G12*(1-2*mu11-mu21-mu22) + G11*2*mu21 + G14*mu12 + G22*2*mu22 + G23*mu12
  g13 <- G13*(1-mu11-mu12-2*mu21) + G11*2*mu11 + G14*mu22 + G23*mu22 + G33*2*mu12
  g14 <- G14*(1-mu11-mu12-mu21-mu22) + G12*mu11 + G13*mu21 + G24*mu22 + G34*mu12
  g22 <- G22*(1-2*mu11-2*mu22) + G12*mu21 + G24*mu12
  g23 <- G23*(1-mu11-mu12-mu21-mu22) + G12*mu11 + G13*mu21 + G24*mu22 + G34*mu12
  g24 <- G24*(1-mu11-mu12-2*mu22)  + G14*mu21 + G22*2*mu11 + G23*mu21 + G44*2*mu12
  g33 <- G33*(1-2*mu12-2*mu21) + G13*mu11 + G34*mu22
  g34 <- G34*(1-2*mu12-mu21-mu22) + G14*mu11 + G23*mu11 + G33*2*mu21 + G44*2*mu22
  g44 <- G44*(1-2*mu12-2*mu22) + G24*mu11 + G34*mu21
  newgeno <- c(g11,g12,g13,g14,g22,g23,g24,g33,g34,g44)
  return(newgeno)
}


###################### REPRODUCTION  #######################################
#' @title reproduction
#' @description Function that return the genotype frequencies after reproduction where selfing is allowed
#' @param geno: vector of genotype frequencies
# VECTEUR de FREQUENCES INITALES ou valeur seectif associe a chacun
# G11 A0A0;B0B0   
# G12 A0A0;B1B0   
# G13 A0A1;B0B0      
# G14 A0A1;B0B1   
# G22 A0A0;B1B1   
# G23 A0A1;B1B0 
# G24 A0A1;B1B1   
# G33 A1A1;B0B0   
# G34 A1A1;B0B1   
# G44 A1A1;B1B1  
#' @param self: selfing rate
#' @param rec: recombination rate
#' @return the genotype frequencies after reproduction

# Equation from Hedrick 1980 (Genetics 94: 791-808)
reproduction <- function(geno,self,rec) {  
  #Before reproduction
  G11 <- geno[1]
  G12 <- geno[2]
  G13 <- geno[3]
  G14 <- geno[4]
  G22 <- geno[5]
  G23 <- geno[6]
  G24 <- geno[7]
  G33 <- geno[8]
  G34 <- geno[9]
  G44 <- geno[10]
  #Les Xi sont les fréquences de gamètes haploïdes après méiose + recombinaison, pas les fréquences alléliques.
  X1 <- G11 + (G12+G13+G14)/2 - rec*(G14-G23)/2 # comprends pas la logique
  X2 <- G22 + (G12+G23+G24)/2 + rec*(G14-G23)/2
  X3 <- G33 + (G13+G23+G34)/2 + rec*(G14-G23)/2
  X4 <- G44 + (G14+G24+G34)/2 - rec*(G14-G23)/2
  # After reproduction # pas de rec sans self ? 
  g11 <- self*(G11+(G12+G13+G23*rec^2+G14*(1-rec)^2)/4) + (1-self)*X1^2
  g22 <- self*(G22+(G12+G24+G14*rec^2+G23*(1-rec)^2)/4) + (1-self)*X2^2
  g33 <- self*(G33+(G13+G34+G14*rec^2+G23*(1-rec)^2)/4) + (1-self)*X3^2
  g44 <- self*(G44+(G24+G34+G23*rec^2+G14*(1-rec)^2)/4) + (1-self)*X4^2
  g12 <- self*(G12+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X1*X2
  g13 <- self*(G13+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X1*X3
  g24 <- self*(G24+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X2*X4
  g34 <- self*(G34+rec*(1-rec)*(G14+G23))/2 + (1-self)*2*X3*X4
  g14 <- self*(G14*(1-rec)^2+G23*rec^2)/2 + (1-self)*2*X1*X4
  g23 <- self*(G23*(1-rec)^2+G14*rec^2)/2 + (1-self)*2*X2*X3
  newgeno <- c(g11,g12,g13,g14,g22,g23,g24,g33,g34,g44)
  return(newgeno)
  }



##################### SELECTION   #########################################
#' @title selection
#' @description Function that return the genotype frequencies after selection
#' @param geno: vector of the three genotype frequencies
#' @param fit: vector of fitness
#' @return the genotype frequencies after selection

selection <- function(geno,fit) {
  geno.sel <- geno * fit
  return(geno.sel/sum(geno.sel))
}

############ Comprends pas cette fucntion ? donne le vecteur fit ? matrix 
# Function to set fitness matrix from selecion,dominance and epistasis coefficients#
# The fitness matrix can also be set manually #
fitness <- function(ha,sa,hb,sb,ehh,eah,ehb,eab){
  c(1,
    1 + hb*sb,
    1 + ha*sa,   
    1 + ha*sa + hb*sb + ehh,
    1 + sb,
    1 + ha*sa + hb*sb + ehh,
    1 + ha*sa + sb + ehb,
    1 + sa,
    1 + sa + hb*sb + eah,
    1 + sa + sb + eab  
  )
}



##################### DRIFT   #########################################
#' @title drift
#' @description Function that return the genotype frequencies after drift
#' We assume that Ne = N.
#' @param geno: vector of the three genotype frequencies #10 ? 
#' @param Npop: poopulation size
#' @return the genotype frequencies after drift
drift <- function(geno,Npop){
  newgeno <- rmultinom(n=1,size=Npop,prob = geno)/Npop
  return(newgeno)
}




############################### #
# Examples of simulations ######
############################# #


######################## ONE GENERATION ##################
#' @title generation
#' @description Function that return the genotype frequencies after one generation
#' @param geno: vector of the three genotype frequencies
#' @param mu11: mutation rate from A0 to A1
#' @param mu12: mutation rate from A1 to A0
#' @param mu21: mutation rate from B0 to B1
#' @param mu22: mutation rate from B1 to B0
#' @param self: selfing rate
#' @param rec: recombination rate
#' @param fit: vector of fitness
#' @param Npop: population size
#' @return the genotype frequencies after one generation
generation <- function(geno,Npop,self,rec,mu11,mu12,mu21,mu22,fit){
  g.mut <- mutation(geno,mu11,mu12,mu21,mu22)
  g.repro <- reproduction(g.mut,self,rec)
  g.sel <- selection(g.repro,fit)
  newgeno <- drift(g.sel,Npop)
  return(newgeno)
}



######################## ONE SIMULATION ##################
# Exemple 1: simulation of one population through time

#' @title simulation
#' @description Function that return vectors of allelic and genotypic frequencies thought time
#' @param geno: vector of the three genotype frequencies
#' @param mu11: mutation rate from A0 to A1
#' @param mu12: mutation rate from A1 to A0
#' @param mu21: mutation rate from B0 to B1
#' @param mu22: mutation rate from B1 to B0
#' @param self: selfing rate
#' @param rec: recombination rate
#' @param fit: vector of fitness
#' @param Npop: population size
#' @param Tmax: the total number of generations to simulate
#' @return list of: frequencies of allele a th
simulation <- function(geno0,Npop,self,rec,mu11,mu12,mu21,mu22,fit,Tmax){
  geno <- geno0
  tab.a <- c()
  tab.b <- c()
  tab.LD <- c()
  for(t in 1:Tmax){
    newgeno <- generation(geno,Npop,self,rec,mu11,mu12,mu21,mu22,fit)
    # Genotype frequencies
    G11 <- newgeno[1]
    G12 <- newgeno[2]
    G13 <- newgeno[3]                # comprends pas ou est l'influence du Fit
    G14 <- newgeno[4]                # no use of fitness ? 
    G22 <- newgeno[5]
    G23 <- newgeno[6]
    G24 <- newgeno[7]
    G33 <- newgeno[8]
    G34 <- newgeno[9]
    G44 <- newgeno[10]
    # Haplotype frequencies
    X1 <- G11 + (G12+G13+G14)/2 - rec*(G14-G23)/2 #a_b
    X2 <- G22 + (G12+G23+G24)/2 + rec*(G14-G23)/2 #A_b
    X3 <- G33 + (G13+G23+G34)/2 + rec*(G14-G23)/2 #a_B
    X4 <- G44 + (G14+G24+G34)/2 - rec*(G14-G23)/2 #A_B
    # Allele frequencies
    fa <- X3 + X4
    fb <- X2 + X4
    LD <- X1*X4 - X2*X3
    tab.a <- c(tab.a,fa)
    tab.b <- c(tab.b,fb)
    tab.LD <- c(tab.LD,LD)
    geno <- newgeno
  }
  return(list("a"=tab.a,"b"=tab.b,"LD"=tab.LD))
}


# Vector of genotypes geno=c(G11,G12,G13,G14,G22,G23,G24,G33,G34,G44)

# G11 A0A0;B0B0   1
# G12 A0A0;B1B0   1 + hb*sb
# G13 A0A1;B0B0   1 + ha*sa   
# G14 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G22 A0A0;B1B1   1 + sb
# G23 A0A1;B0B1   1 + ha*sa + hb*sb + ehh
# G24 A0A1;B1B1   1 + ha*sa + sb + ehb
# G33 A1A1;B0B0   1 + sa
# G34 A1A1;B0B1   1 + sa + hb*sb + eah
# G44 A1A1;B1B1   1 + sa + sb + eab

# Example
# set.seed(123)

# Graphique pour taux de recombinaison non nul


freq <- function(nb,rec) {
  s = c(1:nb)/nb
  final = numeric(nb)
  for (i in (1:nb)) {
    results =logical(20)
    for (repetition in (1:20)) {
      fit.matrix <- fitness(0.5,s[i],0.5,0,0,0,0,0) # Two loci with additive selection
      sim <- simulation(geno0 = c(0.25*0.995,0,0.5*0.995,1*0.005,0,0,0,0.25*0.995,0,0),Npop = 100000,self = 0,rec=rec,mu11 = 0,mu12 = 0,mu21 = 0,mu22=0,fit = fit.matrix,Tmax = 300)
      results[repetition] <- (sim$a[299] >= 0.95)
      }
    final[i] <- mean(results)
  }
  return(final)
}

nb = 40
rec = 0.5

rapport = s/rec  # s not define ? 
final = freq(nb,rec)
plot(x = rapport, y = final)




fit.matrix <- fitness(0.5,0.1,0.5,0,0,0,0,0) # Two loci with additive selection
sim <- simulation(geno0 = c(0.25*0.995,0,0.5*0.995,0.5*0.005,0,0,0,0.25*0.995,0.25*0.005,0),Npop = 100000,self = 0,rec=0,mu11 = 0,mu12 = 0,mu21 = 0,mu22=0,fit = fit.matrix,Tmax = m)
plot(NULL,xlim = c(0,200),ylim=c(-0.1,1))
lines(sim$a,col="lightblue")
lines(sim$b,col="orange")
lines(sim$LD,col="black")
abline(0,0,lty=2)



######### SImulation generation throut time #########################

# HOw it works make a loop with thes 4 functions 
#g.mut <- mutation(geno,mu11,mu12,mu21,mu22)
# g.repro <- reproduction(g.mut,self,rec)
# g.sel <- selection(g.repro,fit)
# newgeno <- drift(g.sel,Npop)

#### crrer une function d' une generation  ####---------------------


generation2= function(g.mut,   #g.mut frequences des allèes des toutes les combinaison posible A a B b 
                      self,   #  selfing rate
                      rec,  # recombination rate
                      fit ){   # vecteur associe a chaque genotype -fitennes function
                     # Npop # pop TOalt pour le drif)
                            
                       # pas de process de mutation ici 
                                g.repro=reproduction(g.mut,self, rec ) # process de reproduction
                                g.sel <- selection(g.repro,fit) # proces de selection
                                #newgeno <- drift(g.sel,Npop) # process drift PAS OBLIGATOIRE dans notre question
                                return(g.sel)
}


#### HYPERPARAMETRISATION DE LA FOCNTION DE GENERATION ####
## FRequences intiale ##
# g.mut_v1 : FREQUENCES DES GENOTYPE INITIAUX DANS NOTRE CAS ETUDES
# on part d 'une pop a l'equilibre d'Hedeinvin berge HW avec p=0.5
# on a deux alles mais dans notre simulation on part d' une pop avec
#pas de diversité sur b, on a donc A en p proportion avant equilibre 
# a en q : (1-p) 
# HW pr un allefe donne -> ( aa;XX 0.25 |  Aa;XX: 0.5 | AA;XX :0.25 ) * 0.99 car toute sont (XX bb)
# sauf un ou un faible % qui est * 0.01 XX;Bb) apparition de l'alle 

#BENEFIQUE sur 1 des 2 GENotype  HOMOZYGOTE car si Alle Benefeique B est associé 

# de maniere equilibre vis à vis du locus A/a alors il y a pas 
# de LINKAGE 


#  allel A et a dans la pop peut importe PP selon heindivein bergs
p=0.5
q=1-p
AA_XX=p^2
Aa_xx=2*p*q
aa_XX=q^2

maj=0.995 # majorit" de la pop ne present pas de diversite sur b
B_app= 1-maj # apparition la la diversité benetique sur b->B


  # freq AA_bb  :  AA_XX*maj OU freq Aa_bb  : Aa_xx*maj  OU  freq aa_bb :   aa_XX*maj 
  # freq AA_Bb : AA_XX*B_app   #  CREATION LINKage  entre A et B
   # tout les autres genotype ne sont pas présent
  
# freq allel ORDIN2 selon la logique de code de Sylvain
g.mut_v1=c (
   aa_XX*maj # G11 A0A0;B0B0   
  , 0# G12 A0A0;B1B0     
  , Aa_xx*maj # G13 A0A1;B0B0    
  , 0# G14 A0A1;B0B1   
  , 0# G22 A0A0;B1B1   
  , 0# G23 A0A1;B0B1   ## G23 A0A1;B1B0
  , 0# G24 A0A1;B1B1   
  ,AA_XX*maj # G33 A1A1;B0B0   
  ,AA_XX*B_app  # G34 A1A1;B0B1   
  , 0# G44 A1A1;B1B1
)



## fitness pour selection ##
fit= fitness(ha=0,sa=0,hb=0.5 # additif pour commencer
             ,sb=0.15 # avantage selectif, comme aucune idee biologique
             
            ,ehh=0,eah=0,ehb=0,eab=0 # Aucune epistasie pour l instant
            )
  

# generation2( g.mut_v1 # freq allel ORDIN2 selon la logique de code de Sylvain
#              , self= 0
#              , rec=0, fit=fit)

##### generations throught time #####

n= 200# nb of generation 
matrix_resultat=matrix(data=0,ncol=length(g.mut_v1), nrow =n )

matrix_resultat[1,]=g.mut_v1
for( i in 2:n)  { 
  matrix_resultat[i,]=generation2( matrix_resultat[i-1,] # freq allel ORDIN2 selon la logique de code de Sylvain
               , self= 0 , rec=0, fit=fit)
   }
df_result= as.data.frame(matrix_resultat)
colnames(df_result)= c("Gaa_bb_G11","Gaa_BbG12","GaA_bbG13","GaA_bBG14",
                       "Gaa_BBG22","GaA_BbG23","GaA-BBG24","GAA_bbG33",
                       "GAA_bBG34","AA_BBG44")
# G23 A0A1;B1B0 
### vilsualtion de la simulation poour parametre fixe ###
### ggg ####
plot(x=1:n, y =df_result[,"Gaa_bb_G11"]+df_result[,"GaA_bbG13"]+df_result[,"GAA_bbG33"], type="l", col="red", ylim=c(0,1))
lines(x=1:n ,y = df_result[,"AA_BBG44"]+ df_result[,"GaA_bBG14"]+ df_result[,"GAA_bBG34"], col="blue")
legend(x = 0 , y=0.6  ,        # Position
       legend = c("XX_bb_freq", "XX_BB_freq"),  # Legend texts
       
       col = c("red", "blue"),           # Line colors
       lwd = 2)                 # Line width


# visulation avec linked desequilibrium dans  -> mesure pa la linked desquilibirum
plot(x=1:n, y =df_result[,"Gaa_bb_G11"]+df_result[,"GaA_bbG13"]+df_result[,"GAA_bbG33"], type="l", col="red", ylim=c(0,1))
lines(x=1:n, y=df_result[,"Gaa_bb_G11"], col="orange")
lines(x=1:n, y=df_result[,"GaA_bbG13"], col="tomato")

legend(x ="topleft"  ,        # Position
       legend = c("XX_bb_freq", "Gaa_bb", "GaA_bb"),  # Legend texts
       cex=0.5 ,
       col = c("red", "orange", "tomato"),           # Line colors
       lwd = 2)                 # Line width
# 
# # compute desequillibre equilibirum  tiral####
# # le double heterozygote il a un forme cis et une form trans mais nous on a ssoit l'un soit l'autre donc pas de problème
# df_result["freA_assob"]=df_result[,"GaA_bbG13"]*0.5+df_result[,"GAA_bbG33"]+df_result[,"GAA_bBG34"]*0.5  + df_result[,"GaA_BbG23"]*0.5
# df_result["frea_assob"]=df_result[,"Gaa_bb_G11"]+df_result[,"Gaa_BbG12"]*0.5 +df_result[,"GaA_bbG13"]*0.5+df_result[,"GaA_bBG14"]*0.5 
# 
# df_result["linkageA_to_b_norm"]=df_result["freA_assob"]/(df_result["freA_assob"]+df_result["frea_assob"])
# df_result["linkagea_to_b_norm"]=df_result["frea_assob"]/(df_result["freA_assob"]+df_result["frea_assob"])
# 
# df_result["freA_assoB"]=df_result[,"GaA_bbG13"]*0.5+df_result[,"AA_BBG44"]+df_result[,"GAA_bBG34"]*0.5  + df_result[,"GaA_bBG14"]*0.5
# df_result["frea_assoB"]=df_result[,"Gaa_BBG22"]+df_result[,"Gaa_BbG12"]*0.5 +df_result[,"GaA-BBG24"]*0.5+df_result[,"GaA_BbG23"]*0.5 
# 
# 
# df_result["linkageA_to_B_norm"]=df_result["freA_assoB"]/(df_result["freA_assoB"]+df_result["frea_assoB"])
# df_result["linkagea_to_B_norm"]=df_result["frea_assoB"]/(df_result["freA_assoB"]+df_result["frea_assoB"])
# 
# # compute desequillibre equilibirum ####
# 
# plot(x=1:n, y =df_result[,"linkageA_to_b_norm"], type="l", col="red", ylim=c(0,1))
# lines(x=1:n ,y = df_result[,"linkagea_to_b_norm"], col="blue")
# legend(x = "topleft" ,        # Position
#        legend = c("linkageA_to_b_norm", "linkagea_to_b_norm"),  # Legend texts
#        
#        col = c("red", "blue"),           # Line colors
#        lwd = 2)                 # Line width

colnames( df_result) =c("G11","G12","G13","G14","G22","G23","G24","G33","G34","G44")
rec=0
df_result$X1 <- df_result$G11 + (df_result$G12+df_result$G13+df_result$G14)/2 - rec*(df_result$G14-df_result$G23)/2
df_result$X2 <- df_result$G22 + (df_result$G12+df_result$G23+df_result$G24)/2 + rec*(df_result$G14-df_result$G23)/2
df_result$X3 <- df_result$G33 + (df_result$G13+df_result$G23+df_result$G34)/2 + rec*(df_result$G14-df_result$G23)/2
df_result$df_result$X4 <- df_result$G44 + (df_result$G14+df_result$G24+df_result$G34)/2 - rec*(df_result$G14-df_result$G23)/2
# Allele frequencies
df_result$fa <- df_result$X3 + df_result$X4
df_result$fb <- df_result$X2 + df_result$X4
df_result$LD <- df_result$X1*df_result$X4 - df_result$X2*df_result$X3

plot(x = 1:n , y=df_result$LD, type = "l", ylim=c(0,0.5), col="red" 
     ,main= " Linkage desequilibirum along generations with 2 locus \n initial apparation of B advanteageous allele"
     ,cex.main=0.8, xlab="generation 0->200", ylab="frequence 0->0.5"
     )
lines(x=1:n , y=df_result$X4 , col="slategray")
abline(v=which(df_result$X4 > 0.5)[1])
legend(x="topright" ,
       legend=c(" linkage desequilibrium", "XX_AB freq in pop", " freq  XX_AB >= 0.5", "=0 recombinaison ")
       , col=c("red","slategray", "black", "black") ,cex=0.8, lty=c(1,1,NA,NA), pch =c(NA, NA, "|","c"))

# X1#a_b
# X2  #A_b
# X3  #a_B
# X4  #A_B

########### generation throught time and s space parameter ###########


#  allel A et a dans la pop peut importe PP selon heindivein bergs
p=0.5 ; q=1-p ; AA_XX=p^2 ; Aa_xx=2*p*q ; aa_XX=q^2
maj=0.995 # majorit" de la pop ne present pas de diversite sur b
B_app= 1-maj # apparition la la diversité benetique sur b->B

# freq AA_bb  :  AA_XX*maj OU freq Aa_bb  : Aa_xx*maj  OU  freq aa_bb :   aa_XX*maj  # freq AA_Bb : AA_XX*B_app   #  CREATION LINKage  entre A et B # tout les autres genotype ne sont pas présent
# freq allel ORDIN2 selon la logique de code de Sylvain
g.mut_v1=c (
  aa_XX*maj # G11 A0A0;B0B0   
  , 0# G12 A0A0;B1B0     
  , Aa_xx*maj # G13 A0A1;B0B0    
  , 0# G14 A0A1;B0B1   
  , 0# G22 A0A0;B1B1   
  , 0# G23 A0A1;B0B1   ## G23 A0A1;B1B0
  , 0# G24 A0A1;B1B1   
  ,AA_XX*maj # G33 A1A1;B0B0   
  ,AA_XX*B_app  # G34 A1A1;B0B1   
  , 0# G44 A1A1;B1B1
)

## fitness pour selection ##
fit= fitness(ha=0,sa=0,hb=0.5 # additif pour commencer
             ,sb=0.15 # avantage selectif, comme aucune idee biologique
             ,ehh=0,eah=0,ehb=0,eab=0 # Aucune epistasie pour l instant
)
s_=seq(from=0.01, to=0.2, by=0.016) #  valeur de s explorer 
s=c(0.0005, 0.005,0.01, s_)
m=400
n= m* length(s)+1 # nb of generatio
matrix_resultat=matrix(data=0,ncol=length(g.mut_v1), nrow =n )
matrix_resultat[1,]=g.mut_v1
for (j in 1:length(s)) { # fait la boucle suivant pr chacune des valeurs de s
  fit= fitness(ha=0,sa=0,hb=0.5 # additif pour commencer
               ,sb=s[j] # avantage selectif, comme aucune idee biologique
               ,ehh=0,eah=0,ehb=0,eab=0 # Aucune epistasie pour l instant
  )
  for( i in 2:m)  { 
    matrix_resultat[i+m*(j-1),]=generation2( matrix_resultat[i-1+m*(j-1),] # freq allel ORDIN2 selon la logique de code de Sylvain
                                     , self= 0 , rec=0, fit=fit)
  if ( i==m) {
   
    matrix_resultat[1+m*(j),]=g.mut_v1
    }
  }
}
df_result= as.data.frame(matrix_resultat)

# compute desequillibre equilibirum 
colnames( df_result) =c("G11","G12","G13","G14","G22","G23","G24","G33","G34","G44")
rec=0
df_result$X1 <- df_result$G11 + (df_result$G12+df_result$G13+df_result$G14)/2 - rec*(df_result$G14-df_result$G23)/2
df_result$X2 <- df_result$G22 + (df_result$G12+df_result$G23+df_result$G24)/2 + rec*(df_result$G14-df_result$G23)/2
df_result$X3 <- df_result$G33 + (df_result$G13+df_result$G23+df_result$G34)/2 + rec*(df_result$G14-df_result$G23)/2
df_result$X4 <- df_result$G44 + (df_result$G14+df_result$G24+df_result$G34)/2 - rec*(df_result$G14-df_result$G23)/2
# Allele frequencies
df_result$fa <- df_result$X3 + df_result$X4
df_result$fb <- df_result$X2 + df_result$X4
df_result$LD <- df_result$X1*df_result$X4 - df_result$X2*df_result$X3
#df_result$log_LD=log(df_result$LD) pas bonne idée

nlines <- length(s)
cols <- colorRampPalette(c("blue", "red3"))(nlines)
plot(x = 1:m , y=df_result[1:m,17], type = "l", ylim=c(0,0.2), col="red" 
     ,main= " Linkage desequilibirum along generations with 2 locus \n initial apparation of B advanteageous allele \n with selection strengh variation"
     ,cex.main=0.8, xlab="generation 0->m", ylab="frequence 0->0.20"
     , xlim=c(-50,400))
for (j in 1:length(s)) {
  lines(x=1:m , y=df_result[(1+m*(j-1)):(m+m*(j-1)),17]
  ,col = cols[j],
  lwd = 1.5)
}

legend(x=-50,y=0.2,
       legend = paste0("s = ", round(s, 4)),  col    = cols, lty    = 1, lwd    = 2,
       cex    = 0.55, inset  = 0.02)

########### generation throught time and rec space parameter ###########


#  allel A et a dans la pop peut importe PP selon heindivein bergs
p=0.5 ; q=1-p ; AA_XX=p^2 ; Aa_xx=2*p*q ; aa_XX=q^2
maj=0.995 # majorit" de la pop ne present pas de diversite sur b
B_app= 1-maj # apparition la la diversité benetique sur b->B

# freq AA_bb  :  AA_XX*maj OU freq Aa_bb  : Aa_xx*maj  OU  freq aa_bb :   aa_XX*maj  # freq AA_Bb : AA_XX*B_app   #  CREATION LINKage  entre A et B # tout les autres genotype ne sont pas présent
# freq allel ORDIN2 selon la logique de code de Sylvain
g.mut_v1=c (
  aa_XX*maj # G11 A0A0;B0B0   
  , 0# G12 A0A0;B1B0     
  , Aa_xx*maj # G13 A0A1;B0B0    
  , 0# G14 A0A1;B0B1   
  , 0# G22 A0A0;B1B1   
  , 0# G23 A0A1;B0B1   ## G23 A0A1;B1B0
  , 0# G24 A0A1;B1B1   
  ,AA_XX*maj # G33 A1A1;B0B0   
  ,AA_XX*B_app  # G34 A1A1;B0B1   
  , 0# G44 A1A1;B1B1
)


s_=seq(from=0.016, to=0.2, by=0.016) #  valeur de s explorer 
s=c(0,0.0005, 0.005,s_)
m=400 # generation per simulation
n= m* length(s)+1 # nb of generatio
matrix_resultat=matrix(data=0,ncol=length(g.mut_v1), nrow =n )
matrix_resultat[1,]=g.mut_v1
fit= fitness(ha=0.5,sa=0,hb=0.5 # additif pour commencer
             ,sb=0.074 # avantage selectif, comme aucune idee biologique
             ,ehh=0,eah=0,ehb=0,eab=0 # Aucune epistasie pour l instant
)
for (j in 1:length(s)) { # fait la boucle suivant pr chacune des valeurs de s
  
  for( i in 2:m)  { 
    matrix_resultat[i+m*(j-1),]=generation2( matrix_resultat[i-1+m*(j-1),] # freq allel ORDIN2 selon la logique de code de Sylvain
                                             , self= 0 , rec=s[j], fit=fit)
    if ( i==m) {
      
      matrix_resultat[1+m*(j),]=g.mut_v1
    }
  }
}
df_result= as.data.frame(matrix_resultat)

# compute desequillibre equilibirum 
colnames( df_result) =c("G11","G12","G13","G14","G22","G23","G24","G33","G34","G44")
# initialisation
df_result[, 11:17] <- NA

for (j in 1:length(s)) {
  
  rec_j <- s[j]
  
  start <- 1 + m*(j-1)
  end   <- min(m + m*(j-1), nrow(df_result))
  idx   <- start:end
  
  # X1
  df_result[idx, 11] <-
    df_result[idx, 1] +
    (df_result[idx, 2] + df_result[idx, 3] + df_result[idx, 4]) / 2 -
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # X2
  df_result[idx, 12] <-
    df_result[idx, 5] +
    (df_result[idx, 2] + df_result[idx, 6] + df_result[idx, 7]) / 2 +
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # X3
  df_result[idx, 13] <-
    df_result[idx, 8] +
    (df_result[idx, 3] + df_result[idx, 6] + df_result[idx, 9]) / 2 +
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # X4
  df_result[idx, 14] <-
    df_result[idx,10] +
    (df_result[idx, 4] + df_result[idx, 7] + df_result[idx, 9]) / 2 -
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # fréquences alléliques
  df_result[idx, 15] <- df_result[idx, 13] + df_result[idx, 14]  # fa
  df_result[idx, 16] <- df_result[idx, 12] + df_result[idx, 14]  # fb
  
  # linkage disequilibrium
  df_result[idx, 17] <-
    df_result[idx, 11] * df_result[idx, 14] -
    df_result[idx, 12] * df_result[idx, 13]
}


nlines <- length(s)
cols <- colorRampPalette(c("green", "red"))(nlines)
plot(x = 1:m , y=df_result[1:m,17], type = "l", ylim=c(0,0.01), col="red" 
     ,main= " Linkage desequilibirum along generations with 2 locus \n initial apparation of B advanteageous allele \n with recombinaison variation"
     ,cex.main=0.8, xlab="generation 0->m", ylab="frequence 0->0.004"
     , xlim=c(-50,400))
for (j in 1:length(s)) {
  lines(x=1:m , y=df_result[(1+m*(j-1)):(m+m*(j-1)),17]
        ,col = cols[j],
        lwd = 1.5)
}

legend(x=-50,y=0.009,
       legend = paste0("rec = ", round(s, 4)),  col    = cols, lty    = 1, lwd    = 2,
       cex    = 0.6, inset  = 0.02)
#----------other visualtion 
nlines <- length(s)
cols <- colorRampPalette(c("purple", "lightblue"))(nlines)
plot(x = 1:m , y=df_result[1:m,14], type = "l", col="red" 
     ,main= "A_B frequences in the pop, B being advantageous \n with recombinaison variation"
     ,cex.main=0.8, xlab="generation 0->m", ylab="frequence 0->0.004"
     , xlim=c(-50,400))
for (j in 1:length(s)) {
  lines(x=1:m , y=df_result[(1+m*(j-1)):(m+m*(j-1)),14]
        ,col = cols[j],
        lwd = 1.5)
}

legend(x=-50,y=0.8,
       legend = paste0("rec = ", round(s, 4)),  col    = cols, lty    = 1, lwd    = 2,
       cex    = 0.55, inset  = 0.02)

# notes ########### ##
#Un individu A0A1 ; B0B1 peut exister sous deux organisations chromosomiques différentes :
# pk G13 et G23 , le cAS des doubles heterozygotes 

# Interprétation de h
# h	Signification biologique
# h = 0	B1 récessif
# h = 0.5	additif
# h = 1	B1 dominant
# h > 1	surdominance
# h < 0	sous-dominance

# interpreation epistasie 
# Paramètre	Question biologique
# ehh	Les deux hétérozygotes interagissent-ils ?
#   eah	A1A1 modifie-t-il l’effet de B0B1 ?
#   ehb	B1B1 modifie-t-il l’effet de A0A1 ?
#   eab	Les deux homozygotes mutants interagissent-ils ? 


# le double heterozygote il a un forme cis et une form trans mais nous on a ssoit l'un soit l'autre donc pas de problème

### modelisation avec B dominant et varation de la recombinaisons ####


#  allel A et a dans la pop peut importe PP selon heindivein bergs
p=0.5 ; q=1-p ; AA_XX=p^2 ; Aa_xx=2*p*q ; aa_XX=q^2
maj=0.995 # majorit" de la pop ne present pas de diversite sur b
B_app= 1-maj # apparition la la diversité benetique sur b->B

# freq AA_bb  :  AA_XX*maj OU freq Aa_bb  : Aa_xx*maj  OU  freq aa_bb :   aa_XX*maj  # freq AA_Bb : AA_XX*B_app   #  CREATION LINKage  entre A et B # tout les autres genotype ne sont pas présent
# freq allel ORDIN2 selon la logique de code de Sylvain
g.mut_v1=c (
  aa_XX*maj # G11 A0A0;B0B0   
  , 0# G12 A0A0;B1B0     
  , Aa_xx*maj # G13 A0A1;B0B0    
  , 0# G14 A0A1;B0B1   
  , 0# G22 A0A0;B1B1   
  , 0# G23 A0A1;B0B1   ## G23 A0A1;B1B0
  , 0# G24 A0A1;B1B1   
  ,AA_XX*maj # G33 A1A1;B0B0   
  ,AA_XX*B_app  # G34 A1A1;B0B1   
  , 0# G44 A1A1;B1B1
)


s_=seq(from=0.016, to=0.2, by=0.016) #  valeur de s explorer 
s=c(0,0.0005, 0.005,s_)
m=400 # generation per simulation
n= m* length(s)+1 # nb of generatio
matrix_resultat=matrix(data=0,ncol=length(g.mut_v1), nrow =n )
matrix_resultat[1,]=g.mut_v1
fit= fitness(ha=0,sa=0,hb=1 # additif pour commencer
             ,sb=0.074 # avantage selectif, comme aucune idee biologique
             ,ehh=0,eah=0,ehb=0,eab=0 # Aucune epistasie pour l instant
)
for (j in 1:length(s)) { # fait la boucle suivant pr chacune des valeurs de s
  
  for( i in 2:m)  { 
    matrix_resultat[i+m*(j-1),]=generation2( matrix_resultat[i-1+m*(j-1),] # freq allel ORDIN2 selon la logique de code de Sylvain
                                             , self= 0 , rec=s[j], fit=fit)
    if ( i==m) {
      
      matrix_resultat[1+m*(j),]=g.mut_v1
    }
  }
}
df_result= as.data.frame(matrix_resultat)

# compute desequillibre equilibirum 
colnames( df_result) =c("G11","G12","G13","G14","G22","G23","G24","G33","G34","G44")
# initialisation
df_result[, 11:17] <- NA

for (j in 1:length(s)) {
  
  rec_j <- s[j]
  
  start <- 1 + m*(j-1)
  end   <- min(m + m*(j-1), nrow(df_result))
  idx   <- start:end
  
  # X1
  df_result[idx, 11] <-
    df_result[idx, 1] +
    (df_result[idx, 2] + df_result[idx, 3] + df_result[idx, 4]) / 2 -
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # X2
  df_result[idx, 12] <-
    df_result[idx, 5] +
    (df_result[idx, 2] + df_result[idx, 6] + df_result[idx, 7]) / 2 +
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # X3
  df_result[idx, 13] <-
    df_result[idx, 8] +
    (df_result[idx, 3] + df_result[idx, 6] + df_result[idx, 9]) / 2 +
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # X4
  df_result[idx, 14] <-
    df_result[idx,10] +
    (df_result[idx, 4] + df_result[idx, 7] + df_result[idx, 9]) / 2 -
    rec_j * (df_result[idx, 4] - df_result[idx, 6]) / 2
  
  # fréquences alléliques
  df_result[idx, 15] <- df_result[idx, 13] + df_result[idx, 14]  # fa
  df_result[idx, 16] <- df_result[idx, 12] + df_result[idx, 14]  # fb
  
  # linkage disequilibrium
  df_result[idx, 17] <-
    df_result[idx, 11] * df_result[idx, 14] -
    df_result[idx, 12] * df_result[idx, 13]
}


nlines <- length(s)
cols <- colorRampPalette(c("green", "red"))(nlines)
plot(x = 1:m , y=df_result[1:m,17], type = "l", ylim=c(0,0.01), col="red" 
     ,main= " Linkage desequilibirum along generations with 2 locus \n initial apparation of B advanteageous allele \n with recombinaison variation"
     ,cex.main=0.8, xlab="generation 0->m", ylab="frequence 0->0.004"
     , xlim=c(-50,400))
for (j in 1:length(s)) {
  lines(x=1:m , y=df_result[(1+m*(j-1)):(m+m*(j-1)),17]
        ,col = cols[j],
        lwd = 1.5)
}

legend(x=-50,y=0.009,
       legend = paste0("rec = ", round(s, 4)),  col    = cols, lty    = 1, lwd    = 2,
       cex    = 0.6, inset  = 0.02)
#----------other visualtion 
nlines <- length(s)
cols <- colorRampPalette(c("purple", "lightblue"))(nlines)
plot(x = 1:m , y=df_result[1:m,14], type = "l", col="red" 
     ,main= "A_B frequences in the pop, B being advantageous \n with recombinaison variation"
     ,cex.main=0.8, xlab="generation 0->m", ylab="frequence 0->0.004"
     , xlim=c(-50,400))
for (j in 1:length(s)) {
  lines(x=1:m , y=df_result[(1+m*(j-1)):(m+m*(j-1)),14]
        ,col = cols[j],
        lwd = 1.5)
}

legend(x=-50,y=0.8,
       legend = paste0("rec = ", round(s, 4)),  col    = cols, lty    = 1, lwd    = 2,
       cex    = 0.55, inset  = 0.02)



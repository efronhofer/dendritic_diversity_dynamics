#####################################################################################
# Dendritic network structure and dispersal affect temporal dynamics of diversity and 
# species persistence 

# Mathew Seymour, Emanuel A. Fronhofer, Florian Altermatt

####################################################################################

# Copyright (C) 2015  Emanuel A. Fronhofer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

####################################################################################
# simulation model for scenario "competition colonization trade-off and variation in 
# carrying capacity"

####################################################################################
# GENERAL PARAMETERS
# number of patches (this value should remain fixed, as the connectivity matrices 
# only exist of n_patch =15)
n_patches <- 15
# number species (if this value is changed the interaction matrix has to be adapted)
n_species <- 15
# number of iterations
n_days <- 50
# number of replicates
n_replicates <- 25

# competition colonization trade-off (zero: no trade-off; the larger the stronger)
tau <- 4

####################################################################################
# COMPETITION COLONIZATION TRADE-OFF
competition_coef <- function(d_focal, d_partner){
  return(exp(-tau*(d_focal-d_partner)))
}

####################################################################################
# SPECIES SPECIFIC PARAMETERS
# carrying capacity of each species when alone
K0 <- sample(c(10,100,1000,10000),n_species,replace=T)
# growth rate
lambda0 <- rep(2,n_species)

###################################################################################
# LOGISITIC GROWTH FUNCTION
# inter-specific competition follows the model described in Kubisch et al. 2014
new_pop_size <- function(old_pop_sizes,spec,K,l,as) {
  # first calculate effect of interspecific competition and antagonistic interactions
  K_real <- K * old_pop_sizes[spec]/K0[spec] / (old_pop_sizes[spec]/K0[spec] + sum(as[-spec]*old_pop_sizes[-spec]/K0[-spec]))
  # calculate new population size now following the Beverton-Holt model
  n <- sum(rpois(old_pop_sizes[spec], l/(1+((l-1)/K_real)*old_pop_sizes[spec])))
  if(is.na(n)){n <- 0}
  # return value
  return(n)
}

#####################################################################
# CONNECTIVITY MATRICES
# dendritic [to,from]
connect_dendritic <- matrix(nrow=n_patches, ncol=n_patches)
connect_dendritic[1,] <- c(0,1,0,0,0,0,0,0,1,0,0,0,0,0,0)
connect_dendritic[2,] <- c(1,0,1,0,0,1,0,0,0,0,0,0,0,0,0)
connect_dendritic[3,] <- c(0,1,0,1,1,0,0,0,0,0,0,0,0,0,0)
connect_dendritic[4,] <- c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
connect_dendritic[5,] <- c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
connect_dendritic[6,] <- c(0,1,0,0,0,0,1,1,0,0,0,0,0,0,0)
connect_dendritic[7,] <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
connect_dendritic[8,] <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0)
connect_dendritic[9,] <- c(1,0,0,0,0,0,0,0,0,1,0,0,1,0,0)
connect_dendritic[10,] <- c(0,0,0,0,0,0,0,0,1,0,1,1,0,0,0)
connect_dendritic[11,] <- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
connect_dendritic[12,] <- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0)
connect_dendritic[13,] <- c(0,0,0,0,0,0,0,0,1,0,0,0,0,1,1)
connect_dendritic[14,] <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)
connect_dendritic[15,] <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0)

# linear [to,from]
connect_linear <- matrix(nrow=n_patches, ncol=n_patches)
connect_linear[1,] <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)
connect_linear[2,] <- c(1,0,1,0,0,0,0,0,0,0,0,0,0,0,0)
connect_linear[3,] <- c(0,1,0,1,0,0,0,0,0,0,0,0,0,0,0)
connect_linear[4,] <- c(0,0,1,0,1,0,0,0,0,0,0,0,0,0,0)
connect_linear[5,] <- c(0,0,0,1,0,1,0,0,0,0,0,0,0,0,0)
connect_linear[6,] <- c(0,0,0,0,1,0,1,0,0,0,0,0,0,0,0)
connect_linear[7,] <- c(0,0,0,0,0,1,0,1,0,0,0,0,0,0,0)
connect_linear[8,] <- c(0,0,0,0,0,0,1,0,1,0,0,0,0,0,0)
connect_linear[9,] <- c(0,0,0,0,0,0,0,1,0,1,0,0,0,0,0)
connect_linear[10,] <- c(0,0,0,0,0,0,0,0,1,0,1,0,0,0,0)
connect_linear[11,] <- c(0,0,0,0,0,0,0,0,0,1,0,1,0,0,0)
connect_linear[12,] <- c(0,0,0,0,0,0,0,0,0,0,1,0,1,0,0)
connect_linear[13,] <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,1,0)
connect_linear[14,] <- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,1)
connect_linear[15,] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0)

##################################################################################
# OUTPUT
# population sizes in the dendritic system
N_dendritic <- array(dim=c(n_replicates, n_patches, n_species, n_days))
# population sizes in the linear system
N_linear <- array(dim=c(n_replicates, n_patches, n_species, n_days))

###################################################################################
# SIMULATION

for (repli in 1:n_replicates){
  
  # INITIALIZATION
  # randomly assign dispersal rates
  d <- runif(n_species,0.05,0.3)
  # strength of inter-specific interactions (0: none; 1: as intra specific; above: inter and larger than intra-specific interactions)
  alpha <- matrix(nrow=n_species, ncol=n_species)
  # determine the interaction matrix alpha[effect on species Y, effect of species X]
  # according to the dispersal rates and the trade-off function
  for (focal in c(1:n_species)){
    for (partner in c(1:n_species)){
      alpha[partner,focal] <- competition_coef(d[focal], d[partner])
    }
  }
  
  # dendritic landscape: initialize with zeros
  N_dendritic[repli,,,] <- 0
  # linear landscape: initialize with zeros
  N_linear[repli,,,] <- 0
  
  # initialize a randomly drawn species in a patch
  # this depends on whether there are as many species as patches or not
  if (n_patches == n_species){
    species_no_in_patch <- sample(1:n_species, replace=F)
  }else{
    species_no_in_patch <- sample(1:n_species, n_patches, replace=T)
  }
  for (patch in 1:n_patches){
    # dendritic landscapes
    N_dendritic[repli,patch,species_no_in_patch[patch],1] <- K0[species_no_in_patch[patch]]
    # linear landscape
    N_linear[repli,patch,species_no_in_patch[patch],1] <- K0[species_no_in_patch[patch]]
  }
  
  ########################################################################
  # ITERATION LOOP
  for (day in c(2:n_days)){
    # new population size is old population size
    N_dendritic[repli,,,day] <- N_dendritic[repli,,,day-1]
    N_linear[repli,,,day] <- N_linear[repli,,,day-1]
    
    # PATCH LOOP
    for (patch in c(1:n_patches)){
      #########################################################################
      # DISPERSAL
      # species loop
      for (species in c(1:n_species)){
        # dendritic landscapes
        # emigration only if individuals are there
        if (N_dendritic[repli,patch,species,day-1] > 0){
          # determine the number of emigrants for this patch
          emigrants <- rbinom(1,N_dendritic[repli,patch,species,day-1],1-(1-d[species])^length(which(connect_dendritic[,patch]==1)))
          # save this information in the present population size
          N_dendritic[repli,patch,species,day] <- N_dendritic[repli,patch,species,day] - emigrants
          # now distribute these emigrants according to the specific connectivity pattern
          N_dendritic[repli,which(connect_dendritic[,patch]==1),species,day] <- N_dendritic[repli,which(connect_dendritic[,patch]==1),species,day] + as.numeric(rmultinom(1,emigrants,prob=rep(1/length(which(connect_dendritic[,patch]==1)),length(which(connect_dendritic[,patch]==1)))))
        }
        
        # linear landscapes
        # emigration only if individuals are there
        if (N_linear[repli,patch,species,day-1] > 0){
          # first determine the number of emigrants for this patch
          emigrants <- rbinom(1,N_linear[repli,patch,species,day-1],1-(1-d[species])^length(which(connect_linear[,patch]==1)))
          # save this information in the present population size
          N_linear[repli,patch,species,day] <- N_linear[repli,patch,species,day] - emigrants
          # now distribute these emigrants according to the specific connectivity pattern
          N_linear[repli,which(connect_linear[,patch]==1),species,day] <- N_linear[repli,which(connect_linear[,patch]==1),species,day] + as.numeric(rmultinom(1,emigrants,prob=rep(1/length(which(connect_linear[,patch]==1)),length(which(connect_linear[,patch]==1)))))
          
        }
      }
    }
    
    ###############################################################################################
    # REPRODUCTION
    # patch loop
    for (patch in c(1:n_patches)){
      
      # help variables
      # these are needed to save the current population sizes in this patch for all species
      # as new population sizes have to be calculated from these values to avoid temporal artefacts
      # from the inter-specific interactions
      old_pop_sizes_dendritic <- N_dendritic[repli,patch,,day]
      old_pop_sizes_linear <- N_linear[repli,patch,,day]
      
      # species loop
      for (species in c(1:n_species)){
        
        # dendritic landscapes
        # reproduction only if individuals are there
        if (N_dendritic[repli,patch,species,day] > 0){
          N_dendritic[repli,patch,species,day] <- new_pop_size(old_pop_sizes_dendritic,species,K0[species],lambda0[species],alpha[species,])
        }
        
        # linear landscapes
        # reproduction only if individuals are there
        if (N_linear[repli,patch,species,day] > 0){
          N_linear[repli,patch,species,day] <- new_pop_size(old_pop_sizes_linear,species,K0[species],lambda0[species],alpha[species,])
        }
      }
    }
  }
}

#############################################################################
# CALCULATE DIVERSITY INDICES
library(simba)

alpha_diversity_dendritic <- array(dim=c(n_replicates, n_days))
beta_diversity_dendritic <- array(dim=c(n_replicates, n_days))
gamma_diversity_dendritic <- array(dim=c(n_replicates, n_days))
alpha_diversity_linear <- array(dim=c(n_replicates, n_days))
beta_diversity_linear <- array(dim=c(n_replicates, n_days))
gamma_diversity_linear <- array(dim=c(n_replicates, n_days))

for (repli in c(1:n_replicates)){
  for (day in c(1:n_days)){
    alpha_diversity_dendritic[repli,day] <- as.numeric(trudi(N_dendritic[repli,,,day])[3])
    beta_diversity_dendritic[repli,day] <- as.numeric(trudi(N_dendritic[repli,,,day])[2])
    gamma_diversity_dendritic[repli,day] <- as.numeric(trudi(N_dendritic[repli,,,day])[1])
    alpha_diversity_linear[repli,day] <- as.numeric(trudi(N_linear[repli,,,day])[3])
    beta_diversity_linear[repli,day] <- as.numeric(trudi(N_linear[repli,,,day])[2])
    gamma_diversity_linear[repli,day] <- as.numeric(trudi(N_linear[repli,,,day])[1])
  }
}

# PLOT RESULTS
x11(width=6,height=8)
par(mar=c(0.5,6.5,0.5,0.5),bty="l",mfrow=c(3,1),oma=c(4.5,0,0,0))

plot(gamma_diversity_dendritic[1,]~c(1:n_days),cex.lab=1.7,ylab=expression(paste(gamma," - diversity")),xlab="",ylim=c(0,n_species),type="n",xaxt="n")
for (i in c(1:n_replicates)){
  points(gamma_diversity_dendritic[i,]~c(1:n_days),type="l",col="red",pch=2)
  points(gamma_diversity_linear[i,]~c(1:n_days),type="l",col="blue",pch=1)
}
# add legend
legend("bottomleft",fill=c("red","blue"),legend=c("dendritic","linear"),bty="n", title="network type",cex=1.7)

# plot mean alpha diversity
plot(alpha_diversity_dendritic[1,]~c(1:n_days),cex.lab=1.7,ylab=expression(paste(bar(alpha)," - diversity")),xlab="time [iterations]",ylim=c(0,n_species),type="n",xaxt="n")
for (i in c(1:n_replicates)){
  points(alpha_diversity_dendritic[i,]~c(1:n_days),type="l",col="red",pch=2)
  points(alpha_diversity_linear[i,]~c(1:n_days),type="l",col="blue",pch=1)
}

# plot beta diversity
plot(beta_diversity_dendritic[1,]~c(1:n_days),cex.lab=1.7,ylab=expression(paste(bar(beta)," - diversity")),xlab="time [iterations]",ylim=c(0,n_species),type="n")
for (i in c(1:n_replicates)){
  points(beta_diversity_dendritic[i,]~c(1:n_days),type="l",col="red",pch=2)
  points(beta_diversity_linear[i,]~c(1:n_days),type="l",col="blue",pch=1)
}
# add x label
mtext(side=1,at=n_days/2,line=3,"time [days]",cex=1.1)

##################################################################################
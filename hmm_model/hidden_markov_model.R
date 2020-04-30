# markov implementation

library("reshape2")
library("ggplot2")

#################
# demonstration #
#################

set.seed(1)

##############################################################################################################################################################
##############################################################################################################################################################

# define the simulation context

# fake data

n_years       <- 20
n_communities <- 20 
n_households  <- rnbinom(n_communities, mu = 50, size = 5)
hist(n_households) # visualize the household numbers per community

rest_waste_per_household  <- 50
vgf_waste_per_household   <- 50

# assume there are 2 unobserved states at the household level
states        <- c(1,2)     # separate waste or not


####################
## SIMULATED DATA ##
####################

#----------------------------#
# initial state probablities #
#----------------------------#

# generate some community specific initial state probabilities (proportion household in state 1 and 2 per community)
init_probs    <- runif(n_communities,0.3,0.8)
init_probs    <- cbind(init_probs, 1 - init_probs)

#-------------------------#
# state transition matrix #
#-------------------------#

# this is similar for all individuals, but covariates can be added to make the transition dynamically
# P_trans      <- cbind(c(0.95,0.05),
#                       c(0.05,0.95))

ilogit <- function(x)exp(x)/(1+exp(x))
getP_trans <- function(covars , pars11, pars22){
  p11 <- ilogit(sum(covars * pars11))
  p12 <- 1 - p11
  p22 <- ilogit(sum(covars * pars22)) 
  p21 <- 1 - p22
  return(rbind(c(p11,p12),c(p21,p22)))
}


# define some fake transition matrix parameters
intercept11                     <- 5 
slope_price_ratio_VGF_ro_REST11 <- 4
slope_collection_frequency11    <- -2
pars11 <- c(intercept11, slope_price_ratio_VGF_ro_REST11, slope_collection_frequency11)

intercept22                     <- -2 
slope_price_ratio_VGF_ro_REST22 <- -4
slope_collection_frequency22    <- 2
pars22 <- c(intercept22, slope_price_ratio_VGF_ro_REST22, slope_collection_frequency22)

# create some fake covariate: price ratio VGF vs REST waste, and VGF waste collection frequency
sd_rw <- 0.025
gft_rest_price <- sapply(rnorm(n_communities, mean = .1, sd = .2),function(x, nyears, sd_rw){
  res <- cumsum(c(x,rep(0, (nyears - 1)))+ c(0, rnorm((nyears-1), 0, sd_rw)))
  res[res < 0] <- res[res < 0] * -1
  return(res)
},nyears = n_years, sd_rw = sd_rw)
matplot(gft_rest_price, main = "price ratio VGF vs REST waste per community over time", xlab = "year")

sd_rw <- 0.3
collection_rate <- sapply(runif(n_communities, min = 0, max = 4),function(x, nyears, sd_rw){
  res <- cumsum(c(x,rep(0, (nyears - 1)))+ c(0, rnorm((nyears-1), 0, sd_rw)))
  res[res < 0]  <- res[res < 0] * -1
  res[res > 4]  <- 4 - (res[res > 4] - 4)
  res[res == 0] <- abs(rnorm(sum(res == 0), 0, sd_rw))
  return(ceiling(res))
},nyears = n_years, sd_rw = sd_rw)
matplot(collection_rate, 
        main = "collection frequency per community per month over time", xlab = "year", type = "l", cex = .5)


df_covars <- data.frame("community" = rep(1:n_communities, n_years),
                        "year" = sort(rep(1:n_years, n_communities)),
                        "intercept" = 1,
                        "price_VGF_REST" = c(t(gft_rest_price)),
                        "collection_rate_month" =  c(t(collection_rate)))
df_covars <- df_covars[df_covars$year > 1,]

# get a list of transition matrixes with each list item corresponding to a community and year
P_trans_list <- lapply(split(df_covars[,c(3:5)],seq(nrow(df_covars[,c(3:5)]))), getP_trans, pars11, pars22)


##############################################################################################################################################################
##############################################################################################################################################################
#------------------#
# the actual model #
#------------------#

# simulate the latent (Hidden) Markov states

# list of matrices (a matrix for each community) with first state specified according community specific initial state probability

# t = 1
community_list <- lapply(seq_len(n_communities),function(x, n_years, n_households, init_probs){
  Z     <- matrix(NA,nrow = n_years, ncol = n_households[x])
  Z[1,] <- sample(states, size = n_households[x], replace = T, prob = init_probs[x,])
  return(Z)
}, n_years, n_households, init_probs)

# t = 1 => n_years - 1
for(t in 1:(n_years-1)){
  community_list <- lapply(seq_len(n_communities), function(x, community_list, t, states, P_trans_list){
    Z <- community_list[[x]]
    P_trans <- P_trans_list[[((t-1)*n_communities + x)]]
    Z[t+1,] <- sapply(1:ncol(Z), function(x, state_vec, states, P_trans){
      return(sample(states,1,prob = P_trans[state_vec[x],]))
      },state_vec = Z[t,],states = states, P_trans = P_trans)
    return(Z)
  }, community_list = community_list, t = t, states = states, P_trans_list = P_trans_list)
}

##############################################################################################################################################################
##############################################################################################################################################################

#--------#
# output #
#--------#

# plot mean state at the community level over time
# matrix with each row corresponding to a community, and the columns to the year
mean_state_by_community <- do.call("rbind",lapply(community_list, function(x){
  return(apply(x,1,mean)-1)
}))

# transform to long format
mean_state_by_community_long <- melt(mean_state_by_community)
colnames(mean_state_by_community_long) <- c("community","year","mean_state")

mean_state_by_community_long <- dplyr::left_join(mean_state_by_community_long, df_covars)

cor(mean_state_by_community_long[mean_state_by_community_long$year > 1,c("mean_state","price_VGF_REST","collection_rate_month")])

ggplot(mean_state_by_community_long, aes(year, mean_state)) +
  geom_line() + 
  facet_wrap(~ community) +
  theme_bw()


ggplot(mean_state_by_community_long, aes(year, scale(mean_state))) +
  geom_line() + 
  geom_line(aes(year, scale(price_VGF_REST)), col = "red") + 
  geom_line(aes(year, scale(collection_rate_month)), col = "green") + 
  facet_wrap(~ community) +
  theme_bw()


###################
## OBSERVED DATA ##
###################
names(n_households) <- 1:length(n_households)
mean_state_by_community_long$household_nr <- n_households[match(as.character(mean_state_by_community_long$community),names(n_households))]

# observed values per community over time
mean_state_by_community_long$observed <- apply(mean_state_by_community_long, 1 , function(x){
  mean(rbinom(x["household_nr"], size = 1, prob =x["mean_state"]))
})

# red lines provides the observed data, and black line the state
ggplot(mean_state_by_community_long, aes(year, mean_state)) +
  geom_line() + 
  geom_line(aes(year, observed), col = "red") +
  facet_wrap(~ community) +
  theme_bw()


mean_state_by_community_long$observed_rest_waste <- mean_state_by_community_long$household_nr * rest_waste_per_household + (1-mean_state_by_community_long$observed) * mean_state_by_community_long$household_nr * vgf_waste_per_household
mean_state_by_community_long$observed_vgf_waste  <- mean_state_by_community_long$observed * mean_state_by_community_long$household_nr * vgf_waste_per_household

mean_state_by_community_long$predicted_rest_waste <- mean_state_by_community_long$household_nr * rest_waste_per_household + (1-mean_state_by_community_long$mean_state) * mean_state_by_community_long$household_nr * vgf_waste_per_household
mean_state_by_community_long$predicted_vgf_waste  <- mean_state_by_community_long$mean_state * mean_state_by_community_long$household_nr * vgf_waste_per_household

# transformed to waste
aggr_data <-aggregate(cbind(predicted_vgf_waste,observed_vgf_waste) ~ year, FUN = sum, data = mean_state_by_community_long)

ggplot(aggr_data, aes(year, predicted_vgf_waste)) +
  geom_line() + 
  geom_line(aes(year, observed_vgf_waste), col = "red") +
  theme_bw()


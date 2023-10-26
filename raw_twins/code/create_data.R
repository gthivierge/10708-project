library(readr)
library(dplyr)
library(MASS)
library(mvnfast)
library(mltools)
library(data.table)

X <- read_csv("twin_pairs_X_3years_samesex.csv")
Y <- read_csv("twin_pairs_Y_3years_samesex.csv")
T_dat <- read_csv("twin_pairs_T_3years_samesex.csv")

# get rid of row labels after import
X <- X[,-c(1)] # why are there 2? idk, leave one as identifier i guess
names(X)[1]<-"ID"
Y <- Y[,-1]
T_dat <- T_dat[,-1]


# The treatment t = 1 is being born the heavier twin whereas, the outcome corresponds to the
# mortality of each of the twins in their first year of life

dat <- bind_cols(X, T_dat, Y)

dat_2kg <- dat %>%
  filter(dbirwt_1 < 2000) %>% # 11984 pairs as stated
dplyr::select(-c(data_year, infant_id_0, infant_id_1, brstate_reg, mplbir_reg, stoccfipb_reg))


# sum(dat_2kg$mort_0)/length(dat_2kg$mort_0) # 18.8 as stated
# sum(dat_2kg$mort_1)/length(dat_2kg$mort_1) # 16.4 as stated


# In order to simulate the
# case of hidden confounding with proxies, we based the treatment assignment on a single variable
# which is highly correlated with the outcome: GESTAT10, the number of gestation weeks prior to
# birth. It is ordinal with values from 0 to 9 indicating birth before 20 weeks gestation, birth after
# 20-27 weeks of gestation and so on 4
# . We then set ti
# |xi
# , zi ∼ Bern
# σ(w
#   >
#     o x + wh(z/10 − 0.1))
# ,
# wo ∼ N (0, 0.1 · I), wh ∼ N (5, 0.1), where z is GESTAT10 and x are the 45 other features


# exclude: data_year, infant_id_0, infant_id_1, 
# brstate_reg, mplbir_reg, stoccfipb_reg: there are already variables for the actual states. 
# this is a guess for what they excluded to end up with 46 covariates

# GESTAT10:
# 01 ... Under 20 weeks
# 02 ... 20 - 27 weeks
# 03 ... 28 - 31 weeks
# 04 ... 32 - 35 weeks
# 05 ... 36 weeks
# 06 ... 37 - 39 weekg
# 07 ... 40 weeks
# 08 ... 41 weekg
# 09 ... 42 weeks and over
# 10 ... Not stated

covars <- dat_2kg %>%
  dplyr::select(-c(ID, mort_0, mort_1, dbirwt_0, dbirwt_1, gestat10)) %>%
  data.frame()

gestat <- dat_2kg$gestat10

sample_tx <- function(w_0, w_h, mat){
  w <- length(mat)
  x <- mat[1:(w-1)]
  z <- mat[w]
  y <- t(w_0) %*% as.numeric(x) + w_h * (z/10 - 0.1) 
  p <- 1/(1 + exp(-y))
  tx <- rbinom(1, 1, p)
  return(tx)
}

d <- ncol(covars)
N <- length(dat_2kg$ID)
Sigma <- 0.1*diag(d)
mu <- rep(0, d)
w_0 <- mvrnorm(1, mu, Sigma)
w_h <- rnorm(1, mean=5, sd=0.1)

dat_complete <- dat_2kg[complete.cases(dat_2kg),]

covars_complete <- dat_complete %>%
  dplyr::select(-c(ID, mort_0, mort_1, dbirwt_0, dbirwt_1)) %>%
  as.matrix()

tx_vec <- apply(X=covars_complete, MARGIN=1, FUN=sample_tx, w_0 = w_0, w_h = w_h)

sum(tx_vec)/length(tx_vec)


out_data <- dat_complete %>%
  as.data.frame() %>%
  bind_cols("Treatment" = tx_vec) %>%
  mutate(gestat10 = as.factor(gestat10))

names(out_data)[which(names(out_data)=="gestat10")]<-"gestat"

out2 <- setDT(out_data)

out_df <- one_hot(out2, cols = "gestat") %>%
  as.data.frame() %>%
  mutate(
    across(
      .cols = contains('gestat_', ignore.case = FALSE),
      .names = '{.col}_2'
    )) %>%
  mutate(
    across(
      .cols = matches('gestat_[1-9]$', ignore.case = FALSE),
      .names = '{.col}_3'
    ))
  
  # df %>% mutate(
  #   across(
  #     .cols = contains('a', ignore.case = FALSE),
  #     .names = '{.col}2'
  #   )

    
    # flipping some of gestation codes with some probability
    # comparison paper used 0.2 so i will do that too
    # there is surely a more elegant way to do this but idgaf rn
    
    flip <- function(x, p){
      s <- rbinom(1, 1, p)
      if(s == 1){
        z <- 1 - x
      } else{
        z <- x
      }
      return(z)
      z <- s * (1-x) + (1 - s) * x
    }
    
    flip(1, 0.2)
    
    sum(flip(out_df$gestat_1, 0.2) == out_df$gestat_1)
    
    gestat_cols <- out_df %>%
      dplyr::select(starts_with("gestat"))
    
    tst <- apply(gestat_cols, c(1,2), flip, p=0.2)

    out_flipped <- as.data.frame(tst)
    
    
    out_flipped2 <- out_df %>%
      dplyr::select(-starts_with("gestat")) %>%
      bind_cols(out_flipped)
    
    # hide data from one twin ("treatment")
    final_out <- out_flipped2 %>%
      rowwise() %>%
      mutate(Y = Treatment * mort_1 + (1 - Treatment) * mort_0,
             dbirwt = Treatment * dbirwt_1 + (1 - Treatment) * dbirwt_0,
             bord = Treatment * bord_1 + (1 - Treatment) * bord_0) 
    
    
    write.csv(final_out, "cleaned_twins_flipped.csv")
    
    
    final_out_2 <- out_df %>%
      rowwise() %>%
      mutate(Y = Treatment * mort_1 + (1 - Treatment) * mort_0,
             dbirwt = Treatment * dbirwt_1 + (1 - Treatment) * dbirwt_0,
             bord = Treatment * bord_1 + (1 - Treatment) * bord_0) 
    
    write.csv(final_out, "cleaned_twins.csv")

# w_o <- mvrnorm(n = N, mu = rep(0,N), Sigma=Sigma)

  # sigma - sigmoid? they say logistic but don't give parameters
  
  

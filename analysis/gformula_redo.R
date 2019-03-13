### Setup
library(data.table);library(boot)
rm(list = ls())

### Paths
data.path <- "data/bone_data.csv"

### Code

## Step 1 - generate person-day data from bone marrow transplant data
# Read in data
in.dt <- fread(data.path)

in.dt[, wait := waitdays / 30.5]
in.dt[, agesq := age**2]
in.dt[, agecurs1 := (age>17.0)*(age-17.0)**3-((age>30.0)*(age-30.0)**3)*(41.4-17.0)/(41.4-30.0)]
in.dt[, agecurs2 := (age>25.4)*(age-25.4)**3-((age>41.4)*(age-41.4)**3)*(41.4-25.4)/(41.4-30.0)]

person.day.data <- function(p.id, data) {
  print(p.id)
  p.dt <- data[id == p.id]
  p.day.dt <- rbindlist(lapply(1:p.dt$t, function(day) {
    yesterday = day-1
    daysq = day**2
    daycu = day**3
    yesterday = day-1
    daysq = day**2
    daycu = day**3
    daycurs1 = ((day>63)*((day-63)/63)**3)+((day>716)*((day-716)/63)**3)*(350.0-63) -((day>350)*((day-350)/63)**3)*(716-63)/(716-350)
    daycurs2 = ((day>168)*((day-168)/63)**3)+((day>716)*((day-716)/63)**3)*(350-168) -((day>350)*((day-350)/63)**3)*(716-168)/(716-350)
    d =   ifelse(day>=p.dt$t, p.dt$d_dea, 0)
    gvhd =   ifelse(day>p.dt$t_gvhd, 1, 0)
    relapse =  ifelse(day>p.dt$t_rel, 1, 0)
    platnorm =  ifelse(day>p.dt$t_pla, 1, 0)
    # lagged variables
    gvhdm1 = ifelse(yesterday>p.dt$t_gvhd, 1, 0)
    relapsem1 =  ifelse(yesterday>p.dt$t_rel, 1, 0)
    platnormm1 =  ifelse(yesterday>p.dt$t_pla, 1, 0)
    censeof = ifelse(day == p.dt$t & d == 0 & day == 1825, 1, 0) 
    censlost= ifelse(day == p.dt$t & d == 0 & day != 1825, 1, 0) 
    add.dt <- data.table(id = p.id, day, yesterday, daysq, daycu, daycurs1, daycurs2, d, gvhd, 
                         relapse, platnorm, gvhdm1, relapsem1, platnormm1, censeof, censlost)
    return(add.dt)
  }))
}

all.p.day.dt <- rbindlist(lapply(in.dt$id, person.day.data, data = in.dt))

all.p.day.dt[, daysnorelapse := cumsum(relapse == 0), by = id] 
all.p.day.dt[, daysnoplatnorm := cumsum(platnorm == 0), by = id]
all.p.day.dt[, daysnogvhd := cumsum(gvhd == 0), by = id]
all.p.day.dt[, daysrelapse := cumsum(relapse == 1), by = id]
all.p.day.dt[, daysplatnorm := cumsum(platnorm == 1), by = id]
all.p.day.dt[, daysgvhd := cumsum(gvhd == 1), by = id]

in.dt[, c("t", "t_rel", "d_dea", "t_gvhd", "d_gvhd", "d_rel", "t_pla", "d_pla") := NULL]
all.p.day.dt[, c("yesterday") := NULL]
dt <- merge(all.p.day.dt, in.dt, by = "id")


## Step 2 - estimate modeling coefficients used to generate probabilities

vectorize.fit <- function(model.fit, dt){
  coefs <-  coef(summary(model.fit))[,1]
  vec <- unlist(lapply(names(dt), function(name) {
    if(name %in% names(coefs)) {
      coefs[name]
    } else {
      0
    }
  }))
  vec <- c(coefs[1], vec)
}

get.dep.var <- function(model.fit) {
  rownames(attr(summary(model.fit)$terms, "factors"))[1]
}

# Add interaction variables: day*gvhd + daysq*gvhd + daycu*gvhd
dt[, day_gvhd := day * gvhd]
dt[, daysq_gvhd := daysq * gvhd]
dt[, daycu_gvhd := daycu * gvhd]


# Model for probability of platnorm=1 at day k
m.platnorm <- glm(platnorm ~ all + cmv + male + age + gvhdm1 + daysgvhd + daysnorelapse + agecurs1 + agecurs2 + wait, data = dt[platnormm1 == 0], 
                  family = binomial(link = "logit"))
v.platnorm <- vectorize.fit(m.platnorm, dt)
pred.vars <- get.dep.var(m.platnorm)

# Model for probablity of relapse = 1 at day k
m.relapse <- glm(relapse ~ all + cmv + male + age + gvhdm1 + daysgvhd  + platnormm1 + daysnoplatnorm + agecurs1 + agecurs2 + day + daysq + wait, data = dt[relapsem1 == 0], 
             family = binomial(link = "logit"))
v.relapse <- vectorize.fit(m.relapse, dt)
pred.vars <- get.dep.var(m.relapse)

# Model for probability of exposure=1 at day k
m.gvhd <- glm(gvhd ~ all + cmv + male + age + platnormm1 + daysnoplatnorm + relapsem1 + daysnorelapse + agecurs1 + agecurs2 + day + daysq + wait, data = dt[gvhdm1 == 0], 
                 family = binomial(link = "logit"))
v.gvhd <- vectorize.fit(m.gvhd, dt)
pred.vars <- c(pred.vars, get.dep.var(m.gvhd))

# Model for probability of censoring=1 at day k
m.censlost <- glm(censlost ~ all + cmv + male + age + daysgvhd + daysnoplatnorm + daysnorelapse + agesq + agecurs2 + day + daysq + daycu + wait, data = dt, 
              family = binomial(link = "logit"))
v.censlost <- vectorize.fit(m.censlost, dt)
pred.vars <- c(pred.vars, get.dep.var(m.censlost))

# Model for probability of outcome=1 at day k
m.d <- glm(d ~ all + cmv + male + age + gvhd + platnorm + daysnoplatnorm + relapse + daysnorelapse + agesq + day + daysq + daycu + wait + day_gvhd + daysq_gvhd + daycu_gvhd, data = dt, 
              family = binomial(link = "logit"))
v.d <- vectorize.fit(m.d, dt)
pred.vars <- c(pred.vars, get.dep.var(m.d))

coef.dt <- data.table(rbind(v.relapse, v.platnorm, v.gvhd, v.censlost, v.d))
colnames(coef.dt)[2:ncol(coef.dt)] <- names(dt)
coef.dt[, variable := c("platnorm", "relapse",  "gvhd", "censlost", "d")]
coef.dt[is.na(coef.dt)] <- 0

## Step 3 - sample with replacement from data
# Predict log.odds, convert to probability, and sample from Bernoulli
gen.draws <- function(var, matrix) {
  var.idx <- which(colnames(M) == var)
  temp.vec <- matrix[, var.idx]
  coef.vec <- as.matrix(coef.dt[variable == var, 1:(ncol(coef.dt) - 1)])[1,]
  temp.matrix <- cbind(rep(1, nrow(matrix)), matrix)
  log.odds <- temp.matrix %*% coef.vec
  prob <- inv.logit(log.odds)
  draws <- rbinom(nrow(prob), 1, prob)
  temp.vec[temp.vec == 0] <- draws[temp.vec == 0]
  return(temp.vec)
}

# Update cumulative values on each iteration
update.cum <- function(var, M) {
  # Update days no platnorm
  M[, paste0("daysno", var)] <- M[, paste0("daysno", var)] + (1 - M[, var])
  # Update platnorm
  M[, paste0("days", var)] <- M[, paste0("days", var)] + M[, var]
  if(var == "gvhd") {
    # Interaction terms
    M[, "day_gvhd"] <- M[, "gvhd"] * M[, var]
    M[, "daysq_gvhd"] <- M[, "gvhd"] * M[, var]
    M[, "daycu_gvhd"] <- M[, "gvhd"] * M[, var]
  }
  return(M)
}

# Update lags on each iteration
lag.list <- c("gvhd", "relapse", "platnorm")
update.lags <- function(var.list, matrix) {
  for(var in var.list) {
    matrix[, paste0(var, "m1")] <- matrix[, var]
  }
  return(matrix)
}

# Update time varying predictors: daysq, daycu
update.time <- function(M) {
  day <- M[,"day"][1] + 1
  M[,"day"] <- M[,"day"] + 1
  M[, "daysq"] <- M[,"day"]**2
  M[, "daycu"] <- M[,"day"]**3
  M[, "daycurs1"] <- ((day>63)*((M[,"day"]-63)/63)**3)+((day>716)*((M[,"day"]-716)/63)**3)*(350.0-63) -((day>350)*((M[,"day"]-350)/63)**3)*(716-63)/(716-350)
  M[, "daycurs2"] <- ((day>168)*((M[,"day"]-168)/63)**3)+((day>716)*((M[,"day"]-716)/63)**3)*(350-168) -((day>350)*((M[,"day"]-350)/63)**3)*(716-168)/(716-350)
  return(M)
}

gen.baseline <- function(n) {
  M.n <- n * nrow(in.dt)
  M.ids <- sample(in.dt$id, M.n, replace = T)
  dt[, c("gvhd", "platnorm", "relapse") := 0]
  dt.matrix <- as.matrix(dt[day == 1])
  M <- cbind(dt.matrix[M.ids,])
  return(M)
}

simulate <- function(matrix, intervene, coef.dt) {
  out.dt <- data.table()
  # Subset variable list for intervention run
  if(intervene) {
    var.list <- setdiff(coef.dt$variable, c("gvhd", "censlost"))
  } else {
    var.list <- coef.dt$variable 
  }
  for(i in 1:(unique(max(dt$day)) - 1)) {
    print(paste0(i, " of ", (max(dt$day) - 1), ": ", nrow(matrix), " alive"))
    for(var in var.list) {
      var.idx <- which(colnames(matrix) == var)
      matrix[, var.idx] <- gen.draws(var, matrix)
    }
    for (var in lag.list) {
      matrix <- update.cum(var, matrix)
    }
    matrix <- update.lags(lag.list, matrix)
    matrix <- update.time(matrix)
    dead.dt <- as.data.table(matrix[matrix[, "d"] == 1 | matrix[, "censlost"] == 1, ,drop = F])
    if(nrow(dead.dt) > 0) {
      out.dt <- rbind(out.dt, dead.dt)
      matrix <- matrix[!(matrix[, "d"] == 1 | matrix[, "censlost"] == 1), ,drop = F]
    }
    if(nrow(matrix) == 0) {
      break
    }
  }
  out.dt <- rbind(out.dt, as.data.table(matrix))
  return(out.dt)
}
  
# Natural course
n.draws <- 1000
M <- gen.baseline(n.draws)
sim.nat <- as.data.table(simulate(M, intervene = F, coef.dt))

# Intervention
M <- gen.baseline(n.draws)
sim.int <- as.data.table(simulate(M, intervene = T, coef.dt))

# Summary stats
nrow(sim.nat[d == 0]) / nrow(sim.nat) * 100
nrow(sim.int[d == 0]) / nrow(sim.int) * 100

## Step 6 - concatentate intervetion data sets and run Cox model

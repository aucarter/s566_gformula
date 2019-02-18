### Setup
library(data.table)

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
    p.day.dt <- rbind(p.day.dt, add.dt)
  }))
}

all.p.day.dt <- rbindlist(lapply(dt$id, person.day.data, data = dt))

all.p.day.dt[, daysnorelapse := cumsum(relapse == 0), by = id] 
all.p.day.dt[, daysnoplatnorm := cumsum(platnorm == 0), by = id]
all.p.day.dt[, daysnogvhd := cumsum(gvhd == 0), by = id]
all.p.day.dt[, daysrelapse := cumsum(relapse == 1), by = id]
all.p.day.dt[, daysplatnorm := cumsum(platnorm == 1), by = id]
all.p.day.dt[, daysgvhd := cumsum(gvhd == 1), by = id]

in.dt[, c("t", "t_rel", "d_dea", "t_gvhd", "d_gvhd", "d_rel", "t_pla", "d_pla") := NULL]

dt <- merge(all.p.day.dt, in.dt, by = "id")


## Step 2 - estimate modeling coefficients used to generate probabilities
# Model for probablity of relapse = 1 at day k
m.relapse <- glm(relapse ~ all + cmv + male + age + gvhdm1 + daysgvhd  + platnorm + daysnoplatnorm + agecurs1 + agecurs2 + day + daysq + wait, data = dt[relapsem1 == 0], 
             family = binomial(link = "logit"))

# Model for probability of platnorm=1 at day k
m.platnorm <- glm(platnorm ~ all + cmv + male + age + gvhdm1 + daysgvhd + daysnorelapse + agecurs1 + agecurs2 + wait, data = dt[platnormm1 == 0], 
       family = binomial(link = "logit"))

# Model for probability of exposure=1 at day k
m.gvhd <- glm(gvhd ~ all + cmv + male + age + platnormm1 + daysnoplatnorm + relapsem1 + daysnorelapse + agecurs1 + agecurs2 + day + daysq + wait, data = dt[gvhdm1 == 0], 
                 family = binomial(link = "logit"))

# Model for probability of censoring=1 at day k
m.censlost <- glm(censlost ~ all + cmv + male + age + daysgvhd + daysnoplatnorm + daysnorelapse + agesq + agecurs2 + day + daysq + daycu + wait, data = dt, 
              family = binomial(link = "logit"))

# Model for probability of outcome=1 at day k
m.d <- glm(d ~ all + cmv + male + age + gvhd + platnorm + daysnoplatnorm + relapse + daysnorelapse + agesq + day + daysq + daycu + wait + day*gvhd + daysq*gvhd + daycu*gvhd, data = dt, 
              family = binomial(link = "logit"))

## Step 3 - sample with replacement from data

## Step 4 and 5 - run Monte Carlo sample for natural course, always and never exposed

## Step 6 - concatentate intervetion data sets and run Cox model

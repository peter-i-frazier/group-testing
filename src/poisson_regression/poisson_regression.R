library(ggplot2)

df <- read.csv('case-hosp-death.csv', header=TRUE)
df$t <- c(1:dim(df)[1])
df$tsquared <- df$t^2
df$CUM_CASE_COUNT <- cumsum(df$CASE_COUNT)


# 1a
m1a <- glm(CUM_CASE_COUNT ~ t + tsquared, family=poisson(link="identity"), start=c(1,1,1), data=df)
summary(m1a)
ypredm1a <- predict(m1a) 

ggplot(df, aes(x = t, y = CUM_CASE_COUNT)) +
  geom_point(color='blue') +
  geom_line(aes(x = t, y = ypredm1a), size = 1, color='green')

# 1b
m1b <- glm(CASE_COUNT ~ t + tsquared, family=poisson(link="identity"), start=c(1,1,1), data=df)
summary(m1b)

ypredm1b <- predict(m1b) 

ggplot(df, aes(x = t, y = CASE_COUNT)) +
  geom_point(color='blue') +
  geom_line(aes(x = t, y = ypredm1b), size = 1, color='green')

# 2a
shift <- function(x, n){
  c(rep(NA, n), x[seq(length(x) - n)])
}

for (k in c(1:7))
{
  colname <- paste0('CASE_MINUS_', k)
  df[colname]<- shift(df$CASE_COUNT, k)
}

m2a <- glm(CASE_COUNT ~ CASE_MINUS_1 + CASE_MINUS_2 + 
             CASE_MINUS_3 + CASE_MINUS_4 + CASE_MINUS_5 + 
             CASE_MINUS_6 + CASE_MINUS_7, family=poisson(link="identity"), 
              start=c(1,1,1,1,1,1,1,1), data=df[8:nrow(df),])
summary(m2a)

ypredm2a <- predict(m2a) 

ggplot(df[8:nrow(df),], aes(x = t, y = CASE_COUNT)) +
  geom_point(color='blue') +
  geom_line(aes(x = t, y = ypredm2a), size = 1, color='green')

m2a5 <- glm(CASE_COUNT ~ CASE_MINUS_1 + CASE_MINUS_2 + 
             CASE_MINUS_3 + CASE_MINUS_4 + CASE_MINUS_5, family=poisson(link="identity"), 
           start=c(1,1,1,1,1,1), data=df[6:nrow(df),])

summary(m2a5)

ypredm2a5 <- predict(m2a5) 

ggplot(df[6:nrow(df),], aes(x = t, y = CASE_COUNT)) +
  geom_point(color='blue') +
  geom_line(aes(x = t, y = ypredm2a5), size = 1, color='green')


m2a3 <- glm(CASE_COUNT ~ CASE_MINUS_1 + CASE_MINUS_2 + 
              CASE_MINUS_3, family=poisson(link="identity"), 
            start=c(1,1,1,1), data=df[4:nrow(df),])

summary(m2a3)

ypredm2a3 <- predict(m2a3) 

ggplot(df[4:nrow(df),], aes(x = t, y = CASE_COUNT)) +
  geom_point(color='blue') +
  geom_line(aes(x = t, y = ypredm2a3), size = 1, color='green')


# 2b


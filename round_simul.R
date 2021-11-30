###simulation with round(beta,3)

###use energy data

lm(Y2 ~ X1+X2+X3,data=energy)
reg2.eng2 <- lm(Y2 ~ X1+,data=energy)$coef[-1]
reg3.eng2 <- lm(Y2 ~ X2,data=energy)$coef[-1]

lm(Y2 ~ X1+X2+X3,data=energy)$coef[-1]
lm(Y2 ~ X1+X2,data=energy)$coef[-1]
lm(Y2 ~ X2+X3,data=energy)$coef[-1]

round.coef1.eng <- round(lm(Y2 ~ X1+X2+X3,data=energy)$coef[-1],1)
round.coef2.eng <- round(lm(Y2 ~ X1+X2,data=energy)$coef[-1],1)
round.coef3.eng <- round(lm(Y2 ~ X2+X3,data=energy)$coef[-1],1)

round.beta.eng <- list(round.coef1.eng,round.coef2.eng,round.coef3.eng)
round.eng.index <- list(c(1,2,3),c(1,2),c(2,3))

round.eng.gensa <- get.gensa.res(round.beta.eng,round.eng.index)
round.eng.gensa$value


###use fish data
head(fish)

round.coef1.fish <- round(lm(Weight ~ Length1 + Length2 + Length3,data=fish)$coef[-1],1)
round.coef2.fish <- round(lm(Weight ~ Length2 + Length3,data=fish)$coef[-1],1)
round.coef3.fish <- round(lm(Weight ~ Length1 + Length3,data=fish)$coef[-1],1)

round.beta.fish <- list(round.coef1.fish,round.coef2.fish,round.coef3.fish)
round.fish.index <- list(c(1,2,3),c(2,3),c(1,3))

round.fish.gensa <- get.gensa.res(round.beta.fish,round.fish.index)
round.fish.gensa$value







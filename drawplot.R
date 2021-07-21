library(ggplot2)
#draw two normal distribution graph


x=seq(0,40,0.01)
normal1 <- data.frame(x, dnorm(x,mean=10,sd=5))

normal2 <- data.frame(x, dnorm(x,mean=25,sd=5))
cols = c("x","y")
colnames(normal1) = cols
colnames(normal2) = cols

p = ggplot() + 
  geom_line(data = normal1, aes(x = x, y = y), color = "blue") +
  geom_line(data = normal2, aes(x = x, y = y), color = "red") +
  xlab('x') +
  ylab('y')

print(p)


s <- seq(0,15,0.01)
plot(s, dnorm(s,5, 2), type="l")
lines(s, dnorm(s,10, 2), col="blue")

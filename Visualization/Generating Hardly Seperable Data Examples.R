x = c(runif(100,0.2,0.3), runif(500,0,1),runif(100,0.7,0.8))
y = rep(0, length(x))
for(i in 1:length(x)){
  y[i] = rnorm(1, mean = x[i], sd = 0.25)
}
plot(x,y)
symbols(x=x, y=y, circles=rep(0.05,length(x)), add=T, inches=F)
require("ggplot2")
x = c(-5,x)
y = c(-1,y)
x = c(5,x)
y = c(-1,y)
df = data.frame(cbind(x,y)); colnames(df) = c("x","y")

commonTheme = list(labs(color="Density",fill="Density",
                        x="runif(100,0.2,0.3), runif(500,0,1),runif(100,0.7,0.8)",
                        y="N(x,0.25)"),
                    theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

m = ggplot(data=df,aes(x,y)) + 
  coord_cartesian(xlim = c(-.1,1.1), ylim = c(-0.5,1.5))

m +  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
  scale_color_hue(l = 20, c = 200) +
  guides(alpha="none") +
  ggtitle("Density-plot of Simulated Observations") +
  geom_point(color = "black") + commonTheme
  
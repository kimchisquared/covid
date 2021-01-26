library(deSolve)
state<-c(S=90,
         E=10,
         IP=0,
         IS=0,
         R=0)
N=100
# define parameter values for equation
parameters<-c(ap=0.6,
              as=0.2,
              ge=0.16,
              gp=0.4,
              rho=0.1)

# define DE function or system of DEs
LV<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dS <- -(ap*IP+as*IS)*S/N
    dE <- (ap*IP+as*IS)*S/N-ge*E
    dIP <- ge*E -gp*IP
    dIS <-gp*IP -rho*IS
    dR <- rho*IS
    # return the rate of change
    list(c(dS,dE, dIP, dIS,dR))
  }) 
}
# provide times (10 years, monthly points)
times<-seq(0,100,length.out=120)
# run and store output
out <- ode(y = state, times = times, func = LV, parms = parameters)
# look at output head(out)
head(out)

# plot output from solver
par(mar = c(5,2,2,5)) # Widen's the margin
plot(out[,1],out[,2],type="b",col="brown",ylab="S",xlab="Time",main=paste("SEIIR state ", state[1],state[2], state[3], state[4],state[5])) #plots S
par(new = T) #allows plot overlay
plot(out[,1],out[,3], type="b",axes=F, xlab=NA, ylab=NA, col="pink") #plots E
par(new = T)
plot(out[,1],out[,4], type="b",axes=F, xlab=NA, ylab=NA, col="orange") #plots I_p
par(new = T)
plot(out[,1],out[,5], type="b",axes=F, xlab=NA, ylab=NA, col="yellow") #plots I_s
par(new = T)
plot(out[,1],out[,6], type="b",axes=F, xlab=NA, ylab=NA, col="darkgreen") #plots R
legend("right", legend=c("S", "E", "I_P", "I_S", "R"),
       col=c("brown", "pink", "orange", "yellow", "darkgreen"), lty=1, cex=0.8)
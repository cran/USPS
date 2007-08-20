require(graphics); require(USPS)
if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
        (.Device %in% c("X11", "GTK", "gnome", "windows","quartz","JavaGD")))

#input the abciximab study data of Kereiakes et al. (2000).
data(lindner)
# outcomes: lifepres & cardbill

cat("\n Press ENTER for UNSUPERVISED analyses...")
scan()
# Define Divisive Cluster Hierarchy for UNSUPERVISED analyses...
UPSxvars <- c("stent","height","female","diabetic","acutemi","ejecfrac","ves1proc")
UPSharch <- UPShclus(lindner, UPSxvars)
UPSharch

# Save UPSpars settings for NN/IV analyses of the "lifepres" outcome
UPSaccum(UPSharch, lindner, abcix, lifepres, faclev=1, accobj="ABClife")

lif001nn <- UPSnnltd(  1)
lif002nn <- UPSnnltd(  2)
lif005nn <- UPSnnltd(  5)
lif010nn <- UPSnnltd( 10)
lif020nn <- UPSnnltd( 20)
plot(lif020nn)  # graphical display
lif030nn <- UPSnnltd( 30)
lif040nn <- UPSnnltd( 40)
plot(lif040nn)  # graphical display
lif050nn <- UPSnnltd( 50)
lif060nn <- UPSnnltd( 60)
summary(lif060nn) # brief console output
lif070nn <- UPSnnltd( 70)

lif030iv <- UPSivadj( 30)
lif040iv <- UPSivadj( 40)
plot(lif040iv)  # graphical display
lif050iv <- UPSivadj( 50)
lif060iv <- UPSivadj( 60)
lif070iv <- UPSivadj( 70)
lif080iv <- UPSivadj( 80)
lif090iv <- UPSivadj( 90)
lif100iv <- UPSivadj(100)
summary(lif100iv) # brief console output
lif200iv <- UPSivadj(200)
lif300iv <- UPSivadj(300)
lif996iv <- UPSivadj(996)

# Overall "Sensitivity Analysis" Summary...
UPSgraph()

# Display contents of UPSdf...
ABClife

# Compare NN/LTD and Artificial LTD distributions for 50 clusters
abcdf <- UPSaltdd(lindner, abcix, lifepres, faclev=1, NNobj=lif050nn)
abcdf
plot(abcdf)

cat("\n Press ENTER to make alternative SUPERVISED analyses...")
scan()

# Define Logit Model for Treatment Choice in SUPERVISED analyses
PStreat <- abcix~stent+height+female+diabetic+acutemi+ejecfrac+ves1proc
logtSPS <- SPSlogit(lindner, PStreat, PSfit, PSrnk, PSbin, appn="lindSPS")
# Store Propensity Score info (default=5 bins) in frame named "lindSPS"
logtSPS

# Testing for Within-Bin Balance on Continuous Covariate...
SPSbalvs <- SPSbalan(lindSPS, abcix, PSbin, ves1proc)
plot(SPSbalvs)
SPSbalvs

# Testing for Within-Bin Balance on Binary (Dichotomous) Covariate...
SPSbalam <- SPSbalan(lindSPS, abcix, PSbin, acutemi)
SPSbalam

# Test for Within-Bin Outcome Differences...
lindlife <- SPSoutco(lindSPS, abcix, PSbin, lifepres, faclev=1)
plot(lindlife)
lindlife

lindcost <- SPSoutco(lindSPS, abcix, PSbin, cardbill)
plot(lindcost) 
lindcost

# Cubic smoothing spline analyses...
SPScbill <- SPSsmoot(lindSPS, abcix, PSfit, cardbill)
plot(SPScbill)

# Loess "symmetric" smoothing analyses...
SPScblss <- SPSloess(lindSPS, abcix, PSfit, cardbill)
plot(SPScblss)

par(opar)
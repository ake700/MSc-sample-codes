/* Inhibition of Cronobacter sakazakii in the Simulator of the Human Intestinal Microbial Ecosystem */

/*Cs population in SHIME*/ 
/*Model = Treatment + population + time + replicate + error */
data cs_shime;
infile '/folders/myfolders/sasuser.v94/[filename].csv' firstobs=2 dlm=",";
length treatment$ 25;
input treatment$ replicate pop time; 
proc print data=cs_shime;
title 'Cs inhibition in SHIME';
run;

proc glimmix data=cs_shime order=data nobound;
class treatment replicate time;
model pop = treatment|time;
covtest "replicate = 0" 0 .;
random replicate; 
random time / subject=treatment type=cs residual;
lsmeans treatment*time / pdiff adjust=tukey lines slicediff=time;
title "Cs inhibition in SHIME";
output out=second predicted=pred residual=resid residual(noblup)=mresid
student=studentresid student(noblup)=smresid;
ods output contrasts=power;
run; 

/*linearity of fixed effects as scatter and box plot*/
proc sgplot data=second;
vbox smresid / group=treatment datalabel;
run;


/* homogeneity of effects */
proc sgscatter data=second;
plot studentresid*(pred treatment time); 
run;

/* normal distribution */ 
proc univariate data=second normal plot;
var studentresid; 
run;

Data power;
set power;
alpha = 0.05;
nc = numdf*fvalue;
fcrit = finv(1-alpha, numdf, dendf, 0);
power = 1-probf(fcrit, numdf, dendf, nc);
Run;
Proc print data=power;
var label numdf dendf alpha nc fcrit power;
title "Power Calculation";
Run;

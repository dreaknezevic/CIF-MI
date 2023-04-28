
*************************************************************;
*** Generate competing risks data with missing event type ;

data data; do ID=1 to 400; output; end; run;

data data; set data;
 *call streaminit(12345678);
 covar = rand('bernoulli',0.5)+1;

 if covar=1 then do;   
  T1 = rand('weibull', 0.7, 11); 
  T2 = rand('weibull', 1.3, 12); 
  C = rand('weibull', 1, 7); 
 end;

 if covar=2 then do; 
  T1 = rand('weibull', 1, 11);
  T2 = rand('weibull', 1, 12);
  C = rand('weibull', 1, 7);
 end;

 years = min(T1, T2, C); 
 if years=T1 then cause_complete=1; *event of interest ;
  else if years=T2 then cause_complete=2; *competing event ;
  else if years=C then cause_complete=0; *censor ;

 * Impose approximately 30% missing onto event type ;
 if cause_complete=0 then cause=0;
  else cause = ifn(rand('Binomial',.3,1),.,cause_complete);
 
 drop T1 T2 C;
run;

* Frequency of event type and missing ;
proc freq data=data; 
 where cause~=0;
 tables cause_complete cause / missing; 
run;

* Formats to use in display ;
proc format;
  value covar
  1="Group 1"
  2="Group 2"
 ;

 value cause 
  1="Dead of disease"
  2="Dead of other causes"
  3="Dead of unknown causes"
 ;
run;

************************************;
*** Calculate CIF from ground truth ; 
************************************;

ods exclude all;
proc lifetest data=data outcif=estimates;
 time years*cause_complete(0) / eventcode=(1,2);
 strata covar;
run;
ods exclude none;

*************************************************************;
*** Display:
    - Plot of CIF cuves 
      (event of interest and competing event for two covariate groups)
	- 5-year CIF estimates ;

proc sgpanel data=estimates;
 panelby failcode / novarname columns=1 headerattrs=(size=11pt weight=bold);
 format failcode cause. covar covar.;
 step x=years y=cif / group=covar lineattrs=(thickness=2) name="step";
 band x=years lower=cif_lcl upper=cif_ucl / group=covar type=step transparency=0.80 fillattrs=(color=gray);
 label cif="Cumulative Incidence" 
       years = "Years";
 rowaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 colaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 keylegend "step" / title="" valueattrs=(size=11pt weight=bold);
run;

title "5-year CIF estimates";
data estimates_5y; set estimates(where=(years>5)); 
 by failcode covar;
 if first.failcode|first.covar;
run; 

proc print data=estimates_5y noobs;
 var failcode covar cif cif_lcl cif_ucl;
 format failcode cause. covar covar. cif cif_lcl cif_ucl percent10.1;
run;

title;

*************************************************************;
*************************************************************;

***************************;
*** Calculate CIF using MI ; 
***************************;

* Multiple imputation for missing event type
   - Covariate data should be complete for events
   - In real data it is desirable to use many covariates for imputation 
     in our toy simulated data we only have one covariate ;
proc mi data=data out=mi nimpute=20 noprint;
 where cause~=0;
 class cause covar;
 var cause covar years;
 fcs discrim(cause / classeffects=include) nbiter=100; 
run;

data mi; set mi(rename=(_imputation_=n));
 keep ID n cause;
run;

* Create long dataset by imputation number ;
data patients; set data(keep=ID years covar); run;
data repeat; do n=1 to 20; output; end; run;
proc sql; create table list as select * from repeat cross join patients; quit;

data mi; merge list mi; 
 by n ID; 
 if cause=. then cause=0;
run;

proc datasets noprint; 
 delete patients repeat list; 
run; quit;

*************************************************************;
* Estimate CIFs for each multiply-imputed data set  ;

%let covar=covar;

ods exclude all;
proc lifetest data=mi outcif=estimates;
 by n;
 time years*cause(0) / eventcode=(1,2);
 strata &covar;
run;
ods exclude none;

data estimates; set estimates;
 t = years;
 s = cif;
 var = stderr**2;
 keep n failcode &covar t s var;
run;

*************************************************************;
* Combine CIF estimates ;

%include "combine_CIF_macros.sas";

%merge(in=estimates(where=(failcode=1 & &covar=1)),out=dod_out1,n_imp=20);
%merge(in=estimates(where=(failcode=2 & &covar=1)),out=doc_out1,n_imp=20);
%merge(in=estimates(where=(failcode=1 & &covar=2)),out=dod_out2,n_imp=20);
%merge(in=estimates(where=(failcode=2 & &covar=2)),out=doc_out2,n_imp=20);

%rubin(in=dod_out1,out=combined_cif_dod1);
%rubin(in=doc_out1,out=combined_cif_doc1);
%rubin(in=dod_out2,out=combined_cif_dod2);
%rubin(in=doc_out2,out=combined_cif_doc2);

data zero; time=0; s_bar=0; lcl=0; ucl=0; run;
data combined_cif_doc1; set zero combined_cif_doc1; run;
data combined_cif_dod1; set zero combined_cif_dod1; run;
data combined_cif_doc2; set zero combined_cif_doc2; run;
data combined_cif_dod2; set zero combined_cif_dod2; run;

data estimates(rename=(time=years s_bar=cif lcl=cif_lcl ucl=cif_ucl)); 
 retain failcode &covar;
 set combined_cif_dod1(in=a) combined_cif_dod2(in=b)
     combined_cif_doc1(in=c) combined_cif_doc2(in=d);
 if a then do; failcode=1; &covar=1; end;
 if b then do; failcode=1; &covar=2; end;
 if c then do; failcode=2; &covar=1; end;
 if d then do; failcode=2; &covar=2; end;
run;

proc datasets noprint; 
 delete zero dod_out1 doc_out1 dod_out2 doc_out2
        combined_cif_dod1 combined_cif_dod2
        combined_cif_doc1 combined_cif_doc2; 
run; quit;

*************************************************************;
*** Display:
    - Plot of CIF cuves 
      (event of interest and competing event for two covariate groups)
	- 5-year CIF estimates ;

proc sgpanel data=estimates;
 panelby failcode / novarname columns=1 headerattrs=(size=11pt weight=bold);
 format failcode cause. covar covar.;
 step x=years y=cif / group=covar lineattrs=(thickness=2) name="step";
 band x=years lower=cif_lcl upper=cif_ucl / group=covar type=step transparency=0.80 fillattrs=(color=gray);
 label cif="Cumulative Incidence" 
       years = "Years";
 rowaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 colaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 keylegend "step" / title="" valueattrs=(size=11pt weight=bold);
run;

title "5-year CIF estimates";
data estimates_5y; set estimates(where=(years>5)); 
 by failcode covar;
 if first.failcode|first.covar;
run; 

proc print data=estimates_5y noobs;
 var failcode covar cif cif_lcl cif_ucl;
 format failcode cause. covar covar. cif cif_lcl cif_ucl percent10.1;
run;

title;

*************************************************************;
*** For comparison: 
	- Complete case analysis
	- Extra cause analysis ;

******************;
*** Complete case ;
******************;

ods exclude all;
proc lifetest data=data outcif=estimates;
 time years*cause(0) / eventcode=(1,2);
 strata covar;
run;
ods exclude none;

proc sgpanel data=estimates;
 panelby failcode / novarname columns=1 headerattrs=(size=11pt weight=bold);
 format failcode cause. covar covar.;
 step x=years y=cif / group=covar lineattrs=(thickness=2) name="step";
 band x=years lower=cif_lcl upper=cif_ucl / group=covar type=step transparency=0.80 fillattrs=(color=gray);
 label cif="Cumulative Incidence" 
       years = "Years";
 rowaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 colaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 keylegend "step" / title="" valueattrs=(size=11pt weight=bold);
run;

title "5-year CIF estimates";
data estimates_5y; set estimates(where=(years>5)); 
 by failcode covar;
 if first.failcode|first.covar;
run; 

proc print data=estimates_5y noobs;
 var failcode covar cif cif_lcl cif_ucl;
 format failcode cause. covar covar. cif cif_lcl cif_ucl percent10.1;
run;

title;

****************;
*** Extra cause ;
****************;

data data; set data;
 cause_ec = cause;
 if cause=. then cause_ec=3;
run;

ods exclude all;
proc lifetest data=data outcif=estimates;
 time years*cause_ec(0) / eventcode=(1,2,3);
 strata covar;
run;
ods exclude none;

proc sgpanel data=estimates;
 panelby failcode / novarname columns=1 headerattrs=(size=11pt weight=bold);
 format failcode cause. covar covar.;
 step x=years y=cif / group=covar lineattrs=(thickness=2) name="step";
 band x=years lower=cif_lcl upper=cif_ucl / group=covar type=step transparency=0.80 fillattrs=(color=gray);
 label cif="Cumulative Incidence" 
       years = "Years";
 rowaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 colaxis labelattrs=(size=11pt weight=bold) valueattrs=(size=11pt weight=bold);
 keylegend "step" / title="" valueattrs=(size=11pt weight=bold);
run;

title "5-year CIF estimates";
data estimates_5y; set estimates(where=(years>5)); 
 by failcode covar;
 if first.failcode|first.covar;
run; 

proc print data=estimates_5y noobs;
 var failcode covar cif cif_lcl cif_ucl;
 format failcode cause. covar covar. cif cif_lcl cif_ucl percent10.1;
run;

title;

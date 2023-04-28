%macro merge(in=,out=,n_imp=);

data out; stop; format t best12.; run;

%do i=1 %to &n_imp;
 data cif; set &in; if n=&i; run; 
 data cif; set cif(rename=(s=s&i var=var&i));
  keep t s&i var&i;
 run;
data out; merge out cif; by t; run;
%end;

* Carry forward survival estimates/variances where 
  imputed DUC observation was no previously included ;
%do i=1 %to &n_imp;
 data out(drop=s var); set out;
  retain s var;
  if s&i~=. then s=s&i;
  if var&i~=. then var=var&i;
  s&i=s;
  var&i=var;
 run;
%end;

* remove first row of zeros ;
%do i=1 %to &n_imp;
 data out; set out; if s&i=0 then delete; run;
%end;

data &out; set out; format t best12.; run;

proc datasets nolist;
 delete cif out;
run; quit;

%mend;

*****************************************************************;
*****************************************************************;

%macro rubin(in=,out=,n_imp=);

proc iml;
 *reset log print;
 use &in; read all var{t} into time; close &in;
 use &in; read all var('s1':'s20') into s; close &in;
 use &in; read all var('var1':'var20') into var; close &in;

 * Calculate survival curve ;
 Q_hat = log(-log(1-s));
 Q_bar = Q_hat[,+]/ncol(s);
 S_bar = 1-exp(-exp(Q_bar));

 * Calculate variance, confidence limits ;
 temp = var/(log(1-s)#(1-s))##2;
 U_bar = temp[,+]/ncol(s);

 temp = (Q_hat-Q_bar)##2;
 B = temp[,+] # (1/(ncol(s)-1));

 T = U_bar+(1+1/ncol(s))#B;

 LCL = 1-exp(-exp(Q_bar-1.96#T##(1/2)));
 UCL = 1-exp(-exp(Q_bar+1.96#T##(1/2)));

 out = time || S_bar || LCL || UCL;
 create &out from out[colname={time S_bar LCL UCL}]; append from out;
quit;

%mend;

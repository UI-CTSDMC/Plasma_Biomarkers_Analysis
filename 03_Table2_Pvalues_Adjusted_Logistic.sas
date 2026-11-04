/*==========================================================
  03_Table2_Pvalues_Adjusted_Logistic.sas
  Adjusted p-values for binary outcomes using logistic regression
  03/06: Adjusted logistic p-values + continuous biomarker OR/CI
==========================================================*/

%macro _var_exists(ds=, var=, mvar=exists);
  %local dsid vnum rc;
  %let &mvar=0;
  %let dsid=%sysfunc(open(&ds));
  %if &dsid %then %do;
    %let vnum=%sysfunc(varnum(&dsid,&var));
    %let rc=%sysfunc(close(&dsid));
    %if &vnum>0 %then %let &mvar=1;
  %end;
%mend;


%macro auto_cat_pvals(
    data=,
    outcome=,
    groupvar=group,
    agevar=age,
    genvar=gen,
    group_ref=1,
    gen_ref=1,
    event_ref=0,
    cont_predictor=,
    main=,
    order=0,
    out=
);

  ods exclude all;

  /* Group model */
  ods output ModelANOVA=__type3_grp
             ContrastTest=__ctr_grp;

  proc logistic data=&data order=internal;
    class &groupvar &genvar(ref="&gen_ref") / param=ref;
    model &outcome(event="&event_ref") = &groupvar &agevar &genvar;

    contrast "1_vs_2" &groupvar  1 -1  0;
    contrast "1_vs_3" &groupvar  1  0 -1;
    contrast "2_vs_3" &groupvar  0  1 -1;
  run;
  quit;

  data __overall_grp;
    set __type3_grp;
    where strip(upcase(Effect)) = strip(upcase("&groupvar"));
    P_123_n = ProbChiSq;
    keep P_123_n;
  run;

  data __pairwise_grp;
    length P_12_n P_13_n P_23_n 8.;
    retain P_12_n P_13_n P_23_n;
    set __ctr_grp end=last;

    _lab = upcase(strip(Contrast));

    if index(_lab,"1_VS_2") then P_12_n = ProbChiSq;
    else if index(_lab,"1_VS_3") then P_13_n = ProbChiSq;
    else if index(_lab,"2_VS_3") then P_23_n = ProbChiSq;

    if last then output;
    keep P_12_n P_13_n P_23_n;
  run;

  /* Continuous biomarker model */
  ods output ModelANOVA=__type3_cont
             OddsRatios=__or_cont;

  proc logistic data=&data order=internal;
    class &genvar(ref="&gen_ref") / param=ref;
    model &outcome(event="&event_ref") = &cont_predictor &agevar &genvar;
  run;
  quit;

  data __overall_cont;
    set __type3_cont;
    where strip(upcase(Effect)) = strip(upcase("&cont_predictor"));
    overall_cont_n = ProbChiSq;
    keep overall_cont_n;
  run;

  data __or_cont2;
    set __or_cont;
    where strip(upcase(Effect)) = strip(upcase("&cont_predictor"));

    length EST_CI $60.;
    EST_CI = cats(
      put(OddsRatioEst, oddsr8.4),
      " (",
      put(LowerCL, oddsr8.4),
      ", ",
      put(UpperCL, oddsr8.4),
      ")"
    );
    keep EST_CI;
  run;

  data &out;
    length P_123 P_12 P_13 P_23 overall_cont $30
           type $20 transform $12 var $64 EST_CI $60;
    merge __overall_grp __pairwise_grp __overall_cont __or_cont2;

    var="&outcome";
    type="adj_cat";
    transform="logistic";
    main=&main;
    order=&order;

    if missing(P_123_n) then P_123="";
    else if P_123_n < 0.0001 then P_123="<0.0001^{super a}";
    else P_123=cats(put(P_123_n, pvalue6.4),"^{super a}");

    if missing(P_12_n) then P_12="";
    else if P_12_n < 0.0001 then P_12="<0.0001^{super a}";
    else P_12=cats(put(P_12_n, pvalue6.4),"^{super a}");

    if missing(P_13_n) then P_13="";
    else if P_13_n < 0.0001 then P_13="<0.0001^{super a}";
    else P_13=cats(put(P_13_n, pvalue6.4),"^{super a}");

    if missing(P_23_n) then P_23="";
    else if P_23_n < 0.0001 then P_23="<0.0001^{super a}";
    else P_23=cats(put(P_23_n, pvalue6.4),"^{super a}");

    if missing(overall_cont_n) then overall_cont="";
    else if overall_cont_n < 0.0001 then overall_cont="<0.0001^{super a}";
    else overall_cont=cats(put(overall_cont_n, pvalue6.4),"^{super a}");

    keep var main order type transform EST_CI overall_cont
         P_123 P_12 P_13 P_23
         P_123_n P_12_n P_13_n P_23_n overall_cont_n;
  run;

  proc datasets lib=work nolist;
    delete __type3_grp __ctr_grp __overall_grp __pairwise_grp
           __type3_cont __or_cont __overall_cont __or_cont2;
  quit;

  ods exclude none;

%mend;

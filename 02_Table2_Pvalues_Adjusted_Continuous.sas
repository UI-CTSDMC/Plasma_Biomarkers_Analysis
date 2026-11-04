/*==========================================================
  02_Table2_Pvalues_Adjusted_Continuous.sas
  Adjusted p-values for continuous outcomes (auto transform + GLM)
==========================================================*/


ods escapechar="^";

/* Superscript markers: b = adjusted GLM; c = log;
   d = square root; e = squared value; f = rank-based. */
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

/*==========================================================
  Macro 1: group-based adjusted p-values
  Returns exactly 1 row per outcome:
    P_123, P_12, P_13, P_23
==========================================================*/
%macro auto_transform_pvals(
    data=,
    outcome=,
    class=,
    predictors=,
    groupvar=group,
    alpha=0.05,
    main=,
    order=0,
    out=
);

  ods exclude all;

  %local p_raw p_log p_sqrt p_sq chosen yvar type stars;
  %let p_raw=.;  %let p_log=.;  %let p_sqrt=.;  %let p_sq=.;
  %let chosen=;

  /*-------------------------------
    RAW
  -------------------------------*/
  proc glm data=&data;
    class &class;
    model &outcome = &predictors / e3;
    output out=__res_raw r=residual;
  run; quit;

  ods output TestsForNormality=__norm_raw;
  proc univariate data=__res_raw normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_raw;
    where Test="Shapiro-Wilk";
    call symputx("p_raw", pValue, "L");
  run;

  /*-------------------------------
    LOG
  -------------------------------*/
  data __log_data;
    set &data;
    if &outcome > 0 then y_log = log(&outcome);
    else y_log = .;
  run;

  proc glm data=__log_data;
    class &class;
    model y_log = &predictors / e3;
    output out=__res_log r=residual;
  run; quit;

  ods output TestsForNormality=__norm_log;
  proc univariate data=__res_log normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_log;
    where Test="Shapiro-Wilk";
    call symputx("p_log", pValue, "L");
  run;

  /*-------------------------------
    SQRT
  -------------------------------*/
  data __sqrt_data;
    set &data;
    if &outcome >= 0 then y_sqrt = sqrt(&outcome);
    else y_sqrt = .;
  run;

  proc glm data=__sqrt_data;
    class &class;
    model y_sqrt = &predictors / e3;
    output out=__res_sqrt r=residual;
  run; quit;

  ods output TestsForNormality=__norm_sqrt;
  proc univariate data=__res_sqrt normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_sqrt;
    where Test="Shapiro-Wilk";
    call symputx("p_sqrt", pValue, "L");
  run;

  /*-------------------------------
    SQUARE
  -------------------------------*/
  data __sq_data;
    set &data;
    y_sq = (&outcome)**2;
  run;

  proc glm data=__sq_data;
    class &class;
    model y_sq = &predictors / e3;
    output out=__res_sq r=residual;
  run; quit;

  ods output TestsForNormality=__norm_sq;
  proc univariate data=__res_sq normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_sq;
    where Test="Shapiro-Wilk";
    call symputx("p_sq", pValue, "L");
  run;

  /*-------------------------------
    Choose transformation
  -------------------------------*/
  %if %sysevalf(&p_raw  > &alpha) %then %let chosen=raw;
  %else %if %sysevalf(&p_log  > &alpha) %then %let chosen=log;
  %else %if %sysevalf(&p_sqrt > &alpha) %then %let chosen=sqrt;
  %else %if %sysevalf(&p_sq   > &alpha) %then %let chosen=sq;
  %else %let chosen=rank;

  %if "&chosen"="raw" %then %do;
    data __final; set &data; run;
    %let yvar=&outcome;
    %let type=adj_cont;
    %let stars=^{super b};
  %end;
  %else %if "&chosen"="log" %then %do;
    data __final; set __log_data; run;
    %let yvar=y_log;
    %let type=log;
    %let stars=^{super c};
  %end;
  %else %if "&chosen"="sqrt" %then %do;
    data __final; set __sqrt_data; run;
    %let yvar=y_sqrt;
    %let type=sqrt;
    %let stars=^{super d};
  %end;
  %else %if "&chosen"="sq" %then %do;
    data __final; set __sq_data; run;
    %let yvar=y_sq;
    %let type=sq;
    %let stars=^{super e};
  %end;
  %else %do;
    proc rank data=&data out=__final ties=mean;
      var &outcome;
      ranks y_rank;
    run;
    %let yvar=y_rank;
    %let type=rank;
    %let stars=^{super f};
  %end;

  /*-------------------------------
    Final group model
  -------------------------------*/
  ods output ModelANOVA=__anova_final
             Contrasts=__ctr;

  proc glm data=__final;
    class &class;
    model &yvar = &predictors / e3;

    contrast "1_vs_2" &groupvar  1 -1  0;
    contrast "1_vs_3" &groupvar  1  0 -1;
    contrast "2_vs_3" &groupvar  0  1 -1;
  run; quit;

  /* IMPORTANT:
     ModelANOVA may contain both Type I and Type III tests.
     Keep only the Type III row for group.
     Depending on SAS version, the variable may be HypothesisType or Type.
  */
  %local has_hyp has_type;
  %_var_exists(ds=work.__anova_final, var=HypothesisType, mvar=has_hyp);
  %_var_exists(ds=work.__anova_final, var=Type,           mvar=has_type);

  data __overall;
    set __anova_final;
    where upcase(Source)=upcase("&groupvar")
      %if &has_hyp=1 %then and HypothesisType=3;
      %else %if &has_type=1 %then and Type=3;
      ;
    P_123_n = ProbF;
    keep P_123_n;
  run;

  %local has_source has_label has_contrast idcol;
  %_var_exists(ds=work.__ctr, var=Source,   mvar=has_source);
  %_var_exists(ds=work.__ctr, var=Label,    mvar=has_label);
  %_var_exists(ds=work.__ctr, var=Contrast, mvar=has_contrast);

  %if &has_source=1 %then %let idcol=Source;
  %else %if &has_label=1 %then %let idcol=Label;
  %else %if &has_contrast=1 %then %let idcol=Contrast;
  %else %let idcol=;

  data __pairwise;
    length P_12_n P_13_n P_23_n 8.;
    retain P_12_n P_13_n P_23_n;
    set __ctr end=last;

    length _lab $200;
    %if %length(&idcol) %then %do;
      _lab = upcase(strip(&idcol));
    %end;
    %else %do;
      _lab = "";
    %end;

    if index(_lab,"1_VS_2") then P_12_n = ProbF;
    else if index(_lab,"1_VS_3") then P_13_n = ProbF;
    else if index(_lab,"2_VS_3") then P_23_n = ProbF;

    if last then output;
    keep P_12_n P_13_n P_23_n;
  run;

  data &out;
    length P_123 P_12 P_13 P_23 $30 type $12 transform $8 var $64;
    merge __overall __pairwise;

    var="&outcome";
    type="&type";
    transform="&chosen";
    main=&main;
    order=&order;

    if missing(P_123_n) then P_123="";
    else if P_123_n < 0.0001 then P_123=cats("<0.0001 ", "&stars");
    else P_123=cats(put(P_123_n, pvalue6.4), " ", "&stars");

    if missing(P_12_n) then P_12="";
    else if P_12_n < 0.0001 then P_12=cats("<0.0001 ", "&stars");
    else P_12=cats(put(P_12_n, pvalue6.4), " ", "&stars");

    if missing(P_13_n) then P_13="";
    else if P_13_n < 0.0001 then P_13=cats("<0.0001 ", "&stars");
    else P_13=cats(put(P_13_n, pvalue6.4), " ", "&stars");

    if missing(P_23_n) then P_23="";
    else if P_23_n < 0.0001 then P_23=cats("<0.0001 ", "&stars");
    else P_23=cats(put(P_23_n, pvalue6.4), " ", "&stars");

    keep var main order type transform P_123 P_12 P_13 P_23
         P_123_n P_12_n P_13_n P_23_n;
  run;

  proc datasets lib=work nolist;
    delete __res_raw __norm_raw __res_log __norm_log __res_sqrt __norm_sqrt
           __res_sq __norm_sq __log_data __sqrt_data __sq_data
           __final __anova_final __ctr __overall __pairwise;
  quit;

  ods exclude none;

%mend auto_transform_pvals;


/*==========================================================
  Macro 2: continuous biomarker association
  Returns exactly 1 row per outcome:
    EST_CI, overall_cont
  Uses the same chosen transformation model and pulls the
  continuous biomarker p-value from ParameterEstimates.
==========================================================*/
%macro auto_cont_biomarker_assoc(
    data=,
    outcome=,
    cont_predictor=,
    class=gen,
    predictors=,
    alpha=0.05,
    main=,
    order=0,
    out=
);

  ods exclude all;

  %local p_raw p_log p_sqrt p_sq chosen yvar type stars;
  %let p_raw=.;  %let p_log=.;  %let p_sqrt=.;  %let p_sq=.;
  %let chosen=;

  /*-------------------------------
    RAW
  -------------------------------*/
  proc glm data=&data;
    class &class;
    model &outcome = &predictors;
    output out=__res_raw r=residual;
  run; quit;

  ods output TestsForNormality=__norm_raw;
  proc univariate data=__res_raw normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_raw;
    where Test="Shapiro-Wilk";
    call symputx("p_raw", pValue, "L");
  run;

  /*-------------------------------
    LOG
  -------------------------------*/
  data __log_data;
    set &data;
    if &outcome > 0 then y_log = log(&outcome);
    else y_log = .;
  run;

  proc glm data=__log_data;
    class &class;
    model y_log = &predictors;
    output out=__res_log r=residual;
  run; quit;

  ods output TestsForNormality=__norm_log;
  proc univariate data=__res_log normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_log;
    where Test="Shapiro-Wilk";
    call symputx("p_log", pValue, "L");
  run;

  /*-------------------------------
    SQRT
  -------------------------------*/
  data __sqrt_data;
    set &data;
    if &outcome >= 0 then y_sqrt = sqrt(&outcome);
    else y_sqrt = .;
  run;

  proc glm data=__sqrt_data;
    class &class;
    model y_sqrt = &predictors;
    output out=__res_sqrt r=residual;
  run; quit;

  ods output TestsForNormality=__norm_sqrt;
  proc univariate data=__res_sqrt normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_sqrt;
    where Test="Shapiro-Wilk";
    call symputx("p_sqrt", pValue, "L");
  run;

  /*-------------------------------
    SQUARE
  -------------------------------*/
  data __sq_data;
    set &data;
    y_sq = (&outcome)**2;
  run;

  proc glm data=__sq_data;
    class &class;
    model y_sq = &predictors;
    output out=__res_sq r=residual;
  run; quit;

  ods output TestsForNormality=__norm_sq;
  proc univariate data=__res_sq normal;
    var residual;
  run; quit;

  data _null_;
    set __norm_sq;
    where Test="Shapiro-Wilk";
    call symputx("p_sq", pValue, "L");
  run;

  /*-------------------------------
    Choose transformation
  -------------------------------*/
  %if %sysevalf(&p_raw  > &alpha) %then %let chosen=raw;
  %else %if %sysevalf(&p_log  > &alpha) %then %let chosen=log;
  %else %if %sysevalf(&p_sqrt > &alpha) %then %let chosen=sqrt;
  %else %if %sysevalf(&p_sq   > &alpha) %then %let chosen=sq;
  %else %let chosen=rank;

  %if "&chosen"="raw" %then %do;
    data __final; set &data; run;
    %let yvar=&outcome;
    %let type=adj_cont_assoc;
    %let stars=^{super b};
  %end;
  %else %if "&chosen"="log" %then %do;
    data __final; set __log_data; run;
    %let yvar=y_log;
    %let type=log_assoc;
    %let stars=^{super c};
  %end;
  %else %if "&chosen"="sqrt" %then %do;
    data __final; set __sqrt_data; run;
    %let yvar=y_sqrt;
    %let type=sqrt_assoc;
    %let stars=^{super d};
  %end;
  %else %if "&chosen"="sq" %then %do;
    data __final; set __sq_data; run;
    %let yvar=y_sq;
    %let type=sq_assoc;
    %let stars=^{super e};
  %end;
  %else %do;
    proc rank data=&data out=__final ties=mean;
      var &outcome;
      ranks y_rank;
    run;
    %let yvar=y_rank;
    %let type=rank_assoc;
    %let stars=^{super f};
  %end;

  /*-------------------------------
    Final continuous biomarker model
  -------------------------------*/
  ods output ParameterEstimates=__parm_cont;

  proc glm data=__final;
    class &class;
    model &yvar = &predictors / solution clparm;
  run; quit;

  data __parm_cont2;
    set __parm_cont;
    where upcase(Parameter)=upcase("&cont_predictor");

    length EST_CI $60. est_txt lcl_txt ucl_txt $15.;

    /* coefficient p-value from chosen model */
    overall_cont_n = Probt;

  /* Estimate */
  if not missing(Estimate) then do;
    if Estimate > 999.9999 then est_txt = ">999.9999";
    else if 0 < Estimate < 0.0001 then est_txt = "<0.0001";
    else if -0.0001 < Estimate < 0 then est_txt = ">-0.0001";
    else if Estimate < -999.9999 then est_txt = "<-999.9999";
    else est_txt = strip(put(round(Estimate, 0.0001), 10.4));
  end;
  else est_txt = "";

  /* Lower CL */
  if not missing(LowerCL) then do;
    if LowerCL > 999.9999 then lcl_txt = ">999.9999";
    else if 0 < LowerCL < 0.0001 then lcl_txt = "<0.0001";
    else if -0.0001 < LowerCL < 0 then lcl_txt = ">-0.0001";
    else if LowerCL < -999.9999 then lcl_txt = "<-999.9999";
    else lcl_txt = strip(put(round(LowerCL, 0.0001), 10.4));
  end;
  else lcl_txt = "";

  /* Upper CL */
  if not missing(UpperCL) then do;
    if UpperCL > 999.9999 then ucl_txt = ">999.9999";
    else if 0 < UpperCL < 0.0001 then ucl_txt = "<0.0001";
    else if -0.0001 < UpperCL < 0 then ucl_txt = ">-0.0001";
    else if UpperCL < -999.9999 then ucl_txt = "<-999.9999";
    else ucl_txt = strip(put(round(UpperCL, 0.0001), 10.4));
  end;
  else ucl_txt = "";

  if est_txt ne "" then
    EST_CI = cats(est_txt, " (", lcl_txt, ", ", ucl_txt, ")");
  else EST_CI = "";

  keep EST_CI overall_cont_n;
run;


  data &out;
    set __parm_cont2;

    length var $64 type $20 transform $12 overall_cont $30;

    var="&outcome";
    type="&type";
    transform="&chosen";
    main=&main;
    order=&order;

    if missing(overall_cont_n) then overall_cont="";
    else if overall_cont_n < 0.0001 then overall_cont=cats("<0.0001 ", "&stars");
    else overall_cont=cats(put(overall_cont_n, pvalue6.4), " ", "&stars");

    keep var main order type transform EST_CI overall_cont overall_cont_n;
  run;

  proc datasets lib=work nolist;
    delete __res_raw __norm_raw __res_log __norm_log __res_sqrt __norm_sqrt
           __res_sq __norm_sq __log_data __sqrt_data __sq_data
           __final __parm_cont __parm_cont2;
  quit;

  ods exclude none;

%mend auto_cont_biomarker_assoc;

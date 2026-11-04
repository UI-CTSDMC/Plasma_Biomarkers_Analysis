*-------------------------------------------------------*;
* PROGRAMMER:        Krittika Chatterjee                *;
*                                                       *;
* DATE CREATED:      DECEMBER 23, 2025                  *;
*                                                       *;
* FILE NAME:         combined biomarker dataset         *;
*                                                       *;
* DESCRIPTION: Combined biomarkers (pTau217p, NFL, GFAP)*;
*              Table 1                                  *;
*-------------------------------------------------------*;

/***************************************************************************
* 0. LIBRARIES AND STANDARD INCLUDES
***************************************************************************/

%let filepath = /*path where data and macro code are located*/;
libname analysis "&filepath.";

%let report_date = 08DEC2025;

%include
    "filename/Table Macros.sas";

/* unadjusted pvalue macro */

/*Perform Two sample tests*/
%Macro TwoSampleTest(data, Var, by, label, formatCI=comma13.2);
	Data temporary;
		set &Data.;
	run;
	proc sort data=temporary;
		by &by.;
	run;
	PROC UNIVARIATE DATA=temporary NORMAL noprint;
		by &by.;
		VAR &var.;
		output out=Normal/*(where=(Test = "Shapiro-Wilk"))*/ PROBN =PROBN ;
	RUN;
	proc transpose data=Normal out=NormalT Prefix=PROBN;
		id &by.;
		Var PROBN;
	run;
	Data Normal2;
		set NormalT;
		by _NAME_;
		if PROBN0 >=.05 and PROBN1 >=.05 then Normal=1;
		else normal = 0;
		ShapiroNormal = strip(put(PROBN0, pvalue6.4))||" - "||strip(put(PROBN1, pvalue6.4));
		keep ShapiroNormal Normal;
	run;
	proc ttest data=temporary ;
		class &by.;
		var &var.;
		ods output Ttests=TTest1 Equality= Equality1 ConfLimits=ConfLimits statistics=statistics;
	run;
	data ConfLimits2;
		set confLimits;
		Length  MeanDIFF95 $50.;
		MeanDIFF95 = strip(put(Mean, &formatCI.))||" ("|| strip(put(LowerCLMean, &formatCI.))||", "||strip(put(UpperCLMean, &formatCI.))||")"; 
		if method in ("Pooled" "Satterthwaite");
	run;
	proc transpose data=ConfLimits2 out=ConfLimitsT Prefix=CI;
		by Variable;
		id Variances;
		Var MeanDIFF95 ;
	run;
	run;
	Data TTest1;
		merge Ttest1 Equality1;
		by Variable;
		if ProbF < .05 then do;
			if Variances = "Unequal";
		end;
		else do;
			if Variances = "Equal";
		end;
	run;
	
	proc npar1way data=temporary  noprint wilcoxon ;
		class &by.;
		var &var.;
		output out=Wilcoxon;
	run;
	proc quantreg data=temporary order=data CI=RESAMPLING(nrep=10000)/*sparsity(BF)*/;
		class &by.;
		model &var. = &by. / quantile=0.5 seed=1234;
   		estimate 'Diff in Medians' &by. 1 -1 / CL seed=1234;
		 ods output estimates=MedianDIFF;
	run;
	data MedianDIFF;
		Set MedianDIFF;
		MedianDIFF = Estimate;
		pmedian = probt;
		SE_median = StdErr;
		MedianCI = strip(put(Estimate, &formatCI.))||" ("|| strip(put(Lower, &formatCI.))||", "||strip(put(Upper, &formatCI.))||")";
		if Lower <= 0 <= Upper then SigMED = 0;
		else SigMEd = 1;
		keep MedianDIFF pmedian SE_median MedianCI SigMED;
	run;
	Data TTEST1;
		merge TTest1 Wilcoxon ConfLimitsT MedianDIFF Normal2;
		Format P2_wil pvalue6.4;
		keep variable probt ProbF Variances P2_WIL PAIR_COMPARED MedianDIFF pmedian SE_median CIEqual CIUnEqual MedianCI ShapiroNormal Normal SigMED   _KW_ DF_KW P_KW;
	run;	
    data TTest; *rename variables for final results;
		length Variable label $200.;
		set TTest1;        
		Variable = "&Var.";
		label  = "&Label.";
		keep Label variable probt ProbF Variances P2_WIL PAIR_COMPARED MedianDIFF pmedian SE_median CIEqual CIUnEqual MedianCI ShapiroNormal Normal SigMED  _KW_ DF_KW P_KW;        
	run; 
	proc datasets lib=work nolist; 
  		delete TTest1 Wilcoxon equality1 MedianDIFF ConfLimits2 Normal2 temporary; 
	quit; 
	run; 
%Mend; 


/*Perform Chi-square or fisher exact tests*/
%Macro ChiSquareFisher(data, Var, by, label, formatCI=comma13.2);
	Data temporary;
		set &Data.;
		&by.2 = strip(put(&by., $8.));
	run;
	proc sort data=temporary;
		by descending &by. descending &var.;
	run;
 	proc freq data=temporary order=data;   
		tables &by.*&var. / out = Counts CHISQ(WARN=OUTPUT) measures riskdiff Relrisk(CL=Wald) OR(CL=Wald);
		ods output RelativeRiskCLs=RelativeRiskCLs(Rename=(LowerCL=RRLowerCL UpperCL=RRUpperCL)) 
		OddsRatioCLs=OddsRatioCLs(Rename=(LowerCL=ORLowerCL UpperCL=ORUpperCL)) RelativeRisks =RelativeRisks RiskDiffCol1=RiskDiffCol1;
		output out=CTout pchi;
	run;
	data RiskDiffCol1;
		set RiskDiffCol1;
		if row = "Difference";
		diffprop = Risk;
		diffpropse = ASE;
		diffpropCI = strip(put(Risk*100, &formatCI.))||" ("|| strip(put(LowerCL*100, &formatCI.))||", "||strip(put(UpperCL*100, &formatCI.))||")"; 
		keep diffprop diffpropse diffpropCI;
	run;
	Data CTOUTb;
		set CTouT;
		if WARN_PCHI = 1 THEN FLAG +1;
		else FLAG + 0;
	run;
	Data _NULL_;
		set CTOUTb end=end1;
		if end1;
		if FLAG > 0 then Call Symput('fisher', 1);
		Else Call Symput('fisher', 0);
	run;

	%put "Fisher &Fisher.";
	proc freq data=Temporary noprint order=data; 
		where not missing(&var.);
		output out=CTout fisher PCHI;                 
		tables &by.*&var. /chisq cellchi2 nopercent;       
		%if &Fisher. = 1 %Then %DO;
			exact fisher;
		%End;   
	run;
   data CTout1;                 
		set CTout;                 
		needFisher = &fisher.;  
	run;            
	data ChiTest2;                 
		set CTout1;  
	run; 
	proc datasets lib=work nolist; 
  		delete Ctout CTout1 CTOUTb; 
	quit;
	run; 
	data ChiTest2; *rename variables for final results;
		merge ChiTest2 OddsRatioCLs RelativeRiskCLs RiskDiffCol1;        
		CHI_STAT=_PCHI_; DF=DF_PCHI; P_VALUE=P_PCHI; 

		ORS = strip(put(OddsRatio, &formatCI.))||" ("|| strip(put(ORlowerCL, &formatCI.))||", "||strip(put(ORUpperCL, &formatCI.))||")"; 
		RRs = strip(put(RelativeRisk, &formatCI.))||" ("|| strip(put(RRlowerCL, &formatCI.))||", "||strip(put(RRUpperCL, &formatCI.))||")"; 

		if ORlowerCL <= 1 <= ORUpperCL then SigOR = 0;
		else SigOR = 1;

		Variable = "&Var.";
		label  = "&Label.";

		if NeedFisher = 1 then P_VALUE = XP2_FISH;
		
		keep Variable label Chi_stat DF P_PCHI p_value XP2_FISH needFisher ORS RRS diffprop diffpropse diffpropCI SIGOR;        
	run; 

	/*Proc Print data=ChiTest2;
	var Pair_Compared numGroups numComparisons FDR p_value CHI XP2_FISH Fisher;
	run;*/
	proc datasets lib=work nolist; 
  		delete ChiTest RiskDiffCol1 OddsRatioCLs RelativeRiskCLs CTOUT: counts; 
	quit; 
	
%mend; 

%macro Teststats(Data, Vars, Types, by, labels, title, adj, transforms, appendresults, num);
	title;
		data Data;
			set &data.;
		run;
	ods select none;
	%Do I = 1 %to %sysfunc(countw(&vars.));
		/*Perform Two sample tests for continous variables*/
		%If "%Scan(&Types., &I.)" = "con" %Then %Do;
		%TwoSampleTest(data=&Data., Var=%Scan(&Vars., &I.) , by=&BY., label=);
			data Stats07;
				set TTest;
			run; 
			proc datasets noprint;
				delete TTEST;
			run;
		%end;
		/*Perform Chi-square Fisher's exact tests for ccatgorical variables*/
		%If "%Scan(&Types., &I.)" = "cat" %Then %Do;
			%ChiSquareFisher(data=&Data., Var=%Scan(&Vars., &I.) , by=&BY., label=);
			data Stats07;
				set chitest2;
			run;
			proc datasets noprint;
				delete ChiTest2 ChiTest3 ChiSqDataSplitT:;
			run;
		%end;
		/*Perform Modeling with adjustments for either continuous or categorical variables*/
		%If "%Scan(&Types., &I.)" = "adj" %Then %Do;
			%Logistics(data=&Data., Var=%Scan(&Vars., &I.) , by=&BY., label=, adj=%Scan(&adj., &I., '$'));
			%RobustPoisson(data=&Data., Var=%Scan(&Vars., &I.) , by=&BY., label=, adj=%Scan(&adj., &I., '$'));
			%LinearReg(data=&Data., Var=%Scan(&Vars., &I.) , by=&BY., label=, adj=%Scan(&adj., &I., '$'),transform=%Scan(&transforms., &I.));
			data Stats07;
				length ADJ odds oddsinv RRs $50.;
				merge PE2 Estimates2 EstimatesLinReg;
				ADJ = "%Scan(&adj., &I., '$')";
				pvalue = put(P_VALUE_LOGIT, pvalue6.4);
				pvaluepois = put(P_VALUE_Poisson, pvalue6.4);
				pvaluelinReg = put(P_VALUE_LinREg, pvalue6.4);
				odds = strip(put(OddsRatioEst, comma13.3))||" ("|| strip(put(OddsRatioEstL, comma13.3))||", "||strip(put(OddsRatioEstU, comma13.3))||")"; 
				oddsinv = strip(put(1/OddsRatioEst, comma13.3))||" ("|| strip(put(1/OddsRatioEstU, comma13.3))||", "||strip(put(1/OddsRatioEstL, comma13.3))||")"; 
				RRs = strip(put(RR, comma13.3))||" ("|| strip(put(RRL, comma13.3))||", "||strip(put(RRU, comma13.3))||")"; 
				MeanDIFFS = strip(put(DIFFMean, comma13.3))||" ("|| strip(put(LinRegLowerCI, comma13.3))||", "||strip(put(LinRegUpperCI, comma13.3))||")"; 
			run;
			proc datasets noprint;
				delete Estimates2 EstimatesLinReg PE2 ChiTest2 ChiTest3 ChiSqDataSplitT:;
			run;
		%end;
		
		Data Analysis&I.;
			length label  var $200. pvalue $20.;
			set Stats07(in=in7);
			Varnum = &i.;
			numvar = &num.;
			Var = "%Scan(&Vars., &I.)";
			Label = "%Scan(&labels., &I., '$')";
			Obs = &I.;
		run;

		proc datasets noprint;
			delete Stats07;
		run;
	%End;
	quit;
	ods select all;
	Data Analysis;
		length Var label $79.;
		set Analysis1-Analysis%sysfunc(countw(&vars.));
	run;
	proc datasets noprint;
		delete Analysis1-Analysis%sysfunc(countw(&vars.));
	run;

	/*Append results into an provided dataset*/
	%If "&AppendResults" ne "" %Then %Do;
		Data &AppendResults.;
			set &AppendResults. Analysis;
		run;
	%End;
%Mend Teststats;


/***************************************************************************
* CREATE UNADJUSTED P-VALUES FOR ONE PAIRWISE COMPARISON
***************************************************************************/

%macro pvalues(data=, label=);

    /*---------------------------------------------------------
      1. Create temporary pairwise dataset

      All required binary variables are already created in DF.
    ---------------------------------------------------------*/

    data temp_&label.;
        set &data.;
    run;


    /*---------------------------------------------------------
      2. Run requested unadjusted tests
    ---------------------------------------------------------*/

    data All_TESTS;
    run;

%Teststats(Data=temp_&label., Vars=gen,               Types=cat, by=group, num=13, AppendResults=All_TESTS);
%Teststats(Data=temp_&label., Vars=age,               Types=con, by=group, num=14, AppendResults=All_TESTS);
%Teststats(Data=temp_&label., Vars=duration_yrs,      Types=con, by=group, num=15, AppendResults=All_TESTS);
%Teststats(Data=temp_&label., Vars=duration_yrs_gt2,  Types=cat, by=group, num=16, AppendResults=All_TESTS);
%Teststats(Data=temp_&label., Vars=race_bin,          Types=cat, by=group, num=27, AppendResults=All_TESTS);
%Teststats(Data=temp_&label., Vars=upsit_pctl,        Types=con, by=group, num=38, AppendResults=All_TESTS);
%Teststats(Data=temp_&label., Vars=upsit_pctl_le15,   Types=cat, by=group, num=39, AppendResults=All_TESTS);
%Teststats(Data=temp_&label., Vars=lowput_ratio_cwm,  Types=con, by=group, num=49, AppendResults=All_TESTS);


    /*---------------------------------------------------------
      3. Select the appropriate unadjusted p-value
    ---------------------------------------------------------*/

    data all_tests_&label.;
        set All_TESTS;

        length pvalue_&label. 8;

        if not missing(P2_WIL) then
            pvalue_&label. = P2_WIL;

        else if not missing(P_KW) then
            pvalue_&label. = P_KW;

        else if not missing(P_PCHI) then
            pvalue_&label. = P_PCHI;

        if needFisher = 1 and not missing(XP2_FISH) then
            pvalue_&label. = XP2_FISH;

        format pvalue_&label. 8.4;


        /* Put categorical p-values on the requested category row */
        if lowcase(strip(var)) in (
            "race_bin",
            "gen",
            "mds1_cog_bin",
            "upsit_pctl_le15"
        ) then order = 1;
        else order = 0;

        main = numvar;

        keep var pvalue_&label. main order;

        if missing(numvar) then delete;
    run;

    proc sort data=all_tests_&label.;
        by main order;
    run;


    /*---------------------------------------------------------
      4. Cleanup
    ---------------------------------------------------------*/

    proc datasets library=work nolist;
        delete All_TESTS temp_&label.;
    quit;

%mend pvalues;


/***************************************************************************
* MERGE ALL PAIRWISE UNADJUSTED P-VALUES
***************************************************************************/

%macro make_pvalues(groups=);

    /*---------------------------------------------------------
      1. Merge numeric p-values and rename them
    ---------------------------------------------------------*/

    data pvalues_noadjust;
        merge all_tests_:;
        by main order;

        %let i = 1;

        %do %while(%length(%scan(&groups., &i.)));

            %let g = %scan(&groups., &i.);

            rename pvalue_cat_&g. = P_&g.;

            %let i = %eval(&i. + 1);

        %end;
    run;


    /*---------------------------------------------------------
      2. Convert numeric p-values to character
    ---------------------------------------------------------*/

    data pvalues_noadjust_char;

        set pvalues_noadjust(
            rename=(
                %let i = 1;

                %do %while(%length(%scan(&groups., &i.)));

                    %let g = %scan(&groups., &i.);

                    P_&g. = _P_&g.

                    %let i = %eval(&i. + 1);

                %end;
            )
        );

        length
            %let i = 1;

            %do %while(%length(%scan(&groups., &i.)));

                %let g = %scan(&groups., &i.);

                P_&g.

                %let i = %eval(&i. + 1);

            %end;
            $20
        ;

        %let i = 1;

        %do %while(%length(%scan(&groups., &i.)));

            %let g = %scan(&groups., &i.);

            P_&g. = put(_P_&g., pvalue8.4);

            %let i = %eval(&i. + 1);

        %end;

        drop _P_:;
    run;

%mend make_pvalues;


/* Automatic transformation macro */

/********************************************
 * Auto-transform + rank fallback GLM macro *
 *------------------------------------------*
 * - Checks residual normality (Shapiro–Wilk)
 *   for: raw, log, sqrt, square
 * - Order: raw -> log -> sqrt -> sq
 * - First with SW p > alpha is chosen
 * - If none pass: use rank-based model
 * - Final dataset: anova_&out
 *   - one row: Source="&groupvar"
 *   - ProbF_2 = p-value with superscript model markers
 ********************************************/

%macro check_transformations(
    data=,
    outcome=,
    class=,
    predictors=,
    groupvar=group,
    alpha=0.05,
    out=
);

  ods exclude all; 

  %local p_raw p_log p_sqrt p_sq chosen stars;
  %let p_raw  = .;
  %let p_log  = .;
  %let p_sqrt = .;
  %let p_sq   = .;
  %let chosen = ;

  /***************
   * RAW MODEL   *
   ***************/
  ods output ModelANOVA=__anova_raw;
  proc glm data=&data;
    class &class;
    model &outcome = &predictors / e3;
    output out=__res_raw r=residual;
  run; quit;

  ods output TestsForNormality=__norm_raw;
  proc univariate data=__res_raw normal;
    var residual;
  run;
  ods output close;

  data _null_;
    set __norm_raw;
    where Test = "Shapiro-Wilk";
    call symputx("p_raw", pValue);
  run;

  /****************
   * LOG MODEL    *
   ****************/
  data __log_data;
    set &data;
    if &outcome > 0 then log_&outcome = log(&outcome);
    else log_&outcome = .;
  run;

  ods output ModelANOVA=__anova_log;
  proc glm data=__log_data;
    class &class;
    model log_&outcome = &predictors / e3;
    output out=__res_log r=residual;
  run; quit;

  ods output TestsForNormality=__norm_log;
  proc univariate data=__res_log normal;
    var residual;
  run;
  ods output close;

  data _null_;
    set __norm_log;
    where Test = "Shapiro-Wilk";
    call symputx("p_log", pValue);
  run;

  /*****************
   * SQRT MODEL    *
   *****************/
  data __sqrt_data;
    set &data;
    if &outcome >= 0 then sqrt_&outcome = sqrt(&outcome);
    else sqrt_&outcome = .;
  run;

  ods output ModelANOVA=__anova_sqrt;
  proc glm data=__sqrt_data;
    class &class;
    model sqrt_&outcome = &predictors / e3;
    output out=__res_sqrt r=residual;
  run; quit;

  ods output TestsForNormality=__norm_sqrt;
  proc univariate data=__res_sqrt normal;
    var residual;
  run;
  ods output close;

  data _null_;
    set __norm_sqrt;
    where Test = "Shapiro-Wilk";
    call symputx("p_sqrt", pValue);
  run;

  /*****************
   * SQUARE MODEL  *
   *****************/
  data __sq_data;
    set &data;
    sq_&outcome = (&outcome)**2;
  run;

  ods output ModelANOVA=__anova_sq;
  proc glm data=__sq_data;
    class &class;
    model sq_&outcome = &predictors / e3;
    output out=__res_sq r=residual;
  run; quit;

  ods output TestsForNormality=__norm_sq;
  proc univariate data=__res_sq normal;
    var residual;
  run;
  ods output close;

  data _null_;
    set __norm_sq;
    where Test = "Shapiro-Wilk";
    call symputx("p_sq", pValue);
  run;

  /****************************************
   * Decide best transform (in this order):
   * raw -> log -> sqrt -> sq
   ****************************************/
  %if %sysevalf(&p_raw  > &alpha) %then %let chosen = raw;
  %else %if %sysevalf(&p_log  > &alpha) %then %let chosen = log;
  %else %if %sysevalf(&p_sqrt > &alpha) %then %let chosen = sqrt;
  %else %if %sysevalf(&p_sq   > &alpha) %then %let chosen = sq;

  /****************************************
   * If none normal -> rank-based fallback
   ****************************************/
  %if "&chosen" = "" %then %do;
    %let chosen = rank;

    data __rank_base;
      set &data;
    run;

    proc rank data=__rank_base out=__rank_data ties=mean;
      var &outcome;
      ranks r_&outcome;
    run;

    ods output ModelANOVA=__anova_rank;
    proc glm data=__rank_data;
      class &class;
      model r_&outcome = &predictors / e3;
    run; quit;

    %let stars = ^{super f};

    data anova_&out;
      set __anova_rank;
      where Source = "&groupvar";
      length ProbF_2 $20 Transform $8;
      Transform = "&chosen";
      ProbF_2   = cats(put(ProbF, pvalue8.4), " &stars");
    run;

  %end;
  %else %do;

    /* Map chosen transform -> superscript model marker */
    %if "&chosen" = "raw"  %then %let stars = ^{super b};
    %else %if "&chosen" = "log"   %then %let stars = ^{super c};
    %else %if "&chosen" = "sqrt"  %then %let stars = ^{super d};
    %else %if "&chosen" = "sq"    %then %let stars = ^{super e};

    data anova_&out;
      set __anova_&chosen;
      where Source = "&groupvar";
      length ProbF_2 $20 Transform $8;
      Transform = "&chosen";
      ProbF_2   = cats(put(ProbF, pvalue8.4), " &stars");
    run;

  %end;

  /***************
   * Clean up    *
   ***************/
  proc datasets lib=work nolist;
    delete __res_raw __norm_raw __anova_raw
           __log_data __res_log __norm_log __anova_log
           __sqrt_data __res_sqrt __norm_sqrt __anova_sqrt
           __sq_data __res_sq __norm_sq __anova_sq
           __rank_base __rank_data __anova_rank;
  quit;

  ods exclude none; 

%mend check_transformations;


/********************************************
 * Helper 1: add_main_order                 *
 * - Adds main/order                        *
 * - Renames ProbF_2/ProbF -> P_<suffix>    *
 * - NEW: Deletes ds if it has no p-value   *
 ********************************************/
%macro add_main_order(ds, mainval, suff);
  %local lib mem dsid rc srcvar hasProbF has_p;

  /* Parse libref if provided; default WORK */
  %if %index(&ds,.) %then %do;
    %let lib=%scan(&ds,1,.);
    %let mem=%scan(&ds,2,.);
  %end;
  %else %do;
    %let lib=work;
    %let mem=&ds;
  %end;

  /* If dataset doesn't exist, skip this pair */
  %if not %sysfunc(exist(&lib..&mem)) %then %do;
    %put NOTE: &lib..&mem does not exist - skipping add_main_order.;
    %return;
  %end;

  /* Figure out which p-value variable exists (ProbF_2 preferred) and if ProbF exists */
  %let srcvar=;
  %let hasProbF=0;
  %let dsid=%sysfunc(open(&lib..&mem,i));
  %if &dsid %then %do;
    %if %sysfunc(varnum(&dsid,ProbF_2))>0 %then %let srcvar=ProbF_2;
    %else %if %sysfunc(varnum(&dsid,ProbF))>0 %then %let srcvar=ProbF;
    %let hasProbF=%sysfunc(varnum(&dsid,ProbF));
    %let rc=%sysfunc(close(&dsid));
  %end;

  /* If there is no ProbF variable at all, treat as "no p-value" and delete */
  %if &hasProbF = 0 %then %do;
    proc datasets lib=&lib nolist;
      delete &mem;
    quit;
    %return;
  %end;

  /* Check if there is at least one non-missing ProbF */
  %let has_p = 0;
  data _null_;
    set &lib..&mem end=last;
    retain has_p 0;
    if not missing(ProbF) then has_p = 1;
    if last then call symputx('has_p', has_p);
  run;

  /* If all ProbF are missing, delete this dataset and exit */
  %if &has_p = 0 %then %do;
    proc datasets lib=&lib nolist;
      delete &mem;
    quit;
    %return;
  %end;

  /* Stamp main/order now that we know this ds has a valid p-value */
  data &lib..&mem;
    set &lib..&mem;
    main  = &mainval.;
    order = 0;
  run;

  /* Rename p-value column to P_<suffix>, if present */
  %if %length(&srcvar) %then %do;
    proc datasets lib=&lib nolist;
      modify &mem;
        rename &srcvar = P_&suff;
    quit;
  %end;
  %else %put NOTE: &lib..&mem has no ProbF_2 or ProbF - rename skipped.;

%mend add_main_order;
/********************************************
 * Helper 2: batch wrapper                  *
 * - For each suffix in &pairs:             *
 *   calls add_main_order(prefix._suffix)   *
 ********************************************/
%macro add_many(prefix=, mainval=, pairs=);
  %local i suff;
  %do i = 1 %to %sysfunc(countw(&pairs));
    %let suff = %scan(&pairs,&i);
    %add_main_order(&prefix._&suff, &mainval, &suff);
  %end;
%mend add_many;


/***************************************************************************
* 1. CREATE ANALYSIS DATASET
***************************************************************************/

data combined_analysis_dataset;
    set analysis.combined_biomarker_with_clinical;

    if not missing(Lilly_ptau217p_date)
       or not missing(IU_Plasma_date)
       or not missing(C2N_AB_date);

    if visit_range = "Early (BL-V02)";
proc contents; 
run;


/***************************************************************************
* 2. ASSIGN GROUPS AND CREATE DERIVED VARIABLES
***************************************************************************/

data df;
    set combined_analysis_dataset;

    if COHORT in ("1", "2");


    /*---------------------------------------------------------
      Analysis groups
    ---------------------------------------------------------*/

    if subgroup_definition = "Healthy Control"
       and NSD = 0 then
        group = 1;

    else if cohort = "1"
        and subgroup_definition = "Sporadic PD"
        and NSD = 1 then
        group = 2;

    if missing(group) then delete;


    /*---------------------------------------------------------
      MoCA category
    ---------------------------------------------------------*/

    if not missing(moca) then do;

        if moca > 25 then
            moca_cat = 1;

        else if 21 <= moca <= 25 then
            moca_cat = 2;

        else
            moca_cat = 3;

    end;
    else moca_cat = .;

/*    /*---------------------------------------------------------*/
/*      Revised UPSIT percentile <= 15*/
/*    ---------------------------------------------------------*/
/*    if upsit_pctl <= 15 and not missing(upsit_pctl) then*/
/*        upsit_pctl_le15 = 1;*/
/*    else if upsit_pctl > 15 then*/
/*        upsit_pctl_le15 = 0;*/
/*    else if missing(upsit_pctl) then*/
/*        upsit_pctl_le15 = .;*/


    /*---------------------------------------------------------
      Race: White versus other race
    ---------------------------------------------------------*/

    if race = 1 then
        race_bin = 0;

    else if not missing(race) and race > 1 then
        race_bin = 1;

    else
        race_bin = .;


    /*---------------------------------------------------------
      Cognitive categorization
    ---------------------------------------------------------*/

    if cogstate = 1 then
        cogstate_bin = 0;

    else if cogstate > 1 then
        cogstate_bin = 1;

    if missing(cogstate) then
        cogstate_bin = .;


    /*---------------------------------------------------------
      NSD stage

      Per supplied instructions:
          stage 1, 2, and 2.5 = 1
          stage 3 through 6   = 0

      Stage 1.5 is not assigned to either category.
    ---------------------------------------------------------*/

    if stage_num in (1 1.5 2 2.5) then
        stage_bin = 1;

    else if stage_num in (3 4 5 6) then
        stage_bin = 0;

    if missing(stage_num) then
        stage_bin = .;


    /*---------------------------------------------------------
      MDS-UPDRS item 1.1
    ---------------------------------------------------------*/

    if mds1_cog > 0 then
        mds1_cog_bin = 1;

    else if not missing(mds1_cog) then
        mds1_cog_bin = 0;

    else
        mds1_cog_bin = .;


    /*---------------------------------------------------------
      MoCA binary category
    ---------------------------------------------------------*/

    if moca_cat in (1 2) then
        moca_bin = 1;

    else if not missing(moca_cat) and moca_cat = 3 then
        moca_bin = 0;

    else
        moca_bin = .;


    /*---------------------------------------------------------
      Hoehn and Yahr
    ---------------------------------------------------------*/

    if hy_orig = 0 then
        hy_bin = 0;

    else if hy_orig in (1 2 3 4 5) then
        hy_bin = 1;

    else
        hy_bin = .;


    /*---------------------------------------------------------
      Disease duration > 2 years
    ---------------------------------------------------------*/

    if duration_yrs > 2 then
        duration_yrs_gt2 = 1;

    else if duration_yrs <= 2 and not missing(duration_yrs) then
        duration_yrs_gt2 = 0;

    else
        duration_yrs_gt2 = .;

run;


proc freq data=df;
    tables group;
run;


/***************************************************************************
* 3. BUILD SUMMARY STATISTICS
***************************************************************************/

%macro build_summary(dataset=, out=);

    /* Biomarkers and cutoff groups */
    %cont2(data=&dataset., var=Lilly_ptau217p,      num=1,  dec=2, dec2=2);
    %cat  (data=&dataset., var=ptau217_cutoff,      num=2);

    %cont2(data=&dataset., var=IU_NFL_plasma_mean,  num=3,  dec=2, dec2=2);
    %cat  (data=&dataset., var=NFL_mean_cutoff,     num=4);

    %cont2(data=&dataset., var=IU_GFAP_plasma_mean, num=5,  dec=2, dec2=2);
    %cat  (data=&dataset., var=GFAP_mean_cutoff,    num=6);

    %cont2(data=&dataset., var=C2N_AB40,            num=7,  dec=2, dec2=2);
    %cat  (data=&dataset., var=C2N_AB40_cutoff,     num=8);

    %cont2(data=&dataset., var=C2N_AB42,            num=9,  dec=2, dec2=2);
    %cat  (data=&dataset., var=C2N_AB42_cutoff,     num=10);

    %cont2(data=&dataset., var=C2N_AB4240,          num=11, dec=2, dec2=2);
    %cat  (data=&dataset., var=C2N_AB4240_cutoff,   num=12);


    /* Demographic variables */
    %cat  (data=&dataset., var=gen,                  num=13);
    %cont2(data=&dataset., var=age,                  num=14, dec=1, dec2=1);
    %cont2(data=&dataset., var=duration_yrs,         num=15, dec=1, dec2=1);
    %cat  (data=&dataset., var=duration_yrs_gt2,     num=16);
    %cat  (data=&dataset., var=race,                 num=27);


    /* UPSIT */
    %cont2(data=&dataset., var=upsit,                num=37, dec=1, dec2=1);
    %cont2(data=&dataset., var=upsit_pctl,           num=38, dec=1, dec2=1);
    %cat  (data=&dataset., var=upsit_pctl_le15,      num=39);

    /* NSD stage */
    %cat(data=&dataset.(where=(NSD=1 or missing(NSD))), var=stage_num, num=40);

    /* Clinical variables */
    %cat  (data=&dataset., var=mds1_cog,             num=41);
    %cont2(data=&dataset., var=updrs1_score,         num=42, dec=1, dec2=0);
    %cont2(data=&dataset., var=updrs2_score,         num=43, dec=1, dec2=0);
    %cont2(data=&dataset., var=updrs3_score,         num=44, dec=1, dec2=0);
    %cont2(data=&dataset., var=updrs_totscore,       num=45, dec=1, dec2=0);
    %cont2(data=&dataset., var=moca,                 num=46, dec=1, dec2=0);
    %cont2(data=&dataset., var=scopa,                num=47, dec=1, dec2=0);
    %cat  (data=&dataset., var=hy_orig,              num=48);

    /* DAT variables */
    %cont2(data=&dataset., var=lowput_ratio_cwm,     num=49, dec=2, dec2=2);
    %cont2(data=&dataset., var=mia_striatum_ref_cwm, num=50, dec=2, dec2=2);
    %cont2(data=&dataset., var=mia_putamen_ref_cwm,  num=51, dec=2, dec2=2);

    /* Cognitive categorization */
    %cat(data=&dataset., var=cogstate, num=52);


    /* Compile results */
    data &out.;
        merge stats4_: freq3_:;
        by main order;
    run;

    proc print data=&out.;
    run;


    /* Cleanup temporary summary datasets */
    proc datasets library=work nolist;
        delete stats: freq:;
    quit;

%mend build_summary;


/* Run summary */
%build_summary(
    dataset=df,
    out=combine1_df
);


data combine_1;
    set combine1_df;
    by main order;
run;


/***************************************************************************
* 4. TABLE LABELS
***************************************************************************/

data labels;

    infile datalines truncover;

    input
        @1  main  best4.
        @6  order best4.
        @10 name  $80.
    ;

datalines;
1    0   ^{style[fontweight=bold]Ptau217}, mean (SD)
1    1   Median (min, max)
1    2   Missing
1.5  0   ^{style[fontweight=bold]Ptau217 tertiles}, n (%)
2    1   < 0.131
2    2   0.131 - 0.174
2    3   > 0.174
2    99  Missing
3    0   ^{style[fontweight=bold]NFL}, mean (SD)
3    1   Median (min, max)
3    2   Missing
3.5  0   ^{style[fontweight=bold]NFL tertiles}, n (%)
4    1   < 9.877
4    2   9.877 - 14.538
4    3   > 14.538
4    99  Missing
5    0   ^{style[fontweight=bold]GFAP}, mean (SD)
5    1   Median (min, max)
5    2   Missing
5.5  0   ^{style[fontweight=bold]GFAP tertiles}, n (%)
6    1   < 64.478
6    2   64.478 - 105.341
6    3   > 105.341
6    99  Missing
7    0   ^{style[fontweight=bold]Aß40}, mean (SD)
7    1   Median (min, max)
7    2   Missing
7.5  0   ^{style[fontweight=bold]Aß40 tertiles}, n (%)
8    1   < 374.886
8    2   374.886 - 422.785
8    3   > 422.785
8    99  Missing
9    0   ^{style[fontweight=bold] Aß42}, mean (SD)
9    1   Median (min, max)
9    2   Missing
9.5  0   ^{style[fontweight=bold]Aß42 tertiles}, n (%)
10   1   < 40.219
10   2   40.219 - 45.870
10   3   > 45.870
10   99  Missing
11   0   ^{style[fontweight=bold]Aß42/40}, mean (SD)
11   1   Median (min, max)
11   2   Missing
11.5 0   ^{style[fontweight=bold]Aß42/40 tertiles}, n (%)
12   1   < 0.104
12   2   0.104 - 0.111
12   3   > 0.111
12   99  Missing
12.5 0   ^{style[fontweight=bold]Sex}, n (%)
13   1   Male
13   2   Female
13   99  Missing
14   0   ^{style[fontweight=bold]Age (years)}, mean (SD)
14   1   Median (min, max)
14   2   Delete
15   0   ^{style[fontweight=bold]Disease Duration (years)}, mean (SD)
15   1   Median (min, max)
15   2   Delete
15.5 0   ^{style[fontweight=bold]Disease duration (years)}, n (%)
16   0   <= 2
16   1   > 2
16   99  Delete
26.5 0   ^{style[fontweight=bold]Race}, n (%)
27   1   White
27   2   Black/African American
27   3   Asian
27   4   Other
27   99  Missing
37   0   ^{style[fontweight=bold]UPSIT Raw Score}, mean (SD)
37   1   Median (min, max)
37   2   Missing
38   0   ^{style[fontweight=bold]Revised UPSIT percentile}, mean (SD)
38   1   Median (min, max)
38   2   Delete
39   0   Delete
39   1   ^{style[fontweight=bold]Revised UPSIT percentile ^{unicode 2264} 15}, n (%)
39   99  Delete
39.5 0   ^{style[fontweight=bold]NSD+}, n (%)
40   0   0
40   1   1a
40   1.5 1b
40   2   2a
40   2.5 2b
40   3   3
40   4   4
40   5   5
40   6   6
40   99  Missing
40.5 0   ^{style[fontweight=bold]MDS-UPDRS item 1.1}, n (%)
41   0   Normal
41   1   Slight
41   2   Mild
41   3   Moderate
41   4   Severe
41   99  Missing
42   0   ^{style[fontweight=bold]MDS-UPDRS part I}, mean (SD)
42   1   Median (min, max)
42   2   Missing
43   0   ^{style[fontweight=bold]MDS-UPDRS part II}, mean (SD)
43   1   Median (min, max)
43   2   Missing
44   0   ^{style[fontweight=bold]MDS-UPDRS part III}, mean (SD)
44   1   Median (min, max)
44   2   Missing
45   0   ^{style[fontweight=bold]MDS-UPDRS total score}, mean (SD)
45   1   Median (min, max)
45   2   Missing
46   0   ^{style[fontweight=bold]MoCA total score}, mean (SD)
46   1   Median (min, max)
46   2   Missing
47   0   ^{style[fontweight=bold]SCOPA-AUT total score}, mean (SD)
47   1   Median (min, max)
47   2   Missing
47.5 0   ^{style[fontweight=bold]Hoehn and Yahr}, n (%)
48   0   Stage 0
48   1   Stage 1
48   2   Stage 2
48   3   Stage 3
48   99  Missing
49   0   ^{style[fontweight=bold]% Age/sex-expected low putamen}, mean (SD)
49   1   Median (min, max)
49   2   Missing
50   0   ^{style[fontweight=bold]Mean striatum CWM}, mean (SD)
50   1   Median (min, max)
50   2   Missing
51   0   ^{style[fontweight=bold]Mean putamen CWM}, mean (SD)
51   1   Median (min, max)
51   2   Missing
51.5 0   ^{style[fontweight=bold]Cognitive Categorization}, n (%)
52   1   Normal cognition
52   2   Mild Cognitive Impairment
52   99  Missing
;
run;


/***************************************************************************
* 5. UNADJUSTED PAIRWISE P-VALUES
***************************************************************************/

%pvalues(
    data=df(where=(group in (1 2))),
    label=cat_12
);

%make_pvalues(
    groups=12 
);


proc print data=pvalues_noadjust_char;
run;


/***************************************************************************
* 6. CONTINUOUS ADJUSTED P-VALUES
***************************************************************************/

%macro run_pairs_for_one_var(outcome=, label=);

    %check_transformations(
        data=df(where=(group in (1 2))),
        outcome=&outcome.,
        class=group gen,
        predictors=group age gen,
        groupvar=group,
        out=&label._12
    );

%mend run_pairs_for_one_var;


/* Biomarkers */
%run_pairs_for_one_var(
    outcome=Lilly_ptau217p,
    label=ptau
);

%run_pairs_for_one_var(
    outcome=IU_NFL_plasma_mean,
    label=nfl
);

%run_pairs_for_one_var(
    outcome=IU_GFAP_plasma_mean,
    label=gfap
);

%run_pairs_for_one_var(
    outcome=C2N_AB40,
    label=AB40
);

%run_pairs_for_one_var(
    outcome=C2N_AB42,
    label=AB42
);

%run_pairs_for_one_var(
    outcome=C2N_AB4240,
    label=AB4240
);

/* Additional continuous variables */
%run_pairs_for_one_var(
    outcome=upsit,
    label=upsit
);

%run_pairs_for_one_var(
    outcome=updrs1_score,
    label=updrs1_score
);

%run_pairs_for_one_var(
    outcome=updrs2_score,
    label=updrs2_score
);

%run_pairs_for_one_var(
    outcome=updrs3_score,
    label=updrs3_score
);

%run_pairs_for_one_var(
    outcome=updrs_totscore,
    label=updrs_totscore
);

%run_pairs_for_one_var(
    outcome=moca,
    label=moca
);

%run_pairs_for_one_var(
    outcome=scopa,
    label=scopa
);

%run_pairs_for_one_var(
    outcome=mia_striatum_ref_cwm,
    label=mia_striatum_ref_cwm
);

%run_pairs_for_one_var(
    outcome=mia_putamen_ref_cwm,
    label=mia_putamen_ref_cwm
);


/***************************************************************************
* 7. ADD CONTINUOUS MODEL P-VALUES TO TABLE LOCATIONS
***************************************************************************/

%let pairs_all = 12;


/* Biomarkers */
%add_many(
    prefix=anova_ptau,
    mainval=1,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_nfl,
    mainval=3,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_gfap,
    mainval=5,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_ab40,
    mainval=7,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_ab42,
    mainval=9,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_ab4240,
    mainval=11,
    pairs=&pairs_all.
);


/* Additional continuous variables */
%add_many(
    prefix=anova_upsit,
    mainval=37,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_updrs1_score,
    mainval=42,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_updrs2_score,
    mainval=43,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_updrs3_score,
    mainval=44,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_updrs_totscore,
    mainval=45,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_moca,
    mainval=46,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_scopa,
    mainval=47,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_mia_striatum_ref_cwm,
    mainval=50,
    pairs=&pairs_all.
);

%add_many(
    prefix=anova_mia_putamen_ref_cwm,
    mainval=51,
    pairs=&pairs_all.
);


/* Merge continuous adjusted p-values */
data pvals_transform;
    merge anova_:;
    by main order;

    where Source = "group";

    keep P_: main order Source Transform;
proc print;
run;


/***************************************************************************
* 8. CATEGORICAL ADJUSTED LOGISTIC MODELS
***************************************************************************/

%macro run_logit_pairs(
    data=,
    var=,
    prefix=,
    main=,
    pairs=
);

    %local
        i
        pair
        g1
        g2
        linkopt
        has_p
    ;


    /* Cumulative logit for three-level cutoff variables */
    %if %upcase(&var.) = NFL_MEAN_CUTOFF
        or %upcase(&var.) = C2N_AB4240_CUTOFF
        or %upcase(&var.) = GFAP_MEAN_CUTOFF
        or %upcase(&var.) = C2N_AB40_CUTOFF
        or %upcase(&var.) = C2N_AB42_CUTOFF
        or %upcase(&var.) = PTAU217_CUTOFF
    %then
        %let linkopt = / link=CLogit;
    %else
        %let linkopt = ;


    %do i = 1 %to %sysfunc(countw(&pairs.));

        %let pair = %scan(&pairs., &i.);


        /* Extract the two comparison groups */
        %if %index(&pair., _) %then %do;

            %let g1 = %scan(&pair., 1, _);
            %let g2 = %scan(&pair., 2, _);

        %end;
        %else %do;

            %let g1 = %substr(&pair., 1, 1);
            %let g2 = %substr(&pair., 2);

        %end;


        ods exclude all;
        ods noresults;

        ods output ModelANOVA=_raw_&prefix._&pair.;

        proc logistic data=&data.;
            where group in (&g1., &g2.);

            class
                group(ref="&g1.")
                gen(ref="1")
                / param=ref
            ;

            model &var. = group age gen &linkopt.;
        run;

        ods results;
        ods exclude none;


        /* Retain only the group-effect p-value */
        %let has_p = 0;

        %if %sysfunc(exist(_raw_&prefix._&pair.)) %then %do;

            data &prefix._&pair.;
                set _raw_&prefix._&pair.;

                where upcase(Effect) = "GROUP"
                      and not missing(ProbChiSq);

                main  = &main.;
                order = 0;

                length P_&pair. $20;

                P_&pair. = cats(
                    put(ProbChiSq, pvalue8.4),
                    " ^{super a}"
                );

                call symputx("has_p", 1, "L");

                keep main order P_&pair.;
            run;

        %end;


        /* Delete empty output */
        %if &has_p. = 0 %then %do;

            %if %sysfunc(exist(&prefix._&pair.)) %then %do;

                proc datasets library=work nolist;
                    delete &prefix._&pair.;
                quit;

            %end;

        %end;


        proc datasets library=work nolist;
            delete _raw_&prefix._&pair.;
        quit;

    %end;

%mend run_logit_pairs;


%macro run_all_logit_sets(var=, prefix=, main=);

    %run_logit_pairs(
        data=df,
        var=&var.,
        prefix=&prefix.,
        main=&main.,
        pairs=12 13 23
    );

%mend run_all_logit_sets;


/* Biomarker cutoff variables */
%run_all_logit_sets(
    var=ptau217_cutoff,
    prefix=cutoffs_anova_ptau,
    main=1.5
);

%run_all_logit_sets(
    var=NFL_mean_cutoff,
    prefix=cutoffs_anova_nfl,
    main=3.5
);

%run_all_logit_sets(
    var=gfap_mean_cutoff,
    prefix=cutoffs_anova_gfap,
    main=5.5
);

%run_all_logit_sets(
    var=C2N_AB40_cutoff,
    prefix=cutoffs_anova_ab40,
    main=7.5
);

%run_all_logit_sets(
    var=C2N_AB42_cutoff,
    prefix=cutoffs_anova_ab42,
    main=9.5
);

%run_all_logit_sets(
    var=C2N_AB4240_cutoff,
    prefix=cutoffs_anova_ab4240,
    main=11.5
);


/* Requested binary clinical variables */
%run_all_logit_sets(
    var=stage_bin,
    prefix=bin_stage,
    main=39.5
);

%run_all_logit_sets(
    var=mds1_cog_bin,
    prefix=bin_mds1,
    main=40.5
);

%run_all_logit_sets(
    var=hy_bin,
    prefix=bin_hy,
    main=47.5
);

%run_all_logit_sets(
    var=cogstate_bin,
    prefix=bin_cog,
    main=51.5
);


/* Merge categorical adjusted p-values */
data pvals_logistic;
    merge cutoffs_anova_:
          bin_:;
    by main order;
run;


proc print data=pvals_logistic;
run;


/***************************************************************************
* 9. COMBINE ADJUSTED P-VALUES
***************************************************************************/

data pvals_adjusted_all;
    set pvals_transform
        pvals_logistic;
run;

proc sort data=pvals_adjusted_all;
    by main order;
run;

proc print data=pvals_adjusted_all;
run;


/***************************************************************************
* 10. COMBINE ADJUSTED AND UNADJUSTED P-VALUES
***************************************************************************/

data pvals_all(drop=var j);

    merge pvals_adjusted_all
          pvalues_noadjust_char;

    by main order;

    array pvals(*) $ P_:;

    do j = 1 to dim(pvals);

        if strip(pvals[j]) = "" then
            pvals[j] = "NA";

        else if strip(scan(pvals[j], 1, "^")) = "." then
            pvals[j] = "NA";

    end;

run;


proc print data=pvals_all;
run;


/***************************************************************************
* 11. GROUP COUNTS
***************************************************************************/

proc sql noprint;

    select count(*)
        into :g1 - :g3
    from df
    group by group
    order by group;

quit;


/***************************************************************************
* 12. COMBINE LABELS, SUMMARY STATISTICS, AND P-VALUES
***************************************************************************/

data combine(
    drop=i j transform source all_na0
);

    merge labels(in=a)
          combine_1
          pvals_all;

    by main order;

    if a;

    if name = "Delete" then delete;


    /* All character summary columns beginning with underscore */
    array col(*) $ _:;


    do i = 1 to dim(col);

        /* Convert zero percentage displays to simple zero */
        if col[i] = "0 (0%)" then
            col[i] = "0";

        else if col[i] = "0 (.%)" then
            col[i] = "0";


        /* Assign NA only to non-heading rows */
        else if strip(col[i]) = ""
            and mod(int(main * 2), 2) = 0 then
            col[i] = "NA";

    end;


    /* Delete rows where every summary column is NA or zero */
    all_na0 = 1;

    do i = 1 to dim(col);

        if upcase(compress(strip(col[i]), ". ")) not in ("NA", "0") then do;

            all_na0 = 0;
            leave;

        end;

    end;

    if all_na0 then delete;

run;


proc sort data=combine;
    by main order;
run;

proc print data=combine;
run;


/* Split the final combined dataset */
data table1_clinical table2_biomarkers;
    set combine;

    /* Sex onwards */
    if main >= 12.5 then output table1_clinical;

    /* Biomarker statistics and cutoffs */
    else if main >= 1 and main <= 12 then
        output table2_biomarkers;
run;


options
    orientation=landscape
    nodate
    nonumber
    missing=" "
;

ods noptitle;
ods graphics / border=off;
ods escapechar="^";


/* Reuse the same report layout for both tables */
%macro display_table(dataset=, table_number=, table_title=);

    title;
    footnote;

    title
        bold
        h=9pt
        j=left
        f="Arial"
        "TABLE &table_number.: &table_title."
    ;

    /* Common model footnotes */
    footnote1 h=9pt j=left f="Arial"
    "^{super a} P-values calculated using PROC LOGISTIC adjusted for age and sex.^n
    ^{super b} P-values calculated using GLM adjusted for age and sex.^n
    ^{super c} P-values calculated using a log-transformed model.^n
    ^{super d} P-values calculated using a square-root-transformed model.^n
    ^{super e} P-values calculated using a squared-value (power-2) transformation.^n
    ^{super f} P-values calculated using a rank-based nonparametric model.";

    footnote2 h=9pt j=left f="Arial"
    "No superscript markers indicate unadjusted p-values. Kruskal-Wallis and chi-square/Fisher's exact tests were used for overall comparisons. Wilcoxon rank-sum and chi-square/Fisher's exact tests were used for pairwise comparisons.";

    /* Footnotes specific to each table */
    %if &table_number. = 1 %then %do;

        footnote3 h=9pt j=left f="Arial"
        "Race p-value is a comparison of non-white vs. white.^n
        For MDS-UPDRS item 1.1, all groups except 'Normal' were collapsed and the resulting p-value is a comparison of normal vs. any pain between two groups.^n
        For NSD Stage, Stages <2B are collapsed and compared against Stages 3<. For Hoehn and Yahr, Stage 0 is compared against Stages >0.";

    %end;
    %else %if &table_number. = 2 %then %do;

        footnote3 h=9pt j=left f="Arial"
        "P-values for P-tau, NFL, GFAP, and A^{unicode 03B2} cutoff variables were calculated using a cumulative logit link function.";

    %end;

    footnote5 h=9pt j=left f="Arial"
    "Report Generated on Data Submitted as of: &report_date..";


    proc report
        data=&dataset.
        spanrows
        split="\"

        style(column)={
            font_size=8pt
            fontfamily=Arial
        }

        style(header)={
            font_size=9pt
            fontfamily=Arial
        }

        style(report)={
            just=left
        }
    ;

        column
            main
            order
            name

            ("Healthy Control"
                _1
            )

            ("Sporadic PD"
                _2
            )

            ("Pairwise p-values"
                P_12
            )
        ;


        define main /
            display
            noprint
        ;

        define order /
            display
            noprint
        ;

        define name /
            ""
            style={
                width=1.75in
            }
        ;


        define _1 /
            "NSD- \^{style[fontsize=9pt fontweight=medium](N = &g1.)}"
            style={
                width=1.10in
                borderleftcolor=DarkSlateGrey
            }
            center
        ;

        define _2 /
            "NSD+ \^{style[fontsize=9pt fontweight=medium](N = &g2.)}"
            style={
                width=1.10in
                borderleftcolor=DarkSlateGrey
            }
            center
        ;


        define P_12 /
            "HC vs sPD"
            style={
                width=1.10in
                borderleftcolor=DarkSlateGrey
            }
            center
        ;


        /* Horizontal line before each new MAIN group */
        break before main /
            style=[
                bordertopcolor=DarkSlateGrey
                bordertopwidth=1
                bordertopstyle=solid
            ]
        ;


        compute name;

            if index(name, "^{style") = 0 then
                call define(
                    _col_,
                    "style",
                    "style=[paddingleft=1em]"
                );

            else
                call define(
                    _row_,
                    "style",
                    "style=[bordertopcolor=DarkSlateGrey]"
                );

        endcomp;


        /* Highlight statistically significant pairwise p-values */
        compute P_12;

            length _numtxt $20;

            _numtxt = compress(
                scan(P_12, 1, "^"),
                "0123456789.+-Ee",
                "k"
            );

            _pnum = input(_numtxt, ?? best32.);

            if not missing(_pnum) and _pnum < 0.05 then
                call define(
                    _col_,
                    "style",
                    "style=[foreground=red fontweight=bold]"
                );

        endcomp;

    run;

    title;
    footnote;

%mend display_table;


/* Table 1: Clinical onwards */
%display_table(
    dataset=table1_clinical,
    table_number=1,
    table_title=DEMOGRAPHIC AND CLINICAL CHARACTERISTICS
);


/* Table 2: Biomarker statistics and cutoffs */
%display_table(
    dataset=table2_biomarkers,
    table_number=2,
    table_title=PLASMA BIOMARKER CHARACTERISTICS
);



/***************************************************************************
* 14. OPTIONAL WORK-LIBRARY CLEANUP
***************************************************************************/

proc datasets library=work kill nolist;
quit;






/*==========================================================
  01_Table2_Pvalues_Unadjusted_Combined.sas
  Unadjusted p-values for Combined Biomarker Table 2

  ONLY unadjusted variables:
    gen
    age
    duration_yrs
    race_bin
    ethn_hispanic
    upsit_pctl
    upsit_pctl_le15
    lowput_ratio_cwm
==========================================================*/

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


%macro pvalues(data=, label=);

  data temp_&label;
    set &data;

    /* race_bin: White vs non-White */
    if race = 1 then race_bin = 0;
    else if not missing(race) and race > 1 then race_bin = 1;
    else race_bin = .;

  run;


  data All_TESTS;
  run;


  /* Mains aligned to current reduced Table 2 layout */
  %Teststats(Data=temp_&label, Vars=gen,               Types=cat, by=group, num=7,  AppendResults=All_TESTS);
  %Teststats(Data=temp_&label, Vars=age,               Types=con, by=group, num=8,  AppendResults=All_TESTS);
  %Teststats(Data=temp_&label, Vars=duration_yrs,      Types=con, by=group, num=9,  AppendResults=All_TESTS);
  %Teststats(Data=temp_&label, Vars=race_bin,          Types=cat, by=group, num=10, AppendResults=All_TESTS);
  %Teststats(Data=temp_&label, Vars=ethn_hispanic,     Types=cat, by=group, num=11, AppendResults=All_TESTS);
  %Teststats(Data=temp_&label, Vars=upsit_pctl,        Types=con, by=group, num=13, AppendResults=All_TESTS);
  %Teststats(Data=temp_&label, Vars=upsit_pctl_le15,   Types=cat, by=group, num=14, AppendResults=All_TESTS);
  %Teststats(Data=temp_&label, Vars=lowput_ratio_cwm,  Types=con, by=group, num=23, AppendResults=All_TESTS);


  data all_tests_&label;
    set All_TESTS;

    if not missing(P2_WIL)                   then pvalue_&label = P2_WIL;
    if missing(P2_WIL) and not missing(P_KW) then pvalue_&label = P_KW;
    else if not missing(P_PCHI)              then pvalue_&label = P_PCHI;
    if needFisher = 1                        then pvalue_&label = XP2_FISH;

    format pvalue_&label 8.4;

    /*
      Header-row placement:
        gen             num=7  -> main=6.5
        race_bin        num=10 -> main=9.5
        ethn_hispanic   num=11 -> main=10.5
        upsit_pctl_le15 num=14 -> main=14, order=1
    */

    if var = "gen" then main = 6.5;
    else if var = "race_bin" then main = 9.5;
    else if var = "ethn_hispanic" then main = 10.5;
    else main = numvar;

    if var = "upsit_pctl_le15" then order = 1;
    else order = 0;

    keep var pvalue_&label main order;

    if missing(numvar) then delete;

    proc sort;
      by main order;
  run;


  proc datasets lib=work nolist;
    delete All_TESTS temp_&label;
  quit;

%mend pvalues;



%macro make_pvalues(groups=);

  data pvalues_noadjust;
    merge all_tests_:;
    by main order;

    %let i = 1;
    %do %while(%length(%scan(&groups,&i)));
      %let g = %scan(&groups,&i);
      rename pvalue_cat_&g = P_&g;
      %let i = %eval(&i+1);
    %end;
  run;


  data pvalues_noadjust_char;
    set pvalues_noadjust(rename=(
      %let i = 1;
      %do %while(%length(%scan(&groups,&i)));
        %let g = %scan(&groups,&i);
        P_&g = _P_&g
        %let i = %eval(&i+1);
      %end;
    ));

    length
      %let i = 1;
      %do %while(%length(%scan(&groups,&i)));
        %let g = %scan(&groups,&i);
        P_&g
        %let i = %eval(&i+1);
      %end;
      $20;

    %let i = 1;
    %do %while(%length(%scan(&groups,&i)));
      %let g = %scan(&groups,&i);

      if missing(_P_&g) then P_&g = "";
      else if _P_&g < 0.0001 then P_&g = "<0.0001";
      else P_&g = put(_P_&g, pvalue8.4);

      %let i = %eval(&i+1);
    %end;

    drop _P_:;
  run;

%mend make_pvalues;


/*==========================================================
  Unadjusted continuous biomarker model
  For continuous outcomes where we want Beta/CI and biomarker p-value
  without adjustment for age and sex

  Formatting:
    >999.999 shown as >999.999
    positive values <0.001 shown as <0.001
==========================================================*/

%macro auto_cont_biomarker_assoc_unadj(
    data=,
    outcome=,
    cont_predictor=,
    alpha=0.05,
    main=,
    order=0,
    out=
);

    ods exclude all;

    ods output ParameterEstimates=__pe_unadj;

    proc glm data=&data;
        model &outcome = &cont_predictor / solution clparm;
    run;
    quit;

    data &out;
        length P_123 P_12 P_13 P_23 overall_cont $30
               type $20 transform $12 var $64 EST_CI $60
               est_fmt lcl_fmt ucl_fmt $20;

        set __pe_unadj;
        where strip(upcase(Parameter)) = strip(upcase("&cont_predictor"));

        var="&outcome";
        type="unadj_cont";
        transform="none";
        main=&main;
        order=&order;

        overall_cont_n = Probt;

        /* Format beta estimate and CI limits */
        if missing(Estimate) then est_fmt = "";
        else if Estimate > 999.999 then est_fmt = ">999.999";
        else if 0 < Estimate < 0.001 then est_fmt = "<0.001";
        else est_fmt = strip(put(Estimate, 8.4));

        if missing(LowerCL) then lcl_fmt = "";
        else if LowerCL > 999.999 then lcl_fmt = ">999.999";
        else if 0 < LowerCL < 0.001 then lcl_fmt = "<0.001";
        else lcl_fmt = strip(put(LowerCL, 8.4));

        if missing(UpperCL) then ucl_fmt = "";
        else if UpperCL > 999.999 then ucl_fmt = ">999.999";
        else if 0 < UpperCL < 0.001 then ucl_fmt = "<0.001";
        else ucl_fmt = strip(put(UpperCL, 8.4));

        EST_CI = cats(est_fmt, " (", lcl_fmt, ", ", ucl_fmt, ")");

        P_123="";
        P_12="";
        P_13="";
        P_23="";

        P_123_n=.;
        P_12_n=.;
        P_13_n=.;
        P_23_n=.;

        if missing(overall_cont_n) then overall_cont="";
        else if overall_cont_n < 0.0001 then overall_cont="<0.0001";
        else overall_cont=put(overall_cont_n, pvalue8.4);

        keep var main order type transform EST_CI overall_cont
             P_123 P_12 P_13 P_23
             P_123_n P_12_n P_13_n P_23_n overall_cont_n;
    run;

    proc datasets lib=work nolist;
        delete __pe_unadj;
    quit;

    ods exclude none;

%mend auto_cont_biomarker_assoc_unadj;


/*==========================================================
  Unadjusted logistic continuous biomarker model
  For binary outcomes where we want OR/CI and biomarker p-value
  without adjustment for age and sex

  Formatting:
    >999.999 shown as >999.999
    positive values <0.001 shown as <0.001
==========================================================*/

%macro auto_cat_pvals_unadj_cont(
    data=,
    outcome=,
    event_ref=1,
    cont_predictor=,
    main=,
    order=0,
    out=
);

    ods exclude all;

    ods output ParameterEstimates=__pe_log_unadj;

    proc logistic data=&data order=internal;
        model &outcome(event="&event_ref") = &cont_predictor;
    run;
    quit;

    data &out;
        length P_123 P_12 P_13 P_23 overall_cont $30
               type $20 transform $12 var $64 EST_CI $60
               or_fmt lcl_fmt ucl_fmt $20;

        set __pe_log_unadj;
        where strip(upcase(Variable)) = strip(upcase("&cont_predictor"));

        var="&outcome";
        type="unadj_cat_cont";
        transform="logistic";
        main=&main;
        order=&order;

        overall_cont_n = ProbChiSq;

        OR_n      = exp(Estimate);
        LowerOR_n = exp(Estimate - 1.96*StdErr);
        UpperOR_n = exp(Estimate + 1.96*StdErr);

        /* Format OR and CI limits */
        if missing(OR_n) then or_fmt = "";
        else if OR_n > 999.999 then or_fmt = ">999.999";
        else if 0 < OR_n < 0.001 then or_fmt = "<0.001";
        else or_fmt = strip(put(OR_n, 8.4));

        if missing(LowerOR_n) then lcl_fmt = "";
        else if LowerOR_n > 999.999 then lcl_fmt = ">999.999";
        else if 0 < LowerOR_n < 0.001 then lcl_fmt = "<0.001";
        else lcl_fmt = strip(put(LowerOR_n, 8.4));

        if missing(UpperOR_n) then ucl_fmt = "";
        else if UpperOR_n > 999.999 then ucl_fmt = ">999.999";
        else if 0 < UpperOR_n < 0.001 then ucl_fmt = "<0.001";
        else ucl_fmt = strip(put(UpperOR_n, 8.4));

        EST_CI = cats(or_fmt, " (", lcl_fmt, ", ", ucl_fmt, ")");

        P_123="";
        P_12="";
        P_13="";
        P_23="";

        P_123_n=.;
        P_12_n=.;
        P_13_n=.;
        P_23_n=.;

        if missing(overall_cont_n) then overall_cont="";
        else if overall_cont_n < 0.0001 then overall_cont="<0.0001";
        else overall_cont=put(overall_cont_n, pvalue8.4);

        keep var main order type transform EST_CI overall_cont
             P_123 P_12 P_13 P_23
             P_123_n P_12_n P_13_n P_23_n overall_cont_n;
    run;

    proc datasets lib=work nolist;
        delete __pe_log_unadj;
    quit;

    ods exclude none;

%mend auto_cat_pvals_unadj_cont;

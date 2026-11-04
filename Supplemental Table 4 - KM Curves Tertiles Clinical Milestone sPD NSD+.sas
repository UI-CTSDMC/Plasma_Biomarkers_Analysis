/*=========================================================
  sPD NSD+: Biomarker Tertiles and Clinical Milestones

  Uses the prepared clinical_analysis_v3 dataset to examine
  associations between plasma biomarker tertiles and time
  to clinical milestones in sPD NSD+ participants.

  Generates Kaplan-Meier curves and survival summary tables
  using age- and sex-adjusted Cox proportional hazards or
  Weibull accelerated failure time models, selected according
  to the proportional hazards assumption test.

=========================================================*/


%let filepath = /*path where analysis data are located*/;
libname analysis "&filepath.";


%let report_date = 08DEC2025;

ods escapechar="^";

/* Make a local copy of the prepared analysis dataset. */
data clinical_analysis_v3;
    set analysis.clinical_analysis_v3;
run;
%let write_individual_rtf = 0;  /* Preserve suppression of individual summary reports. */
%let write_individual_pdf = 1;  /* Preserve individual analysis titles. */

/*---------------------------------------------------------
  Tertile formats
---------------------------------------------------------*/

proc format;

    value ptau
        1 = "Tertile 1"
        2 = "Tertile 2"
        3 = "Tertile 3";

    value nfl
        1 = "Tertile 1"
        2 = "Tertile 2"
        3 = "Tertile 3";

    value gfap
        1 = "Tertile 1"
        2 = "Tertile 2"
        3 = "Tertile 3";

    value ab40f
        1 = "Tertile 1"
        2 = "Tertile 2"
        3 = "Tertile 3";

    value ab42f
        1 = "Tertile 1"
        2 = "Tertile 2"
        3 = "Tertile 3";

    value ab4240f
        1 = "Tertile 1"
        2 = "Tertile 2"
        3 = "Tertile 3";

run;


/*===============================================================
  1. BASE DATASETS
===============================================================*/

/*===============================================================
  2. BASELONG DATASET
===============================================================*/


/*===============================================================
  3. INITIALIZE COMBINED CLINICAL SUMMARY DATASET
===============================================================*/

proc datasets library=work nolist;
    delete Clinical_Tertile_Summary;
quit;

data Clinical_Tertile_Summary;
    length
        Biomarker $100
        Milestone $100
        Method $10
        Parameter $80
        ClassVal $40
        Comparison_Type $35
        Comparison $12
        MilestoneOrder 8
        SortOrder 8
        GroupDF 8
        WaldStat 8
        Beta 8
        SE 8
        Z_stat 8
        P_value 8
        HR 8
        HR_LCL 8
        HR_UCL 8
        Beta_CI $40
        HR_CI $50
        Beta_CI_4 $70
        SE_4 $20
        Z_4 $20
        P_4 $20
        HR_CI_4 $70
        P_Sig 8
    ;
    stop;
run;


/*===============================================================
  4. MACRO TO MAKE ONE CLINICAL SURVIVAL ENDPOINT
     This mirrors your validated clinical endpoint construction.
===============================================================*/

%macro MakeClinicalEndpoint(
    data=,
    event=,
    time=,
    visitvar=Event_ID,
    extra=
);

    proc sort data=&data. out=Survival1;
        where &event. = 1;
        by Patno &time.;
    run;

    proc sort data=Survival1 nodupkey;
        by Patno;
    run;

    proc sort data=&data. out=Survival0;
        where &event. = 0 &extra.;
        by Patno descending &time.;
    run;

    proc sort data=Survival0 nodupkey;
        by Patno;
    run;

    data Survival0b;
        merge Survival0(in=in1)
              Survival1(in=in2 keep=Patno);
        by Patno;

        if in2 then delete;
    run;

    data &event.;
        merge Survival1(in=in1)
              Survival0b(in=in2);
        by Patno;

        if in1 or in2;

        keep Patno &visitvar. &event. &time.;
    run;

    proc datasets library=work nolist;
        delete Survival1 Survival0 Survival0b;
    quit;

%mend MakeClinicalEndpoint;


/*---------------------------------------------------------
  4B. Tertile categorical survival macro

  Key correction:
      - Temp deletes missing group/time/outcome values
      - Temp deletes baseline events where outcome=1 and time=0
      - PH assumption check runs on Temp
      - KM/log-rank/modeling also run on Temp
      - AUTO chooses:
            PHREG   if no ASSESS PH p-value < 0.05
            LIFEREG if any ASSESS PH p-value < 0.05
      - ASSESS PH uses fixed seed for reproducibility
---------------------------------------------------------*/

%macro SurvTertile(
    data=,
    outcome=,
    group=,
    time=,
    label=,
    milestone_order=,
    method=AUTO,
    covars=age gen,
    grouplabel=,
    groupfmt=,
    title_suffix=for Sporadic PD,
    baseline_vars=,
    stack_group=,
    stack_groupfmt=,
    stack_grouplabel=,
    biomarker_label=
);

    %local pval pval2 pval2n dsid hasLCLM rc;
    %local p12 p13 p23;
    %local _modelLabel _dslabel;
    %local _ph_fail _method _ph_pvar;
    %local sigma var_sigma;
    %local cov_g1_sigma cov_g2_sigma;
    %local var_b1 var_b2 cov_b1b2;

    %let sigma=.;
    %let var_sigma=.;
    %let cov_g1_sigma=.;
    %let cov_g2_sigma=.;
    %let var_b1=.;
    %let var_b2=.;
    %let cov_b1b2=.;

    %if %superq(stack_group)=      %then %let stack_group=&group.;
    %if %superq(stack_groupfmt)=   %then %let stack_groupfmt=&groupfmt.;
    %if %superq(stack_grouplabel)= %then %let stack_grouplabel=&grouplabel.;
    %if %superq(baseline_vars)=    %then %let baseline_vars=Patno &group. &outcome. &time.;

    /* Clean label into valid SAS dataset name suffix */
    %let _dslabel = %sysfunc(compress(&label.,
        abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789, k));

    %let _dslabel = %substr(&_dslabel., 1,
        %sysfunc(min(29, %length(&_dslabel.))));


    /*-------------------------------------------------
      1. Create modeling dataset

      Temp is used for KM/log-rank/PHREG/LIFEREG modeling.
      Baseline events are removed here.

      PH assumption check below does NOT use Temp.
      It uses &data. directly, like your manual code.
    -------------------------------------------------*/

    data Temp;
        set &data.;
        if &outcome. = 1 and &time. = 0 then delete;
    run;


/*-------------------------------------------------
  2. Automatic PH testing

  Uses Temp, the exact same dataset used for
  KM/log-rank/PHREG/LIFEREG modeling.

  Rule:
      If any ASSESS PH p-value < 0.05, choose LIFEREG.
      Otherwise, choose PHREG.

  Defensive update:
      Standardizes the ASSESS PH p-value variable
      before applying the decision rule.
-------------------------------------------------*/

%let _ph_fail = 0;
%let _method  = %upcase(&method.);
%let _ph_pvar = ;

%if %upcase(&method.) = AUTO %then %do;

    %put NOTE: method=AUTO - running ASSESS PH check on TEMP data for &label. and &biomarker_label.;

    ods exclude all;
    ods output ProportionalHazardsSupTest = _ph_assess_raw_;

    proc phreg data=Temp  namelen=100;
        class &group.(ref='3') gen(ref='1') / param=ref;
        model &time.*&outcome.(0) = &group. &covars.;
        assess ph / resample seed=1472160691;
    run;

    ods exclude none;

    %if %sysfunc(exist(_ph_assess_raw_)) %then %do;

/*        title "DEBUG: Contents of ASSESS PH Output - &biomarker_label. - &label.";*/
/*        proc contents data=_ph_assess_raw_;*/
/*        run;*/
/**/
/*        title "DEBUG: Raw ASSESS PH Output - &biomarker_label. - &label.";*/
/*        proc print data=_ph_assess_raw_;*/
/*        run;*/
/*        title;*/

        /* Find the p-value variable name defensively */
        proc sql noprint;
            select name into :_ph_pvar trimmed
            from dictionary.columns
            where libname = "WORK"
              and memname = "_PH_ASSESS_RAW_"
              and upcase(name) in ("PVALUE", "PROBCHISQ", "PRCHISQ", "PROBZ", "PRCHISQ");
        quit;

        %if %superq(_ph_pvar)= %then %do;
            %put WARNING: Could not identify p-value variable in _ph_assess_raw_ for &biomarker_label. - &label.. Defaulting to PHREG.;
            %let _ph_fail = 0;
        %end;
        %else %do;

            %put NOTE: Using &_ph_pvar. as ASSESS PH p-value variable for &biomarker_label. - &label.;

            data _ph_assess_;
                set _ph_assess_raw_;

                length _Effect_for_Check $80;
                PH_PValue = .;
                _Effect_for_Check = "";

                if not missing(Variable) then _Effect_for_Check = strip(Variable);
                else if not missing(Effect) then _Effect_for_Check = strip(Effect);
                else if not missing(Parameter) then _Effect_for_Check = strip(Parameter);

                PH_PValue = &_ph_pvar.;

                keep _Effect_for_Check PH_PValue;
            run;

/*            title "DEBUG: Standardized ASSESS PH Output - &biomarker_label. - &label.";*/
/*            proc print data=_ph_assess_;*/
/*            run;*/
/*            title;*/

            data _null_;
                set _ph_assess_ end=last;

                retain ph_fail 0;

                if not missing(PH_PValue) and PH_PValue < 0.05 then ph_fail = 1;

                if last then call symputx('_ph_fail', ph_fail, 'L');
            run;

        %end;

    %end;
    %else %do;

        %put WARNING: _ph_assess_raw_ was not created for &biomarker_label. - &label.. Defaulting to PHREG.;
        %let _ph_fail = 0;

    %end;

    %if %sysevalf(&_ph_fail. = 0) %then %do;
        %let _method = PHREG;
        %put NOTE: PH not violated for &label. and &biomarker_label. - using Cox PH model.;
    %end;
    %else %do;
        %let _method = LIFEREG;
        %put NOTE: PH violated for &label. and &biomarker_label. - using Weibull AFT model.;
    %end;

    proc datasets lib=work nolist;
        delete _ph_assess_raw_ _ph_assess_;
    quit;

%end;
%else %do;

    %put NOTE: method=&method. specified manually - skipping PH test for &label. and &biomarker_label.;
    %let _method = %upcase(&method.);

%end;

    /*-------------------------------------------------
      3. Pairwise log-rank tests
    -------------------------------------------------*/

    ods exclude all;
    ods noresults;

    ods output HomTests = LR12(where=(upcase(Test)="LOG-RANK"));
    proc lifetest data=Temp(where=(&group. in (1,2)));
        time &time.*&outcome.(0);
        strata &group. / test=logrank;
    run;

    ods output HomTests = LR13(where=(upcase(Test)="LOG-RANK"));
    proc lifetest data=Temp(where=(&group. in (1,3)));
        time &time.*&outcome.(0);
        strata &group. / test=logrank;
    run;

    ods output HomTests = LR23(where=(upcase(Test)="LOG-RANK"));
    proc lifetest data=Temp(where=(&group. in (2,3)));
        time &time.*&outcome.(0);
        strata &group. / test=logrank;
    run;

    ods results;
    ods exclude none;


    /*-------------------------------------------------
      4. Overall KM + log-rank
    -------------------------------------------------*/

    ods select none;

    proc lifetest data=Temp
        plots=survival(cb=hw atrisk=0 to 12 by 1)
        confband=hw bandmin=0 bandmax=12;

        time &time.*&outcome.(0);
        strata &group. / test=logrank;
        id patno;
        label &time. = "Year";

        ods output HomTests    = HomTests
                   SurvivalPlot = survsaa3;
    run;

    ods select all;


    /*-------------------------------------------------
      5. Overall model p-value and AFT covariance pieces
    -------------------------------------------------*/

    %if &_method. = PHREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ModelANOVA = _survp_all;

        proc phreg data=Temp;
            class &group.(ref='3') gen(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &group. / type3;
        run;

        ods results;
        ods exclude none;

        data survp;
            set _survp_all;

            if index(upcase("&group."), upcase(strip(Effect))) > 0
               and not missing(Effect);
        run;

        proc datasets lib=work nolist;
            delete _survp_all;
        quit;

    %end;
    %else %if &_method. = LIFEREG %then %do;

        ods exclude all;
        ods noresults;

        ods output Type3Analysis      = type3
                   Estimates          = param_estimates
                   ParameterEstimates = pe_full
                   CovB               = covb_full;

        proc lifereg data=Temp outest=_aft_outest covout;
            class &group. gen;
            model &time.*&outcome.(0) = &covars. &group. / dist=weibull covb;

            estimate "1 vs 2" &group.  1 -1  0;
            estimate "1 vs 3" &group.  1  0 -1;
            estimate "2 vs 3" &group.  0  1 -1;
        run;

        ods results;
        ods exclude none;

        proc datasets lib=work nolist;
            modify type3;
            rename Variable = Effect;
        quit;

        data survp;
            set type3;

            if index(upcase("&group."), upcase(strip(Effect))) > 0
               and not missing(Effect);
        run;

        data _null_;
            set pe_full;

            if upcase(Parameter) = "SCALE" then
                call symputx("sigma", Estimate);
        run;

        proc sql noprint;
            select Scale into :var_sigma trimmed
            from covb_full
            where upcase(Variable) = "SCALE";

            select Scale into :cov_g1_sigma trimmed
            from covb_full
            where upcase(Variable) = upcase("&group.1");

            select Scale into :cov_g2_sigma trimmed
            from covb_full
            where upcase(Variable) = upcase("&group.2");

            select &group.1, &group.2
                into :var_b1 trimmed, :cov_b1b2 trimmed
            from covb_full
            where upcase(Variable) = upcase("&group.1");

            select &group.2
                into :var_b2 trimmed
            from covb_full
            where upcase(Variable) = upcase("&group.2");
        quit;

        %if %superq(var_sigma)=    %then %let var_sigma=.;
        %if %superq(cov_g1_sigma)= %then %let cov_g1_sigma=.;
        %if %superq(cov_g2_sigma)= %then %let cov_g2_sigma=.;
        %if %superq(var_b1)=       %then %let var_b1=.;
        %if %superq(var_b2)=       %then %let var_b2=.;
        %if %superq(cov_b1b2)=     %then %let cov_b1b2=.;

    %end;
    %else %do;

        %put WARNING: method=&_method. not recognized. Falling back to PHREG.;
        %let _method = PHREG;

        ods exclude all;
        ods noresults;

        ods output ModelANOVA = _survp_all;

        proc phreg data=Temp;
            class &group.(ref='3') gen(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &group. / type3;
        run;

        ods results;
        ods exclude none;

        data survp;
            set _survp_all;

            if index(upcase("&group."), upcase(strip(Effect))) > 0
               and not missing(Effect);
        run;

        proc datasets lib=work nolist;
            delete _survp_all;
        quit;

    %end;


    /*-------------------------------------------------
      6. Overall log-rank and model p-value
    -------------------------------------------------*/

    data _null_;
        set HomTests;
        where upcase(Test) = "LOG-RANK";
        call symputx("pval", put(ProbChiSq, pvalue6.4));
    run;

    data _null_;
        set survp;
        call symputx("pval2",  put(ProbChiSq, pvalue6.4));
        call symputx("pval2n", ProbChiSq);
    run;


    /*-------------------------------------------------
      7. Pairwise model comparisons
    -------------------------------------------------*/

    %if &_method. = PHREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ContrastTest = PHCon
                   HazardRatios = PHHR;

        proc phreg data=Temp;
            class &group.(ref='3') gen(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &group.;

            contrast "1 vs 2" &group.  1 -1;
            contrast "1 vs 3" &group.  1  0;
            contrast "2 vs 3" &group.  0  1;

            hazardratio &group. / diff=all cl=wald;
        run;

        ods results;
        ods exclude none;

    %end;
    %else %if &_method. = LIFEREG %then %do;

        data PHCon;
            set param_estimates(rename=(Label=Contrast));

            ProbChiSq = coalesce(ProbChiSq, ProbZ);

            keep Contrast ProbChiSq Estimate StdErr;
        run;

    %end;


    /*-------------------------------------------------
      8. Normalize pairwise p-values
    -------------------------------------------------*/

    %if &_method. = LIFEREG %then %do;

        data temp_phcon;
            set PHCon;

            length Comparison $12;

            _c = upcase(strip(compbl(Contrast)));

            if      _c = "1 VS 2" then Comparison = "T1 vs T2";
            else if _c = "1 VS 3" then Comparison = "T1 vs T3";
            else if _c = "2 VS 3" then Comparison = "T2 vs T3";
            else delete;

            CoxPH_p = ProbChiSq;

            keep Comparison CoxPH_p Estimate StdErr;
        run;

    %end;
    %else %do;

        data temp_phcon;
            set PHCon(keep=Contrast ProbChiSq);

            length Comparison $12;

            _c = upcase(strip(compbl(Contrast)));

            if      _c = "1 VS 2" then Comparison = "T1 vs T2";
            else if _c = "1 VS 3" then Comparison = "T1 vs T3";
            else if _c = "2 VS 3" then Comparison = "T2 vs T3";
            else delete;

            CoxPH_p = ProbChiSq;

            keep Comparison CoxPH_p;
        run;

    %end;


    /*-------------------------------------------------
      9. HR CI:
         PHREG direct
         LIFEREG delta method
    -------------------------------------------------*/

    %if &_method. = PHREG %then %do;

        data temp_phhr;
            set PHHR;

            length Comparison $12 HR_CI $50 _desc $100;

            _desc = upcase(compbl(Description));

            if      index(_desc, "1 VS 2") then Comparison = "T1 vs T2";
            else if index(_desc, "1 VS 3") then Comparison = "T1 vs T3";
            else if index(_desc, "2 VS 3") then Comparison = "T2 vs T3";
            else delete;

            HR     = HazardRatio;
            HR_LCL = WaldLower;
            HR_UCL = WaldUpper;

            HR_CI = catx(" ",
                     put(HR, 8.3),
                     cats("(",
                          put(HR_LCL, 8.3), ", ",
                          put(HR_UCL, 8.3), ")"));

            keep Comparison HR HR_LCL HR_UCL HR_CI;
        run;

    %end;
    %else %if &_method. = LIFEREG %then %do;

        data cov_beta_sigma;
            length Comparison $12;

            Comparison = "T1 vs T2";
            cov_beta_sigma = &cov_g1_sigma. - &cov_g2_sigma.;
            output;

            Comparison = "T1 vs T3";
            cov_beta_sigma = &cov_g1_sigma.;
            output;

            Comparison = "T2 vs T3";
            cov_beta_sigma = &cov_g2_sigma.;
            output;
        run;

        proc sort data=cov_beta_sigma;
            by Comparison;
        run;

        proc sort data=temp_phcon;
            by Comparison;
        run;

        data temp_phhr;
            merge temp_phcon(in=a) cov_beta_sigma;
            by Comparison;

            length HR_CI $50;

            if a;

            HR = .;
            HR_LCL = .;
            HR_UCL = .;
            logHR = .;
            var_logHR = .;
            se_logHR = .;

            if not missing(Estimate)
               and "&sigma." ne "." and "&sigma." ne ""
               and "&var_sigma." ne "." and "&var_sigma." ne "" then do;

                logHR = -Estimate / &sigma.;
                HR    = exp(logHR);

                if Comparison = "T1 vs T2" then do;

                    var_logHR =
                          (1 / (&sigma.**2)) * &var_b1.
                        + (1 / (&sigma.**2)) * &var_b2.
                        + ((Estimate**2) / (&sigma.**4)) * &var_sigma.
                        - (2 / (&sigma.**2)) * &cov_b1b2.
                        - (2 * Estimate / (&sigma.**3)) * &cov_g1_sigma.
                        + (2 * Estimate / (&sigma.**3)) * &cov_g2_sigma.;

                end;
                else if not missing(StdErr) and not missing(cov_beta_sigma) then do;

                    var_logHR =
                          (1 / (&sigma.**2)) * (StdErr**2)
                        + ((Estimate**2) / (&sigma.**4)) * &var_sigma.
                        - (2 * Estimate / (&sigma.**3)) * cov_beta_sigma;

                end;

                if not missing(var_logHR) and var_logHR >= 0 then do;

                    se_logHR = sqrt(var_logHR);
                    HR_LCL   = exp(logHR - 1.96 * se_logHR);
                    HR_UCL   = exp(logHR + 1.96 * se_logHR);

                    HR_CI = catx(" ",
                             put(HR, 8.3),
                             cats("(",
                                  put(HR_LCL, 8.3), ", ",
                                  put(HR_UCL, 8.3), ")"));

                end;
                else HR_CI = "";

            end;
            else HR_CI = "";

            keep Comparison HR HR_LCL HR_UCL HR_CI;
        run;

    %end;


    /*-------------------------------------------------
      10. Pairwise log-rank table
    -------------------------------------------------*/

    data temp_lr;
        length Comparison $12;

        set LR12(in=a keep=ProbChiSq rename=(ProbChiSq=LogRank_p))
            LR13(in=b keep=ProbChiSq rename=(ProbChiSq=LogRank_p))
            LR23(in=c keep=ProbChiSq rename=(ProbChiSq=LogRank_p));

        if      a then Comparison = "T1 vs T2";
        else if b then Comparison = "T1 vs T3";
        else if c then Comparison = "T2 vs T3";

        keep Comparison LogRank_p;
    run;

    proc sort data=temp_lr;
        by Comparison;
    run;

    proc sort data=temp_phcon;
        by Comparison;
    run;

    proc sort data=temp_phhr;
        by Comparison;
    run;

    data PVALS_pairs;
        merge temp_lr(in=l) temp_phcon(in=c) temp_phhr(in=h);
        by Comparison;

        length LogRank CoxPH $8 HR_CI $50;

        if l;

        LogRank = put(LogRank_p, pvalue6.4);
        CoxPH   = put(CoxPH_p, pvalue6.4);

        keep Comparison LogRank CoxPH HR HR_LCL HR_UCL HR_CI;
    run;

    data PVALS_overall;
        length Comparison $12 LogRank $8 CoxPH $8 HR_CI $50;

        Comparison = "Overall";
        LogRank    = "&pval.";
        CoxPH      = "&pval2.";
        HR_CI      = "NA";

        HR     = .;
        HR_LCL = .;
        HR_UCL = .;

        output;
    run;

    data PVALS;
        set PVALS_overall PVALS_pairs;
    run;

    proc sql noprint;
        select put(CoxPH_p, pvalue6.4) into :p12 trimmed
        from temp_phcon
        where Comparison = "T1 vs T2";

        select put(CoxPH_p, pvalue6.4) into :p13 trimmed
        from temp_phcon
        where Comparison = "T1 vs T3";

        select put(CoxPH_p, pvalue6.4) into :p23 trimmed
        from temp_phcon
        where Comparison = "T2 vs T3";
    quit;

    %if %superq(p12)= %then %let p12=.;
    %if %superq(p13)= %then %let p13=.;
    %if %superq(p23)= %then %let p23=.;


    /*-------------------------------------------------
      11. Additional model statistics
    -------------------------------------------------*/

    %let _modelLabel = %sysfunc(ifc(&_method. = LIFEREG, AFT, CoxPH));


    %if &_method. = PHREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ParameterEstimates = ModelPE
                   Estimates          = PairwiseEst;

        proc phreg data=Temp;
            class &group.(ref='3') gen(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &group.;

            estimate "T1 vs T2" &group.  1 -1 / exp;
            estimate "T1 vs T3" &group.  1  0 / exp;
            estimate "T2 vs T3" &group.  0  1 / exp;
        run;

        ods results;
        ods exclude none;

    %end;
    %else %if &_method. = LIFEREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ParameterEstimates = ModelPE
                   Estimates          = PairwiseEst;

        proc lifereg data=Temp;
            class &group. gen;
            model &time.*&outcome.(0) = &covars. &group. / dist=weibull;

            estimate "T1 vs T2" &group.  1 -1  0;
            estimate "T1 vs T3" &group.  1  0 -1;
            estimate "T2 vs T3" &group.  0  1 -1;
        run;

        ods results;
        ods exclude none;

    %end;


    /*-------------------------------------------------
      12. Overall group omnibus row
    -------------------------------------------------*/

    data OverallGroup;
        length Biomarker $100 Milestone $100 Method $10 Parameter $80 ClassVal $40;
        length Comparison_Type $35 Comparison $12 Beta_CI $40 HR_CI $50;

        set survp;

        Biomarker      = "&biomarker_label.";
        Milestone      = "&label.";
        MilestoneOrder = &milestone_order.;
        Method         = "&_modelLabel.";

        Parameter       = "Overall &grouplabel.";
        ClassVal        = "";
        Comparison      = "Overall";
        Comparison_Type = "Overall omnibus";

        GroupDF   = DF;
        WaldStat  = WaldChiSq;

        Beta      = .;
        SE        = .;
        Beta_CI   = "";

        if not missing(WaldChiSq) and WaldChiSq >= 0 then
            Z_stat = sqrt(WaldChiSq);
        else
            Z_stat = .;

        if not missing(ProbChiSq) then P_value = ProbChiSq;
        else P_value = .;

        HR      = .;
        HR_LCL  = .;
        HR_UCL  = .;
        HR_CI   = "";

        keep Biomarker Milestone MilestoneOrder Method Parameter ClassVal
             Comparison_Type Comparison GroupDF WaldStat Beta SE Beta_CI
             Z_stat P_value HR HR_LCL HR_UCL HR_CI;
    run;


    /*-------------------------------------------------
      13. Model parameter estimates
    -------------------------------------------------*/

    data ModelPE_Tidy;
        length Biomarker $100 Milestone $100 Method $10 Parameter $80 ClassVal $40;
        length Comparison_Type $35 Comparison $12 Beta_CI $40 HR_CI $50;

        set ModelPE;

        Biomarker      = "&biomarker_label.";
        Milestone      = "&label.";
        MilestoneOrder = &milestone_order.;
        Method         = "&_modelLabel.";

        Beta = Estimate;
        SE   = StdErr;

        GroupDF  = .;
        WaldStat = .;

        Comparison      = "";
        Comparison_Type = "Model parameter";

        %if &_method. = PHREG %then %do;

            ClassVal = strip(ClassVal0);

            if not missing(ChiSq) and ChiSq >= 0 then
                Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if SE > 0 then
                Z_stat = Beta / SE;
            else
                Z_stat = .;

            P_value = ProbChiSq;

        %end;
        %else %do;

            ClassVal = strip(Level1);

            if not missing(ChiSq) and ChiSq >= 0 then
                Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if SE > 0 then
                Z_stat = Beta / SE;
            else
                Z_stat = .;

            if not missing(Probt) then
                P_value = Probt;
            else if not missing(ProbChiSq) then
                P_value = ProbChiSq;
            else
                P_value = .;

        %end;

        if not missing(Beta) and not missing(SE) and SE > 0 then
            Beta_CI = catx(" ",
                      put(Beta, 8.4),
                      cats("(",
                           put(Beta - 1.96*SE, 8.4), ", ",
                           put(Beta + 1.96*SE, 8.4), ")"));
        else
            Beta_CI = "";

        HR     = .;
        HR_LCL = .;
        HR_UCL = .;
        HR_CI  = "";

        if DF = 0 then delete;
        if upcase(strip(Parameter)) = upcase("&group.") then delete;
        if upcase(strip(Parameter)) = "SCALE" then delete;
        if upcase(strip(Parameter)) = "WEIBULL SHAPE" then delete;

        keep Biomarker Milestone MilestoneOrder Method Parameter ClassVal
             Comparison_Type Comparison GroupDF WaldStat Beta SE Beta_CI
             Z_stat P_value HR HR_LCL HR_UCL HR_CI;
    run;


    /*-------------------------------------------------
      14. Pairwise estimates, with HR CI merged in
    -------------------------------------------------*/

    data PairwiseEst_Tidy;
        length Biomarker $100 Milestone $100 Method $10 Parameter $80 ClassVal $40;
        length Comparison_Type $35 Comparison $12 Beta_CI $40;

        set PairwiseEst;

        Biomarker      = "&biomarker_label.";
        Milestone      = "&label.";
        MilestoneOrder = &milestone_order.;
        Method         = "&_modelLabel.";

        Parameter = Label;
        ClassVal  = "";
        Comparison_Type = "Pairwise tertile comparison";
        Comparison = strip(Label);

        Beta = Estimate;
        SE   = StdErr;

        GroupDF  = .;
        WaldStat = .;

        %if &_method. = PHREG %then %do;

            if not missing(ChiSq) and ChiSq >= 0 then
                Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if SE > 0 then
                Z_stat = Beta / SE;
            else
                Z_stat = .;

            if not missing(ProbZ) then
                P_value = ProbZ;
            else if not missing(ProbChiSq) then
                P_value = ProbChiSq;
            else
                P_value = .;

        %end;
        %else %do;

            if not missing(ZValue) then
                Z_stat = ZValue;
            else if SE > 0 then
                Z_stat = Beta / SE;
            else
                Z_stat = .;

            if not missing(ProbZ) then
                P_value = ProbZ;
            else if not missing(ProbChiSq) then
                P_value = ProbChiSq;
            else
                P_value = .;

        %end;

        if not missing(Beta) and not missing(SE) and SE > 0 then
            Beta_CI = catx(" ",
                      put(Beta, 8.4),
                      cats("(",
                           put(Beta - 1.96*SE, 8.4), ", ",
                           put(Beta + 1.96*SE, 8.4), ")"));
        else
            Beta_CI = "";

        keep Biomarker Milestone MilestoneOrder Method Parameter ClassVal
             Comparison_Type Comparison GroupDF WaldStat Beta SE Beta_CI
             Z_stat P_value;
    run;

    proc sort data=PairwiseEst_Tidy;
        by Comparison;
    run;

    proc sort data=temp_phhr;
        by Comparison;
    run;

    data PairwiseEst_Tidy;
        merge PairwiseEst_Tidy(in=a)
              temp_phhr(keep=Comparison HR HR_LCL HR_UCL HR_CI);
        by Comparison;

        if a;
    run;


    /*-------------------------------------------------
      15. Combine model statistics
    -------------------------------------------------*/

    data ModelStats;
        length Biomarker $100 Milestone $100 Method $10 Parameter $80 ClassVal $40;
        length Comparison_Type $35 Comparison $12 Beta_CI $40 HR_CI $50;

        set OverallGroup ModelPE_Tidy PairwiseEst_Tidy;
    run;

    proc append base=All_ModelStats data=ModelStats force;
    run;

    data MS_&_dslabel.;
        set ModelStats;
    run;


    /*-------------------------------------------------
      16. Detect CI band
    -------------------------------------------------*/

    %let dsid = %sysfunc(open(survsaa3));

    %if &dsid. %then %do;
        %let hasLCLM = %sysfunc(varnum(&dsid., LCLM));
        %let rc = %sysfunc(close(&dsid.));
    %end;
    %else %let hasLCLM = 0;


    /*-------------------------------------------------
      17. KM Plot
    -------------------------------------------------*/

    ods escapechar='^';

    %if %sysevalf(%superq(pval2n)=, boolean) or %sysevalf(&pval2n >= 0.05) %then %do;
      footnote1 height=8pt
        "&_modelLabel p-value | Overall Interaction P-value: &pval2";
    %end;
    %else %do;
      footnote1 height=8pt
        "&_modelLabel p-values | Overall: &pval2 | T1 vs T2: &p12 | T1 vs T3: &p13 | T2 vs T3: &p23";
    %end;

    title "&label. &title_suffix";

    proc sgplot data=survsaa3;
      styleattrs datacolors=(blue red green)
                 datacontrastcolors=(blue red green);
      step    x=Time y=Survival / group=stratumnum name="km";
      scatter x=Time y=Censored / group=stratumnum markerattrs=(symbol=plus);
      band x=Time upper=HW_UCL lower=HW_LCL /
           group=stratumnum transparency=0.80;
      %if &hasLCLM %then %do;
        band x=Time lower=LCLM upper=UCLM /
             group=stratumnum transparency=0.80;
      %end;
      xaxistable Atrisk / location=inside x=tatrisk
                          class=stratumnum colorgroup=stratumnum;
      xaxis values=(0 to 12 by 1);
      format stratumnum &groupfmt.;
      label  stratumnum="&grouplabel";
      keylegend "km" / title="&grouplabel"
                       location=outside position=bottom
                       across=3 opaque;
    run;

    footnote;
    title;


    /*-------------------------------------------------
      18. P-value / HR Summary Table
    -------------------------------------------------*/

    %if &_method. = PHREG %then %do;
        footnote1 j=left h=8pt
            "HR (95% CI) from Cox proportional hazards model. Reference: T3.";
    %end;
    %else %if &_method. = LIFEREG %then %do;
        footnote1 j=left h=8pt
            "Derived HR (95% CI) from Weibull AFT model using HR=exp(-beta/sigma); delta-method CI on log(HR). Reference: T3.";
    %end;

    title2 "P-Value and HR Summary - &label. &title_suffix";

    proc report data=PVALS nowd headline;
        columns Comparison LogRank CoxPH HR_CI;

        define Comparison / 'Comparison'           width=12;
        define LogRank    / 'Log-Rank p-value'     width=16;
        define CoxPH      / "&_modelLabel p-value" width=14;

        %if &_method. = LIFEREG %then %do;
            define HR_CI / 'Derived HR (95% CI)' width=24;
        %end;
        %else %do;
            define HR_CI / 'HR (95% CI)'         width=24;
        %end;
    run;

    footnote;
    title;


    /*-------------------------------------------------
      19. Compact model statistics table
    -------------------------------------------------*/

    data ModelStats_PDF;
        length Parameter_Display $80;
        set ModelStats;

        Parameter_Display = Parameter;

        if index(upcase(Parameter), "OVERALL") then
            Parameter_Display = "Overall tertile effect";
        else if upcase(strip(Parameter)) = "T1 VS T2" then
            Parameter_Display = "T1 vs T2";
        else if upcase(strip(Parameter)) = "T1 VS T3" then
            Parameter_Display = "T1 vs T3";
        else if upcase(strip(Parameter)) = "T2 VS T3" then
            Parameter_Display = "T2 vs T3";

        if index(upcase(Parameter), "OVERALL") or
           upcase(strip(Parameter)) in ("T1 VS T2", "T1 VS T3", "T2 VS T3");
    run;

    title2 "Adjusted Model Statistics - &label. &title_suffix";

    proc report data=ModelStats_PDF nowd headline headskip spacing=1 split='|'
        style(report)=[cellspacing=1 cellpadding=2 just=center]
        style(header)=[font_size=8pt font_weight=bold just=center]
        style(column)=[font_size=8pt just=center];

        columns Method Parameter_Display GroupDF WaldStat Beta_CI SE Z_stat P_value HR_CI;

        define Method / display "Model"
            style(column)=[cellwidth=0.8in just=left];

        define Parameter_Display / display "Parameter / Comparison"
            style(column)=[cellwidth=1.4in just=left];

        define GroupDF / display "DF" format=4.
            style(column)=[cellwidth=0.45in];

        define WaldStat / display "Wald Chi-Sq" format=8.4
            style(column)=[cellwidth=0.85in];

        define Beta_CI / display "Beta (95% CI)"
            style(column)=[cellwidth=1.6in just=left];

        define SE / display "SE" format=8.4
            style(column)=[cellwidth=0.6in];

        define Z_stat / display "Z" format=8.4
            style(column)=[cellwidth=0.6in];

        define P_value / display "p-value" format=pvalue6.4
            style(column)=[cellwidth=0.7in];

        %if &_method. = LIFEREG %then %do;
            define HR_CI / display "Derived HR (95% CI)"
                style(column)=[cellwidth=1.6in just=left];
        %end;
        %else %do;
            define HR_CI / display "HR (95% CI)"
                style(column)=[cellwidth=1.6in just=left];
        %end;

    run;

    title2;


    /*-------------------------------------------------
      20. Baseline event listing
    -------------------------------------------------*/

    title "Participants meeting &label. at Baseline";

    proc print data=&data. noobs;
        where &outcome. = 1 and &time. = 0;
        var &baseline_vars.;
    run;

    title;


    /*-------------------------------------------------
      21. BL vs Post stacked bar
    -------------------------------------------------*/

    proc format;
        value bl
            1 = "BL"
            2 = "Post";
    run;

    data temp2;
      set &data.;
      if &outcome.=1 then do;
        if      &time.=0 then BL=1;
        else if &time.>0 then BL=2;
      end;
    run;

    data _seed_;
      do BL=1 to 2;
        do &stack_group.=1 to 3;
            _wt=0;
            output;
        end;
      end;
    run;

    data temp2_plot;
      set _seed_ temp2(in=real);
      if real then _wt=1;
    run;

    proc sort data=temp2_plot;
        by BL &stack_group.;
    run;

    proc sql noprint;
      select distinct BL into :bl_list separated by ' '
      from temp2
      where BL in (1,2) and not missing(BL);
    quit;

    %if %length(&bl_list) %then %do;
      title "Event Timing (BL vs Post-BL) - &label. &title_suffix";

      proc sgplot data=temp2_plot(where=(BL in (&bl_list)));
        styleattrs datacolors=(cx556B8A cxA05252 cx518A85);
        vbar BL / response=_wt stat=sum group=&stack_group.
                  grouporder=ascending groupdisplay=stack;
        format &stack_group. &stack_groupfmt. BL bl.;
        label  &stack_group.="&stack_grouplabel" BL="BL-Post";
        yaxis label="Frequency";
        xaxis display=(nolabel);
      run;
    %end;
    %else %do;
      %put NOTE: No BL=1 or BL=2 values present - skipping stacked-bar plot.;
    %end;

    title;


    /*-------------------------------------------------
      22. Cleanup
    -------------------------------------------------*/

    proc datasets library=work nolist;
        delete Temp
               LR12 LR13 LR23
               HomTests survsaa3 survp
               PHCon PHHR temp_phcon temp_phhr temp_lr
               PVALS_pairs PVALS_overall PVALS
               ModelPE PairwiseEst ModelPE_Tidy PairwiseEst_Tidy
               OverallGroup ModelStats ModelStats_PDF
               _aft_outest pe_full covb_full cov_beta_sigma
               _seed_ temp2 temp2_plot
               type3 param_estimates;
    quit;

%mend SurvTertile;

/*===============================================================
  5. BIOMARKER-SPECIFIC CLINICAL SUMMARY MACRO
===============================================================*/

%macro ClinicalTertileSummary(
    biomarker_label=,
    title_suffix=for Sporadic PD
);

    proc sql noprint;
        select memname
            into :_mslist separated by ' '
        from dictionary.tables
        where libname = "WORK"
          and upcase(substr(memname, 1, 3)) = "MS_"
        order by crdate;
    quit;

    %if %superq(_mslist)= %then %do;
        %put WARNING: No MS_ datasets found. Run %SurvTertile calls first.;
        %return;
    %end;

    data All_Milestones;
        set &_mslist.;
    run;

    data Summary_Data;
        set All_Milestones;

        length SortOrder 8;

        if index(upcase(Parameter), "OVERALL") then SortOrder = 1;
        else if upcase(strip(Parameter)) = "T1 VS T2" then SortOrder = 2;
        else if upcase(strip(Parameter)) = "T1 VS T3" then SortOrder = 3;
        else if upcase(strip(Parameter)) = "T2 VS T3" then SortOrder = 4;
        else delete;

        length Beta_CI_4 $70 SE_4 $20 Z_4 $20 P_4 $20 HR_CI_4 $70;

        P_Sig = (not missing(P_value) and P_value < 0.05);

        if not missing(Beta) and not missing(SE) then do;
            Beta_LCL = Beta - 1.96 * SE;
            Beta_UCL = Beta + 1.96 * SE;

            Beta_CI_4 = catx(" ",
                         put(Beta, 8.3),
                         cats("(",
                              put(Beta_LCL, 8.3), ", ",
                              put(Beta_UCL, 8.3), ")"));
        end;
        else Beta_CI_4 = "";

        if not missing(SE) then SE_4 = put(SE, 8.3);
        else SE_4 = "";

        if not missing(Z_stat) then Z_4 = put(Z_stat, 8.3);
        else Z_4 = "";

        if not missing(P_value) then do;
            if P_value < 0.001 then P_4 = "<0.001";
            else P_4 = put(P_value, 8.3);
        end;
        else P_4 = "";

        if not missing(HR) and not missing(HR_LCL) and not missing(HR_UCL) then do;
            HR_CI_4 = catx(" ",
                       put(HR, 8.3),
                       cats("(",
                            put(HR_LCL, 8.3), ", ",
                            put(HR_UCL, 8.3), ")"));
        end;
        else if not missing(HR_CI) then HR_CI_4 = HR_CI;
        else HR_CI_4 = "";

        keep Biomarker Milestone MilestoneOrder Method Parameter ClassVal
             Comparison_Type Comparison SortOrder GroupDF WaldStat
             Beta SE Z_stat P_value HR HR_LCL HR_UCL Beta_CI HR_CI
             Beta_CI_4 SE_4 Z_4 P_4 P_Sig HR_CI_4;
    run;

    /* Always show the overall row. Show pairwise tertile rows only
       when the overall p-value is < 0.05 for that biomarker/milestone. */
    proc sql;
        create table Summary_Data_filtered as
        select a.*
        from Summary_Data as a
        left join
        (
            select Biomarker,
                   Milestone,
                   MilestoneOrder,
                   max(case
                           when SortOrder = 1 and
                                not missing(P_value) and
                                P_value < 0.05
                           then 1
                           else 0
                       end) as Overall_Sig
            from Summary_Data
            group by Biomarker, Milestone, MilestoneOrder
        ) as b
          on  a.Biomarker      = b.Biomarker
          and a.Milestone      = b.Milestone
          and a.MilestoneOrder = b.MilestoneOrder
        where a.SortOrder = 1 or b.Overall_Sig = 1
        order by a.Biomarker, a.MilestoneOrder, a.SortOrder;
    quit;

    data Summary_Data;
        set Summary_Data_filtered;
    run;

    proc sort data=Summary_Data nodupkey;
        by Biomarker MilestoneOrder Milestone SortOrder Parameter Method;
    run;

    /* Replace old rows for this biomarker before appending */
    %if %sysfunc(exist(Clinical_Tertile_Summary)) %then %do;

        data Clinical_Tertile_Summary;
            set Clinical_Tertile_Summary;

            if strip(upcase(Biomarker)) = strip(upcase("&biomarker_label.")) then delete;
        run;

    %end;

    proc append base=Clinical_Tertile_Summary
                data=Summary_Data
                force;
    run;

    proc sort data=Summary_Data;
        by Biomarker MilestoneOrder SortOrder;
    run;

    %if &write_individual_rtf = 1 %then %do;

    options
        orientation = landscape
        leftmargin  = 0.25in
        rightmargin = 0.25in
        topmargin   = 0.40in
        bottommargin= 0.40in
        nodate
        nonumber
    ;


    title "Summary of &biomarker_label. Clinical Milestone Tertile Survival Model Statistics - &title_suffix.";

    footnote1 j=left h=7pt
        "Overall row reports the Type III omnibus test for the tertile biomarker.";
    footnote2 j=left h=7pt
        "Pairwise tertile rows are shown only when the overall p-value is < 0.05. Models adjusted for age and sex.";
    footnote3 j=left h=7pt
        "Cox PH models report HR. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Values rounded to 3 decimals.";

    proc report data=Summary_Data nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=1 cellpadding=2]
        style(header)=[font_size=7pt font_weight=bold]
        style(column)=[font_size=7pt];

        columns Biomarker MilestoneOrder Milestone Method SortOrder Parameter
                GroupDF WaldStat Beta_CI_4 SE_4 Z_4 P_4 HR_CI_4;

        define Biomarker      / group 'Biomarker';
        define MilestoneOrder / order noprint;
        define Milestone      / group 'Clinical Milestone';
        define Method         / display 'Model';
        define SortOrder      / order noprint;
        define Parameter      / display 'Parameter / Comparison';
        define GroupDF        / display 'DF' format=4.;
        define WaldStat       / display 'Wald Chi-Sq' format=10.4;
        define Beta_CI_4      / display 'Beta (95% CI)';
        define SE_4           / display 'SE';
        define Z_4            / display 'Z';
        define P_4            / display 'p-value';
        define HR_CI_4        / display 'HR / Derived HR (95% CI)';

        break after Milestone / skip;

    run;

    footnote;
    title;


    %end;

%mend ClinicalTertileSummary;


/*===============================================================
  6. COMBINED CLINICAL REPORT MACRO
  Biomarker order:
      1. pTau217
      2. NfL
      3. GFAP
      4. C2N AB40
      5. C2N AB42
      6. C2N AB42/40
===============================================================*/

%macro WriteCombinedClinicalTertileRTF(
    title_suffix=for Sporadic PD
);

    %if not %sysfunc(exist(Clinical_Tertile_Summary)) %then %do;
        %put WARNING: Clinical_Tertile_Summary does not exist.;
        %return;
    %end;

    proc sql noprint;
        select count(*)
            into :_n_clinical_combined trimmed
        from Clinical_Tertile_Summary;
    quit;

    %if %sysevalf(&_n_clinical_combined. = 0) %then %do;
        %put WARNING: Clinical_Tertile_Summary is empty. No combined clinical report displayed.;
        %return;
    %end;


    /*-----------------------------------------------------
      Apply preferred biomarker order
    -----------------------------------------------------*/

    data Clinical_Tertile_Summary;
        set Clinical_Tertile_Summary;

        length BiomarkerOrder 8;

        if upcase(strip(Biomarker)) in
            ("PTAU217 PLASMA", "LILLY PTAU217", "PTAU217", "PTAU", "P-TAU217")
            then BiomarkerOrder = 1;

        else if upcase(strip(Biomarker)) in
            ("IU NFL", "NFL", "NF-L", "NFL PLASMA")
            then BiomarkerOrder = 2;

        else if upcase(strip(Biomarker)) in
            ("IU GFAP", "GFAP", "GFAP PLASMA")
            then BiomarkerOrder = 3;

        else if upcase(strip(Biomarker)) in
            ("A^{UNICODE 03B2}40 PLASMA", "C2N AB40", "C2N ABETA40", "AB40", "ABETA40")
            then BiomarkerOrder = 4;

        else if upcase(strip(Biomarker)) in
            ("A^{UNICODE 03B2}42 PLASMA", "C2N AB42", "C2N ABETA42", "AB42", "ABETA42")
            then BiomarkerOrder = 5;

        else if upcase(strip(Biomarker)) in
            ("A^{UNICODE 03B2}42/A^{UNICODE 03B2}40 PLASMA",
             "C2N AB42/40", "C2N AB4240", "C2N ABETA42/40",
             "C2N ABETA4240", "AB42/40", "AB4240", "ABETA42/40")
            then BiomarkerOrder = 6;

        else BiomarkerOrder = 99;
    run;


    proc sort data=Clinical_Tertile_Summary nodupkey;
        by BiomarkerOrder Biomarker MilestoneOrder Milestone SortOrder Parameter Method;
    run;

    proc sort data=Clinical_Tertile_Summary;
        by BiomarkerOrder Biomarker MilestoneOrder SortOrder;
    run;


    options
        orientation = landscape
        leftmargin  = 0.25in
        rightmargin = 0.25in
        topmargin   = 0.40in
        bottommargin= 0.40in
        nodate
        nonumber
    ;



    title1 j=left h=9pt
        "^{style [font_weight=bold] Supplement Table.} Association between baseline plasma biomarker tertiles and time-to-event clinical milestones in";
    title2 j=left h=9pt
        "Sporadic PD (NSD+). Cox proportional hazards (CoxPH) or Weibull accelerated failure time (AFT) models were selected automatically by the";
    title3 j=left h=9pt
        "proportional hazards assumption test. All models adjusted for age and sex. Bold red p-values indicate p < 0.05.";

    footnote1 j=left h=7pt
        "Overall row reports the Type III omnibus test for the tertile biomarker.";
    footnote2 j=left h=7pt
        "Pairwise tertile rows are shown only when the overall p-value is < 0.05. Models adjusted for age and sex.";
    footnote3 j=left h=7pt
        "Cox PH models report HR. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Values rounded to 3 decimals.";


    proc report data=Clinical_Tertile_Summary nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=0 cellpadding=2 frame=hsides rules=groups]
        style(header)=[font_size=7pt font_weight=bold background=white just=center]
        style(column)=[font_size=7pt vjust=middle];

        columns BiomarkerOrder Biomarker MilestoneOrder Milestone Method SortOrder Parameter
                GroupDF WaldStat Beta_CI_4 SE_4 Z_4 P_Sig P_4 HR_CI_4;

        define BiomarkerOrder / order noprint;
        define Biomarker      / group 'Biomarker'
            style(column)=[cellwidth=0.95in font_weight=bold just=left vjust=top];
        define MilestoneOrder / order noprint;
        define Milestone      / group 'Clinical Milestone'
            style(column)=[cellwidth=1.45in just=left];
        define Method         / display 'Model'
            style(column)=[cellwidth=0.45in just=center];
        define SortOrder      / order noprint;
        define Parameter      / display 'Parameter / Comparison'
            style(column)=[cellwidth=0.85in just=left];
        define GroupDF        / display 'DF' format=4.
            style(column)=[cellwidth=0.30in just=center];
        define WaldStat       / display 'Wald Chi-Sq' format=8.3
            style(column)=[cellwidth=0.55in just=center];
        define Beta_CI_4      / display 'Beta (95% CI)'
            style(column)=[cellwidth=1.30in just=center];
        define SE_4           / display 'SE'
            style(column)=[cellwidth=0.40in just=center];
        define Z_4            / display 'Z'
            style(column)=[cellwidth=0.40in just=center];
        define P_Sig          / display noprint;
        define P_4            / display 'p'
            style(header)=[font_style=italic]
            style(column)=[cellwidth=0.45in just=center];
        define HR_CI_4        / display 'HR (95% CI)'
            style(column)=[cellwidth=1.30in just=center];

        compute P_4;
            if P_Sig = 1 then
                call define(_col_, "style", "style=[font_weight=bold foreground=red]");
        endcomp;

        break after Biomarker / skip;
        break after Milestone / skip;

    run;


    footnote;
    title1;
    title2;
    title3;



%mend WriteCombinedClinicalTertileRTF;


/*===============================================================
  7. MAIN BIOMARKER CLINICAL WRAPPER
===============================================================*/

%macro RunOneClinicalTertileBiomarker(
    biomarker_short=,
    biomarker_label=,
    datevar=,
    biomarker_var=,
    cutoff_type=,
    tertile_fmt=,

    method_hy=AUTO,
    method_cog=AUTO,
    method_moca21=AUTO,
    method_moca26=AUTO,
    method_updrs=AUTO
);



    %put NOTE: ============================================================;
    %put NOTE: Running clinical milestone tertile survival analysis for &biomarker_label.;
    %put NOTE: ============================================================;


    /*-----------------------------------------------------
      7A. Create biomarker-specific longitudinal dataset
    -----------------------------------------------------*/

    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data ClinicalSurvData_&biomarker_short.;
        set clinical_analysis_v3;
        where Analysis_Source="SP_TC"
          and Analysis_Biomarker="&biomarker_short.";
    run;



    /*-----------------------------------------------------
      7G. Initialize model stats for this biomarker
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete All_ModelStats MS_: All_Milestones Summary_Data;
    quit;

    data All_ModelStats;
        length
            Biomarker $100
            Milestone $100
            Method $10
            Parameter $80
            ClassVal $40
            Comparison_Type $35
            Comparison $12
            MilestoneOrder 8
            GroupDF 8
            WaldStat 8
            Beta 8
            SE 8
            Beta_CI $40
            Z_stat 8
            P_value 8
            HR 8
            HR_LCL 8
            HR_UCL 8
            HR_CI $50
        ;
        stop;
    run;


    /*-----------------------------------------------------
      7H. Run all 5 clinical milestones
    -----------------------------------------------------*/

    options orientation=landscape;
    ods graphics / reset width=6.8in height=4.8in noborder;

    %if &write_individual_pdf = 1 %then %do;

        title "Clinical Milestone Tertile Survival Analysis - &biomarker_label. - Sporadic PD NSD+";
    %end;


    %SurvTertile(
        data=ClinicalSurvData_&biomarker_short.,
        outcome=STATUS_HY_GE3,
        group=Early_Tertile_Group,
        time=TIME_ENRL_HY_GE3,
        label=Time to HY >= 3,
        milestone_order=1,
        method=&method_hy.,
        covars=Age Gen,
        grouplabel=&biomarker_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group Subgroup_definition STATUS_HY_GE3 TIME_ENRL_HY_GE3,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&biomarker_label.,
        biomarker_label=&biomarker_label.
    );


    %SurvTertile(
        data=ClinicalSurvData_&biomarker_short.,
        outcome=STATUS_COGSTATE_CHANGE,
        group=Early_Tertile_Group,
        time=TIME_ENRL_COGSTATE_CHANGE,
        label=Time to Cognitive State Change,
        milestone_order=2,
        method=&method_cog.,
        covars=Age Gen,
        grouplabel=&biomarker_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group Subgroup_definition STATUS_COGSTATE_CHANGE TIME_ENRL_COGSTATE_CHANGE,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&biomarker_label.,
        biomarker_label=&biomarker_label.
    );


    %SurvTertile(
        data=ClinicalSurvData_&biomarker_short.,
        outcome=STATUS_MOCA_LT21,
        group=Early_Tertile_Group,
        time=TIME_ENRL_MOCA_LT_21,
        label=Time to MoCA < 21,
        milestone_order=3,
        method=&method_moca21.,
        covars=Age Gen,
        grouplabel=&biomarker_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group Subgroup_definition STATUS_MOCA_LT21 TIME_ENRL_MOCA_LT_21,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&biomarker_label.,
        biomarker_label=&biomarker_label.
    );


    %SurvTertile(
        data=ClinicalSurvData_&biomarker_short.,
        outcome=STATUS_MOCA_LT26,
        group=Early_Tertile_Group,
        time=TIME_ENRL_MOCA_LT_26,
        label=Time to MoCA < 26,
        milestone_order=4,
        method=&method_moca26.,
        covars=Age Gen,
        grouplabel=&biomarker_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group Subgroup_definition STATUS_MOCA_LT26 TIME_ENRL_MOCA_LT_26,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&biomarker_label.,
        biomarker_label=&biomarker_label.
    );


    %SurvTertile(
        data=ClinicalSurvData_&biomarker_short.,
        outcome=STATUS_UPDRS_Change_GT5,
        group=Early_Tertile_Group,
        time=TIME_ENRL_UPDRS_Change_GT5,
        label=Time to MDS-UPDRS III Change 5+,
        milestone_order=5,
        method=&method_updrs.,
        covars=Age Gen,
        grouplabel=&biomarker_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group Subgroup_definition STATUS_UPDRS_Change_GT5 TIME_ENRL_UPDRS_Change_GT5,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&biomarker_label.,
        biomarker_label=&biomarker_label.
    );


    %if &write_individual_pdf = 1 %then %do;

        title;
    %end;
    ods graphics / reset;


    /*-----------------------------------------------------
      7I. Display biomarker-specific summary
    -----------------------------------------------------*/

    %ClinicalTertileSummary(
        biomarker_label=&biomarker_label.,
        title_suffix=for Sporadic PD
    );


    /*-----------------------------------------------------
      7J. Cleanup biomarker-specific datasets
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete mydata_&biomarker_short.
               EARLY_&biomarker_short.
               datalong_&biomarker_short.
               lastobsvisit_&biomarker_short.
               allvisits_&biomarker_short.
               SurvData_Long_&biomarker_short.
               ClinicalSurvData_&biomarker_short.
               STATUS_HY_GE3
               STATUS_COGSTATE_CHANGE
               STATUS_MOCA_LT21
               STATUS_MOCA_LT26
               STATUS_UPDRS_Change_GT5
               All_ModelStats
               All_Milestones
               Summary_Data
               Summary_Data_filtered
               MS_:;
    quit;

%mend RunOneClinicalTertileBiomarker;


/*===============================================================
  8. RUN CLINICAL MILESTONES FOR EACH BIOMARKER
===============================================================*/


/*-------------------------
  Lilly pTau217
-------------------------*/

%RunOneClinicalTertileBiomarker(
    biomarker_short=Lilly_pTau217,
    biomarker_label=pTau217 plasma,
    datevar=Lilly_ptau217p_date,
    biomarker_var=Lilly_ptau217p,
    cutoff_type=PTAU,
    tertile_fmt=ptau.,
    method_hy=AUTO,
    method_cog=AUTO,
    method_moca21=AUTO,
    method_moca26=AUTO,
    method_updrs=AUTO
);


/*-------------------------
  IU NfL
-------------------------*/

%RunOneClinicalTertileBiomarker(
    biomarker_short=IU_NFL,
    biomarker_label=NFL,
    datevar=IU_Plasma_date,
    biomarker_var=IU_NFL_plasma_mean,
    cutoff_type=NFL,
    tertile_fmt=nfl.,
    method_hy=AUTO,
    method_cog=AUTO,
    method_moca21=AUTO,
    method_moca26=AUTO,
    method_updrs=AUTO
);


/*-------------------------
  IU GFAP
-------------------------*/

%RunOneClinicalTertileBiomarker(
    biomarker_short=IU_GFAP,
    biomarker_label=GFAP,
    datevar=IU_Plasma_date,
    biomarker_var=IU_GFAP_plasma_mean,
    cutoff_type=GFAP,
    tertile_fmt=gfap.,
    method_hy=AUTO,
    method_cog=AUTO,
    method_moca21=AUTO,
    method_moca26=AUTO,
    method_updrs=AUTO
);


/*-------------------------
  C2N AB40
-------------------------*/

%RunOneClinicalTertileBiomarker(
    biomarker_short=C2N_AB40,
    biomarker_label=A^{unicode 03B2}40 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB40,
    cutoff_type=AB40,
    tertile_fmt=ab40f.,
    method_hy=AUTO,
    method_cog=AUTO,
    method_moca21=AUTO,
    method_moca26=AUTO,
    method_updrs=AUTO
);


/*-------------------------
  C2N AB42
-------------------------*/

%RunOneClinicalTertileBiomarker(
    biomarker_short=C2N_AB42,
    biomarker_label=A^{unicode 03B2}42 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB42,
    cutoff_type=AB42,
    tertile_fmt=ab42f.,
    method_hy=AUTO,
    method_cog=AUTO,
    method_moca21=AUTO,
    method_moca26=AUTO,
    method_updrs=AUTO
);


/*-------------------------
  C2N AB42/40
-------------------------*/

%RunOneClinicalTertileBiomarker(
    biomarker_short=C2N_AB4240,
    biomarker_label=A^{unicode 03B2}42/A^{unicode 03B2}40 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB4240,
    cutoff_type=AB4240,
    tertile_fmt=ab4240f.,
    method_hy=AUTO,
    method_cog=AUTO,
    method_moca21=AUTO,
    method_moca26=AUTO,
    method_updrs=AUTO
);

/*===============================================================
  9. DISPLAY FINAL COMBINED CLINICAL REPORT
===============================================================*/

%WriteCombinedClinicalTertileRTF(
    title_suffix=for Sporadic PD
);


/*===============================================================
  10. OPTIONAL FINAL CLEANUP
===============================================================*/

proc datasets library=work nolist;
    delete Clinical_Tertile_Summary;
quit;

/*=========================================================
  sPD NSD+: Biomarker Tertiles and Progression Milestones

  Uses the prepared progression_analysis_v3 dataset to
  examine associations between plasma biomarker tertiles
  and time to progression milestones in sPD NSD+ participants.

  Generates Kaplan-Meier curves and survival summary tables
  using age- and sex-adjusted Cox proportional hazards or
  Weibull accelerated failure time models, selected according
  to the proportional hazards assumption test.

=========================================================*/


%let filepath = /*path where analysis data are located*/;
libname analysis "&filepath.";

%let report_date = 08DEC2025;
%let write_individual_outputs = 0;  /* 0 = suppress biomarker-specific summary reports */

ods escapechar="^";

/* Make a local copy of the prepared analysis dataset. */
data progression_analysis_v3;
    set analysis.progression_analysis_v3;
run;

/*---------------------------------------------------------
  Tertile formats
---------------------------------------------------------*/

proc format;

    value ptau
        1 = "< 0.131"
        2 = "0.131 - 0.174"
        3 = "> 0.174";

    value nfl
        1 = "< 9.877"
        2 = "9.877 - 14.538"
        3 = "> 14.538";

    value gfap
        1 = "< 64.478"
        2 = "64.478 - 105.341"
        3 = "> 105.341";

    value ab40f
        1 = "< 374.886"
        2 = "374.886 - 422.785"
        3 = "> 422.785";

    value ab42f
        1 = "< 40.219"
        2 = "40.219 - 45.870"
        3 = "> 45.870";

    value ab4240f
        1 = "< 0.104"
        2 = "0.104 - 0.111"
        3 = "> 0.111";

run;

/*---------------------------------------------------------
  Initialize combined all-biomarker tertile summary dataset

  SAS dataset names must be <= 32 characters.
---------------------------------------------------------*/

proc datasets library=work nolist;
    delete Combined_Tertile_Summary;
quit;

data Combined_Tertile_Summary;
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
        Beta_CI_3 $60
        SE_3 $20
        Z_3 $20
        P_3 $20
        HR_CI_3 $60
        P_Sig 8
    ;
    stop;
run;

/*---------------------------------------------------------
  Tertile categorical survival macro

  Key correction:
      - PH assumption check runs directly on &data.
      - Temp deletes baseline events only for KM/log-rank/modeling
      - AUTO chooses:
            PHREG   if no ASSESS PH p-value < 0.05
            LIFEREG if any ASSESS PH p-value < 0.05
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
                       across=3 opaque
                       titleattrs=(family="Arial");
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
      19. Compact model statistics table for PDF
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
%macro SurvTertileSummary(
    biomarker_label=,
    title_suffix=for Sporadic PD
);

    proc sql noprint;
        select memname
            into :_mslist separated by ' '
        from dictionary.tables
        where libname = "WORK"
          and upcase(substr(memname,1,3)) = "MS_"
        order by crdate;
    quit;

    %if %superq(_mslist)= %then %do;
        %put WARNING: No MS_ datasets found in WORK. Run %SurvTertile calls first.;
        %return;
    %end;

    data All_Milestones;
        set &_mslist.;
    run;

    data Summary_Data_raw;
        set All_Milestones;

        length SortOrder 8;

        if index(upcase(Parameter), "OVERALL") then SortOrder = 1;
        else if upcase(strip(Parameter)) = "T1 VS T2" then SortOrder = 2;
        else if upcase(strip(Parameter)) = "T1 VS T3" then SortOrder = 3;
        else if upcase(strip(Parameter)) = "T2 VS T3" then SortOrder = 4;
        else delete;
    run;

    /*-----------------------------------------------------
      Display rule for the final RTF:
        - Always retain the overall tertile row.
        - Retain T1 vs T2, T1 vs T3, and T2 vs T3 only
          when the overall p-value is < 0.05 for that
          biomarker and milestone.
    -----------------------------------------------------*/

    proc sql;
        create table Summary_Data_filtered as
        select a.*
        from Summary_Data_raw as a
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
            from Summary_Data_raw
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

        length Beta_CI_3 $60 SE_3 $20 Z_3 $20 P_3 $20 HR_CI_3 $60;

        P_Sig = (SortOrder = 1 and not missing(P_value) and P_value < 0.05);

        if not missing(Beta) and not missing(SE) then do;
            Beta_LCL = Beta - 1.96*SE;
            Beta_UCL = Beta + 1.96*SE;

            Beta_CI_3 = catx(" ",
                         put(Beta, 8.3),
                         cats("(",
                              put(Beta_LCL, 8.3), ", ",
                              put(Beta_UCL, 8.3), ")"));
        end;
        else Beta_CI_3 = "";

        if not missing(SE) then SE_3 = put(SE, 8.3);
        else SE_3 = "";

        if not missing(Z_stat) then Z_3 = put(Z_stat, 8.3);
        else Z_3 = "";

        if not missing(P_value) then do;
            if P_value < 0.001 then P_3 = "<0.001";
            else P_3 = put(P_value, 8.3);
        end;
        else P_3 = "";

        if not missing(HR) and not missing(HR_LCL) and not missing(HR_UCL) then do;
            HR_CI_3 = catx(" ",
                       put(HR, 8.3),
                       cats("(",
                            put(HR_LCL, 8.3), ", ",
                            put(HR_UCL, 8.3), ")"));
        end;
        else if not missing(HR_CI) then HR_CI_3 = HR_CI;
        else HR_CI_3 = "";

        keep Biomarker Milestone MilestoneOrder Method Parameter ClassVal
             Comparison_Type Comparison SortOrder GroupDF WaldStat
             Beta SE Z_stat P_value HR HR_LCL HR_UCL Beta_CI HR_CI
             Beta_CI_3 SE_3 Z_3 P_3 P_Sig HR_CI_3;
    run;


    /*-----------------------------------------------------
      Remove duplicate rows within this biomarker summary
    -----------------------------------------------------*/

    proc sort data=Summary_Data nodupkey;
        by Biomarker MilestoneOrder Milestone SortOrder Parameter Method;
    run;


    /*-----------------------------------------------------
      IMPORTANT FIX:
      Remove old rows for this biomarker before appending.
      This prevents duplicate rows when rerunning a biomarker.
    -----------------------------------------------------*/

    %if %sysfunc(exist(Combined_Tertile_Summary)) %then %do;

        data Combined_Tertile_Summary;
            set Combined_Tertile_Summary;

            if strip(upcase(Biomarker)) = strip(upcase("&biomarker_label.")) then delete;
        run;

    %end;


    /* Append only the fresh version for this biomarker */
    proc append base=Combined_Tertile_Summary
                data=Summary_Data
                force;
    run;


    /*-----------------------------------------------------
      Display biomarker-specific report
    -----------------------------------------------------*/

    proc sort data=Summary_Data;
        by Biomarker MilestoneOrder SortOrder;
    run;

    %if &write_individual_outputs = 1 %then %do;

    options
        orientation = landscape
        leftmargin  = 0.25in
        rightmargin = 0.25in
        topmargin   = 0.40in
        bottommargin= 0.40in
        nodate
        nonumber
    ;


    title "Summary of &biomarker_label. Tertile Survival Model Statistics - &title_suffix.";

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

        columns Biomarker MilestoneOrder Milestone Method SortOrder Parameter GroupDF WaldStat
                Beta_CI_3 SE_3 Z_3 P_3 HR_CI_3;

        define Biomarker      / group 'Biomarker';
        define MilestoneOrder / order noprint;
        define Milestone      / group 'Milestone';
        define Method         / display 'Model';
        define SortOrder      / order noprint;
        define Parameter      / display 'Parameter / Comparison';
        define GroupDF        / display 'DF' format=4.;
        define WaldStat       / display 'Wald Chi-Sq' format=10.4;
        define Beta_CI_3      / display 'Beta (95% CI)';
        define SE_3           / display 'SE';
        define Z_3            / display 'Z';
        define P_3            / display 'p-value';
        define HR_CI_3        / display 'HR / Derived HR (95% CI)';

        break after Milestone / skip;

    run;

    footnote;
    title;



    %end;

%mend SurvTertileSummary;

/*---------------------------------------------------------
  Main biomarker-level macro
---------------------------------------------------------*/

%macro RunOneTertileBiomarker(
    biomarker_short=,
    biomarker_label=,
    graph_label=,
    datevar=,
    biomarker_var=,
    cutoff_type=,
    tertile_fmt=,

    method_any=AUTO,
    method_wb=AUTO,
    method_mc=AUTO,
    method_cog=AUTO,
    method_auto=AUTO,
    method_fd=AUTO,
    method_adl=AUTO
);


    /* Graph-safe display label.
       Use (*ESC*) Unicode syntax for ODS Graphics. */
    %if %superq(graph_label)= %then %let graph_label=&biomarker_label.;



    %put NOTE: ============================================================;
    %put NOTE: Running tertile survival analysis for &biomarker_label.;
    %put NOTE: Date variable      = &datevar.;
    %put NOTE: Biomarker variable = &biomarker_var.;
    %put NOTE: Cutoff type        = &cutoff_type.;
    %put NOTE: ============================================================;


    /*-----------------------------------------------------
      Create biomarker-specific source dataset
    -----------------------------------------------------*/

    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data SurvData_&biomarker_short.;
        set progression_analysis_v3;
        where Analysis_Source="SP_TP"
          and Analysis_Biomarker="&biomarker_short.";
    run;



    /*-----------------------------------------------------
      Distribution check
    -----------------------------------------------------*/

    title "Distribution of Tertile &biomarker_label. and Age";

    proc freq data=SurvData_&biomarker_short.;
        tables Early_Tertile_Group / missing;
        format Early_Tertile_Group &tertile_fmt.;
    run;

    proc means data=SurvData_&biomarker_short. n mean std median q1 q3 min max maxdec=4;
        var Early_Biomarker_Value age;
    run;

    title;


    /*-----------------------------------------------------
      Initialize combined model stats for this biomarker

      Must include all variables, otherwise PROC APPEND with FORCE
      drops Beta, SE, P_value, HR, etc.
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete All_ModelStats MS_: All_Milestones Summary_Data Summary_Data_raw;
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
      Open biomarker-specific PDF
    -----------------------------------------------------*/

    options orientation=landscape;

    title "Tertile &biomarker_label. Survival Analysis for Sporadic PD NSD+";


    /*-----------------------------------------------------
      Run all milestones
    -----------------------------------------------------*/

    %SurvTertile(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_any,
        group=Early_Tertile_Group,
        time=TIME_ENRL_any,
        label=Any progression milestone,
        milestone_order=1,
        method=&method_any.,
        covars=age gen,
        grouplabel=&graph_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group subgroup_definition STATUS_any TIME_ENRL_any,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.
    );

    %SurvTertile(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_wb_any,
        group=Early_Tertile_Group,
        time=TIME_ENRL_wb_any,
        label=Walking and balance domain,
        milestone_order=2,
        method=&method_wb.,
        covars=age gen,
        grouplabel=&graph_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group subgroup_definition STATUS_wb_any TIME_ENRL_wb_any,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.
    );

    %SurvTertile(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_mc_any,
        group=Early_Tertile_Group,
        time=TIME_ENRL_mc_any,
        label=Motor complications domain,
        milestone_order=3,
        method=&method_mc.,
        covars=age gen,
        grouplabel=&graph_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group subgroup_definition STATUS_mc_any TIME_ENRL_mc_any,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.
    );

    %SurvTertile(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_cog_any,
        group=Early_Tertile_Group,
        time=TIME_ENRL_cog_any,
        label=Cognitive domain,
        milestone_order=4,
        method=&method_cog.,
        covars=age gen,
        grouplabel=&graph_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group subgroup_definition STATUS_cog_any TIME_ENRL_cog_any,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.
    );

    %SurvTertile(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_auto_any,
        group=Early_Tertile_Group,
        time=TIME_ENRL_auto_any,
        label=Autonomic dysfunction domain,
        milestone_order=5,
        method=&method_auto.,
        covars=age gen,
        grouplabel=&graph_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group subgroup_definition STATUS_auto_any TIME_ENRL_auto_any,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.
    );

    %SurvTertile(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_fd_any,
        group=Early_Tertile_Group,
        time=TIME_ENRL_fd_any,
        label=Functional dependence domain,
        milestone_order=6,
        method=&method_fd.,
        covars=age gen,
        grouplabel=&graph_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group subgroup_definition STATUS_fd_any TIME_ENRL_fd_any,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.
    );

    %SurvTertile(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_adl_any,
        group=Early_Tertile_Group,
        time=TIME_ENRL_adl_any,
        label=Activities of daily living domain,
        milestone_order=7,
        method=&method_adl.,
        covars=age gen,
        grouplabel=&graph_label.,
        groupfmt=&tertile_fmt.,
        title_suffix=for Sporadic PD,
        baseline_vars=Patno Early_Biomarker_Value Early_Tertile_Group subgroup_definition STATUS_adl_any TIME_ENRL_adl_any,
        stack_group=Early_Tertile_Group,
        stack_groupfmt=&tertile_fmt.,
        stack_grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.
    );

    title;


    /*-----------------------------------------------------
      Display biomarker-specific report summary
    -----------------------------------------------------*/

    %SurvTertileSummary(
        biomarker_label=&biomarker_label.,
        title_suffix=for Sporadic PD
    );


    /*-----------------------------------------------------
      Clean biomarker-specific datasets
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete mydata_&biomarker_short.
               EARLY_&biomarker_short.
               SurvData_&biomarker_short.
               All_ModelStats
               All_Milestones
               Summary_Data
               Summary_Data_filtered
               Summary_Data_raw
               MS_:;
    quit;

%mend RunOneTertileBiomarker;


/*---------------------------------------------------------
  Final combined all-biomarker report
---------------------------------------------------------*/

%macro WriteCombinedTertileBiomarkerRTF(
    title_suffix=for Sporadic PD
);

    proc sql noprint;
        select count(*)
            into :_n_combined trimmed
        from Combined_Tertile_Summary;
    quit;

    %if %sysevalf(&_n_combined. = 0) %then %do;
        %put WARNING: Combined_Tertile_Summary is empty. No combined report displayed.;
        %return;
    %end;

    /*-----------------------------------------------------
      Apply preferred biomarker ordering for combined report

      Order:
          1. pTau217
          2. NfL
          3. GFAP
          4. C2N AB40
          5. C2N AB42
          6. C2N AB42/40 ratio
    -----------------------------------------------------*/

    data Combined_Tertile_Summary;
        set Combined_Tertile_Summary;

        /* Final table includes the six requested progression milestones only. */
        if MilestoneOrder in (1, 2, 3, 4, 5, 6);

        length BiomarkerOrder 8;

        if upcase(strip(Biomarker)) in 
            ("PTAU217 PLASMA", "LILLY PTAU217", "PTAU217", "PTAU") then BiomarkerOrder = 1;

        else if upcase(strip(Biomarker)) in 
            ("IU NFL", "NFL") then BiomarkerOrder = 2;

        else if upcase(strip(Biomarker)) in 
            ("IU GFAP", "GFAP") then BiomarkerOrder = 3;

        else if upcase(strip(Biomarker)) in 
            ("A^{UNICODE 03B2}40 PLASMA", "C2N AB40", "AB40") then BiomarkerOrder = 4;

        else if upcase(strip(Biomarker)) in 
            ("A^{UNICODE 03B2}42 PLASMA", "C2N AB42", "AB42") then BiomarkerOrder = 5;

        else if upcase(strip(Biomarker)) in 
            ("A^{UNICODE 03B2}42/A^{UNICODE 03B2}40 PLASMA", "C2N AB42/40", "C2N AB4240", "AB42/40", "AB4240") then BiomarkerOrder = 6;

        else BiomarkerOrder = 99;
    run;

    proc sort data=Combined_Tertile_Summary nodupkey;
        by BiomarkerOrder Biomarker MilestoneOrder Milestone SortOrder Parameter Method;
    run;

    proc sort data=Combined_Tertile_Summary;
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
        "^{style [font_weight=bold] Supplement Table.} Association between baseline plasma biomarker tertiles and time-to-event progression milestones in";
    title2 j=left h=9pt
        "Sporadic PD (NSD+). Cox proportional hazards (CoxPH) or Weibull accelerated failure time (AFT) models were selected automatically by the";
    title3 j=left h=9pt
        "proportional hazards assumption test. All models adjusted for age and sex. Significant overall p-values (p < 0.05) are shown in bold red.";

    footnote1 j=left h=7pt
        "Overall row reports the Type III omnibus test for the tertile biomarker.";
    footnote2 j=left h=7pt
        "Pairwise rows (T1 vs T2, T1 vs T3, and T2 vs T3) are shown only when the overall tertile p-value is < 0.05. Models adjusted for age and sex.";
    footnote3 j=left h=7pt
        "Cox PH models report HR. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Values rounded to 3 decimals.";

    proc report data=Combined_Tertile_Summary nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=0 cellpadding=2 frame=hsides rules=groups]
        style(header)=[font_size=7pt font_weight=bold background=white just=center]
        style(column)=[font_size=7pt vjust=middle];

        columns BiomarkerOrder Biomarker MilestoneOrder Milestone SortOrder
                Parameter Method GroupDF WaldStat Beta_CI_3 SE_3 Z_3 P_Sig P_3 HR_CI_3;

        define BiomarkerOrder / order noprint order=internal;
        define Biomarker      / order 'Biomarker' order=data
            style(column)=[cellwidth=0.95in font_weight=bold just=left vjust=top];
        define MilestoneOrder / order noprint order=internal;
        define Milestone      / order 'Milestone' order=data
            style(column)=[cellwidth=1.45in just=left];
        define SortOrder      / order noprint order=internal;

        define Parameter      / display 'Parameter / Comparison'
            style(column)=[cellwidth=0.85in just=left];
        define Method         / display 'Model'
            style(column)=[cellwidth=0.45in just=center];
        define GroupDF        / display 'DF' format=4.
            style(column)=[cellwidth=0.30in just=center];
        define WaldStat       / display 'Wald Chi-Sq' format=8.3
            style(column)=[cellwidth=0.55in just=center];
        define Beta_CI_3      / display 'Beta (95% CI)'
            style(column)=[cellwidth=1.30in just=center];
        define SE_3           / display 'SE'
            style(column)=[cellwidth=0.40in just=center];
        define Z_3            / display 'Z'
            style(column)=[cellwidth=0.40in just=center];
        define P_Sig          / display noprint;
        define P_3            / display 'p'
            style(header)=[font_style=italic]
            style(column)=[cellwidth=0.45in just=center];
        define HR_CI_3        / display 'HR (95% CI)'
            style(column)=[cellwidth=1.30in just=center];

        compute P_3;
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



%mend WriteCombinedTertileBiomarkerRTF;


/*---------------------------------------------------------
  Run all biomarkers
---------------------------------------------------------*/

/* Lilly pTau217 */
%RunOneTertileBiomarker(
    biomarker_short=Lilly_pTau217,
    biomarker_label=pTau217 plasma,
    graph_label=pTau217,
    datevar=Lilly_ptau217p_date,
    biomarker_var=Lilly_ptau217p,
    cutoff_type=PTAU,
    tertile_fmt=ptau.,
    method_any=AUTO,
    method_wb=AUTO,
    method_mc=AUTO,
    method_cog=AUTO,
    method_auto=AUTO,
    method_fd=AUTO,
    method_adl=AUTO
);


/* IU NfL */
%RunOneTertileBiomarker(
    biomarker_short=IU_NFL,
    biomarker_label=NFL plasma,
    graph_label=NFL,
    datevar=IU_Plasma_date,
    biomarker_var=IU_NFL_plasma_mean,
    cutoff_type=NFL,
    tertile_fmt=nfl.,
    method_any=AUTO,
    method_wb=AUTO,
    method_mc=AUTO,
    method_cog=AUTO,
    method_auto=AUTO,
    method_fd=AUTO,
    method_adl=AUTO
);


/* IU GFAP */
%RunOneTertileBiomarker(
    biomarker_short=IU_GFAP,
    biomarker_label=GFAP plasma,
    graph_label=GFAP,
    datevar=IU_Plasma_date,
    biomarker_var=IU_GFAP_plasma_mean,
    cutoff_type=GFAP,
    tertile_fmt=gfap.,
    method_any=AUTO,
    method_wb=AUTO,
    method_mc=AUTO,
    method_cog=AUTO,
    method_auto=AUTO,
    method_fd=AUTO,
    method_adl=AUTO
);


/* C2N AB40 */
%RunOneTertileBiomarker(
    biomarker_short=C2N_AB40,
    biomarker_label=A^{unicode 03B2}40 plasma,
    graph_label=A(*ESC*){unicode '03B2'x}40,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB40,
    cutoff_type=AB40,
    tertile_fmt=ab40f.,
    method_any=AUTO,
    method_wb=AUTO,
    method_mc=AUTO,
    method_cog=AUTO,
    method_auto=AUTO,
    method_fd=AUTO,
    method_adl=AUTO
);


/* C2N AB42 */
%RunOneTertileBiomarker(
    biomarker_short=C2N_AB42,
    biomarker_label=A^{unicode 03B2}42 plasma,
    graph_label=A(*ESC*){unicode '03B2'x}42,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB42,
    cutoff_type=AB42,
    tertile_fmt=ab42f.,
    method_any=AUTO,
    method_wb=AUTO,
    method_mc=AUTO,
    method_cog=AUTO,
    method_auto=AUTO,
    method_fd=AUTO,
    method_adl=AUTO
);


/* C2N AB42/40 */
%RunOneTertileBiomarker(
    biomarker_short=C2N_AB4240,
    biomarker_label=A^{unicode 03B2}42/40 plasma,
    graph_label=A(*ESC*){unicode '03B2'x}42/A(*ESC*){unicode '03B2'x}40,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB4240,
    cutoff_type=AB4240,
    tertile_fmt=ab4240f.,
    method_any=AUTO,
    method_wb=AUTO,
    method_mc=AUTO,
    method_cog=AUTO,
    method_auto=AUTO,
    method_fd=AUTO,
    method_adl=AUTO
);


/*---------------------------------------------------------
  Display final combined report
---------------------------------------------------------*/

%WriteCombinedTertileBiomarkerRTF(
    title_suffix=for Sporadic PD
);


/*---------------------------------------------------------
  Final cleanup
---------------------------------------------------------*/
/**/
/*proc datasets library=work nolist;*/
/*    delete BL1 Combined_Tertile_Summary combined:;*/
/*quit;*/





/*---------------------------------------------------------
  Helper macro: Display one KM plot with adjusted p-values
  Uses original SurvivalPlot=survsaa3 style
---------------------------------------------------------*/

%macro ExportOneKMpng(
    data=,
    outcome=,
    time=,
    group=,
    groupfmt=,
    grouplabel=,
    biomarker_label=,
    label=,
    title_suffix=for Sporadic PD,
    method=PHREG,
    covars=age gen,
    imgname=,
    max_time=11,
    width=6.0in,
    height=3.8in
);

    %local _imgname _dslabel _biolabel;
    %local _method _modelLabel;
    %local pval pval2 pval2n p12 p13 p23;
    %local dsid hasLCLM rc;

    %let _method = %upcase(&method.);
    %let _modelLabel = %sysfunc(ifc(&_method. = LIFEREG, AFT, CoxPH));

    %let pval=.;
    %let pval2=.;
    %let pval2n=.;
    %let p12=.;
    %let p13=.;
    %let p23=.;
    %let hasLCLM=0;

    /* Clean image name */
    %if %superq(imgname)= %then %do;

        %let _dslabel = %sysfunc(compress(&label.,
            abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789, k));

        %let _dslabel = %substr(&_dslabel., 1,
            %sysfunc(min(29, %length(&_dslabel.))));

        %let _biolabel = %sysfunc(compress(&biomarker_label.,
            abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789, k));

        %let _biolabel = %substr(&_biolabel., 1,
            %sysfunc(min(20, %length(&_biolabel.))));

        %let _imgname = &_biolabel._&_dslabel._KM;

    %end;
    %else %do;
        %let _imgname = &imgname.;
    %end;

    %put NOTE: ------------------------------------------------------------;
    %put NOTE: Exporting selected KM plot using original KM style: &_imgname.;
    %put NOTE: Outcome=&outcome. Time=&time. Group=&group. Method=&_method.;
    %put NOTE: ------------------------------------------------------------;


    /*-----------------------------------------------------
      1. Analysis dataset
    -----------------------------------------------------*/

    data _km_temp;
        set &data.;

        if missing(&group.) then delete;
        if missing(&time.) then delete;
        if missing(&outcome.) then delete;

        if &outcome. = 1 and &time. = 0 then delete;
    run;


    /*-----------------------------------------------------
      2. Pairwise log-rank tests
    -----------------------------------------------------*/

    ods exclude all;
    ods noresults;

    ods output HomTests = _km_LR12(where=(upcase(Test)="LOG-RANK"));
    proc lifetest data=_km_temp(where=(&group. in (1,2)));
        time &time.*&outcome.(0);
        strata &group. / test=logrank;
    run;

    ods output HomTests = _km_LR13(where=(upcase(Test)="LOG-RANK"));
    proc lifetest data=_km_temp(where=(&group. in (1,3)));
        time &time.*&outcome.(0);
        strata &group. / test=logrank;
    run;

    ods output HomTests = _km_LR23(where=(upcase(Test)="LOG-RANK"));
    proc lifetest data=_km_temp(where=(&group. in (2,3)));
        time &time.*&outcome.(0);
        strata &group. / test=logrank;
    run;

    ods results;
    ods exclude none;


    /*-----------------------------------------------------
      3. Overall KM + log-rank
      This is the original style.
    -----------------------------------------------------*/

    ods graphics on;
    ods select none;

    proc lifetest data=_km_temp
        plots=survival(cb=hw atrisk=0 to &max_time. by 1)
        confband=hw bandmin=0 bandmax=&max_time.;

        time &time.*&outcome.(0);
        strata &group. / test=logrank;
        id patno;
        label &time. = "Year";

        ods output HomTests    = _km_HomTests
                   SurvivalPlot = _km_survsaa3;
    run;

    ods select all;


    /*-----------------------------------------------------
      4. Overall and pairwise adjusted model p-values
    -----------------------------------------------------*/

    %if &_method. = PHREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ModelANOVA = _km_modelanova;

        proc phreg data=_km_temp;
            class &group.(ref='3') gen(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &group. / type3;
        run;

        ods results;
        ods exclude none;

        data _null_;
            set _km_modelanova;

            if index(upcase("&group."), upcase(strip(Effect))) > 0
               and not missing(Effect) then do;
                call symputx("pval2",
             put(rounde(ProbChiSq, 0.001), pvalue6.3));
                call symputx("pval2n", ProbChiSq);
            end;
        run;

        ods exclude all;
        ods noresults;

        ods output ContrastTest = _km_pairwise;

        proc phreg data=_km_temp;
            class &group.(ref='3') gen(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &group.;

            contrast "1 vs 2" &group.  1 -1;
            contrast "1 vs 3" &group.  1  0;
            contrast "2 vs 3" &group.  0  1;
        run;

        ods results;
        ods exclude none;

        data _km_pairwise2;
            set _km_pairwise;

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
    %else %if &_method. = LIFEREG %then %do;

        ods exclude all;
        ods noresults;

        ods output Type3Analysis = _km_type3
                   Estimates     = _km_pairwise_raw;

        proc lifereg data=_km_temp;
            class &group. gen;
            model &time.*&outcome.(0) = &covars. &group. / dist=weibull;

            estimate "1 vs 2" &group.  1 -1  0;
            estimate "1 vs 3" &group.  1  0 -1;
            estimate "2 vs 3" &group.  0  1 -1;
        run;

        ods results;
        ods exclude none;

        proc datasets lib=work nolist;
            modify _km_type3;
            rename Variable = Effect;
        quit;

        data _null_;
            set _km_type3;

            if index(upcase("&group."), upcase(strip(Effect))) > 0
               and not missing(Effect) then do;
                call symputx("pval2",
             put(rounde(ProbChiSq, 0.001), pvalue6.3));
                call symputx("pval2n", ProbChiSq);
            end;
        run;

        data _km_pairwise2;
            set _km_pairwise_raw(rename=(Label=Contrast));

            length Comparison $12;

            _c = upcase(strip(compbl(Contrast)));

            if      _c = "1 VS 2" then Comparison = "T1 vs T2";
            else if _c = "1 VS 3" then Comparison = "T1 vs T3";
            else if _c = "2 VS 3" then Comparison = "T2 vs T3";
            else delete;

            if not missing(ProbChiSq) then CoxPH_p = ProbChiSq;
            else if not missing(ProbZ) then CoxPH_p = ProbZ;
            else CoxPH_p = .;

            keep Comparison CoxPH_p;
        run;

    %end;
    %else %do;

        %put WARNING: method=&method. not recognized. Defaulting p-values to missing.;

        data _km_pairwise2;
            length Comparison $12 CoxPH_p 8;
            Comparison = "T1 vs T2"; CoxPH_p = .; output;
            Comparison = "T1 vs T3"; CoxPH_p = .; output;
            Comparison = "T2 vs T3"; CoxPH_p = .; output;
        run;

    %end;


    /* Overall log-rank p-value */
    data _null_;
        set _km_HomTests;
        where upcase(Test) = "LOG-RANK";
        call symputx("pval", put(rounde(ProbChiSq, 0.001), pvalue6.3));
    run;


    /* Pairwise adjusted p-value macro variables */
    proc sql noprint;
        select put(rounde(CoxPH_p, 0.001), pvalue6.3) into :p12 trimmed
        from _km_pairwise2
        where Comparison = "T1 vs T2";

        select put(rounde(CoxPH_p, 0.001), pvalue6.3) into :p13 trimmed
        from _km_pairwise2
        where Comparison = "T1 vs T3";

        select put(rounde(CoxPH_p, 0.001), pvalue6.3) into :p23 trimmed
        from _km_pairwise2
        where Comparison = "T2 vs T3";
    quit;

    %if %superq(pval)=   %then %let pval=.;
    %if %superq(pval2)=  %then %let pval2=.;
    %if %superq(pval2n)= %then %let pval2n=.;
    %if %superq(p12)=    %then %let p12=.;
    %if %superq(p13)=    %then %let p13=.;
    %if %superq(p23)=    %then %let p23=.;


    /*-----------------------------------------------------
      5. Detect CI band
    -----------------------------------------------------*/

    %if %sysfunc(exist(_km_survsaa3)) %then %do;

        %let dsid = %sysfunc(open(_km_survsaa3));

        %if &dsid. %then %do;
            %let hasLCLM = %sysfunc(varnum(&dsid., LCLM));
            %let rc = %sysfunc(close(&dsid.));
        %end;

    %end;
    %else %do;
        %put ERROR: _km_survsaa3 was not created. Cannot export KM plot &_imgname..;
        %goto cleanup;
    %end;


    /*-----------------------------------------------------
      6. Display plot with original SGPLOT style
    -----------------------------------------------------*/

    ods escapechar='^';

    ods listing image_dpi=300;

ods graphics on / reset=index
                  width=&width.
                  height=&height.
                  border;

    %if %sysevalf(%superq(pval2n)=, boolean) or %sysevalf(&pval2n >= 0.05) %then %do;
        footnote1 height=10pt
            "&_modelLabel p-value | Overall Interaction P-value = &pval2";
    %end;
    %else %do;
        footnote1 height=10pt
            "&_modelLabel p-values | Overall = &pval2 | T1 vs T2 = &p12 | T1 vs T3 = &p13 | T2 vs T3 = &p23";
    %end;

/*    title height=10pt bold "&label. &title_suffix";*/

    proc sgplot data=_km_survsaa3;
        styleattrs datacolors=(blue red green)
                   datacontrastcolors=(blue red green);

        step x=Time y=Survival /
            group=stratumnum
            name="km"
            lineattrs=(thickness=1.2);

        scatter x=Time y=Censored /
            group=stratumnum
            markerattrs=(symbol=plus size=5);

        band x=Time upper=HW_UCL lower=HW_LCL /
            group=stratumnum
            transparency=0.80;

        %if &hasLCLM %then %do;
            band x=Time lower=LCLM upper=UCLM /
                group=stratumnum
                transparency=0.80;
        %end;

        xaxistable Atrisk /
            location=inside
            x=tatrisk
            class=stratumnum
            colorgroup=stratumnum
            valueattrs=(size=12pt);

        xaxis values=(0 to &max_time. by 1)
              label="Year"
              labelattrs=(size=12pt weight=bold)
              valueattrs=(size=12pt);

        yaxis label="Survival Probability"
              labelattrs=(size=12pt weight=bold)
              valueattrs=(size=12pt);

        format stratumnum &groupfmt.;
        label  stratumnum="&grouplabel";

        keylegend "km" /
            title=""
            location=outside
            position=bottom
            across=3
            opaque
            valueattrs=(size=11pt)
            titleattrs=(family="Arial" size=11pt weight=normal);
    run;

    ods graphics / reset;

    footnote;
/*    title;*/


    %cleanup:

    proc datasets library=work nolist;
        delete _km_temp
               _km_LR12
               _km_LR13
               _km_LR23
               _km_HomTests
               _km_survsaa3
               _km_modelanova
               _km_type3
               _km_pairwise
               _km_pairwise_raw
               _km_pairwise2;
    quit;

%mend ExportOneKMpng;
/*---------------------------------------------------------
  Biomarker-level standalone KM image macro

  Purpose:
    - Build the same early tertile SurvData dataset
    - Call %ExportOneKMpng for selected Kaplan-Meier plots

  plot_set options:
    COG      = Cognitive domain only
    SELECTED = Cognitive + Motor complications + Autonomic dysfunction
               + Functional dependence
    ALL      = Any progression + Walking/balance + Motor complications
               + Cognitive + Autonomic dysfunction + Functional dependence
               + ADL
---------------------------------------------------------*/

%macro RunOneTertileKMImages(
    biomarker_short=,
    biomarker_label=,
    graph_label=,
    datevar=,
    biomarker_var=,
    cutoff_type=,
    tertile_fmt=,

    method_any=PHREG,
    method_wb=PHREG,
    method_mc=PHREG,
    method_cog=PHREG,
    method_auto=PHREG,
    method_fd=PHREG,
    method_adl=PHREG,
    plot_set=SELECTED
);

    /* Graph-safe display label for ODS Graphics. */
    %if %superq(graph_label)= %then %let graph_label=&biomarker_label.;

    %put NOTE: ============================================================;
    %put NOTE: Creating selected KM plots for &biomarker_label.;
    %put NOTE: Plot set = &plot_set.;
    %put NOTE: ============================================================;


    /*-----------------------------------------------------
      Create biomarker-specific source dataset
    -----------------------------------------------------*/

    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data SurvData_KM_&biomarker_short.;
        set progression_analysis_v3;
        where Analysis_Source="SP_TP"
          and Analysis_Biomarker="&biomarker_short.";
    run;



    /*-----------------------------------------------------
      Export KM PNGs

      COG:
        Cognitive only

      SELECTED:
        Cognitive
        Motor complications
        Autonomic dysfunction
        Functional dependence

      ALL:
        Any progression
        Walking and balance
        Motor complications
        Cognitive
        Autonomic dysfunction
        Functional dependence
        Activities of daily living
    -----------------------------------------------------*/


    /* Cognitive domain */
%let plot_set = %upcase(&plot_set.);

/* Cognitive domain */
%if &plot_set. = COG or &plot_set. = SELECTED or &plot_set. = ALL %then %do;

    %ExportOneKMpng(
        data=SurvData_KM_&biomarker_short.,
        outcome=STATUS_cog_any,
        time=TIME_ENRL_cog_any,
        group=Early_Tertile_Group,
        groupfmt=&tertile_fmt.,
        grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.,
        label=Cognitive domain,
        method=&method_cog.,
        covars=age gen,
        max_time=11,
        width=6.0in,
        height=3.5in
    );

%end;


/* Selected manuscript domains */
%if &plot_set. = SELECTED or &plot_set. = ALL %then %do;

    %ExportOneKMpng(
        data=SurvData_KM_&biomarker_short.,
        outcome=STATUS_mc_any,
        time=TIME_ENRL_mc_any,
        group=Early_Tertile_Group,
        groupfmt=&tertile_fmt.,
        grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.,
        label=Motor complications domain,
        method=&method_mc.,
        covars=age gen,
        max_time=11,
        width=6.0in,
        height=3.5in
    );

    %ExportOneKMpng(
        data=SurvData_KM_&biomarker_short.,
        outcome=STATUS_auto_any,
        time=TIME_ENRL_auto_any,
        group=Early_Tertile_Group,
        groupfmt=&tertile_fmt.,
        grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.,
        label=Autonomic dysfunction domain,
        method=&method_auto.,
        covars=age gen,
        max_time=11,
        width=6.0in,
        height=3.5in
    );

    %ExportOneKMpng(
        data=SurvData_KM_&biomarker_short.,
        outcome=STATUS_fd_any,
        time=TIME_ENRL_fd_any,
        group=Early_Tertile_Group,
        groupfmt=&tertile_fmt.,
        grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.,
        label=Functional dependence domain,
        method=&method_fd.,
        covars=age gen,
        max_time=11,
        width=6.0in,
        height=3.5in
    );

%end;


/* Additional domains only if plot_set=ALL */
%if &plot_set. = ALL %then %do;

    %ExportOneKMpng(
        data=SurvData_KM_&biomarker_short.,
        outcome=STATUS_any,
        time=TIME_ENRL_any,
        group=Early_Tertile_Group,
        groupfmt=&tertile_fmt.,
        grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.,
        label=Any progression milestone,
        method=&method_any.,
        covars=age gen,
        max_time=11,
        width=6.0in,
        height=3.5in
    );

    %ExportOneKMpng(
        data=SurvData_KM_&biomarker_short.,
        outcome=STATUS_wb_any,
        time=TIME_ENRL_wb_any,
        group=Early_Tertile_Group,
        groupfmt=&tertile_fmt.,
        grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.,
        label=Walking and balance domain,
        method=&method_wb.,
        covars=age gen,
        max_time=11,
        width=6.0in,
        height=3.5in
    );

    %ExportOneKMpng(
        data=SurvData_KM_&biomarker_short.,
        outcome=STATUS_adl_any,
        time=TIME_ENRL_adl_any,
        group=Early_Tertile_Group,
        groupfmt=&tertile_fmt.,
        grouplabel=&graph_label.,
        biomarker_label=&biomarker_label.,
        label=Activities of daily living domain,
        method=&method_adl.,
        covars=age gen,
        max_time=11,
        width=6.0in,
        height=3.5in
    );

%end;

    /*-----------------------------------------------------
      Cleanup
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete mydata_km_&biomarker_short.
               early_km_&biomarker_short.
               SurvData_KM_&biomarker_short.;
    quit;

%mend RunOneTertileKMImages;

/*---------------------------------------------------------
  Selected Kaplan-Meier plots with adjusted p-values
---------------------------------------------------------*/

/* Display the selected Kaplan-Meier plots. */

ods _all_ close;
ods listing;
ods graphics on / reset=all;


/* Lilly pTau217 */
%RunOneTertileKMImages(
    biomarker_short=Lilly_pTau217,
    biomarker_label=pTau217 plasma,
    graph_label=pTau217,
    datevar=Lilly_ptau217p_date,
    biomarker_var=Lilly_ptau217p,
    cutoff_type=PTAU,
    tertile_fmt=ptau.,
    method_cog=PHREG,
    method_mc=PHREG,
    method_auto=PHREG,
    method_fd=LIFEREG,
    plot_set=SELECTED
);


/* IU NfL */
%RunOneTertileKMImages(
    biomarker_short=IU_NFL,
    biomarker_label=NFL plasma,
    graph_label=NFL,
    datevar=IU_Plasma_date,
    biomarker_var=IU_NFL_plasma_mean,
    cutoff_type=NFL,
    tertile_fmt=nfl.,
    method_cog=PHREG,
    method_mc=PHREG,
    method_auto=PHREG,
    method_fd=LIFEREG,
    plot_set=SELECTED
);


/* IU GFAP */
%RunOneTertileKMImages(
    biomarker_short=IU_GFAP,
    biomarker_label=GFAP plasma,
    graph_label=GFAP,
    datevar=IU_Plasma_date,
    biomarker_var=IU_GFAP_plasma_mean,
    cutoff_type=GFAP,
    tertile_fmt=gfap.,
    method_cog=PHREG,
    method_mc=PHREG,
    method_auto=PHREG,
    method_fd=LIFEREG,
    plot_set=SELECTED
);


/* C2N AB42/40 */
%RunOneTertileKMImages(
    biomarker_short=C2N_AB4240,
    biomarker_label=AB42/40 plasma,
    graph_label=A(*ESC*){unicode '03B2'x}42/A(*ESC*){unicode '03B2'x}40,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB4240,
    cutoff_type=AB4240,
    tertile_fmt=ab4240f.,
    method_cog=PHREG,
    method_mc=PHREG,
    method_auto=PHREG,
    method_fd=PHREG,
    plot_set=SELECTED
);


ods _all_ close;
ods listing;



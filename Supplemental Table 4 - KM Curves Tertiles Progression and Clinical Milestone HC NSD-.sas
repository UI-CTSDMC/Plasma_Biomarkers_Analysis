/*=========================================================
  HC NSD-: Biomarker Tertiles and Progression/Clinical Milestones

  Uses the prepared progression_analysis_v3 and
  clinical_analysis_v3 datasets to examine associations
  between pTau217, NFL, and GFAP tertiles and time to
  progression and clinical milestones in HC NSD- participants.

  Generates Kaplan-Meier curves and a combined survival
  summary table using age- and sex-adjusted Cox proportional
  hazards or Weibull accelerated failure time models,
  selected according to the proportional hazards test.
=========================================================*/
%let filepath = /*path where analysis data are located*/;
libname analysis "&filepath.";


/*****************************************************************************************
  Reports and plots:
      - Analyses and plots for all biomarkers and successful milestones
      - One combined summary table with N, Events, and Censored counts

  Output order:
      1. pTau217 plasma
            a. Progression milestones
            b. Clinical milestones
      2. NFL
            a. Progression milestones
            b. Clinical milestones
      3. GFAP
            a. Progression milestones
            b. Clinical milestones

  Model-selection rule:
      - Build Temp for each biomarker/outcome.
      - Temp removes missing group/time/outcome and baseline events:
            outcome=1 and time=0
      - Run ASSESS PH on Temp.
      - Use Cox PH if no ASSESS PH p-value < 0.05.
      - Use Weibull AFT if any ASSESS PH p-value < 0.05.
      - Fixed seed used for reproducible ASSESS PH p-values.

  Skip rule:
      - If a milestone has insufficient usable records, no events, fewer than 2
        tertile groups, or does not produce model-statistics rows, it is skipped
        and is not appended to the final table.

  Required prepared analysis datasets:
      progression_analysis_v3
      clinical_analysis_v3
*****************************************************************************************/

%let report_date = 08DEC2025;

ods escapechar="^";

/* Make local copies of the prepared analysis datasets. */
data progression_analysis_v3;
    set analysis.progression_analysis_v3;
run;

data clinical_analysis_v3;
    set analysis.clinical_analysis_v3;
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


/*===============================================================
  1. INITIALIZE FINAL COMBINED TABLE DATASET
===============================================================*/

proc datasets library=work nolist;
    delete Final_HCNSD_Tertile_Summary
           Final_HCNSD_Tertile_Summary_RTF
           _Overall_Significance;
quit;

data Final_HCNSD_Tertile_Summary;
    length
        Biomarker $100
        AnalysisType $30
        Milestone $120
        Method $10
        Parameter $80
        Comparison $12
        BiomarkerOrder 8
        AnalysisOrder 8
        MilestoneOrder 8
        SortOrder 8
        N_Used 8
        N_Events 8
        N_Censored 8
        GroupDF 8
        WaldStat 8
        Beta 8
        SE 8
        Z_stat 8
        P_value 8
        HR 8
        HR_LCL 8
        HR_UCL 8
        Beta_CI_4 $70
        SE_4 $20
        Z_4 $20
        P_4 $20
        HR_CI_4 $70
    ;
    stop;
run;


/*===============================================================
  2. SMALL UTILITY MACRO: COUNT OBSERVATIONS
===============================================================*/

%macro _get_nobs(ds, mvar);
    %if %symexist(&mvar.) %then %do;
        %let &mvar = 0;
    %end;
    %else %do;
        %global &mvar;
        %let &mvar = 0;
    %end;

    %if %sysfunc(exist(&ds.)) %then %do;
        %local _dsid _rc;
        %let _dsid = %sysfunc(open(&ds.));
        %if &_dsid. %then %do;
            %let &mvar = %sysfunc(attrn(&_dsid., nlobs));
            %let _rc = %sysfunc(close(&_dsid.));
        %end;
    %end;
%mend _get_nobs;


/*===============================================================
  3. CREATE EARLY HC NSD- BIOMARKER DATASET
===============================================================*/

%macro MakeEarlyHCNSDMinusBiomarker(
    biomarker_short=,
    biomarker_label=,
    biomarker_order=,
    datevar=,
    biomarker_var=,
    cutoff_type=
);
    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data Early_&biomarker_short.;
        set progression_analysis_v3;
        where Analysis_Source="HC_T"
          and Analysis_Biomarker="&biomarker_short.";
    run;
    proc sort data=Early_&biomarker_short. nodupkey; by Patno; run;
%mend MakeEarlyHCNSDMinusBiomarker;


/*===============================================================
  4. CREATE PROGRESSION SURVIVAL DATASET
===============================================================*/

%macro MakeProgressionSurvData(
    biomarker_short=
);
    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data ProgressionSurv_&biomarker_short.;
        set progression_analysis_v3;
        where Analysis_Source="HC_T"
          and Analysis_Biomarker="&biomarker_short.";
    run;
%mend MakeProgressionSurvData;


/*===============================================================
  5. CREATE CLINICAL SURVIVAL ENDPOINT DATASETS
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


%macro MakeClinicalSurvData(
    biomarker_short=
);
    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data ClinicalSurv_&biomarker_short.;
        set clinical_analysis_v3;
        where Analysis_Source="HC_T"
          and Analysis_Biomarker="&biomarker_short.";
    run;
%mend MakeClinicalSurvData;


/*===============================================================
  6. CORE MODEL MACRO
===============================================================*/

%macro RunOneTertileModel(
    data=,
    outcome=,
    time=,
    group=Early_Tertile_Group,
    covars=age gen,
    biomarker=,
    biomarker_order=,
    analysis_type=,
    analysis_order=,
    milestone=,
    milestone_order=,
    group_label=,
    groupfmt=
);

    %local _n_used _n_events _n_censored _n_groups _method _ph_fail _ph_pvar _n_type3 _n_est _n_final _bl_list;
    %local pval pval2 pval2n p12 p13 p23 _modelLabel;
    %local dsid hasLCLM rc;
    %local sigma;

    %let _method     = AUTO;
    %let _ph_fail    = 0;
    %let _ph_pvar    = ;
    %let _n_censored = 0;
    %let sigma       = .;

    proc datasets library=work nolist;
        delete Temp
               _ph_assess_raw_ _ph_assess_
               _type3 _type3b _est _est2 _overall _pairs _model_params _final_one
               _pdf_modelstats _pdf_modelstats_fmt _pe _covb _aft_outest
               HomTests LR12 LR13 LR23 temp_lr _pair_model PVALS_pairs PVALS_overall PVALS;
    quit;

    /*-----------------------------------------------------
      Build Temp:
        - same dataset for PH check and modeling
        - removes missing group/time/outcome
        - removes baseline events
    -----------------------------------------------------*/

    data Temp;
        set &data.;
        if &outcome. = 1 and &time. = 0 then delete;
    run;

    proc sql noprint;
        select count(*),
               sum(&outcome. = 1),
               count(distinct &group.)
            into :_n_used trimmed,
                 :_n_events trimmed,
                 :_n_groups trimmed
        from Temp;
    quit;

    %if %superq(_n_used)= %then %let _n_used=0;
    %if %superq(_n_events)= %then %let _n_events=0;
    %if %superq(_n_groups)= %then %let _n_groups=0;

    %let _n_censored = %eval(&_n_used. - &_n_events.);

    %if %sysevalf(&_n_used. < 5 or &_n_events. < 1 or &_n_groups. < 2) %then %do;
        %put WARNING: Skipping &biomarker. - &analysis_type. - &milestone.: insufficient data. N=&_n_used Events=&_n_events Censored=&_n_censored Groups=&_n_groups.;
        %goto cleanup;
    %end;

    /*-----------------------------------------------------
      AUTO model selection using ASSESS PH on Temp
    -----------------------------------------------------*/

    ods exclude all;
    ods noresults;

    ods output ProportionalHazardsSupTest = _ph_assess_raw_;

    proc phreg data=Temp  namelen=100;
        class &group.(ref='3') gen(ref='1') / param=ref;
        model &time.*&outcome.(0) = &group. &covars. ;
        assess ph / resample seed=1472160691;
    run;

    ods results;
    ods exclude none;

    %if %sysfunc(exist(_ph_assess_raw_)) %then %do;

        proc sql noprint;
            select name into :_ph_pvar trimmed
            from dictionary.columns
            where libname = "WORK"
              and memname = "_PH_ASSESS_RAW_"
              and upcase(name) in ("PVALUE", "PROBCHISQ", "PRCHISQ", "PROBZ");
        quit;

        %if %superq(_ph_pvar)= %then %do;
            %put WARNING: Could not identify ASSESS PH p-value variable for &biomarker. - &milestone.. Defaulting to PHREG.;
            %let _ph_fail = 0;
        %end;
        %else %do;

            data _ph_assess_;
                set _ph_assess_raw_;

                length _Effect_for_Check $80;
                PH_PValue = &_ph_pvar.;

                if not missing(Variable) then _Effect_for_Check = strip(Variable);
                else _Effect_for_Check = "";

                if not missing(PH_PValue) and PH_PValue < 0.05 then PH_Fail_Row = 1;
                else PH_Fail_Row = 0;

                keep _Effect_for_Check PH_PValue PH_Fail_Row;
            run;

            proc sql noprint;
                select max(PH_Fail_Row)
                    into :_ph_fail trimmed
                from _ph_assess_;
            quit;

            %if %superq(_ph_fail)= %then %let _ph_fail = 0;

        %end;

    %end;
    %else %do;
        %put WARNING: ASSESS PH output not created for &biomarker. - &milestone.. Defaulting to PHREG.;
        %let _ph_fail = 0;
    %end;

    %if %sysevalf(&_ph_fail. = 1) %then %let _method = LIFEREG;
    %else %let _method = PHREG;

    %if &_method. = LIFEREG %then %let _modelLabel = AFT;
    %else %let _modelLabel = CoxPH;

    %put NOTE: Selected &_method. for &biomarker. - &analysis_type. - &milestone.;


    /*-----------------------------------------------------
      Run final model and collect overall/pairwise rows
    -----------------------------------------------------*/

    %if &_method. = PHREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ModelANOVA        = _type3
                   ParameterEstimates = _pe
                   Estimates          = _est;

        proc phreg data=Temp;
            class &group.(ref='3') gen(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &group. / type3;

            estimate "T1 vs T2" &group.  1 -1 / exp;
            estimate "T1 vs T3" &group.  1  0 / exp;
            estimate "T2 vs T3" &group.  0  1 / exp;
        run;

        ods results;
        ods exclude none;

        %_get_nobs(_type3, _n_type3);
        %_get_nobs(_est,   _n_est);

        %if %sysevalf(&_n_type3. = 0 or &_n_est. = 0) %then %do;
            %put WARNING: Skipping &biomarker. - &analysis_type. - &milestone.: PHREG model output not created.;
            %goto cleanup;
        %end;

        data _overall;
            length Biomarker $100 AnalysisType $30 Milestone $120 Method $10 Parameter $80 Comparison $12;
            set _type3;

            if upcase(strip(Effect)) ne upcase("&group.") then delete;

            Biomarker      = "&biomarker.";
            BiomarkerOrder = &biomarker_order.;
            AnalysisType   = "&analysis_type.";
            AnalysisOrder  = &analysis_order.;
            Milestone      = "&milestone.";
            MilestoneOrder = &milestone_order.;
            Method         = "CoxPH";

            Parameter  = "Overall &group_label.";
            Comparison = "Overall";
            SortOrder  = 1;

            N_Used     = &_n_used.;
            N_Events   = &_n_events.;
            N_Censored = &_n_censored.;

            GroupDF  = DF;
            WaldStat = WaldChiSq;
            Beta     = .;
            SE       = .;
            Z_stat   = sqrt(WaldChiSq);
            P_value  = ProbChiSq;

            HR     = .;
            HR_LCL = .;
            HR_UCL = .;

            keep BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder
                 Biomarker AnalysisType Milestone Method Parameter Comparison
                 N_Used N_Events N_Censored GroupDF WaldStat Beta SE Z_stat P_value
                 HR HR_LCL HR_UCL;
        run;

        data _pairs;
            length Biomarker $100 AnalysisType $30 Milestone $120 Method $10 Parameter $80 Comparison $12;
            set _est;

            _label = upcase(strip(coalescec(Label, Contrast)));

            if      _label = "T1 VS T2" then do; Comparison = "T1 vs T2"; SortOrder = 2; end;
            else if _label = "T1 VS T3" then do; Comparison = "T1 vs T3"; SortOrder = 3; end;
            else if _label = "T2 VS T3" then do; Comparison = "T2 vs T3"; SortOrder = 4; end;
            else delete;

            Biomarker      = "&biomarker.";
            BiomarkerOrder = &biomarker_order.;
            AnalysisType   = "&analysis_type.";
            AnalysisOrder  = &analysis_order.;
            Milestone      = "&milestone.";
            MilestoneOrder = &milestone_order.;
            Method         = "CoxPH";

            Parameter = Comparison;

            N_Used     = &_n_used.;
            N_Events   = &_n_events.;
            N_Censored = &_n_censored.;

            GroupDF  = .;
            WaldStat = .;
            Beta     = Estimate;
            SE       = StdErr;

            if not missing(ChiSq) and ChiSq >= 0 then Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if not missing(ZValue) then Z_stat = ZValue;
            else if not missing(SE) and SE > 0 then Z_stat = Estimate / SE;
            else Z_stat = .;

            if not missing(ProbChiSq) then P_value = ProbChiSq;
            else if not missing(ProbZ) then P_value = ProbZ;
            else if not missing(Probt) then P_value = Probt;
            else P_value = .;

            HR     = exp(Estimate);
            HR_LCL = exp(Estimate - 1.96*StdErr);
            HR_UCL = exp(Estimate + 1.96*StdErr);

            keep BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder
                 Biomarker AnalysisType Milestone Method Parameter Comparison
                 N_Used N_Events N_Censored GroupDF WaldStat Beta SE Z_stat P_value
                 HR HR_LCL HR_UCL;
        run;

    %end;
    %else %if &_method. = LIFEREG %then %do;

        ods exclude all;
        ods noresults;

        ods output Type3Analysis      = _type3
                   Estimates          = _est
                   ParameterEstimates = _pe
                   CovB               = _covb;

        proc lifereg data=Temp outest=_aft_outest covout;
            class &group. gen;
            model &time.*&outcome.(0) = &covars. &group. / dist=weibull covb;

            estimate "T1 vs T2" &group.  1 -1  0;
            estimate "T1 vs T3" &group.  1  0 -1;
            estimate "T2 vs T3" &group.  0  1 -1;
        run;

        ods results;
        ods exclude none;

        %_get_nobs(_type3, _n_type3);
        %_get_nobs(_est,   _n_est);

        %if %sysevalf(&_n_type3. = 0 or &_n_est. = 0) %then %do;
            %put WARNING: Skipping &biomarker. - &analysis_type. - &milestone.: LIFEREG model output not created.;
            %goto cleanup;
        %end;

        data _null_;
            set _pe;
            if upcase(Parameter) = "SCALE" then call symputx("sigma", Estimate, "L");
        run;

        %if %superq(sigma)= %then %let sigma=.;

        data _type3b;
            set _type3;

            length Effect $80;

            if not missing(Variable) then Effect = strip(Variable);
        run;

        data _overall;
            length Biomarker $100 AnalysisType $30 Milestone $120 Method $10 Parameter $80 Comparison $12;
            set _type3b;

            if upcase(strip(Effect)) ne upcase("&group.") then delete;

            Biomarker      = "&biomarker.";
            BiomarkerOrder = &biomarker_order.;
            AnalysisType   = "&analysis_type.";
            AnalysisOrder  = &analysis_order.;
            Milestone      = "&milestone.";
            MilestoneOrder = &milestone_order.;
            Method         = "AFT";

            Parameter  = "Overall &group_label.";
            Comparison = "Overall";
            SortOrder  = 1;

            N_Used     = &_n_used.;
            N_Events   = &_n_events.;
            N_Censored = &_n_censored.;

            GroupDF  = DF;
            WaldStat = WaldChiSq;
            Beta     = .;
            SE       = .;
            Z_stat   = sqrt(WaldChiSq);
            P_value  = ProbChiSq;

            HR     = .;
            HR_LCL = .;
            HR_UCL = .;

            keep BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder
                 Biomarker AnalysisType Milestone Method Parameter Comparison
                 N_Used N_Events N_Censored GroupDF WaldStat Beta SE Z_stat P_value
                 HR HR_LCL HR_UCL;
        run;

        data _pairs;
            length Biomarker $100 AnalysisType $30 Milestone $120 Method $10 Parameter $80 Comparison $12;
            set _est;

            _label = upcase(strip(coalescec(Label, Contrast)));

            if      _label = "T1 VS T2" then do; Comparison = "T1 vs T2"; SortOrder = 2; end;
            else if _label = "T1 VS T3" then do; Comparison = "T1 vs T3"; SortOrder = 3; end;
            else if _label = "T2 VS T3" then do; Comparison = "T2 vs T3"; SortOrder = 4; end;
            else delete;

            Biomarker      = "&biomarker.";
            BiomarkerOrder = &biomarker_order.;
            AnalysisType   = "&analysis_type.";
            AnalysisOrder  = &analysis_order.;
            Milestone      = "&milestone.";
            MilestoneOrder = &milestone_order.;
            Method         = "AFT";

            Parameter = Comparison;

            N_Used     = &_n_used.;
            N_Events   = &_n_events.;
            N_Censored = &_n_censored.;

            GroupDF  = .;
            WaldStat = .;
            Beta     = Estimate;
            SE       = StdErr;

            if not missing(ZValue) then Z_stat = ZValue;
            else if not missing(ChiSq) and ChiSq >= 0 then Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if not missing(SE) and SE > 0 then Z_stat = Estimate / SE;
            else Z_stat = .;

            if not missing(ProbChiSq) then P_value = ProbChiSq;
            else if not missing(ProbZ) then P_value = ProbZ;
            else if not missing(Probt) then P_value = Probt;
            else P_value = .;

            HR     = .;
            HR_LCL = .;
            HR_UCL = .;

            if not missing(Estimate)
               and not missing(StdErr)
               and "&sigma." ne "."
               and "&sigma." ne "" then do;

                HR     = exp(-Estimate / &sigma.);
                HR_LCL = exp(-(Estimate + 1.96*StdErr) / &sigma.);
                HR_UCL = exp(-(Estimate - 1.96*StdErr) / &sigma.);
            end;

            keep BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder
                 Biomarker AnalysisType Milestone Method Parameter Comparison
                 N_Used N_Events N_Censored GroupDF WaldStat Beta SE Z_stat P_value
                 HR HR_LCL HR_UCL;
        run;

    %end;

    /*-----------------------------------------------------
      Covariate rows for the adjusted model statistics table.
      These are kept in the displayed report to match the original analyses,
      but the long combined table still uses the overall/pairwise rows.
    -----------------------------------------------------*/

    %if &_method. = PHREG %then %do;

        data _model_params;
            length Biomarker $100 AnalysisType $30 Milestone $120 Method $10 Parameter $80 Comparison $12 ClassVal $40;
            set _pe;

            Biomarker      = "&biomarker.";
            BiomarkerOrder = &biomarker_order.;
            AnalysisType   = "&analysis_type.";
            AnalysisOrder  = &analysis_order.;
            Milestone      = "&milestone.";
            MilestoneOrder = &milestone_order.;
            Method         = "CoxPH";

            Parameter = strip(Parameter);
            ClassVal  = strip(ClassVal0);
            Comparison = "";
            SortOrder  = 1.5;

            N_Used     = &_n_used.;
            N_Events   = &_n_events.;
            N_Censored = &_n_censored.;

            GroupDF  = .;
            WaldStat = .;
            Beta     = Estimate;
            SE       = StdErr;

            if not missing(ChiSq) and ChiSq >= 0 then Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if not missing(SE) and SE > 0 then Z_stat = Estimate / SE;
            else Z_stat = .;

            P_value = ProbChiSq;

            HR     = .;
            HR_LCL = .;
            HR_UCL = .;

            if DF = 0 then delete;
            if upcase(strip(Parameter)) = upcase("&group.") then delete;

            keep BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder
                 Biomarker AnalysisType Milestone Method Parameter Comparison ClassVal
                 N_Used N_Events N_Censored GroupDF WaldStat Beta SE Z_stat P_value
                 HR HR_LCL HR_UCL;
        run;

    %end;
    %else %if &_method. = LIFEREG %then %do;

        data _model_params;
            length Biomarker $100 AnalysisType $30 Milestone $120 Method $10 Parameter $80 Comparison $12 ClassVal $40;
            set _pe;

            Biomarker      = "&biomarker.";
            BiomarkerOrder = &biomarker_order.;
            AnalysisType   = "&analysis_type.";
            AnalysisOrder  = &analysis_order.;
            Milestone      = "&milestone.";
            MilestoneOrder = &milestone_order.;
            Method         = "AFT";

            Parameter = strip(Parameter);
            ClassVal  = strip(Level1);
            Comparison = "";
            SortOrder  = 1.5;

            N_Used     = &_n_used.;
            N_Events   = &_n_events.;
            N_Censored = &_n_censored.;

            GroupDF  = .;
            WaldStat = .;
            Beta     = Estimate;
            SE       = StdErr;

            if not missing(ZValue) then Z_stat = ZValue;
            else if not missing(ChiSq) and ChiSq >= 0 then Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if not missing(SE) and SE > 0 then Z_stat = Estimate / SE;
            else Z_stat = .;

            if not missing(ProbChiSq) then P_value = ProbChiSq;
            else if not missing(ProbZ) then P_value = ProbZ;
            else if not missing(Probt) then P_value = Probt;
            else P_value = .;

            HR     = .;
            HR_LCL = .;
            HR_UCL = .;

            if DF = 0 then delete;
            if upcase(strip(Parameter)) = upcase("&group.") then delete;
            if upcase(strip(Parameter)) = "SCALE" then delete;
            if upcase(strip(Parameter)) = "WEIBULL SHAPE" then delete;

            keep BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder
                 Biomarker AnalysisType Milestone Method Parameter Comparison ClassVal
                 N_Used N_Events N_Censored GroupDF WaldStat Beta SE Z_stat P_value
                 HR HR_LCL HR_UCL;
        run;

    %end;

    data _final_one;
        length Beta_CI_4 $70 SE_4 $20 Z_4 $20 P_4 $20 HR_CI_4 $70;
        set _overall _pairs;

        if not missing(Beta) and not missing(SE) then do;
            Beta_CI_4 = catx(" ",
                         put(Beta, 8.3),
                         cats("(",
                              put(Beta - 1.96*SE, 8.3), ",",
                              put(Beta + 1.96*SE, 8.3), ")"));
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
                            put(HR_LCL, 8.3), ",",
                            put(HR_UCL, 8.3), ")"));
        end;
        else HR_CI_4 = "";
    run;

    %_get_nobs(_final_one, _n_final);

    %if %sysevalf(&_n_final. > 0) %then %do;

        /*-------------------------------------------------
          Displayed report should match the original style:
            1) KM SGPLOT
            2) P-value / HR summary table
            3) Adjusted model statistics table
            4) BL vs Post-BL SGPLOT
        -------------------------------------------------*/

        ods graphics on / reset width=6.8in height=4.6in noborder;
        ods escapechar='^';

        /* Pairwise log-rank tests for the P-value / HR summary table */
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

        /* Overall KM + log-rank, then draw with PROC SGPLOT */
        ods select none;
        proc lifetest data=Temp
            plots=survival(cb=hw atrisk=0 to 12 by 1)
            confband=hw bandmin=0 bandmax=12;
            time &time.*&outcome.(0);
            strata &group. / test=logrank;
            id patno;
            label &time. = "Year";
            ods output HomTests=HomTests SurvivalPlot=survsaa3;
        run;
        ods select all;

        data _null_;
            set HomTests;
            where upcase(Test) = "LOG-RANK";
            call symputx("pval", put(ProbChiSq, pvalue6.4), "L");
        run;

        data _null_;
            set _final_one;
            where Comparison = "Overall";
            call symputx("pval2",  P_4, "L");
            call symputx("pval2n", P_value, "L");
        run;

        %if %superq(pval)=  %then %let pval=.;
        %if %superq(pval2)= %then %let pval2=.;
        %if %superq(pval2n)= %then %let pval2n=. ;

        data temp_lr;
            length Comparison $12;
            set LR12(in=a keep=ProbChiSq rename=(ProbChiSq=LogRank_p))
                LR13(in=b keep=ProbChiSq rename=(ProbChiSq=LogRank_p))
                LR23(in=c keep=ProbChiSq rename=(ProbChiSq=LogRank_p));

            if      a then Comparison = "1 vs 2";
            else if b then Comparison = "1 vs 3";
            else if c then Comparison = "2 vs 3";

            keep Comparison LogRank_p;
        run;

        data _pair_model;
            length Comparison $12 CoxPH $8 HR_CI $70;
            set _final_one;

            if      Comparison = "T1 vs T2" then Comparison = "1 vs 2";
            else if Comparison = "T1 vs T3" then Comparison = "1 vs 3";
            else if Comparison = "T2 vs T3" then Comparison = "2 vs 3";
            else delete;

            CoxPH = P_4;
            HR_CI = HR_CI_4;

            keep Comparison CoxPH HR_CI P_value;
        run;

        proc sort data=temp_lr;     by Comparison; run;
        proc sort data=_pair_model; by Comparison; run;

        data PVALS_pairs;
            merge temp_lr(in=l) _pair_model(in=m);
            by Comparison;

            length LogRank $8;
            if l or m;

            if not missing(LogRank_p) then LogRank = put(LogRank_p, pvalue6.4);
            else LogRank = ".";

            keep Comparison LogRank CoxPH HR_CI;
        run;

        data PVALS_overall;
            length Comparison $12 LogRank $8 CoxPH $8 HR_CI $70;
            Comparison = "Overall";
            LogRank = "&pval.";
            CoxPH = "&pval2.";
            HR_CI = "NA";
            output;
        run;

        data PVALS;
            set PVALS_overall PVALS_pairs;
        run;

        proc sql noprint;
            select CoxPH into :p12 trimmed from PVALS where Comparison = "1 vs 2";
            select CoxPH into :p13 trimmed from PVALS where Comparison = "1 vs 3";
            select CoxPH into :p23 trimmed from PVALS where Comparison = "2 vs 3";
        quit;

        %if %superq(p12)= %then %let p12=.;
        %if %superq(p13)= %then %let p13=.;
        %if %superq(p23)= %then %let p23=. ;

        %let dsid = %sysfunc(open(survsaa3));
        %if &dsid. %then %do;
            %let hasLCLM = %sysfunc(varnum(&dsid., LCLM));
            %let rc = %sysfunc(close(&dsid.));
        %end;
        %else %let hasLCLM = 0;

        %if %sysevalf(%superq(pval2n)=, boolean) or %sysevalf(&pval2n. >= 0.05) %then %do;
            footnote1 height=8pt
                "&_modelLabel. p-value | Overall Interaction P-value: &pval2.";
        %end;
        %else %do;
            footnote1 height=8pt
                "&_modelLabel. p-values | Overall: &pval2. | T1 vs T2: &p12. | T1 vs T3: &p13. | T2 vs T3: &p23.";
        %end;

        title "&milestone. for Healthy Control";

        proc sgplot data=survsaa3;
            styleattrs datacolors=(blue red green)
                       datacontrastcolors=(blue red green);
            step    x=Time y=Survival / group=stratumnum name="km";
            scatter x=Time y=Censored / group=stratumnum markerattrs=(symbol=plus);
            band x=Time upper=HW_UCL lower=HW_LCL /
                 group=stratumnum transparency=0.80;
            %if &hasLCLM. %then %do;
                band x=Time lower=LCLM upper=UCLM /
                     group=stratumnum transparency=0.80;
            %end;
            xaxistable Atrisk / location=inside x=tatrisk
                                class=stratumnum colorgroup=stratumnum;
            xaxis values=(0 to 12 by 1);
            format stratumnum &groupfmt.;
            label  stratumnum="&group_label.";
            keylegend "km" / title="&group_label."
                             location=outside position=bottom
                             across=3 opaque;
        run;

        footnote;
        title;

        /* P-value / HR summary table */
        %if &_method. = PHREG %then %do;
            footnote1 j=left h=8pt
                "HR (95% CI) from Cox proportional hazards model. Reference: T3.";
        %end;
        %else %if &_method. = LIFEREG %then %do;
            footnote1 j=left h=8pt
                "Derived HR (95% CI) from Weibull AFT model using HR=exp(-beta/sigma); delta-method CI on log(HR). Reference: T3.";
        %end;

        title2 "P-Value and HR Summary - &milestone. for Healthy Control";

        proc report data=PVALS nowd headline;
            columns Comparison LogRank CoxPH HR_CI;
            define Comparison / 'Comparison'           width=12;
            define LogRank    / 'Log-Rank|p-value'     width=16;
            define CoxPH      / "&_modelLabel.|p-value"    width=14;
            %if &_method. = LIFEREG %then %do;
                define HR_CI / 'Derived HR (95% CI)' width=24;
            %end;
            %else %do;
                define HR_CI / 'HR (95% CI)'         width=24;
            %end;
        run;

        footnote;
        title;
        title2;

        /* Adjusted model statistics table */
        data _pdf_modelstats;
            length ClassVal $40;
            set _overall _model_params _pairs;
        run;

        proc sort data=_pdf_modelstats;
            by SortOrder Parameter;
        run;

        data _pdf_modelstats_fmt;
            length Beta_CI_4 $70 SE_4 $20 Z_4 $20 P_4 $20 HR_CI_4 $70;
            set _pdf_modelstats;

            if not missing(Beta) and not missing(SE) then do;
                Beta_CI_4 = catx(" ",
                             put(Beta, 8.3),
                             cats("(",
                                  put(Beta - 1.96*SE, 8.3), ",",
                                  put(Beta + 1.96*SE, 8.3), ")"));
            end;
            else Beta_CI_4 = "";

            if not missing(SE) then SE_4 = put(SE, 8.3);
            else SE_4 = "-";

            if not missing(Z_stat) then Z_4 = put(Z_stat, 8.3);
            else Z_4 = "";

            if not missing(P_value) then do;
                if P_value < 0.001 then P_4 = "<0.001";
                else P_4 = put(P_value, 8.3);
            end;
            else P_4 = "";

            if not missing(HR) and not missing(HR_LCL) and not missing(HR_UCL) then do;
                HR_CI_4 = cats(put(HR, 8.3), "(",
                               put(HR_LCL, 8.3), ",",
                               put(HR_UCL, 8.3), ")");
            end;
            else HR_CI_4 = "";
        run;

        title2 "Adjusted Model Statistics - &milestone. for Healthy Control";

        proc report data=_pdf_modelstats_fmt nowd headline;
            columns Milestone Method Parameter ClassVal Beta_CI_4 SE_4 Z_4 P_4 HR_CI_4;
            define Milestone  / 'Milestone';
            define Method     / 'Model';
            define Parameter  / 'Parameter|Comparison';
            define ClassVal   / 'Level';
            define Beta_CI_4  / 'Beta (95% CI)' width=28;
            define SE_4       / 'Standard|Error';
            define Z_4        / 'Z-Statistic';
            define P_4        / 'p-value';
            %if &_method. = LIFEREG %then %do;
                define HR_CI_4 / 'Derived HR (95% CI)' width=24;
            %end;
            %else %do;
                define HR_CI_4 / 'HR (95% CI)' width=24;
            %end;
        run;

        title2;

        /* BL vs Post-BL stacked bar plot */
        proc format;
            value bl
                1 = "BL"
                2 = "Post";
        run;

        data _event_timing;
            set &data.;

            if &outcome. = 1 then do;
                if &time. = 0 then BL = 1;
                else if &time. > 0 then BL = 2;
            end;
        run;

        data _seed_bar;
            do BL = 1 to 2;
                do &group. = 1 to 3;
                    _wt = 0;
                    output;
                end;
            end;
        run;

        data _event_timing_plot;
            set _seed_bar _event_timing(in=real);
            if real then _wt = 1;
        run;

        proc sort data=_event_timing_plot;
            by BL &group.;
        run;

        %let _bl_list=;

        proc sql noprint;
            select distinct BL
                into :_bl_list separated by ' '
            from _event_timing
            where BL in (1, 2) and not missing(BL);
        quit;

        %if %length(&_bl_list.) %then %do;

            title "Event Timing (BL vs Post-BL) - &milestone. for Healthy Control";

            proc sgplot data=_event_timing_plot(where=(BL in (&_bl_list.)));
                styleattrs datacolors=(cx556B8A cxA05252 cx518A85);
                vbar BL / response=_wt stat=sum group=&group.
                          grouporder=ascending groupdisplay=stack;
                format &group. &groupfmt. BL bl.;
                label  &group.="&group_label." BL="BL-Post";
                yaxis label="Frequency";
                xaxis display=(nolabel);
            run;

            title;

        %end;
        %else %do;
            %put NOTE: No BL or post-baseline events for &biomarker. - &analysis_type. - &milestone. - skipping stacked-bar plot.;
        %end;

        footnote;
        title;
        title2;

        proc append base=Final_HCNSD_Tertile_Summary
                    data=_final_one
                    force;
        run;

        %put NOTE: Appended &biomarker. - &analysis_type. - &milestone. to final table.;
    %end;
    %else %do;
        %put WARNING: Skipping &biomarker. - &analysis_type. - &milestone.: no final rows created.;
    %end;

    %cleanup:

    proc datasets library=work nolist;
        delete Temp
               _ph_assess_raw_ _ph_assess_
               _type3 _type3b _est _est2 _overall _pairs _model_params _final_one
               _pdf_modelstats _pdf_modelstats_fmt _pe _covb _aft_outest
               HomTests LR12 LR13 LR23 temp_lr _pair_model PVALS_pairs PVALS_overall PVALS
               survsaa3 _event_timing _seed_bar _event_timing_plot;
    quit;

%mend RunOneTertileModel;


/*===============================================================
  7. RUN PROGRESSION + CLINICAL MILESTONES FOR ONE BIOMARKER
===============================================================*/

%macro RunAllForOneBiomarker(
    biomarker_short=,
    biomarker_label=,
    biomarker_order=,
    datevar=,
    biomarker_var=,
    cutoff_type=,
    groupfmt=,
    group_label=
);

    options orientation=landscape nodate nonumber;
    ods graphics on / reset=all;

    title "HC NSD- Tertile Survival Analysis - &biomarker_label.";

    %MakeEarlyHCNSDMinusBiomarker(
        biomarker_short=&biomarker_short.,
        biomarker_label=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        datevar=&datevar.,
        biomarker_var=&biomarker_var.,
        cutoff_type=&cutoff_type.
    );

    %MakeProgressionSurvData(
        biomarker_short=&biomarker_short.
    );

    /*-------------------------
      Progression milestones
    -------------------------*/

    %RunOneTertileModel(
        data=ProgressionSurv_&biomarker_short.,
        outcome=STATUS_any,
        time=TIME_ENRL_any,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Progression,
        analysis_order=1,
        milestone=Any progression milestone,
        milestone_order=1,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ProgressionSurv_&biomarker_short.,
        outcome=STATUS_wb_any,
        time=TIME_ENRL_wb_any,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Progression,
        analysis_order=1,
        milestone=Walking and balance domain,
        milestone_order=2,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ProgressionSurv_&biomarker_short.,
        outcome=STATUS_mc_any,
        time=TIME_ENRL_mc_any,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Progression,
        analysis_order=1,
        milestone=Motor complications domain,
        milestone_order=3,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ProgressionSurv_&biomarker_short.,
        outcome=STATUS_cog_any,
        time=TIME_ENRL_cog_any,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Progression,
        analysis_order=1,
        milestone=Cognitive domain,
        milestone_order=4,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ProgressionSurv_&biomarker_short.,
        outcome=STATUS_auto_any,
        time=TIME_ENRL_auto_any,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Progression,
        analysis_order=1,
        milestone=Autonomic dysfunction domain,
        milestone_order=5,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ProgressionSurv_&biomarker_short.,
        outcome=STATUS_fd_any,
        time=TIME_ENRL_fd_any,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Progression,
        analysis_order=1,
        milestone=Functional dependence domain,
        milestone_order=6,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    /* Activities of daily living domain intentionally excluded from the
       requested progression milestone output. */


    /*-------------------------
      Clinical milestones
    -------------------------*/

    %MakeClinicalSurvData(
        biomarker_short=&biomarker_short.
    );

    %RunOneTertileModel(
        data=ClinicalSurv_&biomarker_short.,
        outcome=STATUS_HY_GE3,
        time=TIME_ENRL_HY_GE3,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Clinical,
        analysis_order=2,
        milestone=Time to HY >= 3,
        milestone_order=1,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ClinicalSurv_&biomarker_short.,
        outcome=STATUS_COGSTATE_CHANGE,
        time=TIME_ENRL_COGSTATE_CHANGE,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Clinical,
        analysis_order=2,
        milestone=Time to Cognitive State Change,
        milestone_order=2,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ClinicalSurv_&biomarker_short.,
        outcome=STATUS_MOCA_LT21,
        time=TIME_ENRL_MOCA_LT_21,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Clinical,
        analysis_order=2,
        milestone=Time to MoCA < 21,
        milestone_order=3,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ClinicalSurv_&biomarker_short.,
        outcome=STATUS_MOCA_LT26,
        time=TIME_ENRL_MOCA_LT_26,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Clinical,
        analysis_order=2,
        milestone=Time to MoCA < 26,
        milestone_order=4,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    %RunOneTertileModel(
        data=ClinicalSurv_&biomarker_short.,
        outcome=STATUS_UPDRS_Change_GT5,
        time=TIME_ENRL_UPDRS_Change_GT5,
        biomarker=&biomarker_label.,
        biomarker_order=&biomarker_order.,
        analysis_type=Clinical,
        analysis_order=2,
        milestone=Time to MDS-UPDRS III Change 5+,
        milestone_order=5,
        group_label=&group_label.,
        groupfmt=&groupfmt.
    );

    ods graphics / reset;

    proc datasets library=work nolist;
        delete Early_&biomarker_short.
               ProgressionSurv_&biomarker_short.
               ClinicalSurv_&biomarker_short.;
    quit;

%mend RunAllForOneBiomarker;


/*===============================================================
  8. RUN ALL BIOMARKERS IN REQUESTED ORDER
===============================================================*/

options orientation=landscape nodate nonumber;
ods graphics on / reset=all;


%RunAllForOneBiomarker(
    biomarker_short=PTAU217,
    biomarker_label=pTau217 plasma,
    biomarker_order=1,
    datevar=Lilly_ptau217p_date,
    biomarker_var=Lilly_ptau217p,
    cutoff_type=PTAU,
    groupfmt=ptau.,
    group_label=pTau217 plasma
);

%RunAllForOneBiomarker(
    biomarker_short=NFL,
    biomarker_label=NFL,
    biomarker_order=2,
    datevar=IU_Plasma_date,
    biomarker_var=IU_NFL_plasma_mean,
    cutoff_type=NFL,
    groupfmt=nfl.,
    group_label=NFL
);

%RunAllForOneBiomarker(
    biomarker_short=GFAP,
    biomarker_label=GFAP,
    biomarker_order=3,
    datevar=IU_Plasma_date,
    biomarker_var=IU_GFAP_plasma_mean,
    cutoff_type=GFAP,
    groupfmt=gfap.,
    group_label=GFAP
);

ods graphics / reset;



/*===============================================================
  9. FINAL COMBINED SUMMARY TABLE
===============================================================*/

proc sort data=Final_HCNSD_Tertile_Summary;
    by BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder;
run;

/* Always retain the overall row. Retain T1 vs T2, T1 vs T3, and
   T2 vs T3 only when the corresponding overall p-value is < 0.05. */
proc sql;
    create table _Overall_Significance as
    select BiomarkerOrder,
           AnalysisOrder,
           MilestoneOrder,
           max(case
                   when SortOrder = 1
                    and not missing(P_value)
                    and P_value < 0.05
                   then 1
                   else 0
               end) as Overall_Sig
    from Final_HCNSD_Tertile_Summary
    group by BiomarkerOrder, AnalysisOrder, MilestoneOrder;

    create table Final_HCNSD_Tertile_Summary_RTF as
    select a.*
    from Final_HCNSD_Tertile_Summary as a
    left join _Overall_Significance as b
      on a.BiomarkerOrder = b.BiomarkerOrder
     and a.AnalysisOrder  = b.AnalysisOrder
     and a.MilestoneOrder = b.MilestoneOrder
    where a.SortOrder = 1
       or coalesce(b.Overall_Sig, 0) = 1
    order by a.BiomarkerOrder,
             a.AnalysisOrder,
             a.MilestoneOrder,
             a.SortOrder;
quit;

options
    orientation = landscape
    leftmargin  = 0.25in
    rightmargin = 0.25in
    topmargin   = 0.40in
    bottommargin= 0.40in
    nodate
    nonumber
;


title "Supplement Table. Association between baseline plasma biomarker tertiles and time-to-event progression and clinical milestones in HC (NSD-)";

footnote1 j=left h=7pt
    "Rows are ordered by biomarker, then progression milestones, then clinical milestones.";
footnote2 j=left h=7pt
    "Overall row reports the Type III omnibus test for the tertile biomarker. Pairwise tertile rows are shown only when the overall p-value is < 0.05.";
footnote3 j=left h=7pt
    "Cox PH models report HR. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Models adjusted for age and sex. Values rounded to 3 decimals.";

proc report data=Final_HCNSD_Tertile_Summary_RTF nowd split='|'
    style(report)=[outputwidth=10.75in cellspacing=1 cellpadding=2]
    style(header)=[font_size=7pt font_weight=bold]
    style(column)=[font_size=7pt];

    columns BiomarkerOrder AnalysisOrder MilestoneOrder SortOrder
            Biomarker AnalysisType Milestone Method Parameter
            N_Used N_Events N_Censored GroupDF WaldStat Beta_CI_4 SE_4 Z_4 P_value P_4 HR_CI_4;

    define BiomarkerOrder / order noprint;
    define AnalysisOrder  / order noprint;
    define MilestoneOrder / order noprint;
    define SortOrder      / order noprint;

    define Biomarker      / group 'Biomarker';
    define AnalysisType   / group 'Analysis|Type';
    define Milestone      / group 'Milestone';
    define Method         / display 'Model';
    define Parameter      / display 'Parameter /|Comparison';

    define N_Used         / display 'N' format=8.;
    define N_Events       / display 'Events' format=8.;
    define N_Censored     / display 'Censored' format=8.;
    define GroupDF        / display 'DF' format=4.;
    define WaldStat       / display 'Wald Chi-Sq' format=8.3;
    define Beta_CI_4      / display 'Coefficient (95% CI)';
    define SE_4           / display 'SE';
    define Z_4            / display 'Z';
    define P_value        / noprint;
    define P_4            / display 'p';
    define HR_CI_4        / display 'HR (95% CI)';

    compute P_4;
        if not missing(P_value) and P_value < 0.05 then
            call define(_col_, 'style', 'style=[foreground=red font_weight=bold]');
    endcomp;

    break after Biomarker / skip;
    break after AnalysisType / skip;
    break after Milestone / skip;

run;

footnote;
title;



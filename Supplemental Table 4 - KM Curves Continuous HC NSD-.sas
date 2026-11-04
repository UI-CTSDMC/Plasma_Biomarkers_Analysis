/*=========================================================
  HC NSD-: Continuous Biomarkers and Survival Milestones

  Uses the prepared progression_analysis_v3 and
  clinical_analysis_v3 datasets to examine associations
  between continuous plasma biomarker levels and time to
  progression and clinical milestones in HC NSD- participants.

  Fits age and sex-adjusted Cox proportional hazards or
  Weibull accelerated failure time models, selected according
  to the proportional hazards assumption test, and generates
  separate progression and clinical survival summary tables.
  The original rules for skipping unsuccessful models are retained.

=========================================================*/


%let filepath = /*path where analysis data are located*/;
libname analysis "&filepath.";

/*==========================================================================
  COMBINED HC NSD- PROGRESSION AND CLINICAL MILESTONE SURVIVAL PROGRAMS

  Execution order is intentional:
    1. Progression milestone program
    2. Clinical milestone program
==========================================================================*/

%let report_date = 08DEC2025;

/* Make local copies of the prepared analysis datasets. */
data progression_analysis_v3;
    set analysis.progression_analysis_v3;
run;

data clinical_analysis_v3;
    set analysis.clinical_analysis_v3;
run;
%let write_individual_outputs = 0;

ods escapechar="^";

proc datasets library=work nolist;
    delete Combined_HC_Progression_Summary
           Combined_HC_Clinical_Summary;
quit;


/*=========================================================
  Continuous survival macro

  Important update:
      If overall biomarker p-value is missing or model output fails,
      the milestone is skipped from displayed reports.
=========================================================*/

%macro SurvCont(
    data=,
    outcome=,
    time=,
    biomarker=,
    label=,
    method=AUTO,
    covars=age gen,
    classvars=gen,
    biomarker_label=,
    title_suffix=for Healthy Control NSD-,
    baseline_vars=
);

    %local _method _ph_fail _modelLabel _dslabel;
    %local sigma var_sigma cov_bio_sigma var_bio_beta;
    %local _nobs _nevents _ncensors _overall_ok _model_ok;

    %let sigma=.;
    %let var_sigma=.;
    %let cov_bio_sigma=.;
    %let var_bio_beta=.;
    %let _overall_ok=0;
    %let _model_ok=0;

    %if %superq(baseline_vars)= %then 
        %let baseline_vars=Patno &biomarker. &outcome. &time.;

    %let _dslabel = %sysfunc(compress(&label.,
        abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789, k));

    %let _dslabel = %substr(&_dslabel., 1,
        %sysfunc(min(29, %length(&_dslabel.))));


    /*-------------------------------------------------
      1. Remove missing values and baseline events
    -------------------------------------------------*/

    data Temp;
        set &data.;

        if missing(&biomarker.) then delete;
        if missing(&time.) then delete;
        if missing(&outcome.) then delete;

        if &outcome. = 1 and &time. = 0 then delete;
    run;

    proc sql noprint;
        select count(*),
               sum(&outcome. = 1),
               sum(&outcome. = 0)
            into :_nobs trimmed,
                 :_nevents trimmed,
                 :_ncensors trimmed
        from Temp;
    quit;

    %if %sysevalf(&_nobs. = 0) or %sysevalf(&_nevents. = 0) %then %do;
        %put WARNING: Skipping &label. for &biomarker_label. - no usable observations or no events after exclusions.;
        proc datasets library=work nolist;
            delete Temp;
        quit;
        %return;
    %end;


    /*-------------------------------------------------
      2. Automatic PH testing
    -------------------------------------------------*/

    %let _ph_fail = 0;
    %let _method  = %upcase(&method.);

    %if %upcase(&method.) = AUTO %then %do;

        %put NOTE: method=AUTO - running PH assessment for &label. and &biomarker_label.;

        ods exclude all;
        ods noresults;

        ods output ProportionalHazardsSupTest = _ph_assess_;

        proc phreg data=Temp;
            class &classvars.(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &biomarker. / maxiter=100;
            assess ph / resample seed=1472160691;;
        run;

        ods results;
        ods exclude none;

        %if %sysfunc(exist(_ph_assess_)) %then %do;

            data _null_;
                set _ph_assess_ end=last;
                retain ph_fail 0;

                if not missing(pvalue) and pvalue < 0.05 then ph_fail = 1;

                if last then call symputx('_ph_fail', ph_fail, 'L');
            run;

            %if &_ph_fail. = 0 %then %do;
                %let _method = PHREG;
                %put NOTE: PH not violated for &label. and &biomarker_label. - using Cox PH model.;
            %end;
            %else %do;
                %let _method = LIFEREG;
                %put NOTE: PH violated for &label. and &biomarker_label. - using Weibull AFT model.;
            %end;

        %end;
        %else %do;
            %put WARNING: PH assessment output not created for &label. and &biomarker_label.. Trying Cox PH model.;
            %let _method = PHREG;
        %end;

        proc datasets lib=work nolist;
            delete _ph_assess_;
        quit;

    %end;
    %else %do;
        %put NOTE: method=&method. specified manually - skipping PH test for &label. and &biomarker_label.;
        %let _method = %upcase(&method.);
    %end;


    /*-------------------------------------------------
      3. Fit selected model
    -------------------------------------------------*/

    %if &_method. = PHREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ParameterEstimates = ModelPE
                   ModelANOVA         = Type3PE
                   HazardRatios       = HRPE;

        proc phreg data=Temp;
            class &classvars.(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &biomarker. / type3;
            hazardratio &biomarker. / units=1 cl=wald;
        run;

        ods results;
        ods exclude none;

    %end;
    %else %if &_method. = LIFEREG %then %do;

        ods exclude all;
        ods noresults;

        ods output ParameterEstimates = ModelPE
                   Type3Analysis      = Type3PE
                   CovB               = CovBPE;

        proc lifereg data=Temp covout;
            class &classvars.;
            model &time.*&outcome.(0) = &covars. &biomarker. / dist=weibull covb;
        run;

        ods results;
        ods exclude none;

        %if %sysfunc(exist(Type3PE)) %then %do;
            proc datasets lib=work nolist;
                modify Type3PE;
                rename Variable = Effect;
            quit;
        %end;

        %if %sysfunc(exist(ModelPE)) %then %do;
            data _null_;
                set ModelPE;

                if upcase(Parameter) = "SCALE" then 
                    call symputx('sigma', Estimate);
            run;
        %end;

        %if %sysfunc(exist(CovBPE)) %then %do;
            proc sql noprint;
                select Scale into :var_sigma trimmed
                from CovBPE
                where upcase(Variable) = "SCALE";

                select Scale into :cov_bio_sigma trimmed
                from CovBPE
                where upcase(Variable) = upcase("&biomarker.");

                select &biomarker. into :var_bio_beta trimmed
                from CovBPE
                where upcase(Variable) = upcase("&biomarker.");
            quit;
        %end;

        %if %superq(var_sigma)= %then %let var_sigma=.;
        %if %superq(cov_bio_sigma)= %then %let cov_bio_sigma=.;
        %if %superq(var_bio_beta)= %then %let var_bio_beta=.;

    %end;
    %else %do;

        %put WARNING: method=&_method. not recognized. Trying Cox PH model.;
        %let _method = PHREG;

        ods exclude all;
        ods noresults;

        ods output ParameterEstimates = ModelPE
                   ModelANOVA         = Type3PE
                   HazardRatios       = HRPE;

        proc phreg data=Temp;
            class &classvars.(ref='1') / param=ref;
            model &time.*&outcome.(0) = &covars. &biomarker. / type3;
            hazardratio &biomarker. / units=1 cl=wald;
        run;

        ods results;
        ods exclude none;

    %end;


    /*-------------------------------------------------
      4. Skip milestone if model output or overall p-value is absent
    -------------------------------------------------*/

    %if not %sysfunc(exist(Type3PE)) or not %sysfunc(exist(ModelPE)) %then %do;
        %put WARNING: Skipping &label. for &biomarker_label. - required model output dataset missing.;
        proc datasets library=work nolist;
            delete Temp ModelPE Type3PE HRPE CovBPE;
        quit;
        %return;
    %end;

    data OverallBiomarker;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set Type3PE;

        where upcase(strip(Effect)) = upcase("&biomarker.");

        Biomarker = "&biomarker_label.";
        Milestone = "&label.";
        Method    = "%sysfunc(ifc(&_method.=LIFEREG,AFT,CoxPH))";
        Parameter = "Overall &biomarker_label.";
        ClassVal  = "";

        GroupDF   = DF;
        WaldStat  = WaldChiSq;

        Beta      = .;
        SE        = .;
        Beta_CI   = "";

        if not missing(WaldChiSq) and WaldChiSq >= 0 then 
            Z_stat = sqrt(WaldChiSq);
        else 
            Z_stat = .;

        P_value = ProbChiSq;

        HR      = .;
        HR_LCL  = .;
        HR_UCL  = .;
        HR_CI   = "";

        keep Biomarker Milestone Method Parameter ClassVal GroupDF WaldStat
             Beta SE Beta_CI Z_stat P_value HR HR_LCL HR_UCL HR_CI;
    run;

    proc sql noprint;
        select count(*)
            into :_overall_ok trimmed
        from OverallBiomarker
        where not missing(P_value);
    quit;

    %if %sysevalf(&_overall_ok. = 0) %then %do;
        %put WARNING: Skipping &label. for &biomarker_label. - overall biomarker p-value missing.;
        proc datasets library=work nolist;
            delete Temp OverallBiomarker ModelPE Type3PE HRPE CovBPE;
        quit;
        %return;
    %end;


    /*-------------------------------------------------
      5. Model label
    -------------------------------------------------*/

    %let _modelLabel = %sysfunc(ifc(&_method.=LIFEREG,AFT,CoxPH));


    /*-------------------------------------------------
      6. Parameter estimates table
    -------------------------------------------------*/

    data ModelPE_Tidy;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set ModelPE;

        Biomarker = "&biomarker_label.";
        Milestone = "&label.";
        Method    = "&_modelLabel.";

        Beta = Estimate;
        SE   = StdErr;

        GroupDF  = .;
        WaldStat = .;

        %if &_method. = PHREG %then %do;

            ClassVal = strip(ClassVal0);

            if upcase(strip(Parameter)) = upcase("&biomarker.") then 
                Parameter = "&biomarker_label. continuous";

            if not missing(ChiSq) and ChiSq >= 0 then 
                Z_stat = sign(Estimate) * sqrt(ChiSq);
            else if SE > 0 then 
                Z_stat = Beta / SE;
            else 
                Z_stat = .;

            P_value = ProbChiSq;

            if not missing(Beta) and not missing(SE) and SE > 0 then
                Beta_CI = catx(" ",
                               put(Beta,8.4),
                               cats("(",
                                    put(Beta - 1.96*SE,8.4), ", ",
                                    put(Beta + 1.96*SE,8.4), ")"));
            else 
                Beta_CI = "";

            if upcase(strip(Parameter)) = upcase("&biomarker_label. CONTINUOUS") then do;
                HR     = exp(Beta);
                HR_LCL = exp(Beta - 1.96*SE);
                HR_UCL = exp(Beta + 1.96*SE);

                HR_CI = catx(" ",
                             put(HR,8.4),
                             cats("(",
                                  put(HR_LCL,8.4), ", ",
                                  put(HR_UCL,8.4), ")"));
            end;
            else do;
                HR     = .;
                HR_LCL = .;
                HR_UCL = .;
                HR_CI  = "";
            end;

        %end;
        %else %do;

            ClassVal = strip(Level1);

            if upcase(strip(Parameter)) = upcase("&biomarker.") then 
                Parameter = "&biomarker_label. continuous";

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

            if not missing(Beta) and not missing(SE) and SE > 0 then
                Beta_CI = catx(" ",
                               put(Beta,8.4),
                               cats("(",
                                    put(Beta - 1.96*SE,8.4), ", ",
                                    put(Beta + 1.96*SE,8.4), ")"));
            else 
                Beta_CI = "";

            HR     = .;
            HR_LCL = .;
            HR_UCL = .;
            HR_CI  = "";

            if upcase(strip(Parameter)) = upcase("&biomarker_label. CONTINUOUS") then do;

                if not missing(Beta) and 
                   "&sigma." ne "." and "&sigma." ne "" and 
                   "&var_sigma." ne "." and "&var_sigma." ne "" and 
                   "&var_bio_beta." ne "." and "&var_bio_beta." ne "" then do;

                    logHR = -Beta / &sigma.;
                    HR    = exp(logHR);

                    var_logHR =
                        (1 / (&sigma.**2)) * &var_bio_beta.
                      + ((Beta**2) / (&sigma.**4)) * &var_sigma.
                      - (2 * Beta / (&sigma.**3)) * &cov_bio_sigma.;

                    if not missing(var_logHR) and var_logHR >= 0 then do;
                        se_logHR = sqrt(var_logHR);
                        HR_LCL   = exp(logHR - 1.96*se_logHR);
                        HR_UCL   = exp(logHR + 1.96*se_logHR);

                        HR_CI = catx(" ",
                                     put(HR,8.4),
                                     cats("(",
                                          put(HR_LCL,8.4), ", ",
                                          put(HR_UCL,8.4), ")"));
                    end;
                end;
            end;

        %end;

        if DF = 0 then delete;
        if upcase(strip(Parameter)) = "SCALE" then delete;
        if upcase(strip(Parameter)) = "WEIBULL SHAPE" then delete;

        keep Biomarker Milestone Method Parameter ClassVal GroupDF WaldStat
             Beta SE Beta_CI Z_stat P_value HR HR_LCL HR_UCL HR_CI;
    run;


    /*-------------------------------------------------
      7. Combine output
    -------------------------------------------------*/

    data ModelStats;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set OverallBiomarker ModelPE_Tidy;
    run;

    proc append base=All_ModelStats data=ModelStats;
    run;

    data MS_&_dslabel.;
        set ModelStats;
    run;


    /*-------------------------------------------------
      8. Print compact model statistics
    -------------------------------------------------*/

    data ModelStats_PDF;
        length Parameter_Display $100;
        set ModelStats;

        Parameter_Display = Parameter;

        if index(upcase(Parameter), "CONTINUOUS") then 
            Parameter_Display = "Continuous biomarker";
    run;

    title "Adjusted Survival Model Statistics - &label. &title_suffix";
    title2 "Continuous Biomarker: &biomarker_label.";

    %if &_method. = PHREG %then %do;
        footnote1 j=left h=8pt
            "Cox proportional hazards model adjusted for age and sex. HR is per 1-unit increase in &biomarker_label..";
    %end;
    %else %if &_method. = LIFEREG %then %do;
        footnote1 j=left h=8pt
            "Weibull AFT model adjusted for age and sex. Derived HR uses HR=exp(-beta/sigma), where sigma is the Weibull AFT scale parameter.";
    %end;

    proc report data=ModelStats_PDF nowd split='|'
        style(report)=[frame=box rules=groups cellspacing=0 cellpadding=3]
        style(header)=[font_size=8pt font_weight=bold background=white just=center]
        style(column)=[font_size=8pt];

        columns Method Parameter_Display ClassVal GroupDF WaldStat Beta_CI SE Z_stat P_value HR_CI;

        define Method / display 'Model'
            style(column)=[just=center];

        define Parameter_Display / display 'Parameter'
            style(column)=[just=left];

        define ClassVal / display 'Level'
            style(column)=[just=center];

        define GroupDF / display 'DF' format=4.
            style(column)=[just=center];

        define WaldStat / display 'Wald|Chi-Sq' format=8.4
            style(column)=[just=center];

        define Beta_CI / display 'Beta|(95% CI)'
            style(column)=[just=center];

        define SE / display 'SE' format=8.4
            style(column)=[just=center];

        define Z_stat / display 'Z' format=8.4
            style(column)=[just=center];

        define P_value / display 'p-value' format=pvalue6.4
            style(column)=[just=center];

        %if &_method. = PHREG %then %do;
            define HR_CI / display 'HR|(95% CI)'
                style(column)=[just=center];
        %end;
        %else %do;
            define HR_CI / display 'Derived HR|(95% CI)'
                style(column)=[just=center];
        %end;

    run;

    footnote;
    title;
    title2;


    /*-------------------------------------------------
      9. Baseline event listing
    -------------------------------------------------*/

    title "Participants Meeting &label. at Baseline";
    title2 "Continuous Biomarker: &biomarker_label.";

    proc print data=&data. noobs;
        where &outcome. = 1 and &time. = 0;
        var &baseline_vars.;
    run;

    title;
    title2;


    /*-------------------------------------------------
      10. BL vs Post-BL SGPLOT
    -------------------------------------------------*/

    proc format;
        value bl 
            1 = "BL"
            2 = "Post-BL";
    run;

    data temp2;
        set &data.;

        if &outcome. = 1 then do;
            if &time. = 0 then BL = 1;
            else if &time. > 0 then BL = 2;
        end;

        if BL in (1, 2);
    run;

    proc sql noprint;
        select distinct BL 
            into :bl_list separated by ' '
        from temp2
        where BL in (1, 2);
    quit;

    %if %length(&bl_list) %then %do;

        title "Event Timing (BL vs Post-BL) - &label. &title_suffix";
        title2 "Continuous Biomarker: &biomarker_label.";

        proc sgplot data=temp2(where=(BL in (&bl_list)));
            vbar BL / stat=freq;
            format BL bl.;
            xaxis display=(nolabel);
            yaxis label="Frequency";
        run;

    %end;
    %else %do;
        %put NOTE: No BL=1 or BL=2 values present - skipping BL vs Post-BL plot.;
    %end;

    title;
    title2;


    /*-------------------------------------------------
      11. Cleanup temporary datasets
    -------------------------------------------------*/

    proc datasets library=work nolist;
        delete Temp 
               OverallBiomarker 
               ModelPE 
               ModelPE_Tidy 
               ModelStats
               ModelStats_PDF
               Type3PE 
               HRPE 
               CovBPE
               temp2;
    quit;

%mend SurvCont;


/*=========================================================
  Summary report macro

  Important update:
      keeps only successful continuous biomarker rows where
      P_value is not missing.
=========================================================*/

%macro SurvContSummary(
    biomarker_label=,
    combined_out=,
    title_suffix=for Healthy Control NSD-
);

    proc sql noprint;
        select memname
            into :_mslist separated by ' '
        from dictionary.tables
        where libname = 'WORK'
          and upcase(substr(memname,1,3)) = 'MS_'
        order by crdate;
    quit;

    %if %superq(_mslist)= %then %do;
        %put WARNING: No successful MS_ datasets found in WORK. No summary report displayed for &biomarker_label..;
        %return;
    %end;

    data All_Milestones;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set &_mslist.;
    run;

    data Summary_Data;
        set All_Milestones;

        length SortOrder 8;

        if index(upcase(Parameter), "CONTINUOUS") and not missing(P_value) then SortOrder = 1;
        else delete;
    run;

    proc sql noprint;
        select count(*) into :_summary_n trimmed
        from Summary_Data;
    quit;

    %if %sysevalf(&_summary_n. = 0) %then %do;
        %put WARNING: No non-missing continuous biomarker p-values found. No summary report displayed for &biomarker_label..;
        %return;
    %end;

    proc sort data=Summary_Data;
        by Biomarker Milestone SortOrder;
    run;

    %if not %sysfunc(exist(&combined_out.)) %then %do;
        data &combined_out.;
            set Summary_Data(obs=0);
        run;
    %end;

    proc append base=&combined_out. data=Summary_Data force;
    run;

    %if &write_individual_outputs = 1 %then %do;

    options 
        orientation=landscape
        leftmargin=0.25in
        rightmargin=0.25in
        topmargin=0.40in
        bottommargin=0.40in
        nodate
        nonumber
    ;


    title "Summary of Continuous &biomarker_label. Survival Model Statistics - &title_suffix";

    footnote1 j=left h=7pt
        "Only milestones with available non-missing overall biomarker p-values are included.";
    footnote2 j=left h=7pt
        "Cox PH models report HR per 1-unit increase in biomarker. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Models adjusted for age and sex.";

    proc report data=Summary_Data nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=1 cellpadding=1]
        style(header)=[font_size=6pt font_weight=bold]
        style(column)=[font_size=6pt];

        columns Biomarker Milestone Method Beta_CI SE Z_stat P_value HR_CI;

        define Biomarker / group 'Biomarker'
            style(column)=[cellwidth=1.05in];

        define Milestone / group 'Milestone'
            style(column)=[cellwidth=2.00in];

        define Method / group 'Model'
            style(column)=[cellwidth=0.45in];

        define Beta_CI / display 'Beta|(95% CI)'
            style(column)=[cellwidth=1.55in];

        define SE / display 'SE' format=8.4
            style(column)=[cellwidth=0.55in];

        define Z_stat / display 'Z' format=8.4
            style(column)=[cellwidth=0.55in];

        define P_value / display 'p-value' format=pvalue6.4
            style(column)=[cellwidth=0.60in];

        define HR_CI / display 'HR / Derived HR|(95% CI)'
            style(column)=[cellwidth=1.70in];

        break after Milestone / skip;
    run;

    footnote;
    title;



    %end;

%mend SurvContSummary;


/*=========================================================
  Combined HC NSD- summary report
=========================================================*/

%macro WriteCombinedHCRTF(
    data=,
    title_text=
);

    %if not %sysfunc(exist(&data.)) %then %do;
        %put WARNING: &data. does not exist. No combined report displayed.;
        %return;
    %end;

    proc sql noprint;
        select count(*) into :_combined_n trimmed
        from &data.;
    quit;

    %if %sysevalf(&_combined_n. = 0) %then %do;
        %put WARNING: &data. is empty. No combined report displayed.;
        %return;
    %end;

    data _Combined_HC_Report;
        set &data.;

        if strip(Biomarker) = "pTau217 plasma" then BiomarkerOrder = 1;
        else if strip(Biomarker) = "NFL" then BiomarkerOrder = 2;
        else if strip(Biomarker) = "GFAP" then BiomarkerOrder = 3;
        else BiomarkerOrder = 99;

        if strip(Milestone) = "Any progression milestone" then MilestoneOrder = 1;
        else if strip(Milestone) = "Walking and balance domain" then MilestoneOrder = 2;
        else if strip(Milestone) = "Motor complications domain" then MilestoneOrder = 3;
        else if strip(Milestone) = "Cognitive domain" then MilestoneOrder = 4;
        else if strip(Milestone) = "Autonomic dysfunction domain" then MilestoneOrder = 5;
        else if strip(Milestone) = "Functional dependence domain" then MilestoneOrder = 6;
        else if strip(Milestone) = "Activities of daily living domain" then MilestoneOrder = 7;
        else if strip(Milestone) = "Time to HY >= 3" then MilestoneOrder = 1;
        else if strip(Milestone) = "Time to Cognitive State Change" then MilestoneOrder = 2;
        else if strip(Milestone) = "Time to MoCA < 21" then MilestoneOrder = 3;
        else if strip(Milestone) = "Time to MoCA < 26" then MilestoneOrder = 4;
        else if strip(Milestone) = "Time to MDS-UPDRS III Change 5+" then MilestoneOrder = 5;
        else MilestoneOrder = 99;
    run;

    proc sort data=_Combined_HC_Report;
        by BiomarkerOrder MilestoneOrder SortOrder;
    run;

    options
        orientation=landscape
        leftmargin=0.25in
        rightmargin=0.25in
        topmargin=0.40in
        bottommargin=0.40in
        nodate
        nonumber;


    title "&title_text.";

    footnote1 j=left h=7pt
        "Only milestones with available non-missing overall biomarker p-values are included.";
    footnote2 j=left h=7pt
        "Cox PH models report HR per 1-unit increase in biomarker. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Models adjusted for age and sex.";

    proc report data=_Combined_HC_Report nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=1 cellpadding=1]
        style(header)=[font_size=6pt font_weight=bold]
        style(column)=[font_size=6pt];

        columns Biomarker Milestone Method Beta_CI SE Z_stat P_value HR_CI;

        define Biomarker / group order=data 'Biomarker'
            style(column)=[cellwidth=1.05in];
        define Milestone / group order=data 'Milestone'
            style(column)=[cellwidth=2.00in];
        define Method / group 'Model'
            style(column)=[cellwidth=0.45in];
        define Beta_CI / display 'Beta|(95% CI)'
            style(column)=[cellwidth=1.55in];
        define SE / display 'SE' format=8.4
            style(column)=[cellwidth=0.55in];
        define Z_stat / display 'Z' format=8.4
            style(column)=[cellwidth=0.55in];
        define P_value / display 'p-value' format=pvalue6.4
            style(column)=[cellwidth=0.60in];
        define HR_CI / display 'HR / Derived HR|(95% CI)'
            style(column)=[cellwidth=1.70in];

        compute P_value;
            if not missing(P_value) and P_value < 0.05 then
                call define(_col_, "style", "style=[font_weight=bold foreground=red]");
        endcomp;

        break after Biomarker / skip;
    run;

    footnote;
    title;


    proc datasets library=work nolist;
        delete _Combined_HC_Report;
    quit;

%mend WriteCombinedHCRTF;


/*=========================================================
  Main HC NSD- biomarker-level macro
=========================================================*/

%macro RunOneBiomarker_HC_NSDminus(
    biomarker_short=,
    biomarker_label=,
    datevar=,
    biomarker_var=
);

    %put NOTE: ============================================================;
    %put NOTE: Running HC NSD- continuous survival analysis for &biomarker_label.;
    %put NOTE: Date variable      = &datevar.;
    %put NOTE: Biomarker variable = &biomarker_var.;
    %put NOTE: ============================================================;

    /*-----------------------------------------------------
      Create biomarker-specific source dataset for HC
    -----------------------------------------------------*/

    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data SurvData_&biomarker_short.;
        set progression_analysis_v3;
        where Analysis_Source="HC_C"
          and Analysis_Biomarker="&biomarker_short.";
    run;



    proc means data=SurvData_&biomarker_short. n mean std median q1 q3 min max maxdec=4;
        var Early_Biomarker age;
        title "Distribution of Continuous &biomarker_label. and Age - Healthy Control NSD-";
    run;

    title;


    /*-----------------------------------------------------
      Initialize combined model stats for this biomarker
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete All_ModelStats MS_: All_Milestones Summary_Data;
    quit;

    data All_ModelStats;
        length 
            Biomarker $100
            Milestone $100
            Method $10
            Parameter $100
            ClassVal $40
            Beta_CI $40
            HR_CI $50
        ;

        call missing(
            Biomarker,
            Milestone,
            Method,
            Parameter,
            ClassVal,
            GroupDF,
            WaldStat,
            Beta,
            SE,
            Beta_CI,
            Z_stat,
            P_value,
            HR,
            HR_LCL,
            HR_UCL,
            HR_CI
        );

        stop;
    run;


    /*-----------------------------------------------------
      Prepare biomarker-specific analysis title
    -----------------------------------------------------*/

    options 
        orientation=landscape
        leftmargin=0.35in
        rightmargin=0.35in
        topmargin=0.40in
        bottommargin=0.40in
        nodate
        nonumber
    ;

    ods graphics / reset width=9in height=5in;

    %if &write_individual_outputs = 1 %then %do;

        title "Continuous &biomarker_label. Survival Analysis for Healthy Control NSD-";
    %end;


    /*-----------------------------------------------------
      Run all progression milestones
      Skipped automatically if overall biomarker p-value is missing.
    -----------------------------------------------------*/

    %SurvCont(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_any,
        time=TIME_ENRL_any,
        biomarker=Early_Biomarker,
        label=Any progression milestone,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_any TIME_ENRL_any
    );

    %SurvCont(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_wb_any,
        time=TIME_ENRL_wb_any,
        biomarker=Early_Biomarker,
        label=Walking and balance domain,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_wb_any TIME_ENRL_wb_any
    );

    %SurvCont(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_mc_any,
        time=TIME_ENRL_mc_any,
        biomarker=Early_Biomarker,
        label=Motor complications domain,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_mc_any TIME_ENRL_mc_any
    );

    %SurvCont(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_cog_any,
        time=TIME_ENRL_cog_any,
        biomarker=Early_Biomarker,
        label=Cognitive domain,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_cog_any TIME_ENRL_cog_any
    );

    %SurvCont(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_auto_any,
        time=TIME_ENRL_auto_any,
        biomarker=Early_Biomarker,
        label=Autonomic dysfunction domain,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_auto_any TIME_ENRL_auto_any
    );

    %SurvCont(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_fd_any,
        time=TIME_ENRL_fd_any,
        biomarker=Early_Biomarker,
        label=Functional dependence domain,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_fd_any TIME_ENRL_fd_any
    );

    %SurvCont(
        data=SurvData_&biomarker_short.,
        outcome=STATUS_adl_any,
        time=TIME_ENRL_adl_any,
        biomarker=Early_Biomarker,
        label=Activities of daily living domain,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_adl_any TIME_ENRL_adl_any
    );

    %if &write_individual_outputs = 1 %then %do;

        title;
    %end;


    /*-----------------------------------------------------
      Display biomarker-specific summary report
    -----------------------------------------------------*/

    %SurvContSummary(
        biomarker_label=&biomarker_label.,
        combined_out=Combined_HC_Progression_Summary,
        title_suffix=for Healthy Control NSD-
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
               MS_:;
    quit;

%mend RunOneBiomarker_HC_NSDminus;


/*=========================================================
  Run all biomarkers for HC NSD-
=========================================================*/

%RunOneBiomarker_HC_NSDminus(
    biomarker_short=Lilly_pTau217,
    biomarker_label=pTau217 plasma,
    datevar=Lilly_ptau217p_date,
    biomarker_var=Lilly_ptau217p
);

%RunOneBiomarker_HC_NSDminus(
    biomarker_short=IU_GFAP,
    biomarker_label=GFAP,
    datevar=IU_Plasma_date,
    biomarker_var=IU_GFAP_plasma_mean
);

%RunOneBiomarker_HC_NSDminus(
    biomarker_short=IU_NFL,
    biomarker_label=NFL,
    datevar=IU_Plasma_date,
    biomarker_var=IU_NFL_plasma_mean
);

%WriteCombinedHCRTF(
    data=Combined_HC_Progression_Summary,
    title_text=Combined Continuous Biomarker Progression Milestone Survival Model Summary for Healthy Control NSD-
);



/*==========================================================================
  HC NSD- CLINICAL MILESTONE PROGRAM
  This section follows the progression program because it uses the macros
  and work datasets established above.
==========================================================================*/

/*===============================================================
  HC NSD- CLINICAL MILESTONES - CONTINUOUS BIOMARKER AUTOMATION

  This code runs clinical milestone survival models for:
      1. Lilly pTau217
      2. IU GFAP
      3. IU NfL
      4. C2N AB40
      5. C2N AB42
      6. C2N AB42/40

  Cohort:
      Healthy Control NSD-

  Rule:
      If any milestone has missing/absent overall biomarker p-value,
      that milestone is skipped from displayed reports by %SurvCont.
================================================================*/

%let report_date = 08DEC2025;



/*===============================================================
  Helper macro:
  Select first event or last non-event for each clinical outcome
================================================================*/

%macro MakeClinicalTTE(data=, event=, time=, visitvar=Event_ID, extra=);

    proc sort data=&data. out=survival1;
        where &event. = 1;
        by Patno &time.;
    run;

    proc sort data=survival1 nodupkey;
        by patno;
    run;

    proc sort data=&data. out=survival0;
        where &event. = 0 &extra.;
        by Patno descending &time.;
    run;

    proc sort data=survival0 nodupkey;
        by patno;
    run;

    data survival0b;
        merge survival0(in=in1) survival1(in=in2 keep=patno);
        by patno;

        if in2 then delete;
    run;

    data &event.;
        merge survival1(in=in1) survival0b(in=in2);
        by patno;

        if in1 or in2;

        keep patno &visitvar. &event. &time.;
    run;

    proc datasets library=work nolist;
        delete survival1 survival0 survival0b;
    quit;

%mend MakeClinicalTTE;


/*===============================================================
  Main macro:
  Build HC NSD- clinical survival dataset for one biomarker,
  then run the clinical milestones.
================================================================*/

%macro RunOneClinicalBiomarker(
    biomarker_short=,
    biomarker_label=,
    datevar=,
    biomarker_var=
);


    
    

    %put NOTE: ============================================================;
    %put NOTE: Running HC NSD- clinical milestone survival analysis for &biomarker_label.;
    %put NOTE: Date variable      = &datevar.;
    %put NOTE: Biomarker variable = &biomarker_var.;
    %put NOTE: ============================================================;


    /*-----------------------------------------------------
      Biomarker-specific source dataset

      Standardized continuous biomarker:
          Early_Biomarker
    -----------------------------------------------------*/

    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data SurvDataClinical_&biomarker_short.;
        set clinical_analysis_v3;
        where Analysis_Source="HC_C"
          and Analysis_Biomarker="&biomarker_short.";
    run;



    /*-----------------------------------------------------
      Distribution check
    -----------------------------------------------------*/

    proc means data=SurvDataClinical_&biomarker_short. n mean std median q1 q3 min max maxdec=4;
        var Early_Biomarker age;
        title "Distribution of Continuous &biomarker_label. and Age - HC NSD- Clinical Milestone Dataset";
    run;

    title;


    /*-----------------------------------------------------
      Initialize model stats for this biomarker
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete All_ModelStats MS_: All_Milestones Summary_Data;
    quit;

    data All_ModelStats;
        length
            Biomarker $100
            Milestone $100
            Method $10
            Parameter $100
            ClassVal $40
            Beta_CI $40
            HR_CI $50
        ;

        call missing(
            Biomarker,
            Milestone,
            Method,
            Parameter,
            ClassVal,
            GroupDF,
            WaldStat,
            Beta,
            SE,
            Beta_CI,
            Z_stat,
            P_value,
            HR,
            HR_LCL,
            HR_UCL,
            HR_CI
        );

        stop;
    run;


    /*-----------------------------------------------------
      Prepare biomarker-specific analysis title
    -----------------------------------------------------*/

    options
        orientation=landscape
        leftmargin=0.35in
        rightmargin=0.35in
        topmargin=0.40in
        bottommargin=0.40in
        nodate
        nonumber;

    ods graphics / reset width=9in height=5in;

    %if &write_individual_outputs = 1 %then %do;

        title "Continuous &biomarker_label. Clinical Milestone Survival Analysis for Healthy Control NSD-";
    %end;


    /*-----------------------------------------------------
      Run clinical milestones

      The %SurvCont macro should skip any milestone with:
          - no events,
          - missing model output,
          - missing overall biomarker p-value.
    -----------------------------------------------------*/

    %SurvCont(
        data=SurvDataClinical_&biomarker_short.,
        outcome=STATUS_HY_GE3,
        time=TIME_ENRL_HY_GE3,
        biomarker=Early_Biomarker,
        label=Time to HY >= 3,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_HY_GE3 TIME_ENRL_HY_GE3
    );

    %SurvCont(
        data=SurvDataClinical_&biomarker_short.,
        outcome=STATUS_COGSTATE_CHANGE,
        time=TIME_ENRL_COGSTATE_CHANGE,
        biomarker=Early_Biomarker,
        label=Time to Cognitive State Change,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_COGSTATE_CHANGE TIME_ENRL_COGSTATE_CHANGE
    );

    %SurvCont(
        data=SurvDataClinical_&biomarker_short.,
        outcome=STATUS_MOCA_LT21,
        time=TIME_ENRL_MOCA_LT_21,
        biomarker=Early_Biomarker,
        label=Time to MoCA < 21,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_MOCA_LT21 TIME_ENRL_MOCA_LT_21
    );

    %SurvCont(
        data=SurvDataClinical_&biomarker_short.,
        outcome=STATUS_MOCA_LT26,
        time=TIME_ENRL_MOCA_LT_26,
        biomarker=Early_Biomarker,
        label=Time to MoCA < 26,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_MOCA_LT26 TIME_ENRL_MOCA_LT_26
    );

    %SurvCont(
        data=SurvDataClinical_&biomarker_short.,
        outcome=STATUS_UPDRS_Change_GT5,
        time=TIME_ENRL_UPDRS_Change_GT5,
        biomarker=Early_Biomarker,
        label=Time to MDS-UPDRS III Change 5+,
        method=AUTO,
        covars=age gen,
        classvars=gen,
        biomarker_label=&biomarker_label.,
        title_suffix=for Healthy Control NSD-,
        baseline_vars=Patno Early_Biomarker subgroup_definition STATUS_UPDRS_Change_GT5 TIME_ENRL_UPDRS_Change_GT5
    );


    %if &write_individual_outputs = 1 %then %do;

        title;
    %end;


    /*-----------------------------------------------------
      Summary report for this biomarker
    -----------------------------------------------------*/

    %SurvContSummary(
        biomarker_label=&biomarker_label.,
        combined_out=Combined_HC_Clinical_Summary,
        title_suffix=for Healthy Control NSD-
    );


    /*-----------------------------------------------------
      Biomarker-specific cleanup
    -----------------------------------------------------*/

    proc datasets library=work nolist;
        delete mydata_&biomarker_short.
               EARLY_&biomarker_short.
               datalong_&biomarker_short.
               lastobsvisit_&biomarker_short.
               allvisits_&biomarker_short.
               ClinicalBase_&biomarker_short.
               start_led_&biomarker_short.
               SurvDataClinical_&biomarker_short.
               STATUS_HY_GE3
               STATUS_TOFF_FUNCFLUC_GE3
               STATUS_MOCA_LT21
               STATUS_MOCA_LT26
               STATUS_UPDRS_Change_GT5
               STATUS_DAT_75
               STATUS_DAT2_75
               STATUS_STAGE_CHANGE
               STATUS_STAGE_CHANGE_3
               STATUS_PRIMDIAG_CHANGE
               STATUS_COGSTATE_CHANGE
               All_ModelStats
               All_Milestones
               Summary_Data
               MS_:;
    quit;

%mend RunOneClinicalBiomarker;


/*===============================================================
  Run HC NSD- clinical milestone models for all biomarkers
================================================================*/

%RunOneClinicalBiomarker(
    biomarker_short=Lilly_pTau217,
    biomarker_label=pTau217 plasma,
    datevar=Lilly_ptau217p_date,
    biomarker_var=Lilly_ptau217p
);

%RunOneClinicalBiomarker(
    biomarker_short=IU_GFAP,
    biomarker_label=GFAP,
    datevar=IU_Plasma_date,
    biomarker_var=IU_GFAP_plasma_mean
);

%RunOneClinicalBiomarker(
    biomarker_short=IU_NFL,
    biomarker_label=NFL,
    datevar=IU_Plasma_date,
    biomarker_var=IU_NFL_plasma_mean
);

%WriteCombinedHCRTF(
    data=Combined_HC_Clinical_Summary,
    title_text=Combined Continuous Biomarker Clinical Milestone Survival Model Summary for Healthy Control NSD-
);

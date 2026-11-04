/*=========================================================
  sPD NSD+: Continuous Biomarkers and Survival Milestones

  Uses the prepared progression_analysis_v3 and
  clinical_analysis_v3 datasets to examine associations
  between continuous plasma biomarker levels and time to
  progression and clinical milestones in sPD NSD+ participants.

  Fits age- and sex-adjusted Cox proportional hazards or
  Weibull accelerated failure time models, selected according
  to the proportional hazards assumption test, and generates
  separate progression and clinical survival summary tables.
=========================================================*/


%let filepath = /*path where analysis data are located*/;
libname analysis "&filepath.";

/*==========================================================================
  COMBINED sPD PROGRESSION AND CLINICAL MILESTONE SURVIVAL PROGRAMS

  Execution order is intentional:
    1. Progression milestone program
    2. Clinical milestone program
==========================================================================*/

%let report_date = 08DEC2025;

data progression_analysis_v3;
    set analysis.progression_analysis_v3;
run;

data clinical_analysis_v3;
    set analysis.clinical_analysis_v3;
run;
%let write_individual_outputs = 0;

ods escapechar="^";


/*---------------------------------------------------------
  Initialize combined all-biomarker summary dataset
---------------------------------------------------------*/

proc datasets library=work nolist;
    delete Combined_All_Biomarker_Summary;
quit;

data Combined_All_Biomarker_Summary;
    length 
        Biomarker $100
        Milestone $100
        Method $10
        Beta_CI_3 $50
        SE_3 $20
        Z_3 $20
        P_3 $20
        HR_CI_3 $60
    ;

    length P_Sig 8;

    call missing(
        Biomarker,
        Milestone,
        Method,
        Beta_CI_3,
        SE_3,
        Z_3,
        P_3,
        HR_CI_3,
        SortOrder,
        P_Sig
    );

    stop;
run;


/*---------------------------------------------------------
  Continuous survival macro

  This macro runs ONE milestone for ONE biomarker.
---------------------------------------------------------*/

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
    title_suffix=for Sporadic PD,
    baseline_vars=
);

    %local _method _ph_fail _modelLabel _dslabel;
    %local sigma var_sigma cov_bio_sigma var_bio_beta;

    %let sigma=.;
    %let var_sigma=.;
    %let cov_bio_sigma=.;
    %let var_bio_beta=.;

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
            assess ph / resample seed=1472160691;
        run;

        ods results;
        ods exclude none;

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

        proc datasets lib=work nolist;
            modify Type3PE;
            rename Variable = Effect;
        quit;


        data _null_;
            set ModelPE;

            if upcase(Parameter) = "SCALE" then 
                call symputx('sigma', Estimate);
        run;


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

        %if %superq(var_sigma)= %then %let var_sigma=.;
        %if %superq(cov_bio_sigma)= %then %let cov_bio_sigma=.;
        %if %superq(var_bio_beta)= %then %let var_bio_beta=.;

    %end;
    %else %do;

        %put WARNING: method=&_method. not recognized. Falling back to PHREG.;
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
      4. Model label
    -------------------------------------------------*/

    %let _modelLabel = %sysfunc(ifc(&_method.=LIFEREG,AFT,CoxPH));


    /*-------------------------------------------------
      5. Overall continuous biomarker test row
    -------------------------------------------------*/

    data OverallBiomarker;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set Type3PE;

        where upcase(strip(Effect)) = upcase("&biomarker.");

        Biomarker = "&biomarker_label.";
        Milestone = "&label.";
        Method    = "&_modelLabel";
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


    /*-------------------------------------------------
      6. Parameter estimates table
    -------------------------------------------------*/

    data ModelPE_Tidy;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set ModelPE;

        Biomarker = "&biomarker_label.";
        Milestone = "&label.";
        Method    = "&_modelLabel";

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
                Beta_CI = cats(put(Beta,8.4), " (",
                               put(Beta - 1.96*SE,8.4), ", ",
                               put(Beta + 1.96*SE,8.4), ")");
            else 
                Beta_CI = "";

            if upcase(strip(Parameter)) = upcase("&biomarker_label. CONTINUOUS") then do;
                HR     = exp(Beta);
                HR_LCL = exp(Beta - 1.96*SE);
                HR_UCL = exp(Beta + 1.96*SE);

                HR_CI = cats(put(HR,8.4), " (",
                             put(HR_LCL,8.4), ", ",
                             put(HR_UCL,8.4), ")");
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
                Beta_CI = cats(put(Beta,8.4), " (",
                               put(Beta - 1.96*SE,8.4), ", ",
                               put(Beta + 1.96*SE,8.4), ")");
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

                        HR_CI = cats(put(HR,8.4), " (",
                                     put(HR_LCL,8.4), ", ",
                                     put(HR_UCL,8.4), ")");
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

  Keeps:
      - Intercept
      - age
      - gen
      - continuous biomarker row

  Removes from report display:
      - Biomarker column
      - Milestone column
      - Overall biomarker row

  This prevents PROC REPORT from splitting the table.
-------------------------------------------------*/

data ModelStats_PDF;
    length Parameter_Display $100;
    set ModelStats;

    /* Remove the repeated milestone column from the table body;
       milestone is already in the title. */
    Parameter_Display = Parameter;

    /* Shorten the continuous biomarker row label */
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
%macro SurvContSummary(
    biomarker_label=,
    title_suffix=for Sporadic PD
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
        %put WARNING: No MS_ datasets found in WORK. Run %SurvCont calls first.;
        %return;
    %end;

    data All_Milestones;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set &_mslist.;
    run;

    data Summary_Data_raw;
        set All_Milestones;

        length SortOrder 8;

        /* Keep only continuous biomarker estimate row */
        if index(upcase(Parameter), "CONTINUOUS") = 0 then delete;

        if Milestone = "Any progression milestone" then SortOrder = 1;
        else if Milestone = "Walking and balance domain" then SortOrder = 2;
        else if Milestone = "Motor complications domain" then SortOrder = 3;
        else if Milestone = "Cognitive domain" then SortOrder = 4;
        else if Milestone = "Autonomic dysfunction domain" then SortOrder = 5;
        else if Milestone = "Functional dependence domain" then SortOrder = 6;
        else delete;
    run;

    data Summary_Data;
        set Summary_Data_raw;

        length Beta_CI_3 $50 SE_3 $20 Z_3 $20 P_3 $20 HR_CI_3 $60;

        P_Sig = (not missing(P_value) and P_value < 0.05);

        /* Beta CI rounded to 3 decimals */
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

        /* SE, Z, p-value rounded to 3 decimals */
        if not missing(SE) then SE_3 = put(SE, 8.3);
        else SE_3 = "";

        if not missing(Z_stat) then Z_3 = put(Z_stat, 8.3);
        else Z_3 = "";

        if not missing(P_value) then do;
            if P_value < 0.001 then P_3 = "<0.001";
            else P_3 = put(P_value, 8.3);
        end;
        else P_3 = "";

        /* HR CI rounded to 3 decimals */
        if not missing(HR) and not missing(HR_LCL) and not missing(HR_UCL) then do;
            HR_CI_3 = catx(" ",
                       put(HR, 8.3),
                       cats("(",
                            put(HR_LCL, 8.3), ", ",
                            put(HR_UCL, 8.3), ")"));
        end;
        else HR_CI_3 = "";

        keep Biomarker Milestone Method SortOrder Beta_CI_3 SE_3 Z_3 P_3 P_Sig HR_CI_3;
    run;

    proc sort data=Summary_Data;
        by Biomarker SortOrder;
    run;

    /* Append to combined all-biomarker summary */
    proc append base=Combined_All_Biomarker_Summary data=Summary_Data force;
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
        "Only continuous biomarker estimate rows are shown. Values are rounded to 3 decimals.";
    footnote2 j=left h=7pt
        "Cox PH models report HR per 1-unit increase in biomarker. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Models adjusted for age and sex.";

    proc report data=Summary_Data nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=1 cellpadding=2]
        style(header)=[font_size=7pt font_weight=bold]
        style(column)=[font_size=7pt];

        columns Biomarker Milestone Method Beta_CI_3 SE_3 Z_3 P_3 HR_CI_3;

        define Biomarker / group 'Biomarker';
        define Milestone / display 'Milestone';
        define Method    / display 'Model';
        define Beta_CI_3 / display 'Beta (95% CI)';
        define SE_3      / display 'SE';
        define Z_3       / display 'Z';
        define P_3       / display 'p-value';
        define HR_CI_3   / display 'HR / Derived HR (95% CI)';

    run;

    footnote;
    title;



    %end;

%mend SurvContSummary;
/*---------------------------------------------------------
  Main biomarker-level macro

  This macro:
      1. creates biomarker-specific SurvData
      2. initializes All_ModelStats
      3. prepares the analysis title
      4. runs all milestones
      5. clears the analysis title
      6. prepares one summary report
      7. cleans biomarker-specific work datasets
---------------------------------------------------------*/

%macro RunOneBiomarker(
    biomarker_short=,
    biomarker_label=,
    datevar=,
    biomarker_var=
);


    
    


    %put NOTE: ============================================================;
    %put NOTE: Running continuous survival analysis for &biomarker_label.;
    %put NOTE: Date variable      = &datevar.;
    %put NOTE: Biomarker variable = &biomarker_var.;
    %put NOTE: ============================================================;


    /*-----------------------------------------------------
      Create biomarker-specific source dataset
    -----------------------------------------------------*/

    /* Read the saved analysis input. Existing model calls follow unchanged. */
    data SurvData_&biomarker_short.;
        set progression_analysis_v3;
        where Analysis_Source="SP_C"
          and Analysis_Biomarker="&biomarker_short.";
    run;



    /*-----------------------------------------------------
      Distribution check
    -----------------------------------------------------*/

    proc means data=SurvData_&biomarker_short. n mean std median q1 q3 min max maxdec=4;
        var Early_Biomarker age;
        title "Distribution of Continuous &biomarker_label. and Age";
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

    %if &write_individual_outputs = 1 %then %do;

        title "Continuous &biomarker_label. Survival Analysis for Sporadic PD NSD+";
    %end;


    /*-----------------------------------------------------
      Run all milestones

      I set method=AUTO for all biomarkers because PH behavior can
      differ across biomarkers. If PH is violated, the macro uses
      Weibull AFT; otherwise, it uses Cox PH.
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
           Summary_Data_raw
           MS_:;
quit;

%mend RunOneBiomarker;

%macro WriteCombinedBiomarkerRTF(
    title_suffix=for Sporadic PD
);

    proc sql noprint;
        select count(*)
            into :_n_combined trimmed
        from Combined_All_Biomarker_Summary;
    quit;

    %if %sysevalf(&_n_combined. = 0) %then %do;
        %put WARNING: Combined_All_Biomarker_Summary is empty. No combined report displayed.;
        %return;
    %end;

    data Combined_All_Biomarker_Summary;
        set Combined_All_Biomarker_Summary;

        if strip(Biomarker) = "pTau217 plasma" then BiomarkerOrder = 1;
        else if strip(Biomarker) = "NFL" then BiomarkerOrder = 2;
        else if strip(Biomarker) = "GFAP" then BiomarkerOrder = 3;
        else if strip(Biomarker) = "A^{unicode 03B2}40 plasma" then BiomarkerOrder = 4;
        else if strip(Biomarker) = "A^{unicode 03B2}42 plasma" then BiomarkerOrder = 5;
        else if strip(Biomarker) = "A^{unicode 03B2}42/A^{unicode 03B2}40 plasma" then BiomarkerOrder = 6;
        else BiomarkerOrder = 99;
    run;

    proc sort data=Combined_All_Biomarker_Summary;
        by BiomarkerOrder SortOrder;
    run;

    options 
        orientation=landscape
        leftmargin=0.25in
        rightmargin=0.25in
        topmargin=0.40in
        bottommargin=0.40in
        nodate
        nonumber
    ;


    title "Combined Continuous Biomarker Survival Model Summary - &title_suffix";

    footnote1 j=left h=7pt
        "Only continuous biomarker estimate rows are shown. Values are rounded to 3 decimals.";
    footnote2 j=left h=7pt
        "Milestones are ordered as: Any, Walking/balance, Motor complications, Cognitive, Autonomic dysfunction, Functional dependence.";
    footnote3 j=left h=7pt
        "Cox PH models report HR per 1-unit increase in biomarker. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Models adjusted for age and sex.";

    proc report data=Combined_All_Biomarker_Summary nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=1 cellpadding=2]
        style(header)=[font_size=7pt font_weight=bold]
        style(column)=[font_size=7pt];

        columns Biomarker Milestone Method Beta_CI_3 SE_3 Z_3 P_Sig P_3 HR_CI_3;

        define Biomarker / group order=data 'Biomarker';
        define Milestone / display 'Milestone';
        define Method    / display 'Model';
        define Beta_CI_3 / display 'Beta (95% CI)';
        define SE_3      / display 'SE';
        define Z_3       / display 'Z';
        define P_Sig     / display noprint;
        define P_3       / display 'p-value';
        define HR_CI_3   / display 'HR / Derived HR (95% CI)';

        compute P_3;
            if P_Sig = 1 then
                call define(_col_, "style", "style=[font_weight=bold foreground=red]");
        endcomp;

        break after Biomarker / skip;

    run;

    footnote;
    title;



%mend WriteCombinedBiomarkerRTF;
/*---------------------------------------------------------
  Run all biomarkers

  Date variable rules:
      Lilly pTau217: Lilly_ptau217p_date
      IU GFAP/NfL:  IU_Plasma_date
      C2N markers:  C2N_AB_date
---------------------------------------------------------*/

%RunOneBiomarker(
    biomarker_short=Lilly_pTau217,
    biomarker_label=pTau217 plasma,
    datevar=Lilly_ptau217p_date,
    biomarker_var=Lilly_ptau217p
);


%RunOneBiomarker(
    biomarker_short=IU_GFAP,
    biomarker_label=GFAP,
    datevar=IU_Plasma_date,
    biomarker_var=IU_GFAP_plasma_mean
);


%RunOneBiomarker(
    biomarker_short=IU_NFL,
    biomarker_label=NFL,
    datevar=IU_Plasma_date,
    biomarker_var=IU_NFL_plasma_mean
);


%RunOneBiomarker(
    biomarker_short=C2N_AB40,
    biomarker_label=A^{unicode 03B2}40 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB40
);


%RunOneBiomarker(
    biomarker_short=C2N_AB42,
    biomarker_label=A^{unicode 03B2}42 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB42
);


%RunOneBiomarker(
    biomarker_short=C2N_AB4240,
    biomarker_label=A^{unicode 03B2}42/A^{unicode 03B2}40 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB4240
);

%WriteCombinedBiomarkerRTF(
    title_suffix=for Sporadic PD
);
/*---------------------------------------------------------
  Final cleanup
---------------------------------------------------------*/

proc datasets library=work nolist;
    delete BL1 combined:;
quit;


/*==========================================================================
  CLINICAL MILESTONE PROGRAM
  This section follows the progression program because it uses the macros
  and work datasets established above.
==========================================================================*/

/*===============================================================
  CLINICAL MILESTONES - CONTINUOUS BIOMARKER AUTOMATION

  This code runs clinical milestone survival models for:
      1. Lilly pTau217
      2. IU GFAP
      3. IU NfL
      4. C2N AB40
      5. C2N AB42
      6. C2N AB42/40

  Biomarker-specific report code is retained with the original display switch.

  IMPORTANT:
      This assumes %SurvCont and %SurvContSummary are already defined
      from the continuous progression milestone code.
================================================================*/


%let report_date = 08DEC2025;
%let write_individual_outputs = 0;


ods escapechar="^";

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

/*---------------------------------------------------------
  Initialize combined all-biomarker clinical summary dataset
---------------------------------------------------------*/

proc datasets library=work nolist;
    delete Comb_Clin_Summary;
quit;

data Comb_Clin_Summary;
    length 
        Biomarker $100
        Milestone $100
        Method $10
        Beta_CI_3 $50
        SE_3 $20
        Z_3 $20
        P_3 $20
        HR_CI_3 $60
    ;

    length P_Sig 8;

    call missing(
        Biomarker,
        Milestone,
        Method,
        Beta_CI_3,
        SE_3,
        Z_3,
        P_3,
        HR_CI_3,
        SortOrder,
        P_Sig
    );

    stop;
run;

%macro SurvContSummary(
    biomarker_label=,
    title_suffix=for Sporadic PD
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
        %put WARNING: No MS_ datasets found in WORK. Run %SurvCont calls first.;
        %return;
    %end;

    data All_Milestones;
        length Biomarker $100 Milestone $100 Method $10 Parameter $100 ClassVal $40;
        length Beta_CI $40 HR_CI $50;

        set &_mslist.;
    run;

    data Summary_Data_raw;
        set All_Milestones;

        length SortOrder 8;

        /* Keep only continuous biomarker estimate row */
        if index(upcase(Parameter), "CONTINUOUS") = 0 then delete;

        /* Clinical milestone order */
        if Milestone = "Time to HY >= 3" then SortOrder = 1;
        else if Milestone = "Time to Cognitive State Change" then SortOrder = 2;
        else if Milestone = "Time to MoCA < 21" then SortOrder = 3;
        else if Milestone = "Time to MoCA < 26" then SortOrder = 4;
        else if Milestone = "Time to MDS-UPDRS III Change 5+" then SortOrder = 5;
        else delete;
    run;

    data Summary_Data;
        set Summary_Data_raw;

        length Beta_CI_3 $50 SE_3 $20 Z_3 $20 P_3 $20 HR_CI_3 $60;

        P_Sig = (not missing(P_value) and P_value < 0.05);

        /* Beta CI rounded to 3 decimals */
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

        /* HR CI rounded to 3 decimals */
        if not missing(HR) and not missing(HR_LCL) and not missing(HR_UCL) then do;
            HR_CI_3 = catx(" ",
                       put(HR, 8.3),
                       cats("(",
                            put(HR_LCL, 8.3), ", ",
                            put(HR_UCL, 8.3), ")"));
        end;
        else HR_CI_3 = "";

        keep Biomarker Milestone Method SortOrder Beta_CI_3 SE_3 Z_3 P_3 P_Sig HR_CI_3;
    run;

    proc sort data=Summary_Data;
        by Biomarker SortOrder;
    run;

    /* Append to combined all-biomarker clinical summary */
proc append base=Comb_Clin_Summary data=Summary_Data force;
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


    title "Summary of Continuous &biomarker_label. Clinical Milestone Survival Model Statistics - &title_suffix";

    footnote1 j=left h=7pt
        "Only continuous biomarker estimate rows are shown. Values are rounded to 3 decimals.";
    footnote2 j=left h=7pt
        "Cox PH models report HR per 1-unit increase in biomarker. Weibull AFT models report derived HR using HR=exp(-beta/sigma). Models adjusted for age and sex.";

    proc report data=Summary_Data nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=1 cellpadding=2]
        style(header)=[font_size=7pt font_weight=bold]
        style(column)=[font_size=7pt];

        columns Biomarker Milestone Method Beta_CI_3 SE_3 Z_3 P_3 HR_CI_3;

        define Biomarker / group 'Biomarker';
        define Milestone / display 'Milestone';
        define Method    / display 'Model';
        define Beta_CI_3 / display 'Beta (95% CI)';
        define SE_3      / display 'SE';
        define Z_3       / display 'Z';
        define P_3       / display 'p-value';
        define HR_CI_3   / display 'HR / Derived HR (95% CI)';

    run;

    footnote;
    title;



    %end;

%mend SurvContSummary;



/*===============================================================
  Main macro:
  Build clinical survival dataset for one biomarker,
  then run all clinical milestones.
================================================================*/

%macro RunOneClinicalBiomarker(
    biomarker_short=,
    biomarker_label=,
    datevar=,
    biomarker_var=
);


    %put NOTE: ============================================================;
    %put NOTE: Running clinical milestone survival analysis for &biomarker_label.;
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
        where Analysis_Source="SP_C"
          and Analysis_Biomarker="&biomarker_short.";
    run;



    /*-----------------------------------------------------
      Distribution check
    -----------------------------------------------------*/

    proc means data=SurvDataClinical_&biomarker_short. n mean std median q1 q3 min max maxdec=4;
        var Early_Biomarker age;
        title "Distribution of Continuous &biomarker_label. and Age - Clinical Milestone Dataset";
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

        title "Continuous &biomarker_label. Clinical Milestone Survival Analysis for Sporadic PD NSD+";
    %end;


    /*-----------------------------------------------------
      Run clinical milestones

      Using method=AUTO so each biomarker/outcome decides:
          PHREG if PH is okay
          LIFEREG if PH is violated
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD,
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
        title_suffix=for Sporadic PD
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
               Summary_Data_raw
               MS_:;    
quit;

%mend RunOneClinicalBiomarker;
%macro CombinedBiomarker(
    title_suffix=for Sporadic PD
);

    proc sql noprint;
        select count(*)
            into :_n_combined trimmed
        from Comb_Clin_Summary;
    quit;

    %if %sysevalf(&_n_combined. = 0) %then %do;
        %put WARNING: Comb_Clin_Summary is empty. No combined clinical report displayed.;
        %return;
    %end;

    data Comb_Clin_Summary;
        set Comb_Clin_Summary;

        if strip(Biomarker) = "pTau217 plasma" then BiomarkerOrder = 1;
        else if strip(Biomarker) = "NFL" then BiomarkerOrder = 2;
        else if strip(Biomarker) = "GFAP" then BiomarkerOrder = 3;
        else if strip(Biomarker) = "A^{unicode 03B2}40 plasma" then BiomarkerOrder = 4;
        else if strip(Biomarker) = "A^{unicode 03B2}42 plasma" then BiomarkerOrder = 5;
        else if strip(Biomarker) = "A^{unicode 03B2}42/A^{unicode 03B2}40 plasma" then BiomarkerOrder = 6;
        else BiomarkerOrder = 99;
    run;

    proc sort data=Comb_Clin_Summary;
        by BiomarkerOrder SortOrder;
    run;

    options 
        orientation=landscape
        leftmargin=0.25in
        rightmargin=0.25in
        topmargin=0.40in
        bottommargin=0.40in
        nodate
        nonumber
    ;


    title1 j=left h=9pt
        "^{style [font_weight=bold] Supplement Table.} Association between baseline plasma biomarker levels (continuous) and time-to-event clinical milestones in";
    title2 j=left h=9pt
        "Sporadic PD (NSD+). Cox proportional hazards (CoxPH) or Weibull accelerated failure time (AFT) models were selected automatically by the";
    title3 j=left h=9pt
        "proportional hazards assumption test. All models adjusted for age and sex. Bold red p-values indicate p < 0.05.";

    proc report data=Comb_Clin_Summary nowd split='|'
        style(report)=[outputwidth=10.25in cellspacing=0 cellpadding=4 frame=hsides rules=groups]
        style(header)=[font_size=8pt font_weight=bold background=white just=center]
        style(column)=[font_size=8pt vjust=middle];

        columns Biomarker Milestone Method P_Sig Beta_CI_3 SE_3 Z_3 P_3 HR_CI_3;

        define Biomarker / group order=data 'Biomarker'
            style(column)=[cellwidth=1.10in font_weight=bold just=left vjust=top];
        define Milestone / display 'Milestone'
            style(column)=[cellwidth=2.30in just=left];
        define Method    / display 'Model'
            style(column)=[cellwidth=0.55in just=center];
        define P_Sig     / display noprint;
        define Beta_CI_3 / display 'Coefficient (95% CI)'
            style(column)=[cellwidth=1.55in just=center];
        define SE_3      / display 'SE'
            style(column)=[cellwidth=0.50in just=center];
        define Z_3       / display 'Z'
            style(column)=[cellwidth=0.50in just=center];
        define P_3       / display 'p'
            style(header)=[font_style=italic]
            style(column)=[cellwidth=0.55in just=center];
        define HR_CI_3   / display 'HR (95% CI)'
            style(column)=[cellwidth=1.75in just=center];

        compute Beta_CI_3;
            if P_Sig = 1 then
                call define(_col_, "style", "style=[font_weight=bold]");
        endcomp;

        compute P_3;
            if P_Sig = 1 then
                call define(_col_, "style", "style=[font_weight=bold foreground=red]");
        endcomp;

        compute HR_CI_3;
            if P_Sig = 1 then
                call define(_col_, "style", "style=[font_weight=bold]");
        endcomp;

        break after Biomarker / skip;

    run;

    footnote;
    title1;
    title2;
    title3;



%mend CombinedBiomarker;


/*===============================================================
  Run clinical milestone models for all biomarkers
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


%RunOneClinicalBiomarker(
    biomarker_short=C2N_AB40,
    biomarker_label=A^{unicode 03B2}40 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB40
);


%RunOneClinicalBiomarker(
    biomarker_short=C2N_AB42,
    biomarker_label=A^{unicode 03B2}42 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB42
);


%RunOneClinicalBiomarker(
    biomarker_short=C2N_AB4240,
    biomarker_label=A^{unicode 03B2}42/A^{unicode 03B2}40 plasma,
    datevar=C2N_AB_date,
    biomarker_var=C2N_AB4240
);


%CombinedBiomarker(
    title_suffix=for Sporadic PD
);

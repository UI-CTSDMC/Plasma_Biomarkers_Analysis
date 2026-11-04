/*==========================================================
  04_Table2_Master_and_Run_Combined.sas

  Combined Biomarker Table 2
  Reduced variable set
  sPD NSD+ and HC NSD-

  Includes:
    - Fixed cutoff grouping
    - Tertile summaries
    - Unadjusted p-values for selected variables only
    - Adjusted group p-values
    - Continuous biomarker Beta/OR (95% CI)
    - Continuous biomarker p-value
    - Fixed cutoff labels in column headers

  Update:
    - C2N biomarker statistic rows are skipped for HC NSD- tables
      when skip_c2n_stats=1.
    - UPSIT percentile, UPSIT percentile <=15, and lowput ratio
      continuous biomarker models are unadjusted.
    - Adjusted logistic OR/continuous p-values are preserved.
==========================================================*/


/***************************************************************************
* 0) LIBRARIES & STANDARD INCLUDES
***************************************************************************/
%let filepath = /*path where data and macro code are located*/;
%let report_date = 08DEC2025;

libname analysis    "&filepath.";

/* Summary table macros */
%include "&filepath.\Table Macros.sas";
/* P-value macros */
%include "&filepath.\01_Table2_Pvalues_Unadjusted_Combined.sas";
%include "&filepath.\02_Table2_Pvalues_Adjusted_Continuous.sas";
%include "&filepath.\03_Table2_Pvalues_Adjusted_Logistic.sas";


/***************************************************************************
* 0B) UNADJUSTED CONTINUOUS-BIOMARKER MODEL MACROS
*     Use these only for UPSIT percentile, UPSIT <=15, and lowput ratio.
***************************************************************************/

/*==========================================================
  Unadjusted continuous biomarker model
  For continuous outcomes where we want Beta/CI and biomarker p-value
  without adjustment for age and sex
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
               type $20 transform $12 var $64 EST_CI $60;

        set __pe_unadj;
        where strip(upcase(Parameter)) = strip(upcase("&cont_predictor"));

        var="&outcome";
        type="unadj_cont";
        transform="none";
        main=&main;
        order=&order;

        overall_cont_n = Probt;

        EST_CI = cats(
            put(Estimate, 8.4),
            " (",
            put(LowerCL, 8.4),
            ", ",
            put(UpperCL, 8.4),
            ")"
        );

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

%mend;


/*==========================================================
  Unadjusted logistic continuous biomarker model
  For binary outcomes where we want OR/CI and biomarker p-value
  without adjustment for age and sex
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
               type $20 transform $12 var $64 EST_CI $60;

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

        EST_CI = cats(
            put(OR_n, 8.4),
            " (",
            put(LowerOR_n, 8.4),
            ", ",
            put(UpperOR_n, 8.4),
            ")"
        );

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

%mend;


/***************************************************************************
* 1) LABELS
***************************************************************************/

data labels;
    infile datalines truncover;
    input @1 main best4. @6 order best4. @10 name $100.;
datalines;
1    0   ^{style[fontweight=bold]pTau217 plasma}, mean (SD)
1    1   Median (min, max)
1    2   Missing
2    0   ^{style[fontweight=bold]NFL}, mean (SD)
2    1   Median (min, max)
2    2   Missing
3    0   ^{style[fontweight=bold]GFAP}, mean (SD)
3    1   Median (min, max)
3    2   Missing
4    0   ^{style[fontweight=bold]A^{unicode 03B2}40 plasma}, mean (SD)
4    1   Median (min, max)
4    2   Missing
5    0   ^{style[fontweight=bold]A^{unicode 03B2}42 plasma}, mean (SD)
5    1   Median (min, max)
5    2   Missing
6    0   ^{style[fontweight=bold]A^{unicode 03B2}42/A^{unicode 03B2}40 plasma}, mean (SD)
6    1   Median (min, max)
6    2   Missing
6.5  0   ^{style[fontweight=bold]Sex}, n (%)
7    1   Male
7    2   Female
7    99  Missing
8    0   ^{style[fontweight=bold]Age (years)}, mean (SD)
8    1   Median (min, max)
8    2   Delete
9    0   ^{style[fontweight=bold]Duration (years)}, mean (SD)
9    1   Median (min, max)
9    2   Missing
9.5  0   ^{style[fontweight=bold]Race}, n (%)
10   1   White
10   2   Black/African American
10   3   Asian
10   4   Other
10   99  Missing
10.5 0   ^{style[fontweight=bold]Ethnicity}, n (%)
11   0   Not Hispanic/Latino
11   1   Hispanic/Latino
11   99  Missing
12   0   ^{style[fontweight=bold]UPSIT Raw Score}, mean (SD)
12   1   Median (min, max)
12   2   Delete
13   0   ^{style[fontweight=bold]Revised UPSIT percentile}, mean (SD)
13   1   Median (min, max)
13   2   Delete
14   0   Delete
14   1   ^{style[fontweight=bold]Revised UPSIT percentile ^{unicode 2264} 15}, n (%)
14   99  Missing
14.5 0   ^{style[fontweight=bold]NSD+}, n (%)
15   1   1a
15   1.5 1b
15   2   2a
15   2.5 2b
15   3   3
15   4   4
15   5   5
15   6   6
15   99  Missing
15.5 0   ^{style[fontweight=bold]MDS-UPDRS item 1.1}, n (%)
16   0   Normal
16   1   Slight
16   2   Mild
16   3   Moderate
16   4   Severe
16   99  Missing
17   0   ^{style[fontweight=bold]MDS-UPDRS part I}, mean (SD)
17   1   Median (min, max)
17   2   Missing
18   0   ^{style[fontweight=bold]MDS-UPDRS part II}, mean (SD)
18   1   Median (min, max)
18   2   Missing
19   0   ^{style[fontweight=bold]MDS-UPDRS part III}, mean (SD)
19   1   Median (min, max)
19   2   Missing
20   0   ^{style[fontweight=bold]MDS-UPDRS total score}, mean (SD)
20   1   Median (min, max)
20   2   Missing
21   0   ^{style[fontweight=bold]MoCA total score}, mean (SD)
21   1   Median (min, max)
21   2   Missing
22   0   ^{style[fontweight=bold]SCOPA-AUT total score}, mean (SD)
22   1   Median (min, max)
22   2   Missing
23   0   ^{style[fontweight=bold]% Age/sex-expected low putamen}, mean (SD)
23   1   Median (min, max)
23   2   Missing
24   0   ^{style[fontweight=bold]Mean striatum CWM}, mean (SD)
24   1   Median (min, max)
24   2   Missing
25   0   ^{style[fontweight=bold]Mean putamen CWM}, mean (SD)
25   1   Median (min, max)
25   2   Missing
25.5 0   ^{style[fontweight=bold]Cognitive categorization}, n (%)
26   1   Normal cognition
26   2   Mild cognitive impairment (MCI)
26   3   Dementia
26   99  Missing
;
run;

proc sort data=labels;
    by main order;
run;


/***************************************************************************
* 2) BUILD SUMMARY MACRO
***************************************************************************/

%macro build_summary(dataset=, out=, skip_c2n_stats=0);

    proc datasets lib=work nolist;
        delete stats: freq:;
    quit;

    %cont2(data=&dataset, var=Lilly_ptau217p,       num=1,  dec=2, dec2=2);
    %cont2(data=&dataset, var=IU_NFL_plasma_mean,   num=2,  dec=2, dec2=2);
    %cont2(data=&dataset, var=IU_GFAP_plasma_mean,  num=3,  dec=2, dec2=2);

    %if &skip_c2n_stats = 0 %then %do;
        %cont2(data=&dataset, var=C2N_AB40,         num=4,  dec=2, dec2=2);
        %cont2(data=&dataset, var=C2N_AB42,         num=5,  dec=2, dec2=2);
        %cont2(data=&dataset, var=C2N_AB4240,       num=6,  dec=2, dec2=2);
    %end;

    %cat  (data=&dataset, var=gen,                  num=7);
    %cont2(data=&dataset, var=age,                  num=8,  dec=1, dec2=1);
    %cont2(data=&dataset, var=duration_yrs,         num=9,  dec=1, dec2=1);
    %cat  (data=&dataset, var=race,                 num=10);
    %cat  (data=&dataset, var=ethn_hispanic,        num=11);

    %cont2(data=&dataset, var=UPSIT,                num=12, dec=1, dec2=1);
    %cont2(data=&dataset, var=upsit_pctl,           num=13, dec=1, dec2=1);
    %cat  (data=&dataset, var=upsit_pctl_le15,      num=14);

    %cat  (data=&dataset, var=stage_num,            num=15);
    %cat  (data=&dataset, var=mds1_cog,             num=16);

    %cont2(data=&dataset, var=updrs1_score,         num=17, dec=1, dec2=0);
    %cont2(data=&dataset, var=updrs2_score,         num=18, dec=1, dec2=0);
    %cont2(data=&dataset, var=updrs3_score,         num=19, dec=1, dec2=0);
    %cont2(data=&dataset, var=updrs_totscore,       num=20, dec=1, dec2=0);
    %cont2(data=&dataset, var=moca,                 num=21, dec=1, dec2=0);
    %cont2(data=&dataset, var=scopa,                num=22, dec=1, dec2=0);

    %cont2(data=&dataset, var=lowput_ratio_cwm,     num=23, dec=2, dec2=2);
    %cont2(data=&dataset, var=mia_striatum_ref_cwm, num=24, dec=2, dec2=2);
    %cont2(data=&dataset, var=mia_putamen_ref_cwm,  num=25, dec=2, dec2=2);

    %cat  (data=&dataset, var=cogstate,             num=26);

    data &out;
        merge stats4_: freq3_:;
        by main order;
    run;

    proc sort data=&out;
        by main order;
    run;

    proc datasets lib=work nolist;
        delete stats: freq:;
    quit;

%mend build_summary;


/***************************************************************************
* 3) MAIN TABLE MACRO
***************************************************************************/

%macro table2(
    var=,
    cutoff_type=,
    date=,
    cohort_filter=,
    cohort_label=,
    table_label=,
    outfile_label=,
    t1_label=,
    t2_label=,
    t3_label=,
    skip_stage=0,
    skip_c2n_stats=0
);

    %local biom;
    %let biom=&var;

    /***********************************************************************
    * 3A) Prepare analysis dataset
    ***********************************************************************/

    data combined_analysis_dataset;
        set analysis.combined_biomarker_with_clinical;

        if not missing(&date);
        if visit_range = "Early (BL-V02)";
        if &cohort_filter;

        if not missing(moca) then do;
            if      moca ge 26 then moca_cat = 1;
            else if moca ge 21 then moca_cat = 2;
            else                    moca_cat = 3;
        end;
    run;


    data df;
        set combined_analysis_dataset;

        /*======================================================
          Fixed biomarker cutoff groups
        ======================================================*/

%if %upcase(&cutoff_type.) = PTAU %then %do;

    if not missing(Lilly_ptau217p) then do;
        if Lilly_ptau217p < 0.131 then group = 1;
        else if Lilly_ptau217p >= 0.131
             and Lilly_ptau217p <= 0.174 then group = 2;
        else if Lilly_ptau217p > 0.174 then group = 3;
        else group = .;
    end;
    else if Lilly_ptau217p_LLOD = 1 then group = 1;
    else if Lilly_ptau217p_ULOD = 1 then group = 3;
    else group = .;

%end;

%else %if %upcase(&cutoff_type.) = NFL %then %do;

    if not missing(IU_NFL_plasma_mean) then do;
        if IU_NFL_plasma_mean < 9.877088349 then group = 1;
        else if IU_NFL_plasma_mean >= 9.877088349
             and IU_NFL_plasma_mean <= 14.53844902 then group = 2;
        else if IU_NFL_plasma_mean > 14.53844902 then group = 3;
        else group = .;
    end;
    else group = .;

%end;

%else %if %upcase(&cutoff_type.) = GFAP %then %do;

    if not missing(IU_GFAP_plasma_mean) then do;
        if IU_GFAP_plasma_mean < 64.47757262 then group = 1;
        else if IU_GFAP_plasma_mean >= 64.47757262
             and IU_GFAP_plasma_mean <= 105.3411683 then group = 2;
        else if IU_GFAP_plasma_mean > 105.3411683 then group = 3;
        else group = .;
    end;
    else group = .;

%end;

%else %if %upcase(&cutoff_type.) = AB40 %then %do;

    if not missing(C2N_AB40) then do;
        if C2N_AB40 < 374.886 then group = 1;
        else if C2N_AB40 >= 374.886
             and C2N_AB40 <= 422.785 then group = 2;
        else if C2N_AB40 > 422.785 then group = 3;
        else group = .;
    end;
    else group = .;

%end;

%else %if %upcase(&cutoff_type.) = AB42 %then %do;

    if not missing(C2N_AB42) then do;
        if C2N_AB42 < 40.219 then group = 1;
        else if C2N_AB42 >= 40.219
             and C2N_AB42 <= 45.87 then group = 2;
        else if C2N_AB42 > 45.87 then group = 3;
        else group = .;
    end;
    else group = .;

%end;

%else %if %upcase(&cutoff_type.) = AB4240 %then %do;

    if not missing(C2N_AB4240) then do;
        if C2N_AB4240 < 0.10402 then group = 1;
        else if C2N_AB4240 >= 0.10402
             and C2N_AB4240 <= 0.110585 then group = 2;
        else if C2N_AB4240 > 0.110585 then group = 3;
        else group = .;
    end;
    else group = .;

%end;

%else %do;
    group = .;
%end;

        if missing(group) then delete;

        /* Binary variables for adjusted logistic models */
        if cogstate = 1 then cogstate_bin = 0;
        else if cogstate > 1 then cogstate_bin = 1;
        else cogstate_bin = .;

        if stage_num in (1, 1.5, 2, 2.5) then stage_bin = 0;
        else if stage_num in (3, 4, 5, 6) then stage_bin = 1;
        else stage_bin = .;

        if mds1_cog > 0 then mds1_cog_bin = 1;
        else if not missing(mds1_cog) then mds1_cog_bin = 0;
        else mds1_cog_bin = .;

        if race = 1 then race_bin = 0;
        else if not missing(race) and race > 1 then race_bin = 1;
        else race_bin = .;

    run;


    /***********************************************************************
    * 3B) Summary statistics
    ***********************************************************************/

    %build_summary(
        dataset=df,
        out=combine_1,
        skip_c2n_stats=&skip_c2n_stats
    );


    /***********************************************************************
    * 3C) Unadjusted p-values
    ***********************************************************************/

    %pvalues(data=df(where=(group in (1 2 3))), label=cat_123);
    %pvalues(data=df(where=(group in (1 2))),   label=cat_12);
    %pvalues(data=df(where=(group in (2 3))),   label=cat_23);
    %pvalues(data=df(where=(group in (1 3))),   label=cat_13);

    %make_pvalues(groups=123 12 23 13);

    proc sort data=pvalues_noadjust_char;
        by main order;
    run;


    /***********************************************************************
    * 3D) Adjusted continuous group p-values
    ***********************************************************************/

    %auto_transform_pvals(data=df, outcome=Lilly_ptau217p,       class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=1,  order=0, out=p_adj_ptau217);
    %auto_transform_pvals(data=df, outcome=IU_NFL_plasma_mean,   class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=2,  order=0, out=p_adj_nfl);
    %auto_transform_pvals(data=df, outcome=IU_GFAP_plasma_mean,  class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=3,  order=0, out=p_adj_gfap);

    %if &skip_c2n_stats = 0 %then %do;
        %auto_transform_pvals(data=df, outcome=C2N_AB40,         class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=4,  order=0, out=p_adj_ab40);
        %auto_transform_pvals(data=df, outcome=C2N_AB42,         class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=5,  order=0, out=p_adj_ab42);
        %auto_transform_pvals(data=df, outcome=C2N_AB4240,       class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=6,  order=0, out=p_adj_ab4240);
    %end;

    %auto_transform_pvals(data=df, outcome=UPSIT,                class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=12, order=0, out=p_adj_upsit);
    %auto_transform_pvals(data=df, outcome=updrs1_score,         class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=17, order=0, out=p_adj_updrs1);
    %auto_transform_pvals(data=df, outcome=updrs2_score,         class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=18, order=0, out=p_adj_updrs2);
    %auto_transform_pvals(data=df, outcome=updrs3_score,         class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=19, order=0, out=p_adj_updrs3);
    %auto_transform_pvals(data=df, outcome=updrs_totscore,       class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=20, order=0, out=p_adj_updrst);
    %auto_transform_pvals(data=df, outcome=moca,                 class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=21, order=0, out=p_adj_moca);
    %auto_transform_pvals(data=df, outcome=scopa,                class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=22, order=0, out=p_adj_scopa);
    %auto_transform_pvals(data=df, outcome=mia_striatum_ref_cwm, class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=24, order=0, out=p_adj_striat);
    %auto_transform_pvals(data=df, outcome=mia_putamen_ref_cwm,  class=group gen, predictors=group age gen, groupvar=group, alpha=0.05, main=25, order=0, out=p_adj_putam);

    data p_adj_all(drop=P_123_n P_12_n P_13_n P_23_n EST_CI overall_cont overall_cont_n);
        set
            p_adj_ptau217
            p_adj_nfl
            p_adj_gfap

            %if &skip_c2n_stats = 0 %then %do;
                p_adj_ab40
                p_adj_ab42
                p_adj_ab4240
            %end;

            p_adj_upsit
            p_adj_updrs1
            p_adj_updrs2
            p_adj_updrs3
            p_adj_updrst
            p_adj_moca
            p_adj_scopa
            p_adj_striat
            p_adj_putam
        ;
    run;

    proc sort data=p_adj_all;
        by main order;
    run;


    /***********************************************************************
    * 3E) Continuous biomarker predictor models
    *     Adds EST_CI and overall_cont columns
    ***********************************************************************/

    %auto_cont_biomarker_assoc(data=df, outcome=Lilly_ptau217p,       cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=1,  order=0, out=biom_assoc_ptau217);
    %auto_cont_biomarker_assoc(data=df, outcome=IU_NFL_plasma_mean,   cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=2,  order=0, out=biom_assoc_nfl);
    %auto_cont_biomarker_assoc(data=df, outcome=IU_GFAP_plasma_mean,  cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=3,  order=0, out=biom_assoc_gfap);

    %if &skip_c2n_stats = 0 %then %do;
        %auto_cont_biomarker_assoc(data=df, outcome=C2N_AB40,         cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=4,  order=0, out=biom_assoc_ab40);
        %auto_cont_biomarker_assoc(data=df, outcome=C2N_AB42,         cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=5,  order=0, out=biom_assoc_ab42);
        %auto_cont_biomarker_assoc(data=df, outcome=C2N_AB4240,       cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=6,  order=0, out=biom_assoc_ab4240);
    %end;

    %auto_cont_biomarker_assoc(data=df, outcome=UPSIT,                cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=12, order=0, out=biom_assoc_upsit);
    %auto_cont_biomarker_assoc_unadj(data=df, outcome=upsit_pctl,     cont_predictor=&biom, alpha=0.05, main=13, order=0, out=biom_assoc_upsitpctl);
    %auto_cont_biomarker_assoc(data=df, outcome=updrs1_score,         cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=17, order=0, out=biom_assoc_updrs1);
    %auto_cont_biomarker_assoc(data=df, outcome=updrs2_score,         cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=18, order=0, out=biom_assoc_updrs2);
    %auto_cont_biomarker_assoc(data=df, outcome=updrs3_score,         cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=19, order=0, out=biom_assoc_updrs3);
    %auto_cont_biomarker_assoc(data=df, outcome=updrs_totscore,       cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=20, order=0, out=biom_assoc_updrst);
    %auto_cont_biomarker_assoc(data=df, outcome=moca,                 cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=21, order=0, out=biom_assoc_moca);
    %auto_cont_biomarker_assoc(data=df, outcome=scopa,                cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=22, order=0, out=biom_assoc_scopa);
    %auto_cont_biomarker_assoc_unadj(data=df, outcome=lowput_ratio_cwm, cont_predictor=&biom, alpha=0.05, main=23, order=0, out=biom_assoc_lowput);
    %auto_cont_biomarker_assoc(data=df, outcome=mia_striatum_ref_cwm, cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=24, order=0, out=biom_assoc_striat);
    %auto_cont_biomarker_assoc(data=df, outcome=mia_putamen_ref_cwm,  cont_predictor=&biom, class=gen, predictors=&biom age gen, alpha=0.05, main=25, order=0, out=biom_assoc_putam);

    data biom_all;
        set biom_assoc_:;
    run;

    proc sort data=biom_all;
        by main order;
    run;

    proc datasets lib=work nolist;
        delete biom_assoc_:;
    quit;


    /***********************************************************************
    * 3F) Adjusted categorical p-values
    *     Also keeps OR (95% CI) + continuous biomarker p-value
    ***********************************************************************/

    %if &skip_stage = 0 %then %do;

        %auto_cat_pvals(
            data=df,
            outcome=stage_bin,
            groupvar=group,
            agevar=age,
            genvar=gen,
            gen_ref=1,
            event_ref=1,
            cont_predictor=&biom,
            main=14.5,
            order=0,
            out=pc_adj_stage
        );

    %end;

    %auto_cat_pvals_unadj_cont(data=df, outcome=upsit_pctl_le15, event_ref=0, cont_predictor=&biom, main=14, order=1, out=pc_adj_upsit15);

    %auto_cat_pvals(
        data=df,
        outcome=mds1_cog_bin,
        groupvar=group,
        agevar=age,
        genvar=gen,
        gen_ref=1,
        event_ref=0,
        cont_predictor=&biom,
        main=15.5,
        order=0,
        out=pc_adj_mds1
    );

    %auto_cat_pvals(
        data=df,
        outcome=cogstate_bin,
        groupvar=group,
        agevar=age,
        genvar=gen,
        gen_ref=1,
        event_ref=1,
        cont_predictor=&biom,
        main=25.5,
        order=0,
        out=pc_adj_cog
    );


/*----------------------------------------------------------
  Keep adjusted categorical GROUP p-values only.
  Drop EST_CI/overall_cont here so those are handled separately.
----------------------------------------------------------*/

data pc_adj_group_pvals(drop=P_123_n P_12_n P_13_n P_23_n overall_cont_n EST_CI overall_cont);
    set
        %if &skip_stage = 0 %then %do;
            pc_adj_stage
        %end;
        pc_adj_mds1
        pc_adj_cog
    ;
run;

proc sort data=pc_adj_group_pvals;
    by main order;
run;


/*----------------------------------------------------------
  Adjusted logistic OR + continuous biomarker p-value.
  Rename columns so they CANNOT be overwritten by other datasets.
----------------------------------------------------------*/

data pc_adj_cont;
    set
        %if &skip_stage = 0 %then %do;
            pc_adj_stage
        %end;
        pc_adj_mds1
        pc_adj_cog
    ;

    length EST_CI_cat overall_cont_cat $60;

    EST_CI_cat = EST_CI;
    overall_cont_cat = overall_cont;

    keep main order EST_CI_cat overall_cont_cat;
run;

proc sort data=pc_adj_cont;
    by main order;
run;


/*----------------------------------------------------------
  Unadjusted logistic OR + continuous biomarker p-value
  for UPSIT percentile <=15 only.
----------------------------------------------------------*/

data pc_upsit15_cont;
    set pc_adj_upsit15;

    length EST_CI_upsit15 overall_cont_upsit15 $60;

    EST_CI_upsit15 = EST_CI;
    overall_cont_upsit15 = overall_cont;

    keep main order EST_CI_upsit15 overall_cont_upsit15;
run;

proc sort data=pc_upsit15_cont;
    by main order;
run;


/*----------------------------------------------------------
  Continuous outcome biomarker association outputs.
  Rename columns so they CANNOT overwrite categorical ORs.
----------------------------------------------------------*/

data biom_all_clean;
    set biom_all;

    length EST_CI_cont overall_cont_cont $60;

    EST_CI_cont = EST_CI;
    overall_cont_cont = overall_cont;

    keep main order EST_CI_cont overall_cont_cont;
run;

proc sort data=biom_all_clean;
    by main order;
run;


/***********************************************************************
* 3G) Merge p-values and continuous biomarker columns
***********************************************************************/

data pvals_all;
    merge
        pvalues_noadjust_char
        p_adj_all
        pc_adj_group_pvals
        biom_all_clean
        pc_adj_cont
        pc_upsit15_cont
    ;
    by main order;

    length EST_CI overall_cont $60;

    /*
      Priority:
        1. Adjusted categorical logistic rows:
           NSD stage, MDS-UPDRS item 1.1, COGSTATE
        2. Unadjusted UPSIT <=15 logistic row
        3. Continuous outcome rows
    */

    EST_CI = coalescec(EST_CI_cat, EST_CI_upsit15, EST_CI_cont);
    overall_cont = coalescec(overall_cont_cat, overall_cont_upsit15, overall_cont_cont);

    drop
        EST_CI_cat
        EST_CI_upsit15
        EST_CI_cont
        overall_cont_cat
        overall_cont_upsit15
        overall_cont_cont
    ;
run;

proc sort data=pvals_all;
    by main order;
run;


    /***********************************************************************
    * 3H) Final combine
    ***********************************************************************/

    data combine(drop=i all_na0 overall_num transform type var);
        merge labels(in=a)
              combine_1
              pvals_all;
        by main order;
        if a;

        if name = "Delete" then delete;

        /* Remove pTau217, NFL, GFAP, Aß40, Aß42, and Aß42/Aß40 rows  - based on Tom's modification (Sep 19th)*/
        if main in (1, 2, 3, 4, 5, 6) then delete;
        /* Remove C2N biomarker rows for HC NSD- tables */
        if &skip_c2n_stats = 1 and main in (4, 5, 6) then delete;

        /* Remove NSD stage rows for HC NSD- */
        if &skip_stage = 1 and main in (14.5, 15) then delete;

        array col(*) $ _:;

        do i = 1 to dim(col);

            if col[i] = "0 (0%)" then col[i] = "0";
            else if col[i] = "0 (.%)" then col[i] = "0";

            else if strip(col[i]) = ""
                 and mod(int(main*2), 2) = 0
            then col[i] = "NA";

        end;

        all_na0 = 1;
        do i = 1 to dim(col);
            if upcase(compress(strip(col[i]), ". ")) not in ("NA", "0") then do;
                all_na0 = 0;
                leave;
            end;
        end;

        if all_na0 then delete;

        /* Convert the overall p-value from character to numeric */
        if substr(strip(P_123), 1, 1) = "<" then
            overall_num = 0.0001;
        else
            overall_num = input(
                compress(scan(P_123, 1, "^"), "* "),
                ?? best32.
            );
        /* Hide pairwise p-values when the overall p-value is not significant */
        if not missing(overall_num) and overall_num >= 0.05 then do;
            P_12 = "";
            P_13 = "";
            P_23 = "";
        end;

    run;

    proc sort data=combine;
        by main order;
    run;


    /***********************************************************************
    * 3I) Group Ns
    ***********************************************************************/

    proc sql noprint;
        select count(distinct patno) into :g1 trimmed
        from df
        where group = 1;

        select count(distinct patno) into :g2 trimmed
        from df
        where group = 2;

        select count(distinct patno) into :g3 trimmed
        from df
        where group = 3;
    quit;


    /***********************************************************************
    * 3J) Output report
    ***********************************************************************/

    options orientation=landscape nodate nonumber missing=" ";
    ods escapechar="^";
    ods graphics / border=off;
    ods noptitle;

    title bold h=10pt j=left f="Arial"
        "TABLE 2: &table_label. &cohort_label. Participants by Tertiles";

    footnote1 h=9pt j=left f="Arial"
        "Report Generated on Data Submitted as of: &report_date..";

    footnote2 h=9pt j=left f="Arial"
        "*  P-values calculated using PROC LOGISTIC adjusted for age and sex.^n"
        "**  P-values calculated using GLM adjusted for age and sex.^n"
        "*** P-values calculated using log-transformed model.^n"
        "**** P-values calculated using square-root-transformed model.^n"
        "***** P-values calculated using squared-value (power-2) transformation.^n"
        "****** P-values calculated using rank-based (nonparametric) model.";

    footnote3 h=9pt j=left f="Arial"
        "No asterisks indicate unadjusted p-values. Kruskal-Wallis and Chi-Square/Fisher's Exact Test were used for overall p-values. Wilcoxon Rank Sum tests and Chi-Square/Fisher's Exact Test were used for pairwise comparisons.";

    footnote4 h=9pt j=left f="Arial"
        "In p-value calculations, the race variable collapsed all non-white groups together and the resulting p-value is a comparison of non-white vs. white between groups.";

    footnote5 h=9pt j=left f="Arial"
        "For MDS-UPDRS item 1.1 and COGSTATE, all groups except 'Normal' were collapsed and the resulting p-value compares normal vs. abnormal categories.";

    footnote6 h=9pt j=left f="Arial"
        "For NSD Stage, stages 1-2.5 are collapsed and compared against stages 3-6, where applicable.";

    footnote7 h=9pt j=left f="Arial"
        "Beta / OR (95% CI) and continuous biomarker p-values are from models using the table biomarker as a continuous predictor. Unless otherwise specified, models are adjusted for age and sex. UPSIT percentile, UPSIT percentile <=15, and % age/sex-expected low putamen models are unadjusted.";

proc report data=combine spanrows split="\"
    style(report)={just=l width=100% cellpadding=1 cellspacing=0}
    style(column)={font_size=7pt fontfamily=Arial}
    style(header)={font_size=7pt fontfamily=Arial};

    column main order name
           ("&table_label. tertiles" _1 _2 _3)
           P_123 P_12 P_13 P_23 EST_CI overall_cont;

    define main  / order order=internal noprint;
    define order / order order=internal noprint;

    define name / display "Variable"
        style(column)={width=1.45in}
        style(header)={width=1.45in};

    define _1 / display center
        "&t1_label.^n^{style[fontsize=7pt fontweight=MEDIUM](N = &g1)}"
        style(column)={width=1.05in borderleftcolor=DarkSlateGrey}
        style(header)={width=1.05in};

    define _2 / display center
        "&t2_label.^n^{style[fontsize=7pt fontweight=MEDIUM](N = &g2)}"
        style(column)={width=1.20in}
        style(header)={width=1.20in};

    define _3 / display center
        "&t3_label.^n^{style[fontsize=7pt fontweight=MEDIUM](N = &g3)}"
        style(column)={width=1.05in borderrightcolor=DarkSlateGrey}
        style(header)={width=1.05in};

    define P_123 / display center "Overall\P-value"
        style(column)={width=.65in}
        style(header)={width=.65in};

    define P_12 / display center "T1 vs. T2\P-value"
        style(column)={width=.65in}
        style(header)={width=.65in};

    define P_13 / display center "T1 vs. T3\P-value"
        style(column)={width=.65in}
        style(header)={width=.65in};

    define P_23 / display center "T2 vs. T3\P-value"
        style(column)={width=.65in}
        style(header)={width=.65in};

    define EST_CI / display center
        "Beta / OR\(95% CI)"
        style(column)={width=1.05in}
        style(header)={width=1.05in};

    define overall_cont / display center
        "Continuous\biomarker\P-value"
        style(column)={width=.75in}
        style(header)={width=.75in};

        compute name;
            if index(name, "^{style[fontweight=bold]") = 0 then
                call define(_col_, "style", "style=[paddingleft=3em]");
            else
                call define(_row_, "style", "style=[bordertopcolor=DarkSlateGrey bordertopwidth=1 bordertopstyle=solid]");
        endcomp;

        compute P_123;
            length _numtxt $20;
            _numtxt = compress(scan(P_123, 1, "^"), "0123456789.+-Ee", "k");
            _pnum = input(_numtxt, ?? best32.);
            if not missing(_pnum) and _pnum < 0.05 then
                call define(_col_, "style", "style=[foreground=red fontweight=bold]");
        endcomp;

        compute P_12;
            length _numtxt $20;
            _numtxt = compress(scan(P_12, 1, "^"), "0123456789.+-Ee", "k");
            _pnum = input(_numtxt, ?? best32.);
            if not missing(_pnum) and _pnum < 0.05 then
                call define(_col_, "style", "style=[foreground=red fontweight=bold]");
        endcomp;

        compute P_13;
            length _numtxt $20;
            _numtxt = compress(scan(P_13, 1, "^"), "0123456789.+-Ee", "k");
            _pnum = input(_numtxt, ?? best32.);
            if not missing(_pnum) and _pnum < 0.05 then
                call define(_col_, "style", "style=[foreground=red fontweight=bold]");
        endcomp;

        compute P_23;
            length _numtxt $20;
            _numtxt = compress(scan(P_23, 1, "^"), "0123456789.+-Ee", "k");
            _pnum = input(_numtxt, ?? best32.);
            if not missing(_pnum) and _pnum < 0.05 then
                call define(_col_, "style", "style=[foreground=red fontweight=bold]");
        endcomp;

        compute overall_cont;
            length _numtxt $20;
            _numtxt = compress(scan(overall_cont, 1, "^"), "0123456789.+-Ee", "k");
            _pnum = input(_numtxt, ?? best32.);
            if not missing(_pnum) and _pnum < 0.05 then
                call define(_col_, "style", "style=[foreground=red fontweight=bold]");
        endcomp;

    run;

    title;
    footnote;


    /***********************************************************************
    * 3K) Cleanup
    ***********************************************************************/

proc datasets lib=work nolist;
    delete
        combined_analysis_dataset
        df
        combine_1
        combine
        pvalues_noadjust
        pvalues_noadjust_char
        all_tests_:
        p_adj_:
        pc_adj_:
        pc_adj_group_pvals
        pc_adj_cont
        pc_upsit15_cont
        pvals_all
        biom_all
        biom_all_clean
    ;
quit;


%mend table2;


/***************************************************************************
* 4) BATCH CALLS
***************************************************************************/
/*==========================================================
  sPD NSD+
  C2N statistics included
==========================================================*/

%table2(
    var=Lilly_ptau217p,
    cutoff_type=PTAU,
    date=Lilly_ptau217p_date,
    cohort_filter=subgroup_definition="Sporadic PD" and NSD=1,
    cohort_label=sPD NSD+,
    table_label=pTau217 plasma,
    outfile_label=Lilly_pTau217,
    t1_label=%str(<0.131),
    t2_label=%str(0.131-0.174),
    t3_label=%str(>0.174),
    skip_stage=0,
    skip_c2n_stats=0
);

%table2(
    var=IU_NFL_plasma_mean,
    cutoff_type=NFL,
    date=IU_Plasma_date,
    cohort_filter=subgroup_definition="Sporadic PD" and NSD=1,
    cohort_label=sPD NSD+,
    table_label=NFL,
    outfile_label=IU_NfL,
    t1_label=%str(<9.877088349),
    t2_label=%str(9.877088349-14.53844902),
    t3_label=%str(>14.53844902),    
    skip_stage=0,
    skip_c2n_stats=0
);

%table2(
    var=IU_GFAP_plasma_mean,
    cutoff_type=GFAP,
    date=IU_Plasma_date,
    cohort_filter=subgroup_definition="Sporadic PD" and NSD=1,
    cohort_label=sPD NSD+,
    table_label=GFAP,
    outfile_label=IU_GFAP,
    t1_label=%str(<64.47757262),
    t2_label=%str(64.47757262-105.3411683),
    t3_label=%str(>105.3411683),    
    skip_stage=0,
    skip_c2n_stats=0
);

%table2(
    var=C2N_AB40,
    cutoff_type=AB40,
    date=C2N_AB_date,
    cohort_filter=subgroup_definition="Sporadic PD" and NSD=1,
    cohort_label=sPD NSD+,
    table_label=A^{unicode 03B2}40 plasma,
    outfile_label=C2N_AB40,
    t1_label=%str(<374.886),
    t2_label=%str(374.886-422.785),
    t3_label=%str(>422.785),
    skip_stage=0,
    skip_c2n_stats=0
);

%table2(
    var=C2N_AB42,
    cutoff_type=AB42,
    date=C2N_AB_date,
    cohort_filter=subgroup_definition="Sporadic PD" and NSD=1,
    cohort_label=sPD NSD+,
    table_label=A^{unicode 03B2}42 plasma,
    outfile_label=C2N_AB42,
    t1_label=%str(<40.219),
    t2_label=%str(40.219-45.87),
    t3_label=%str(>45.87),
    skip_stage=0,
    skip_c2n_stats=0
);

%table2(
    var=C2N_AB4240,
    cutoff_type=AB4240,
    date=C2N_AB_date,
    cohort_filter=subgroup_definition="Sporadic PD" and NSD=1,
    cohort_label=sPD NSD+,
    table_label=A^{unicode 03B2}42/A^{unicode 03B2}40 plasma,
    outfile_label=C2N_AB4240,
    t1_label=%str(<0.10402),
    t2_label=%str(0.10402-0.110585),
    t3_label=%str(>0.110585),
    skip_stage=0,
    skip_c2n_stats=0
);


/*==========================================================
  HC NSD-
  C2N statistic rows skipped inside these tables
==========================================================*/

%table2(
    var=Lilly_ptau217p,
    cutoff_type=PTAU,
    date=Lilly_ptau217p_date,
    cohort_filter=subgroup_definition="Healthy Control" and NSD=0,
    cohort_label=HC NSD-,
    table_label=pTau217 plasma,
    outfile_label=Lilly_pTau217,
    t1_label=%str(<0.131),
    t2_label=%str(0.131-0.174),
    t3_label=%str(>0.174),
    skip_stage=1,
    skip_c2n_stats=1
);

%table2(
    var=IU_NFL_plasma_mean,
    cutoff_type=NFL,
    date=IU_Plasma_date,
    cohort_filter=subgroup_definition="Healthy Control" and NSD=0,
    cohort_label=HC NSD-,
    table_label=NFL,
    outfile_label=IU_NfL,
    t1_label=%str(<9.877088349),
    t2_label=%str(9.877088349-14.53844902),
    t3_label=%str(>14.53844902),
    skip_stage=1,
    skip_c2n_stats=1
);

%table2(
    var=IU_GFAP_plasma_mean,
    cutoff_type=GFAP,
    date=IU_Plasma_date,
    cohort_filter=subgroup_definition="Healthy Control" and NSD=0,
    cohort_label=HC NSD-,
    table_label=GFAP,
    outfile_label=IU_GFAP,
    t1_label=%str(<64.47757262),
    t2_label=%str(64.47757262-105.3411683),
    t3_label=%str(>105.3411683),
    skip_stage=1,
    skip_c2n_stats=1
);

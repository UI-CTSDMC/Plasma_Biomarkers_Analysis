/*--------------------------------------------------------------------------
  Supplemental Table 1
  Plasma biomarker distributions by diagnostic group
  Descriptive statistics only; no p-values
--------------------------------------------------------------------------*/
%let filepath = /*path where data and macro code are located*/;

%let report_date = 08DEC2025;

libname analysis
    "&filepath.";

%include
    "&filepath\Table Macros.sas";


/*--------------------------------------------------------------------------
  1. Select the early-window analysis records
--------------------------------------------------------------------------*/

data combined_analysis_dataset;
    set analysis.combined_biomarker_with_clinical;

    if not missing(Lilly_ptau217p_date)
       or not missing(IU_Plasma_date)
       or not missing(C2N_AB_date);

    if visit_range = "Early (BL-V02)";
run;


/*--------------------------------------------------------------------------
  2. Create the five table columns

     1 = Overall
     2 = Sporadic PD
     3 = Healthy Control
     4 = LRRK2 PD
     5 = GBA PD

  Each participant is output twice: once to Overall and once to the relevant
  diagnostic group. NSD status is not used to define these columns.
--------------------------------------------------------------------------*/

data df;
    set combined_analysis_dataset;

    cohort_group = .;

    if cohort = "1" and subgroup_definition = "Sporadic PD" then
        cohort_group = 2;
    else if cohort = "2" and subgroup_definition = "Healthy Control" then
        cohort_group = 3;
    else if cohort = "1" and subgroup_definition = "LRRK2" then
        cohort_group = 4;
    else if cohort = "1" and subgroup_definition = "GBA" then
        cohort_group = 5;

    if missing(cohort_group) then delete;

    group = 1;
    output;

    group = cohort_group;
    output;
run;


/* Confirm that each source record contributes to Overall and one subgroup. */
proc freq data=df;
    tables group / missing;
run;

proc sql;
    select group,
           count(*) as Number_of_rows,
           count(distinct PATNO) as Number_of_participants
    from df
    group by group
    order by group;
quit;


/*--------------------------------------------------------------------------
  3. Build biomarker summaries
--------------------------------------------------------------------------*/

%macro build_summary(dataset=, out=);

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

    data &out.;
        merge stats4_: freq3_:;
        by main order;
    run;

    proc sort data=&out.;
        by main order;
    run;

    proc datasets library=work nolist;
        delete stats: freq:;
    quit;

%mend build_summary;

%build_summary(dataset=df, out=combine1_df);

proc sort data=combine1_df;
    by main order;
run;


/*--------------------------------------------------------------------------
  4. Row labels
--------------------------------------------------------------------------*/

data labels;
    infile datalines truncover;
    input
        @1  main  best4.
        @6  order best4.
        @10 name  $100.;

datalines;
1    0   ^{style[fontweight=bold]pTau217 plasma (Lilly)}, mean (SD)
1    1   Median (min, max)
1    2   Missing
1.5  0   ^{style[fontweight=bold]pTau217 tertiles (Lilly)}, n (%)
2    1   < 0.131
2    2   0.131 - 0.174
2    3   > 0.174
2    99  Missing
3    0   ^{style[fontweight=bold]NfL}, mean (SD)
3    1   Median (min, max)
3    2   Missing
3.5  0   ^{style[fontweight=bold]NfL tertiles}, n (%)
4    1   < 9.8771
4    2   9.8771 - 14.5384
4    3   > 14.5384
4    99  Missing
5    0   ^{style[fontweight=bold]GFAP}, mean (SD)
5    1   Median (min, max)
5    2   Missing
5.5  0   ^{style[fontweight=bold]GFAP tertiles}, n (%)
6    1   < 64.4776
6    2   64.4776 - 105.3412
6    3   > 105.3412
6    99  Missing
7    0   ^{style[fontweight=bold]A^{unicode 03B2}40 plasma}, mean (SD)
7    1   Median (min, max)
7    2   Missing
7.5  0   ^{style[fontweight=bold]A^{unicode 03B2}40 plasma tertiles}, n (%)
8    1   < 374.886
8    2   374.886 - 422.785
8    3   > 422.785
9    0   ^{style[fontweight=bold]A^{unicode 03B2}42 plasma}, mean (SD)
9    1   Median (min, max)
9    2   Missing
9.5  0   ^{style[fontweight=bold]A^{unicode 03B2}42 plasma tertiles}, n (%)
10   1   < 40.219
10   2   40.219 - 45.870
10   3   > 45.870
11   0   ^{style[fontweight=bold]A^{unicode 03B2}42/A^{unicode 03B2}40 plasma}, mean (SD)
11   1   Median (min, max)
11   2   Missing
11.5 0   ^{style[fontweight=bold]A^{unicode 03B2}42/A^{unicode 03B2}40 plasma tertiles}, n (%)
12   1   < 0.1040
12   2   0.1040 - 0.1106
12   3   > 0.1106
;
run;

proc sort data=labels;
    by main order;
run;


/*--------------------------------------------------------------------------
  5. Assemble the display dataset
--------------------------------------------------------------------------*/

data combine(drop=i all_na0);
    merge labels(in=in_labels)
          combine1_df;
    by main order;

    if in_labels;

    array summary_columns(*) $ _1-_5;

    do i = 1 to dim(summary_columns);
        if summary_columns[i] in ("0 (0%)", "0 (.%)") then
            summary_columns[i] = "0";
        else if strip(summary_columns[i]) = ""
            and mod(int(main * 2), 2) = 0 then
            summary_columns[i] = "NA";
    end;

    /* Remove a detail row only when all five columns are NA or zero. */
    if mod(int(main * 2), 2) = 0 then do;
        all_na0 = 1;

        do i = 1 to dim(summary_columns);
            if upcase(compress(strip(summary_columns[i]), ". ")) not in
                ("NA", "0") then do;
                all_na0 = 0;
                leave;
            end;
        end;

        if all_na0 then delete;
    end;
run;

proc sort data=combine;
    by main order;
run;

proc print data=combine;
run;


/*--------------------------------------------------------------------------
  6. Column denominators
--------------------------------------------------------------------------*/

proc sql noprint;
    select count(*) into :g1 trimmed from df where group = 1;
    select count(*) into :g2 trimmed from df where group = 2;
    select count(*) into :g3 trimmed from df where group = 3;
    select count(*) into :g4 trimmed from df where group = 4;
    select count(*) into :g5 trimmed from df where group = 5;
quit;


/*--------------------------------------------------------------------------
  7. Produce the RTF table
--------------------------------------------------------------------------*/

options
    orientation=landscape
    papersize=letter
    leftmargin=0.50in
    rightmargin=0.50in
    topmargin=0.35in
    bottommargin=0.35in
    nodate
    nonumber
    missing=" ";

ods escapechar="^";
ods noptitle;
ods graphics / border=off;

title
    bold
    h=9pt
    j=left
    f="Arial"
    "Supplemental Table 1. Plasma Biomarker Distributions Across Sporadic PD, Healthy Control, LRRK2 PD, and GBA PD Groups";

footnote1
    h=7pt
    j=left
    f="Arial"
    "Values are presented as mean (SD), median (minimum, maximum), or n (%), as indicated.";

footnote2
    h=7pt
    j=left
    f="Arial"
    "PD=Parkinson's disease; NfL=neurofilament light chain; GFAP=glial fibrillary acidic protein.";

footnote3
    h=7pt
    j=left
    f="Arial"
    "Report generated using data submitted as of &report_date..";

proc report
    data=combine
    nowd
    spanrows
    split="\"

    style(report)=[
        outputwidth=9.0in
        just=center
        cellspacing=0
        cellpadding=1
    ]

    style(column)=[
        font_size=7pt
        fontfamily=Arial
        vjust=middle
    ]

    style(header)=[
        font_size=7pt
        fontfamily=Arial
        font_weight=bold
        just=center
        vjust=middle
    ];

    column main order name _1 _2 _3 _4 _5;

    define main /
        order order=internal noprint;

    define order /
        order order=internal noprint;

    define name /
        display
        "Variable"
        style(column)=[width=2.25in]
        style(header)=[width=2.25in just=left];

    define _1 /
        display center
        "Overall\^{style[fontsize=6.5pt fontweight=medium](N=&g1.)}"
        style(column)=[width=1.30in borderleftcolor=DarkSlateGrey]
        style(header)=[width=1.30in borderleftcolor=DarkSlateGrey];

    define _2 /
        display center
        "Sporadic PD\^{style[fontsize=6.5pt fontweight=medium](N=&g2.)}"
        style(column)=[width=1.30in]
        style(header)=[width=1.30in];

    define _3 /
        display center
        "Healthy Control\^{style[fontsize=6.5pt fontweight=medium](N=&g3.)}"
        style(column)=[width=1.30in]
        style(header)=[width=1.30in];

    define _4 /
        display center
        "LRRK2 PD\^{style[fontsize=6.5pt fontweight=medium](N=&g4.)}"
        style(column)=[width=1.30in]
        style(header)=[width=1.30in];

    define _5 /
        display center
        "GBA PD\^{style[fontsize=6.5pt fontweight=medium](N=&g5.)}"
        style(column)=[width=1.30in borderrightcolor=DarkSlateGrey]
        style(header)=[width=1.30in borderrightcolor=DarkSlateGrey];

    break before main /
        style=[
            bordertopcolor=DarkSlateGrey
            bordertopwidth=1
            bordertopstyle=solid
        ];

    compute name;
        if index(name, "^{style") = 0 then
            call define(_col_, "style", "style=[paddingleft=1em]");
        else
            call define(
                _row_,
                "style",
                "style=[bordertopcolor=DarkSlateGrey]"
            );
    endcomp;
run;

title;
footnote;

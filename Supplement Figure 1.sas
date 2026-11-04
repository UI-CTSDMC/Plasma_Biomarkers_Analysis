%let filepath = /*path where data are located*/;
libname analysis "&filepath.";

/*============================================================
  1) Biomarker cohort PATNOs
============================================================*/
data df;
    set analysis.combined_biomarker_with_clinical;
    if not missing(Lilly_ptau217p_date)
    or not missing(IU_Plasma_date)
    or not missing(C2N_AB_date);
run;

proc sort data=df(keep=patno) out=biomarker_patnos nodupkey;
    by patno;
run;

/*============================================================
  2) Pull all visits for those PATNOs
============================================================*/
proc sort data=analysis.combined_biomarker_with_clinical out=all_visits;
    by patno event_id;
run;

data cohort_all_visits;
    merge biomarker_patnos(in=a) all_visits(in=b);
    by patno;
    if a and b;
run;

/*============================================================
  3) Create month2 and year
============================================================*/
data cohort_all_visits2;
    set cohort_all_visits;

    month2 = month;

    if missing(month2) then do;
        select (strip(upcase(event_id)));
            when ('BL')  month2=0;
            when ('V01') month2=3;
            when ('V02') month2=6;
            when ('V03') month2=9;
            when ('V04') month2=12;
            when ('V05') month2=18;
            when ('V06') month2=24;
            when ('V07') month2=30;
            when ('V08') month2=36;
            when ('V09') month2=42;
            when ('V10') month2=48;
            when ('V11') month2=54;
            when ('V12') month2=60;
            when ('V13') month2=72;
            when ('V14') month2=84;
            when ('V15') month2=96;
            when ('V16') month2=108;
            when ('V17') month2=120;
            when ('V18') month2=132;
            when ('V19') month2=144;
            when ('V20') month2=156;
            when ('V21') month2=168;
            when ('V22') month2=180;
            when ('V23') month2=192;
            when ('V24') month2=204;
            when ('V25') month2=216;
            otherwise;
        end;
    end;

    year = month2/12;
run;

data cohort_all_visits2; set cohort_all_visits2; where year <= 11; run;

/*============================================================
  6) Number followed at each yearly interval
============================================================*/
proc sort data=cohort_all_visits2(keep=patno cohort subgroup_definition year nsd)
          out=yearly_intervals nodupkey;
    by patno year;
run;

/* sPD */
proc sql;
    create table year_counts_sPD as
    select year, count(distinct patno) as participants
    from yearly_intervals
    where cohort='1' and subgroup_definition='Sporadic PD' 
    group by year
    order by year;
quit;

data year_counts_sPD;
    set year_counts_sPD;
    group = "sPD";
run;

/* HC */
proc sql;
    create table year_counts_HC as
    select year, count(distinct patno) as participants
    from yearly_intervals
    where cohort='2' and subgroup_definition='Healthy Control'
    group by year
    order by year;
quit;

data year_counts_HC;
    set year_counts_HC;
    group = "HC";
run;

/* Combine */
data year_counts_plot;
    set  year_counts_sPD year_counts_HC;
run;

/* Bar plot */
ods graphics / reset width=10in height=5in;

title "Number of Participants Followed at Each Yearly Interval";
proc sgplot data=year_counts_plot;
    vbar year / response=participants group=group groupdisplay=cluster;
    xaxis label="Years Since Baseline" integer;
    yaxis label="Number of Participants";
    keylegend / title="Group";
run;
title;

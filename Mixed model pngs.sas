*******************************;
* Author: Noah Johnson        *;
* File: Mixed model pngs      *;
*						      *;
* Description:				  *;
* Generates mixed model plots *;
* in Figures 1, 2, 3 in       *;
* manuscript                  *;
*******************************;


%let data_filepath= /*path to data*/;
%let output_path= /*path for desired output*/;

*obtain necessary datasets;
libname analysis "&data_filepath.";


/*Macro generates Mixed Model plot images for Figures 1,2,3 in manuscript*/

*inputs are;

*tertile - tertile biomarker;
*biom    - numeric biomarker;
*date    - date biomarker was obtained;
*Var     - longitudinal outcome of interest;
*Label   - Title for longitudinal outcome in output file;
*data    - longitudinal data specified;
*Label   - Label of variable for image;
*label2  - Label of variable in on y-axis;
*label3  - biomarker image name output; 
*age_gen - determines whether to adjust by age and sex;
*miny    - minimum on y-axis;
*maxy    - maximum on x-axis;
*by      - interval between each tick on y axis;
*len     - x-axis length;


%macro estimates(tertile=,biom=,date=,Var=, Label=,data=,num=,label2=,age_gen=age gen,miny=20,maxy=30,by=5,len=11, label3=);


/*Load in dataset*/
data datalong2;
	set analysis.biomarker_mixed_model_data;
	*subset;
	if subgroup_definition="Sporadic PD" and BL_NSD=1;

	*early vars;
	Early_biom=Early_&biom.;
	Early_tertile=Early_&tertile.;

	proc sort; by Patno event_ID;
run;


/*Create formats to be used for figures*/
proc format;
value $event_ID
"BL" ="BL"
	"V02" = "M6"
	"V04" = "Y1"
	"V06" = "Y2"
	"V08" = "Y3"
	"V10" = "Y4"
	"V12" = "Y5"
	"V13" = "Y6"
	"V14" = "Y7"
	"V15" = "Y8"
	"V16" = "Y9"
	"V17" = "Y10"
	"V18" = "Y11"
	"V19" = "Y12"
	"V20" = "Y13"
	"V22" = "Y15"
;
Value cutoffs
1 = "T1" 
2 = "T2"
3 = "T3";
run;





*Dataset to limit visits to 11 years;
Data datalong3;
	set datalong2;


	if event_ID in ("V04" "V06" "V08" "V10" "V12" "V13" "V14" "V15" "V16" "V17" "V18");


run;



*Dataset to limit visits for DAT;
Data datalong4;
	set datalong2;


	if event_ID in ("V04" "V06" "V10");

run;


proc sort data=datalong3;
	by Early_Tertile;
run;

/*Obtain subjects with at least 1 followup visit*/
proc sort data=datalong3 out=w_followup(keep=Patno) nodupkey;
	by patno;
run;



/*Obtain subjects with at least 1 followup visit*/
proc sort data=&data out=w_followup(keep=Patno) nodupkey;
	by patno;
run;

/*keep BL overvations for those with at least 1 follow-up visit*/
Data BL_W_Followup;
	merge Datalong2(where=(Event_ID = "BL")) w_followup(in=in1);
	by patno;
	if in1;
run;

proc sort data=&data;
	by patno event_ID;
run;

/*Obtain a dataset with those at Baseline and with follow-up visit*/
Data Analysisdata;
	merge BL_W_Followup &data;
	by Patno event_ID;
run;


/*Obtain N at baseline and all follow-up visits*/
Proc means data=Analysisdata noprint;
	class Early_Tertile year;
	Var &var.;
	output out=NUMS(where=(_Type_ = 3)) N=Num;
run;

Data temp;
	set &data;
	if not missing(&var.);
run;

ods select none;


/*Mixed model with Variable ~ Variable at BL+ age + sex + Early_tertile + time + Early_tertile*time*/
/*Outputs datasets with model estimates and contrasts*/

proc mixed data=TEMP namelen=100 empirical;
Class patno gen Early_Tertile(ref="3") ;
Model &var. = BL_&var. &age_gen Early_Tertile time Early_Tertile*time / covb S cl chisq outp=predicted_data1;
Random intercept time / subject = patno type=un; 
Estimate "T1 vs. T2"  Early_Tertile*time  1 -1  0  / cl;
Estimate "T1 vs. T3"  Early_Tertile*time  1  0 -1  / cl;
Estimate "T2 vs. T3"  Early_Tertile*time  0  1 -1  / cl;
/*Store the Mixed Model for future predictions*/
store out=MixedModel;
ods output Estimates=results1  tests3=t1 SolutionF=Solutions;
run;
ods select all;


/*Obtain Average continuous covariates for those with at least year 1 follow-up*/
proc means data = Temp(where=(event_ID="V04")) stackods noprint;
	var BL_&var. age ;
	output out= score0_bl mean = BL_&var. age ;
run;

/*Expand this for each category (CUTOFFS) also expand for sex*/
data score0_bl2;
	set score0_bl(in=ina) score0_bl(in=inb) score0_bl(in=inc);
	
	/*One dataset for each tertile*/
	if ina then Early_Tertile = 1;
	else if inb then Early_Tertile = 2;
	else if inc then Early_Tertile = 3;
	
	/*Obtain all combinations of Sex and Time*/
	do gen=1 to 2 by 1; 
		time = 0;output;
		time = 1;output;
		time = 2;output;
		time = 3;output;
		time = 4;output;
		time = 5;output;
		time = 6;output;
		time = 7;output;
		time = 8;output;
		time = 9;output;
		time = 10;output;
		time = 11;output;
	end;
	format Early_Tertile age BL_&var.;
run;

ods select none;
/*Obtain estimates using our stored mixed model*/
proc plm restore = MixedModel;
	score data = score0_bl2 out = preds3 predicted lclm uclm;
run;

/*Take average across sexes*/
proc means data=preds3 noprint;
	class Early_Tertile time;
	Var Predicted LCLM UCLM;
	output out=preds4(where=(_TYPE_ = 3)) mean=Predicted LCLM UCLM;
run;

/*merge with the numbers at each visit*/
Data Preds4b;
	merge preds4 Nums(rename=(Year=Time));
	by Early_Tertile time;
	if Early_Tertile = 1 then UCLM1 = UCLM;
	if Early_Tertile = 2 then UCLM2 = UCLM;
	if Early_Tertile = 3 then UCLM3 = UCLM;

	if Early_Tertile = 1 then LCLM1 = LCLM;
	if Early_Tertile = 2 then LCLM2 = LCLM;
	if Early_Tertile = 3 then LCLM3 = LCLM;
run;

proc sql noprint;
	select ProbF into :p_val1
	from t1
	where Effect="time*Early_tertile";
quit;

proc sql noprint;
	select Probt into :p_val2
	from results1
	where Label="T1 vs. T2";
quit;

proc sql noprint;
	select Probt into :p_val3
	from results1
	where Label="T1 vs. T3";
quit;

proc sql noprint;
	select Probt into :p_val4
	from results1
	where Label="T2 vs. T3";
quit;


%if %sysevalf(&&p_val1 < 0.05) %then %do;
    %let line1 = Interaction p-values | Overall %str(=) &&p_val1 |;
	%let line2 =  T1 vs T2 %str(=) &&p_val2 |;
	%let line3 =  T1 vs T3 %str(=) &&p_val3 |;
	%let line4 =  T2 vs T3 %str(=) &&p_val4;
%end;
%else %do;
        %let line1 = Overall Interaction p-value %str(=) &&p_val1;
		%let line2 =;
		%let line3 =;
		%let line4 =;
%end;


ods select all;

/*output images to filepath*/
ods listing gpath="&output_path." image_dpi=600;

ods graphics / reset
    imagename="&label3. &label."
    imagefmt=png
    width=6in
    height=4in;
/*Plot creation*/
proc sgplot data=preds4b noautolegend;

    series x=time y=Predicted / 
        group=Early_tertile
        name="series";

    styleattrs 
        datacolors=(blue red green) 
        datacontrastcolors=(blue red green);

    band lower=LCLM1 upper=UCLM1 x=time /
        group=Early_tertile
        fillattrs=(transparency=.8 color=blue);

    band lower=LCLM2 upper=UCLM2 x=time /
        group=Early_tertile 
        fillattrs=(transparency=.8 color=red);

    band lower=LCLM3 upper=UCLM3 x=time /
        group=Early_tertile
        fillattrs=(transparency=.8 color=green);

    format Early_tertile cutoffs.;

xaxis values=(0 to &len by 1)
    labelattrs=(size=14pt weight=bold)
    valueattrs=(size=12pt);

yaxis values=(&miny to &maxy by &by)
    labelattrs=(size=14pt weight=bold)
    valueattrs=(size=12pt);


    label 
        Predicted="Mean Predicted &label2."
        time="Year";

xaxistable Num /
    location=inside
    x=time
    class=Early_tertile
    colorgroup=Early_tertile
    valueattrs=(size=12pt);

    *title "Clinical Outcome model for &label.";
	*Footnote "Line1. Line2. Line3. Line4.";
run;



ods listing close;

%mend;




/*MOCA calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=MOCA, Label = MOCA,data=datalong3,num=1, label2=MOCA,label3=P-Tau217);

%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=MOCA, Label = MOCA,data=datalong3,num=2, label2=MOCA,label3=NFL);

%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=MOCA, Label = MOCA,data=datalong3,num=3, label2=MOCA,label3=GFAP);

%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=MOCA, Label = MOCA,data=datalong3,num=6, label2=MOCA,label3=Aß42/40);

/*MDS I score calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=updrs1_score, Label = MDS ,data=datalong3,num=1, label2=MDS-UPDRS I,miny=0,maxy=20,by=5,label3=P-Tau217);

%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=updrs1_score, Label = MDS ,data=datalong3,num=2, label2=MDS-UPDRS I,miny=0,maxy=20,by=5,label3=NFL);

%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=updrs1_score, Label = MDS ,data=datalong3,num=3, label2=MDS-UPDRS I,miny=0,maxy=20,by=5,label3=GFAP);


%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=updrs1_score, Label = MDS ,data=datalong3,num=6, label2=MDS-UPDRS I,miny=0,maxy=20,by=5,label3=Aß42/40);

/*MDS III score calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=updrs3_score, Label = MDS OFF ,data=datalong3,num=1, label2=MDS-UPDRS III OFF,miny=0,maxy=60,by=10,label3=P-Tau217);

%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=updrs3_score, Label = MDS OFF ,data=datalong3,num=2, label2=MDS-UPDRS III OFF,miny=0,maxy=60,by=10,label3=NFL);

%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=updrs3_score, Label = MDS OFF ,data=datalong3,num=3, label2=MDS-UPDRS III OFF,miny=0,maxy=60,by=10,label3=GFAP);

%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=updrs3_score, Label = MDS OFF ,data=datalong3,num=6, label2=MDS-UPDRS III OFF,miny=0,maxy=60,by=10,label3=Aß42/40);

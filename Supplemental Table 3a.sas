*******************************;
* Author: Noah Johnson        *;
* File:   				      *;
* Supplemental Table 3a       *;
* 							  *;
* Description:				  *;
* Generates continuous sPD    *;
* supplemental mixed model    *;
* table          			  *;
*******************************;



%let data_filepath= /*path to data*/;
*obtain necessary datasets;
libname analysis "&data_filepath.";




/*Macro generates Mixed Model estimates for given biomarker and longitudinal outcome*/

*inputs are;

*tertile - tertile biomarker;
*biom    - numeric biomarker;
*date    - date biomarker was obtained;
*Var     - longitudinal outcome of interest;
*Label   - Title for longitudinal outcome in output file;
*data    - longitudinal data specified;
*Label   - longitudinal data specified;
*label2  - Label of biomarker in legend;
*age_gen - determines whether to adjust by age and sex;
*num     - used to sort estimates in dataframe at later step by biomarker;
*num2    - used to sort estimates in dataframe at later step by variable;

%macro estimates(tertile=,biom=,date=,Var=, Label=,data=,num=,label2=,num2=,age_gen=age gen);

/*Load in dataset*/
data datalong2;
	set analysis.biomarker_mixed_model_data;
	*subset;
	if subgroup_definition="Sporadic PD" and BL_NSD=1;

	*early biomarkers;
	Early_biom=Early_&biom.;
	Early_tertile=Early_&tertile.;

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
	by Early_biom;
run;
/*Obtain subjects with at least 1 followup visit*/
proc sort data=datalong3 out=w_followup(keep=Patno) nodupkey;
	by patno;
run;


proc sort data=&data;
	by Early_Tertile;
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




Data temp;
	set &data;
	if not missing(&var.);
run;

*ods select none;


/*Mixed model with Variable ~ Variable at BL+ age + sex + Early biomarker + time + Early biomarker*time*/
/*Outputs datasets with model estimates and contrasts*/


proc mixed data=TEMP empirical namelen=100  plots(MAXPOINTS=none)=ALL;
Class patno gen  ;
Model &var. = BL_&var. &age_gen Early_biom time Early_biom*time / covb S cl chisq outp=predicted_data1;
Random intercept time / subject = patno type=un; 

store out=MixedModel;

ods output Estimates=results1  tests3=t1 SolutionF=Solutions;
run;


ods select all;


/*format results*/
	data tc_&num._&num2;
		set Solutions;
		length Comp Biomarker Outcome $50.;
		order=&num2;
		main=&num;
		if effect="Early_biom*time";
		Outcome="&label.";
		Comp="Overall";
		Biomarker="&label2.";
		Est_CI = cat(
		    strip(ifc(abs(Estimate) < 0.0001,
		              put(Estimate, E12.4),
		              put(Estimate, 12.4))),
		    " (",
		    strip(ifc(abs(Lower) < 0.0001,
		              put(Lower, E12.4),
		              put(Lower, 12.4))),
		    ",",
		    strip(ifc(abs(Upper) < 0.0001,
		              put(Upper, E12.4),
		              put(Upper, 12.4))),
		    ")"
		);
		*if Label in ("1 vs. 2" "1 vs. 3" "2 vs. 3");
		SE=put(StdErr,6.4);
		Pval=put(Probt,6.4);
		keep Biomarker Outcome Comp Est_CI SE tValue Pval main order;
	run;


%mend;


/*MDS I calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=updrs1_score, Label = MDS-UPDRS I Score,data=datalong3,num=1, label2=pTau217,num2=1,age_gen= age gen);
%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=updrs1_score, Label = MDS-UPDRS I Score,data=datalong3,num=2, label2=NFL,num2=1,age_gen= age gen);
%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=updrs1_score, Label = MDS-UPDRS I Score,data=datalong3,num=3, label2=GFAP,num2=1,age_gen= age gen);
%estimates(tertile=C2N_AB40_cutoff,biom=C2N_AB40,Var=updrs1_score, Label = MDS-UPDRS I Score,data=datalong3,num=4, label2=Aß40 ,num2=1,age_gen= age gen);
%estimates(tertile=C2N_AB42_cutoff,biom=C2N_AB42,Var=updrs1_score, Label = MDS-UPDRS I Score,data=datalong3,num=5, label2=Aß42 ,num2=1,age_gen= age gen);
%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=updrs1_score, Label = MDS-UPDRS I Score,data=datalong3,num=6, label2=Aß42/Aß40,num2=1,age_gen= age gen);

/*MDS II calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=updrs2_score, Label = MDS-UPDRS II Score,data=datalong3,num=1, label2=pTau217,num2=2,age_gen= age gen);
%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=updrs2_score, Label = MDS-UPDRS II Score,data=datalong3,num=2, label2=NFL,num2=2,age_gen= age gen);
%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=updrs2_score, Label = MDS-UPDRS II Score,data=datalong3,num=3, label2=GFAP,num2=2,age_gen= age gen);
%estimates(tertile=C2N_AB40_cutoff,biom=C2N_AB40,Var=updrs2_score, Label = MDS-UPDRS II Score,data=datalong3,num=4, label2=Aß40 ,num2=2,age_gen= age gen);
%estimates(tertile=C2N_AB42_cutoff,biom=C2N_AB42,Var=updrs2_score, Label = MDS-UPDRS II Score,data=datalong3,num=5, label2=Aß42 ,num2=2,age_gen= age gen);
%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=updrs2_score, Label = MDS-UPDRS II Score,data=datalong3,num=6, label2=Aß42/Aß40,num2=2,age_gen= age gen);

/*MDS III calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=updrs3_score, Label = MDS-UPDRS III Score OFF,data=datalong3,num=1, label2=pTau217,num2=3,age_gen= age gen);
%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=updrs3_score, Label = MDS-UPDRS III Score OFF,data=datalong3,num=2, label2=NFL,num2=3,age_gen= age gen);
%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=updrs3_score, Label = MDS-UPDRS III Score OFF,data=datalong3,num=3, label2=GFAP,num2=3,age_gen= age gen);
%estimates(tertile=C2N_AB40_cutoff,biom=C2N_AB40,Var=updrs3_score, Label = MDS-UPDRS III Score OFF,data=datalong3,num=4, label2=Aß40 ,num2=3,age_gen= age gen);
%estimates(tertile=C2N_AB42_cutoff,biom=C2N_AB42,Var=updrs3_score, Label = MDS-UPDRS III Score OFF,data=datalong3,num=5, label2=Aß42 ,num2=3,age_gen= age gen);
%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=updrs3_score, Label = MDS-UPDRS III Score OFF,data=datalong3,num=6, label2=Aß42/Aß40,num2=3,age_gen= age gen);


/*MOCA calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=MOCA, Label = MOCA,data=datalong3,num=1, label2=pTau217,num2=4,age_gen= age gen);
%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=MOCA, Label = MOCA,data=datalong3,num=2, label2=NFL,num2=4,age_gen= age gen);
%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=MOCA, Label = MOCA,data=datalong3,num=3, label2=GFAP,num2=4,age_gen= age gen);
%estimates(tertile=C2N_AB40_cutoff,biom=C2N_AB40,Var=MOCA, Label = MOCA,data=datalong3,num=4, label2=Aß40 ,num2=4,age_gen= age gen);
%estimates(tertile=C2N_AB42_cutoff,biom=C2N_AB42,Var=MOCA, Label = MOCA,data=datalong3,num=5, label2=Aß42 ,num2=4,age_gen= age gen);
%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=MOCA, Label = MOCA,data=datalong3,num=6, label2=Aß42/Aß40,num2=4,age_gen= age gen);

/*DAT calls*/
%estimates(tertile=Ptau_cutoff,biom=Lilly_ptau217p,Var=lowput_ratio_cwm, Label = Lowest Putamen Ratio,data=datalong4,num=1, label2=pTau217,num2=5,age_gen=);
%estimates(tertile=NFL_mean_cutoff,biom=IU_NFL_plasma_mean,Var=lowput_ratio_cwm, Label = Lowest Putamen Ratio,data=datalong4,num=2, label2=NFL,num2=5,age_gen=);
%estimates(tertile=GFAP_mean_cutoff,biom=IU_GFAP_plasma_mean,Var=lowput_ratio_cwm, Label = Lowest Putamen Ratio,data=datalong4,num=3, label2=GFAP,num2=5,age_gen=);
%estimates(tertile=C2N_AB40_cutoff,biom=C2N_AB40,Var=lowput_ratio_cwm, Label = Lowest Putamen Ratio,data=datalong4,num=4, label2=Aß40 ,num2=5,age_gen=);
%estimates(tertile=C2N_AB42_cutoff,biom=C2N_AB42,Var=lowput_ratio_cwm, Label = Lowest Putamen Ratio,data=datalong4,num=5, label2=Aß42 ,num2=5,age_gen=);
%estimates(tertile=C2N_AB4240_cutoff,biom=C2N_AB4240,Var=lowput_ratio_cwm, Label = Lowest Putamen Ratio,data=datalong4,num=6, label2=Aß42/Aß40,num2=5,age_gen=);


*compile all data together and format;
data output;
	set tc_1: tc_2: tc_3: tc_4: tc_5: tc_6:;
		if order^=1 then Biomarker="";

	retain overall_p;
    /* Save overall p-value */
    if Comp = "Overall" then overall_p = pval;

    /* Remove tertile comparisons if overall p-value > 0.05 */
    if overall_p > 0.05 and Comp in ("T1 vs. T2", "T1 vs. T3", "T2 vs. T3") then delete;

	proc sort; by main order  descending outcome Comp;
run;


/*Reports data in desired format*/


ods escapechar="^";
Title "^S={foreground=black}Supplemental Table 3a: All Continuous sPD Mixed Models Summary";
Footnote1 H=8PT J=LEFT F="Arial" "^S={foreground=black}All Mixed models were adjusted for age and sex except for Lowest Putamen ratio.";	
Footnote2 H=8PT J=LEFT F="Arial" "^S={foreground=black}Estimates and p-values are based on time and continuous biomarker interaction."  ;



proc report data=output nowd split="|";

    column Biomarker Outcome  Est_CI SE tValue Pval;

	define Outcome /  style(column)={width=1.5in};

    define Est_CI / "Coefficient (95% CI)";
    define tValue / "T-stat" style(column)={width=.9in};
    define Pval   / "p-value";

    /* Shade + bold rows where Biomarker is present */
    compute Biomarker;
        if not missing(Biomarker) then do;
            call define(_col_, "style",
                "style={font_weight=bold}");

        end;
    endcomp;

    compute Pval;
        if Pval<0.05 then do;
            call define(_col_, "style",
                "style={font_weight=bold}");

        end;
    endcomp;


run;



*categorical variable descriptive summary;
%macro cat(data = , var = , num = , dec = 0);

	proc freq data=&data;
		table group*&var;
		ods output CrossTabFreqs = freq_&num;
	run;

	proc freq data=&data;
		table group*&var / missing;
		ods output CrossTabFreqs = freq_m_&num;
	run;

	data freq2_&num (keep=group main order n_per);
		*set freq_&num (where=(not missing(rowpercent)))
			freq_m_&num (where=(missing(&var) AND not missing(rowpercent)));
		set freq_&num (where=(not missing(group) AND not missing(&var)))
			freq_m_&num (where=(missing(&var) AND not missing(rowpercent)));


		if not missing(&var) then do;
			if 0.00 < RowPercent < 0.50 then n_per = strip(put(frequency, 5.)) || ' (<1%)';
			else 					  n_per = strip(put(frequency, 5.)) || ' (' || strip(put(rowpercent, 7.&dec.)) || '%)';
		end;
		else n_per = strip(put(frequency, 5.));
		*address scenarios where all participants in a group had missing values (no % calculated);
		if n_per = '0 (%)' then n_per = 'NA';
		if n_per = '0 ( %)' then n_per = 'NA';

		main = &num;
		if not missing(&var) then order = &var;
		else if missing(&var) then order = 99;
				new = order*1;
		drop order;
		rename new=order;

	run;

	proc sort data=freq2_&num; by main order group; run;

	proc transpose data=freq2_&num out=freq3_&num (drop=_name_);
		by main order; id group;
		var n_per;
	run;

%mend cat;

*percent only (categorical variable descriptive summary);
%macro pct(data = , var = , num = , dec = 1);

	proc freq data=&data;
		table group*&var;
		ods output CrossTabFreqs = freq_&num;
	run;

	proc freq data=&data;
		table group*&var / missing;
		ods output CrossTabFreqs = freq_m_&num;
	run;

	data freq2_&num (keep=group main order per);
		set freq_&num (where=(not missing(rowpercent)))
			freq_m_&num (where=(missing(&var) AND not missing(rowpercent)));

		if not missing(&var) then per = strip(put(rowpercent, 7.&dec.)) || '%';
		else per = strip(put(frequency, 5.));

		main = &num;
		if not missing(&var) then order = &var;
		else if missing(&var) then order = 99;
				new = order*1;
		drop order;
		rename new=order;
	run;

	proc sort data=freq2_&num; by main order group; run;

	proc transpose data=freq2_&num out=freq3_&num (drop=_name_);
		by main order; id group;
		var per;
	run;

%mend pct;

*continuous variable descriptive summary;
%macro cont(data = , var = , dec = , dec2 = , num = );

	proc means data=&data;
		class group;
		var &var;
		output out=stats_&num mean = mean  stddev = std
							  median = median q1 = q1 q3 = q3 min = minimum  max = maximum
							  nmiss = missing  n = n;
	run;

	data stats2_&num (keep=main ave med rge miss group);
		set stats_&num (where=( not missing(group) ));

		if not missing(mean) then do;
			ave = strip(put(mean, 9.&dec.)) || " (" || compress(put(std, 9.&dec.))|| ")";
			med = strip(put(median, 9.&dec2.)) || " (" || compress(put(q1, 9.&dec2.))|| "-" || compress(put(q3, 9.&dec2.)) || ")";
			rge = strip(put(minimum, 9.&dec2.)) || "-" || compress(put(maximum, 9.&dec2.));
		end;
		miss = strip(put(missing, 5.));

		main = &num;
	run;

	proc sort data=stats2_&num; by main group; run;

	proc transpose data=stats2_&num out=stats3_&num (drop=_name_);
		by main; id group;
		var ave med rge miss;
	run;

	data stats4_&num;
		set stats3_&num;
		order = _N_ - 1;
	run;

%mend cont;

*continuous variable descriptive summary including mean (SD) and median (min, max) and missing;
%macro cont2(data = , var = , dec = , dec2 = , num = );

	proc means data=&data;
		class group;
		var &var;
		output out=stats_&num mean = mean  stddev = std
							  median = median min = minimum  max = maximum
							  nmiss = missing  n = n;
	run;

	data stats2_&num (keep=main ave med miss group);
		set stats_&num (where=( not missing(group) ));

		*address n=1 missing SD scenarios;
		std_char = compress(put(std, 9.&dec.));
		if n = 1 AND missing(std) then std_char = "NA";

		if not missing(mean) then do;
			ave = strip(put(mean, 9.&dec.)) || " (" || compress(std_char) || ")";
			med = strip(put(median, 9.&dec2.)) || " (" || compress(put(minimum, 9.&dec2.))|| ", " || compress(put(maximum, 9.&dec2.)) || ")";
		end;
		miss = strip(put(missing, 5.));

		*address scenarios where all participants in a group had missing values (no % calculated);
		if _FREQ_ = missing then do;
			ave = "NA";
			med = "NA";	
		end;

		main = &num;
	run;

	proc sort data=stats2_&num; by main group; run;

	proc transpose data=stats2_&num out=stats3_&num (drop=_name_);
		by main; id group;
		var ave med miss;
	run;

	data stats4_&num;
		set stats3_&num;
		order = _N_ - 1;
	run;

%mend cont2;



*Chi-square test (frequencies, percents, p value);
%macro chisq(data = , var = , num = , dec = 0);

	proc freq data=&data;
		table group*&var / chisq;
		ods output CrossTabFreqs = freq_&num
				   ods output ChiSq = chi1_&num;
	run;

	proc freq data=&data;
		table group*&var / missing;
		ods output CrossTabFreqs = freq_m_&num;
	run;

	data freq2_&num (keep=group main order n_per);
		set freq_&num (where=(not missing(rowpercent)))
			freq_m_&num (where=(missing(&var) AND not missing(rowpercent)));

		if not missing(&var) then
			n_per = strip(put(frequency, 5.)) || ' (' || strip(put(rowpercent, 7.&dec)) || '%)';
		else n_per = strip(put(frequency, 5.));

		main = &num;
		if not missing(&var) then order = &var;
		else if missing(&var) then order = 99;
				new = order*1;
		drop order;
		rename new=order;
	run;

	proc sort data=freq2_&num; by main order group; run;

	proc transpose data=freq2_&num out=freq3_&num (drop=_name_);
		by main order; id group;
		var n_per;
	run;

	data chi_&num (keep=main order pval );
		set chi1_&num (where=(Statistic = 'Chi-Square'));
		length pval $8.;
		pval = strip(put(Prob, PVALUE6.4));
		main = &num - 0.5; order = 0;
	run;

	data freq4_&num;
		set chi_&num freq3_&num;
		if order ne 99 then order = _N_ - 1;
	run;

%mend chisq;

*customized chi-square macro for table1;
%MACRO chisq2(var = , num = , main = , order = 0);

	*females;
	proc freq data=popn(where=(group IN (3 4)));
		table group*&var / chisq;
		ods output ChiSq = chi1_f_&num(rename=(Prob = Prob_F)
									   where=(Statistic = 'Chi-Square'));
	run;

	*males;
	proc freq data=popn(where=(group IN (5 6)));
		table group*&var / chisq;
		ods output ChiSq = chi1_m_&num(rename=(Prob = Prob_M)
									   where=(Statistic = 'Chi-Square'));
	run;
	
	*merge;
	data chi2_&num (keep=main order pval: );
		merge chi1_f_&num(keep=Prob_F)
			  chi1_m_&num(keep=Prob_M);
		main   = &main;
		order  = &order;
		pval_f = strip(put(Prob_F, PVALUE6.4));
		pval_m = strip(put(Prob_M, PVALUE6.4));
	run;

%MEND;

*customized t-test macro for table1;
%MACRO ttest2(var = , num = , main = , order = 0);

	*females;
	proc ttest data=popn(where=(group IN (3 4))) plots=none;
		class group;
		var &var;
		ods output ttests = ttest1_f_&var
			       equality = equalvar_f_&var;
	run;

	data ttest_f_&var(keep=pval_f);
		merge ttest1_f_&var
			  equalvar_f_&var(keep=Variable ProbF);
		by Variable;
		if ProbF < 0.05 AND Variances = "Equal" then delete;
		if ProbF ge 0.05 AND Variances = "Unequal" then delete;
		pval_f = strip(put(probt, PVALUE6.4));	
	run;

	*males;
	proc ttest data=popn(where=(group IN (5 6))) plots=none;
		class group;
		var &var;
		ods output ttests = ttest1_m_&var
			       equality = equalvar_m_&var;
	run;

	data ttest_m_&var(keep=pval_m);
		merge ttest1_m_&var
			  equalvar_m_&var(keep=Variable ProbF);
		by Variable;
		if ProbF < 0.05 AND Variances = "Equal" then delete;
		if ProbF ge 0.05 AND Variances = "Unequal" then delete;
		pval_m = strip(put(probt, PVALUE6.4));	
	run;

	*merge;
	data ttest2_&num;
		merge ttest_f_&var
			  ttest_m_&var;
		main   = &main;
		order  = &order;
	run;

%MEND;


*Kruskal-Wallis test;
%macro KruskalWallis(data = , var = , dec = , dec2 = , num = , test = KruskalWallisTest, order = 1);

	proc means data=&data;
		class group;
		var &var;
		output out=stats_&num mean = mean  stddev = std
							  median = median min = minimum  max = maximum
							  nmiss = missing  n = n;
	run;

	data stats2_&num (keep=main ave med miss group);
		set stats_&num (where=( not missing(group) ));

		if not missing(mean) then do;
			ave = strip(put(mean, 9.&dec.)) || " (" || compress(put(std, 9.&dec.))|| ")";
			med = strip(put(median, 9.&dec2.)) || " (" || compress(put(minimum, 9.&dec2.))|| ", " || compress(put(maximum, 9.&dec2.)) || ")";
		end;
		miss = strip(put(missing, 5.));

		main = &num;
	run;

	proc sort data=stats2_&num; by main group; run;

	proc transpose data=stats2_&num out=stats3_&num (drop=_name_);
		by main; id group;
		var ave med miss;
	run;

	data stats3b_&num;
		set stats3_&num;
		order = _N_;
	run;

	*pvalue;
	proc npar1way data=&data wilcoxon /*plots=none*/;
		class group;
		var &var;
		ods output &test = pval1_&var;
	run;

	data pval_&var(keep=main pval order);
		set pval1_&var/*(where=(Name1 = "P2_WIL"))*/;

		pval = strip(put(Prob, PVALUE6.4));	
		main = &num;
		order = &order; 
	run;

	data stats4_&num;
		merge pval_&var stats3b_&num;
		by main order;
		order = _N_;
	run;

%mend KruskalWallis;

*Fisher's exact test (frequencies, percents, p value);
%macro cat_exact(data = , var = , num = , dec = 0, order = 0);

	proc freq data=&data;
		table group*&var / chisq exact;
		ods output CrossTabFreqs = freq_&num
				   FishersExact = chi1_&num;
	run;

	proc freq data=&data;
		table group*&var / missing;
		ods output CrossTabFreqs = freq_m_&num;
	run;

	data freq2_&num (keep=group main order n_per);
		set freq_&num (where=(not missing(rowpercent)))
			freq_m_&num (where=(missing(&var) AND not missing(rowpercent)));

		if not missing(&var) then
			n_per = strip(put(frequency, 5.)) || ' (' || strip(put(rowpercent, 7.&dec.)) || '%)';
		else n_per = strip(put(frequency, 5.));

		main = &num;
		if not missing(&var) then order = &var;
		else if missing(&var) then order = 99;
				new = order*1;
		drop order;
		rename new=order;
	run;

	proc sort data=freq2_&num; by main order group; run;

	proc transpose data=freq2_&num out=freq3_&num (drop=_name_);
		by main order; id group;
		var n_per;
	run;

	data chi_&num (keep=main order pval );
		*change the following to work with both 2-group and 3-group comparisons;
		set chi1_&num (where=(Name1 = 'XP2_FISH'));
		*set chi1_&num (where=(label1 = 'Pr <= P'));
		length pval $8.;

		pval = strip(put(nValue1, PVALUE6.4));
		*if pval = '0.0000' then pval = '<.0001';

		/*
		if nValue1 < .0001 then pval = '<0.0001';
		else if nValue1 < .001 then pval = '<0.001';
		else if nValue1 < 0.1 then pval = strip(put(nValue1, pvalue5.3));
		else if nValue1 ge 0.1 then pval = strip(put(nValue1, 4.2));
		*/

		main = &num; order = &order;
	run;

	data freq4_&num;
		merge chi_&num freq3_&num;
		by main order;
		*if order ne 99 then order = _N_ - 1;
	run;

%mend cat_exact;




/** Codes for Pheno Conversion **/

%macro CHAR_STATS(inData=, var=, group=, rowID=, incMiss=1);
   /* Step 1: Sort the input dataset by patno and the grouping variable */
   proc sort data=&inData. out=_sortedData;
      by patno &group.;
   run;

   ods select none;

   /* Step 2: Use PROC FREQ with different options based on incMiss */
   %if &incMiss = 1 %then %do;
      proc freq data=_sortedData;
         tables &group.*&var. / missing;
         ods output CrossTabFreqs = _freq (where=(_TYPE_ in ("11")));
      run;

      data _freqOut;
         set _freq;
         /* Format the row percentage */
         row_pct = put(RowPercent, 5.);
         if Frequency ne 0 and strip(row_pct) = "0" then row_pct = "<1";
         /* Create a combined frequency/percentage string */
         n_per = strip(put(Frequency, 5.)) || " (" || strip(row_pct) || "%)";
         if n_per = "0 (0%)" then n_per = "0";
         /* For nonmissing values, use &var; missing gets recoded as 99 */
         if not missing(&var.) then i = &var.;
         else i = 99;
         rowID = &rowID.;
         keep &group. n_per i rowID;
      run;
   %end;
   %else %do;
      /* When missing values are not to be included in row percentages */
      proc freq data=_sortedData;
         tables &group.*&var.;
         ods output CrossTabFreqs = _freq (where=(_TYPE_ in ("11")));
      run;

      proc freq data=_sortedData;
         tables &group.*&var. / missing;
         ods output CrossTabFreqs = _freqMiss (where=(_TYPE_ in ("11") and missing(&var.)) drop=RowPercent);
      run;

      data _freqOut;
         set _freq _freqMiss;
         row_pct = put(RowPercent, 5.);
         if Frequency ne 0 and strip(row_pct) = "0" then row_pct = "<1";
         n_per = strip(put(Frequency, 5.)) || " (" || strip(row_pct) || "%)";
         if n_per = "0 (0%)" then n_per = "0";
         if Frequency = 0 and missing(RowPercent) then n_per = "0";
         if not missing(&var.) then i = &var.;
         else do;
            i = 99;
            n_per = strip(put(Frequency, 5.));
         end;
         rowID = &rowID.;
         keep &group. n_per i rowID;
      run;
   %end;

   proc sort data=_freqOut;
      by i;
   run;

   /* Step 3: Transpose the formatted frequency results so that each level of the group appears as a column */
   proc transpose data=_freqOut out=charFreq_&var.(drop=_NAME_);
      by rowID i;
      id &group.;
      var n_per;
   run;

   ods select all;
%mend CHAR_STATS;

/*-------------------------*/
/* Numeric Summary Macro   */
/*-------------------------*/
%macro NUM_STATS(indata=, var=, rowID=);
   proc sort data=&indata. out=sortedNum;
      by group1;
   run;
   
   proc means data=sortedNum mean std median p25 p75 min max nmiss noprint;
      by group1;
      var &var.;
      output out=summary mean=mean std=sd median=median p25=p25 p75=p75 min=min max=max nmiss=nmis;
   run;
   
   data summary2 (keep=group1 ave med range miss);
      set summary;
	ave = strip(put(mean, 5.1)) || " (" || strip(put(sd,5.1)) || ")"; 
	if missing(sd) then ave = strip(put(mean, 5.1)) || " (" || "NA" || ")"; 
	iqr = strip(put(p25,3.)) || ", " || strip(put(p75,3.));
	med = strip(put(median, 5.1)) || " (" || strip(iqr) || ")";
	range = " (" || strip(put(min,3.)) || ", " || strip(put(max,3.)) || ")";
	miss = strip(put(nmis, 3.)); 
      if missing(miss) then miss = "0";
      if missing(mean) then do;
         ave = "-";
         med = "-";
         range = "-";
      end;
   run;
   
   proc transpose data=summary2 out=summary3 (drop=_NAME_);
      id group1;
      var ave med range miss;
   run;
   
   data summary_&var.;
      set summary3;
      i = _N_;
      rowID = &rowID;
   run;
   
   proc sort data=summary_&var.;
      by rowID i;
   run;
%mend NUM_STATS;

This is a guide to reproduce figures from the paper.

Code was written and executed in Matlab 7.12.0 (R2011a).

Figure 1:

-Execute,'Fig1a','Fig1b','Fig1c','Fig1e' on the command line (do not include quotes).

Figure 2: 

-These are freehand figures generated in PowerPoint

Figure 3a: 

-For Fraction Random Parameter Sets Reaching a Steady State:
	+Use the 'TimeDomain_SS_Fraction.m' function.
	+On the command line, type TimeDomain_SS_Fraction('[mat-file name]')
	+Use the mat-file name in question according to the following:

Gluc to isoprene = 'GlucoseToIsoprene'
Gluc-to-isoprene(with regulation) = 'GlucoseToIsopreneRegG6P'
Chimeric glycolysis = 'Chimeric_Glycolysis'
Gluc to poly-hydroxybutyrate = 'Purge_Valve'
MCC (Fpk/Xpk) 1:3 = 'MCC_Xpk_Fpk'
MCC (Xpk only) = 'MCC_Xpk'
MCC (Fpk only) = 'MCC_Fpk'

Results will be displayed as a 3-tuple containing the number of steady states observed out of 1000 for each of three repeats.

-For Fraction Random Parameter Sets Reaching a Steady State:
	+Use the 'Stable_Fraction.m' function.
	+On the command line, type Stable_Fraction('[mat-file name]')
	+Use the mat-file name in question according to the same list as above.

Results will be displayed as a fraction of stable steady states observed out of 1000 for one repeat.

Note: leave out the .mat from the file name i.e.: TimeDomain_SS_Fraction('Purge_Valve')

Figure 3b: Freehand figure generated from PowerPoint

Figure 4: 

-Execute,'Fig4a','Fig4b','Fig4c','Fig4d','Fig4e','Fig4f' on the command line (do not include quotes).

Figs 5b, 6b, 7b: 

- For EMRA profile (Fig 5b, 6b, 7b), Execute Main('[mat-file name]') for the appropriate model.

- Then run Plot_Robustness after pasting the result file name ('Model[mat-file name] Results 1000') into the curly brackets on line 3.

- The double plot in 7b can be generated by running Plot_Robustness twice for the two different result sets, and deleting 'figure' on line 14.  The line color change is accomplished by changing 'b' to 'r' on line 31.

Fig 6d,e:

-Execote 'Fig6de' on the command line.
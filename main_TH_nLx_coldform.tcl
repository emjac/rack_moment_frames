#################################################################################################################
# Procedure to run a series of time history analyses
#
# EJ
#################################################################################################################

# Set rack properties
set h1 [expr 1./10000*round(10000*68*25.4/1000)]; # First story height [m]
set hx [expr 1./10000*round(10000*68*25.4/1000)]; # Upper story heights [m]
set L [expr 1./10000*round(10000*99.0*25.4/1000)]; # Beam length c/c dimensions are used [m]
set beamPrfl "poly-cold-form-beam";
set colPrfl "poly-cold-form-col";
set Fyc 345e6;
set Fyb 345e6;
set boxed_levels 0;
set PL 13345.; # Metric weight of one 3000 lb pallet [N]
set nB 6; # Number of bays [unitless]
set zeta 0.02; # Percent of critical damping (for Rayleigh damping). Ex. 0.02 means 2% of critical.
set muS 0.1; 
set muK 0;
set sprOption "RBS44";
set bplOption "Parallel";
set kc 1.;
set kctop 1.;
set NL_geo "LargePdelta";
set g_tol 0.;
set l_tol 0.;

# Directory where ground motions are kept
set directory [eval pwd]
set nL_list { 3 4 5 6 };
# set kbpl_list { 1.55 1.77 1.98 2.21 };
set kbpl_list { 4.23 4.45 4.67 4.89 }

#################################################################################################################

# Include procedures for a number of different analyses
lappend ::auto_path [eval pwd]
package require rack_pack
namespace import rack_pack::* 

set zeta_list [list 0.01 0.02 0.03 0.04 0.05]

for { set i 0 } { $i <= 3 } { incr i } {
	# for { set j 0 } { $j <= 0 } { incr j } {
		# for { set k 0 } { $k <= 4 } { incr k } {
			
			# Parameters for TH analysis
			
			set directory "Accelerograms/";
			set accRecord "Re_Scaled_F19.6_A112.2_BRUT_Trial1_s.acc"; set dt 0.002;
			# set accRecord "Re_Scaled_NGAKKiK_Major_H2_623_s.acc"; set dt 0.005;
			set nL [lindex $nL_list $i];
			set kbpl [lindex $kbpl_list $i];
			set SF 2.5;
			# set zeta [lindex $zeta_list $k];
		 
			# Create rack model
			rackMRF $h1 $hx $L $beamPrfl $Fyb $colPrfl $Fyc $boxed_levels $PL $nB $nL $zeta $muS $muK $sprOption $bplOption $kc $kctop $kbpl $NL_geo $g_tol $l_tol "$accRecord,muS$muS,SF$SF,tp2";
			
			# DisplayModel2D
			# Apply gravity loads
			record_rack "basenodes" "reaction";
			rack_gravity_static [expr $PL/6] "on_beams"; # Apply gravity loading to include P-delta effects
			set message  "Total vertical reaction is : "
			set dof 2;
			set step 10;
			check_total_reaction $dof $step; # Sum the data in the last line and show result
			
			set out_step ""
			record_rack "ext_col_nodes" "disp" $out_step;
			record_rack "ext_col_nodes" "acc" $out_step;
			record_rack "connectors" "deformation" $out_step;
			record_rack "columns" "force" $out_step;
			record_rack "basenodes" "reaction" $out_step;
			record_rack "palletnodes" "disp" $out_step;
			record_rack "beamnodes" "disp" $out_step;
			record_rack "beamnodes" "acc" $out_step;
			
			rack_TH $dt $directory $accRecord $SF;
		# }
	# }
}
#################################################################################################################
# Procedure to run a series of time history analyses
#
# EJ
#################################################################################################################

# Set rack properties
set h1 1.7; # First story height [m]
set hx 0.; # Upper story heights [m]
set L [expr 1./10000*round(10000*99.0*25.4/1000)]; # Beam length c/c dimensions are used [m]
set beamPrfl "poly-cold-form-beam";
set colPrfl "poly-cold-form-col";
set Fyc 345e6;
set Fyb 345e6;
set boxed_levels 0;
set PL 13345.; # Metric weight of one 3000 lb pallet [N]
set nB 1; # Number of bays [unitless]
set nL 1; # Number of levels [unitless]
set sprOption "RBS44";
set kc 1.;
set kctop 0.;
set NL_geo "LargePdelta";
set l_tol 0.;
set g_tol 0.;
set zeta 0.02;
set muS 0.17;
set muK 0.8

# Directory where ground motions are kept
set directory [eval pwd]
# Accelerations as output by the shake table
set accl_names { "P3_realtest.acc" "P4_realtest.acc" "P5_realtest.acc" "F4_realtest.acc" "F5_realtest.acc" "F6_realtest.acc" };
# set d0_list [list -0.011 0.0057 0.0015 0. 0. 0.]

#################################################################################################################

# Include procedures for a number of different analyses
lappend ::auto_path [eval pwd]
package require rack_pack
namespace import rack_pack::* 

# set zeta_list [list 0.01 0.02 0.03 0.04 0.05 0.06 ]
# set kc_list [list 35428.7 32886.48 35428.7 49448.92	49145.97 49600.4]
# set kbpl_list [list 0. 0. 0. 49586.44 51812.15 48473.59]

for { set i 0 } { $i <= 6 } { incr i } {
		 
	if { $i > 2 } { 
		set bplOption "Parallel"; set kbpl 0.9; 
		# set sprOption "L"; set bplOption "L"; 
	} else {
		set bplOption "MinMax"; set kbpl 1.0; 
		#set bplOption "L"; set kbpl 0;
	}
 
	# Parameters for TH analysis
	set directory "Accelerograms/";
	set accRecord [lindex $accl_names $i];
	#set muS [lindex $muS_list $j];
	
	set d0 0.;
	set dt 0.005;
	set SF 1.0;
 
	# Create rack model
	rackMRF $h1 $hx $L $beamPrfl $Fyb $colPrfl $Fyc $boxed_levels $PL $nB $nL $zeta $muS $muK $sprOption $bplOption $kc $kctop $kbpl $NL_geo $d0 0 "zeta$zeta,muS$muS,muK$muK,acc$accRecord";
	
	# Apply gravity loads
	record_rack "basenodes" "reaction";
	rack_gravity_static [expr $PL/6] "on_beams"; # Apply gravity loading to include P-delta effects
	set message  "Total vertical reaction is : "
	set dof 2;
	set step 10;
	check_total_reaction $dof $step; # Sum the data in the last line and show result
	
	#DisplayModel2D
	set out_time ""
	record_rack "ext_col_nodes" "disp" $out_time;
	record_rack "palletnodes" "disp" $out_time;
	record_rack "beamnodes" "disp" $out_time;
	record_rack "columns"  "force" $out_time;
	record_rack "basenodes" "reaction" $out_time;
	
	rack_TH $dt $directory $accRecord $SF;
}

########################################################################
# Package of commands for building a model of a rack semi-rigid moment
# frame, for running analysis on that frame and recording results
#
# rack_MRF is the principal function that calls other functions in order to
# construct a model of a rack
# 
########################################################################

# Create the rack_pack package

package provide rack_pack 1.0;
package require Tcl 8.5;
model BasicBuilder -ndm 2 -ndf 3

namespace eval ::rack_pack:: {
	
	# Creates the rack_pack namespace with some common functions and variables	
	
	set version 1.0;
	
	# Export functions such that they can be called by scripts outside the rack_pack file
	namespace export rack_period;
	namespace export triangular_pattern;
	namespace export rack_gravity_static;
	namespace export add_rayleigh_damping
	namespace export rack_TH;
	namespace export rackMRF;
	namespace export check_total_reaction;
	namespace export record_rack;
	namespace export use_NL_spring;
	
	# Default output directory
	variable output_dir ""; # Set to null each time the package is called
	
	# See get_tag
	variable tagA 0; # Acceleration tag initialisation
	variable tagT 0; # Time series tag initialisation
	variable tagP 0; # Load pattern tag initialisation
	variable tagM 0; # Material tag initialisation
	variable tagS 0; # Section tag initialisation
	variable tagG 0; # Geometric tags initialisation
	variable tagF 0; # Friction tags initialisation
	
	# Categorized lists of elements and nodes to make adding recorders and loads to groups of nodes easier
	# Nodes are added to these lists when rackMRF is called and then they are available to all procedures in
	# rack_pack
	variable base_node_list {}; # List of nodes at the base of the structure
	variable beam_node_list {}; # List of beam nodes not including the one coincident with the column
	variable col_mid_node_list {}; # List of nodes in the middle of the columns (when fiber columns are used)
	variable bpl_elem_list {}; # List of base-plate connector elements
	variable col_elem_list {}; # List of column nodes
	variable extcol_node_list {}; # List of nodes at each level where exterior column nodes are attached to beam-to-column connectors
	variable intcol_node_list {}; # List of nodes at each level where interior column nodes are attached to beam-to-column connectors
	variable connector_elem_list {}; # List of beam-to-column connector elements
	variable beam_elem_list {}; # List of beam elements
	variable beam_end_list {}; # List of fiber beam ends
	variable pallet_elem_list {}; # List of elements between beams and pallet nodes
	variable pallet_node_list {}; # List of nodes to which pallet masses are assigned
	
	# Some modelling constants
	variable pi [expr 2.0*asin(1.0)]; # Definition of pi
	variable g 9.80665; # Gravity [m/s2]
	variable E [expr 200.0e9]; # Steel elastic modulus [N/m2]
	variable mSteel [expr 77000./$g]; # Unit mass of steel [kg/m3]
	variable Fy_bpl 345e6; # !!!!! Careful : Changing this value will change beam and column plasticity !!!!!
	
}

proc ::rack_pack::get_Tag { tagType } {
	
 #######################################################################
 # Creates time, acceleration and pattern tags such that the same tag is
 # not accidentally used twice by analysis procedures
 #
 # Example of creating a material tag and saving it to a variable:
 # set sprMat [get_Tag "Material"];
 # 
 # 11 july 2016, EJ
 #######################################################################
	
	variable tagA; 
	variable tagT;
	variable tagP;
	variable tagM;
	variable tagS;
	variable tagG;
	variable tagF;
	
	if { $tagType == "Section"} { incr tagS; return $tagS; };
	if { $tagType == "Accel" } { incr tagA; return $tagA; };
	if { $tagType == "Time" } { incr tagT; return $tagT; };
	if { $tagType == "Pattern" } { incr tagP; return $tagP; };
	if { $tagType == "Material" } { incr tagM; return $tagM; };
	if { $tagType == "Geo" } { incr tagG; return $tagG; };
	if { $tagType == "Fric" } { incr tagF; return $tagF; };
	
}

proc ::rack_pack::rackMRF { h1 hx L beamPrfl Fyb colPrfl Fyc boxed_levels PL nB nL zeta muSlow muFast sprOption bplOption { kc 1. } { kctop 1. } { kbpl 1. } {NL_geo "LargePdelta"} {g_tol 0.} {l_tol 0.} {output_dir_append ""} } {

	#################################################################################################################
	# Creates a 2D multi-story MRF with concentrated plasticity at column bases and beam to column joints, concentrated plasticity beams 
	# and distributed plasticity (fibre cross-section discretisation) of the columns allowing in-plane buckling
	# First story column height can be specified differently than other column heights
	# ALL UNITS SI
	#
	# INPUTS
	# h1 := First level height, double [m]
	# hx := Other level heights, double [m]
	# L := Length of beams, double [m]
	# beamPrfl := CISC designation for beam profile (other custom sections are available see proc secProp), string
	# Fyb := Elastic limit for beams, double [N/m2] (if a value less than Fy_bpl = 345 MPa is used, plastic beams will be used)
	# colPrfl := CISC designation for column profile (other custom sections are available see proc secProp), string
	# Fyc: = Elastic limit for columns, double [N/m2] (if a value less than Fy_bpl = 345 MPa is used, plastic columns will be used)
	# boxed_levels := number of boxed levels, integer [unitless]
	# PL := Weight of an individual pallet [N], the procedure assumes that there are two pallets per bay per level and that their weight is supported by 2 parallel MRF
	# nB := Number of bays, integer [unitless]
	# nL := Number of levels, interger [unitless]
	# zeta := percent of critical damping to be used, default is 0. (Rayleigh damping) [%]
	# muSlow := pallet friction coefficient [unitless], set to zero for no sliding
	# muFast := kinetic friction coefficient [unitless]
	# sprOption := "L" for linear springs, name of the non-linear spring if non-linear behaviour is wanted see use_NL_spring for some predefined options
	# bplOption := "L" for linear springs, name of the non-linear spring if non-linear behaviour is wanted see use_NL_spring for some predefined options
	# kc := When linear connections are used, kc represents the beam-column connector stiffness, double [Nm/rad] 
	# 		OR
	# 	 := When non-linear connections are used, kc represents the strength scaling factor to be applied to the non-linear material, double [unitless]
	# kctop := Same as kc, but for the top interior connectors which may need to be a different scaling for strong column capacity design 
	# kbpl := When linear connections are used, kbpl represents the base-plate connector stiffness, double [Nm/rad] 
	# 		OR
	# 	   := When non-linear connections are used, kbpl represents the a scaling factor to be applied to the non-linear material, double [unitless]
	#		  In the special case that fiber sections are used, kbpl represents the thickness of the plate, double [m]
	# NL_geo := defines the kind of geometric non-linearity to be applied (default is to include P-delta effects), string [unitless]
	# g_tol := initial global out of plomb condition (default is zero), double [mm/mm]
	# l_tol := initial local column cambrure (default is zero), double [mm/mm]
	# output_dir_append := a string to append to the default output directory
	#
	# Example :
	#
	# EJ
	#################################################################################################################
	
	# Remove existing models
	wipe; 
	
	# Sets the name of the directory to save results in
	set default_dir "rack_nL[expr $nL]nB[expr $nB]_"
	
	# Creates output directory in which to save results
	make_output_dir $default_dir$output_dir_append
	
	# Saves and outputs modelling info
	puts_and_save "h1: $h1"; puts_and_save "hx: $hx"; puts_and_save "L: $L"; puts_and_save "beamPrfl: $beamPrfl"; puts_and_save "Fyb: $Fyb";
	puts_and_save "colPrfl: $colPrfl"; puts_and_save "Fyc: $Fyc"; puts_and_save "boxed_levels: $boxed_levels";	puts_and_save "PL: $PL";	puts_and_save "nB: $nB";	puts_and_save "nL: $nL";
	puts_and_save "zeta: $zeta"; puts_and_save "muSlow: $muSlow"; puts_and_save "muFast: $muFast"; puts_and_save "sprOption: $sprOption"; puts_and_save "bplOption: $bplOption";	
	puts_and_save "kc : $kc"; puts_and_save "kctop : $kctop"; puts_and_save "kbpl: $kbpl"; puts_and_save "NL_geo: $NL_geo"; puts_and_save "g_tol: $g_tol"; puts_and_save "l_tol: $l_tol";
	
	# Loads some variables
	variable g;
	variable pi;
	variable mSteel;
	variable E;	
	variable Fy_bpl;
	variable output_dir;
	
	variable tagA; # Acceleration tag initialisation
	variable tagT; # Time series tag initialisation
	variable tagP; # Load pattern tag initialisation
	variable tagM; # Material tag initialisation
	variable tagS; # Section tag initialisation
	variable tagG; # Geometric tags initialisation
	variable tagF; # Friction ta initialisation
	
	# Calculates some parameters
	set mL [expr $PL/$g]; # Pallet mass
	set nc [expr ($nB+1)]; # Number of upright base plate connections
	set nb [expr $nL*$nc]; # Number of beam-connector connections
	
	# Define geometric transformations
	# For beam elements
	set BeamTransf [get_Tag "Geo"];
	geomTransf Linear $BeamTransf;
	# For column elements 
	set ColTransf [get_Tag "Geo"];
	if { $NL_geo == "noPdelta" } { geomTransf Linear $ColTransf };
	if { $NL_geo == "Pdelta" } { geomTransf PDelta $ColTransf };
	if { $NL_geo == "LargePdelta" } { geomTransf Corotational $ColTransf };

	# Create array "h" of level heights and an array of elevations "hn"
	set h(1) $h1;
	set hn(1) 0;
	set hn(2) $h1;
	for { set i 2 } { $i <= $nL+1 } { incr i } {
			set h($i) $hx;
			set hn([expr $i+1]) [expr $hn($i)+$hx];
	}

	# Categorized lists of elements and nodes to make adding recorders and loads to groups of nodes easier
	# Elements and nodes are appended to the appropriate lists as they are created 
	variable base_node_list {}; # List of nodes at the base of the structure
	variable beam_end_list {}; # List of the beam sub elements connected directly to connectors
	variable bpl_elem_list {}; # List of base-plate connector elements
	variable col_elem_list {}; # List of column nodes
	variable extcol_node_list {}; # List of nodes at each level where exterior column nodes are attached to beam-to-column connectors
	variable intcol_node_list {}; # List of nodes at each level where interior column nodes are attached to beam-to-column connectors
	variable connector_elem_list {}; # List of beam-to-column connector elements
	variable beam_elem_list {}; # List of beam elements
	variable pallet_elem_list {}; # List of elements between beams and pallet nodes
	variable pallet_node_list {}; # List of nodes to which pallet masses are assigned
	variable col_mid_node_list {}; # List of nodes in the middle of the columns (when fiber columns are used)
	variable beam_node_list {}; # List of nodes along each beam
	
	# Define nodes
	for { set i 1 } { $i <= [expr $nL+1] } { incr i } {
		for { set j 1 } { $j <= $nc } { incr j } {
			# Column nodes
			node $i$j [expr $L*($j-1) + $g_tol*$hn($i)]  $hn($i);
			# Column nodes are added to lists
			if { $i == 1} { lappend base_node_list $i$j }
			if { $i != 1} { 
				if { $j==1 || $j==$nc } { 
					lappend extcol_node_list $i$j;
				} else { 
					lappend intcol_node_list $i$j;
				} 
			}
			if {$i == 1} {	
				# do nothing
			} elseif { $j < $nc } {
				# Beam nodes
				# the array pos contains the relative positions of the beam divisions
				array set pos {
					1 0 
					2 0.03 
					3 0.225 
					4 0.48  
					5 0.52 
					6 0.745 
					7 0.97 
					8 1 
				}
				set nb_div [array size pos]
				for { set k 1 } { $k <= $nb_div } { incr k } {
					if { $k == 1 || $k == $nb_div } {
						# Beam nodes coincident to columns 
						node $i$j$k [expr $L*($j-1) + $g_tol*$hn($i) + $pos($k)*$L]  $hn($i) -mass 0.0 0.0 0;
					} else {
						# Beam nodes not coincident with columns (will be given coincident pallet nodes)
						set node_mass [expr $mL/($nb_div-2)];
						node $i$j$k [expr $L*($j-1) + $g_tol*$hn($i) + $pos($k)*$L]  $hn($i) -mass 1. 1. 0.;
						lappend beam_node_list $i$j$k; 
						# Pallet nodes are given each an equal share of pallet mass
						set l 1; # Note numbering scheme : If beam node is 212 then corresponding pallet node is 2121
						node $i$j$k$l [expr $L*($j-1) + $g_tol*$hn($i) + $pos($k)*$L]  $hn($i) -mass $node_mass $node_mass 0.;
						lappend pallet_node_list $i$j$k$l;
					}
				}
			}
		}
	}

	### Define Materials
	# Define connector material
	set sprMat [get_Tag "Material"];
	set sprMatTop [get_Tag "Material"];
	if { $sprOption == "L" } {
			use_linear_springs $sprMat $kc;
			use_linear_springs $sprMatTop $kctop;
	} else {
			use_NL_spring $sprMat $sprOption $kc;
			use_NL_spring $sprMatTop $sprOption $kctop;
	}
	## Define base-plate material
	set bplMat [get_Tag "Material"];
	if { $bplOption == "L" } {
			use_linear_springs $bplMat $kbpl;
	} else {
			use_NL_spring $bplMat $bplOption $kbpl; 
	}
	

	## Pallet Material
	## Friction Model 
	if { [expr $muSlow] > 0 } { 
		set kInit 100.e6;
		set frnTag [get_Tag "Fric"];
		frictionModel VelDependent $frnTag $muSlow $muFast 0.6
		set matPyTag [get_Tag "Material"];  uniaxialMaterial Elastic $matPyTag 100.e6;
	}	

	## Create column material and section if we are using fibre columns or concentrated plasticity columns
	if { (($l_tol != 0.) || ($Fyc < $Fy_bpl)) } {
		
		set colMat [get_Tag "Material"];
		
		# Kinematic hardening parameters
		set b_k 0.01; set R_0 7.5; set r_1 0.9; set r_2 0.15;
		# Ultimate strength parameters
		set f_u [expr $Fyc/$Fy_bpl*450e6]; set R_u 20;
		uniaxialMaterial Steel4 $colMat $Fyc $E -kin $b_k $R_0 $r_1 $r_2 -ult $f_u $R_u;
		set NELEM 8;
		
		set dc [secProp $colPrfl "D"];
		set bc [secProp $colPrfl "B"];
		set tc [secProp $colPrfl "T"];
		set wc [secProp $colPrfl "W"];
		set xc [secProp $colPrfl "X"];
		set xoc 0.; # Doesnt matter for planar
		set Ac [secProp $colPrfl "A"];
		
		set mDc [expr $mSteel*$Ac];  #Mass per unit length
		set mDboxed [expr $mSteel*2*$Ac];
		
		# Both boxed section and the channel section are made and then when making the individual elements the appropriate tag is used
		set columnChannelTag [make_fibre_channel_section $colMat $dc $bc $tc $wc $xc [expr -$xoc]]
		set boxTag [make_fibre_boxed_channel_section $colMat $dc $bc $tc $wc $xc [expr -$xoc]]
	}
	
	## Create beam material and section if we are using concentrated plasticity beams
	if { $Fyb < $Fy_bpl } {
		
		set beamMat [get_Tag "Material"];
		
		# Kinematic hardening parameters
		set b_k 0.01; set R_0 7.5; set r_1 0.9; set r_2 0.15;
		# Ultimate strength parameters
		set f_u [expr $Fyb/$Fy_bpl*450e6]; set R_u 20;
		uniaxialMaterial Steel4 $beamMat $Fyb $E -kin $b_k $R_0 $r_1 $r_2 -ult $f_u $R_u;
		
		set db [secProp $beamPrfl "D"];
		set bb [secProp $beamPrfl "B"];
		set tb [secProp $beamPrfl "T"];
		set wb [secProp $beamPrfl "W"];
		set xb [secProp $beamPrfl "X"];
		set xob 0.; # Doesnt matter for planar
		set Ab [secProp $beamPrfl "A"];
		
		set mDb [expr $mSteel*$Ab];  #Mass per unit length
		
		set beamChannelTag [make_fibre_channel_section $beamMat $db $bb $tb $wb $xb [expr -$xob]]
	}
	

	
	# Base-plate material
	variable Fy_bpl;
	if { $bplOption == "fiber"} {
	
		## Compression support elements
		set matSup [get_Tag "Material"];
		uniaxialMaterial ENT $matSup [expr $E*100];
	
		set d_col [secProp $colPrfl "D"]; # Column section height
		set wp [expr 8.5*25.4/1000]; # Base-plate width
		set a [expr 2*25.4/1000]; # Base-plate bolt-to-column distance
		set bp [expr 2*$a+$d_col]; # Bolt-to-bolt distance
		set tp $kbpl; # Base-plate thickness
		set nfh 1; set nfw 6;
		set plateTag [get_Tag "Section"]; 
		set plateMat [get_Tag "Material"];
		# Kinematic hardening parameters
		set b_k 0.01; set R_0 7.5; set r_1 0.9; set r_2 0.15;
		# Ultimate strength parameters
		set f_u 450e6; set R_u 20;
		uniaxialMaterial Steel4 $plateMat $Fy_bpl $E -kin $b_k $R_0 $r_1 $r_2 -ult $f_u $R_u;
		section fiberSec $plateTag {
		patch rect $plateMat $nfw $nfh [expr -$tp/2] [expr -$wp/2] [expr $tp/2] [expr $wp/2];
		}
	}

	# Define Elements
	# Columns: named according to their start node and end node
	for { set i 1 } { $i <= $nL } { incr i } {
		for { set j 1 } { $j <= $nc } { incr j } {
			set nodeI $i$j;
			set nodeJ [expr $i+1]$j;
			# Check to see if fibre or other columns should be used	
			if { $l_tol != 0. } {
				# Check to see if the column is a simple channel or a boxed channel	
				if { $i <= $boxed_levels } { 
						use_fiber_columns $nodeI $nodeJ $NELEM $boxTag $mDboxed $ColTransf $l_tol; # Create sub-elements and add them to list of column elements
					} else { 
						use_fiber_columns $nodeI $nodeJ $NELEM $columnChannelTag $mDc $ColTransf $l_tol; # Create sub-elements and add them to list of column elements
					}
			} else {
				if { $Fyc  <  $Fy_bpl} {
					# If Fyc is less than 345 MPa then columns are meant to be checked for plastification at their ends and a concentrated plasticity element is used
					# Check to see if the column is a simple channel or a boxed channel	
					set l_pl [expr 0.05*$hx]; # Length of plastic hinge
					if { $i <= $boxed_levels } {  
							lappend col_elem_list [use_concentrated_plasticity_beamcolumns $nodeI $nodeJ $colPrfl $ColTransf $l_pl $boxTag $mDboxed ];
						} else {
							lappend col_elem_list [use_concentrated_plasticity_beamcolumns $nodeI $nodeJ $colPrfl $ColTransf $l_pl $columnChannelTag $mDc ]; 
						}	
				} else {
					# If Fyc = 345 MPa then plastification is not being checked and elastic elements can be used
					# Check to see if the column is a simple channel or a boxed channel	
					if { $i <= $boxed_levels } { set boxed 2. } else { set boxed 1. }
					lappend col_elem_list [use_L_beamcolumns $nodeI $nodeJ $colPrfl $boxed $ColTransf]; # Create an elastic element and add it to list of column elements
				}
			}
		}
	}
	
	# Beams: named according to their start node and end node
	for { set i 1 } { $i <= $nL } { incr i } {
		for { set j 1 } { $j < $nc } { incr j } {
			for { set k 1 } { $k < $nb_div } { incr k } { 
				set nodeI [expr $i+1]$j$k;
				set nodeJ [expr $i+1]$j[expr $k+1];
				if { $k == 1 || $k == $nb_div } {
					# Beam elements that are attached to connectors can be elastic or plastic
					if { $Fyb < $Fy_bpl } {
						# If Fyb is less than 345 MPa then beams are meant to be checked for plastification at their ends and a fiber element is used between the column and the pallet
						lappend beam_end_list [use_fiber_beam_ends $nodeI $nodeJ $beamChannelTag $mDb $BeamTransf];
					} else {
						# If Fyb = 345 MPa then plastification is not being checked and elastic elements can be used
						lappend beam_elem_list [use_L_beamcolumns $nodeI $nodeJ $beamPrfl 1. $BeamTransf]; # Create element and add it to list of beam elements
					}
				} else {
					# Add beam elements that are not directly connected to connectors
					lappend beam_elem_list [use_L_beamcolumns $nodeI $nodeJ $beamPrfl 1. $BeamTransf]; # Create element and add it to list of beam elements
				}
				#Define pallet elements
				if { $k == 1 || $k == $nb_div } {
					#do nothing
				} else {
				
					set nodeK [expr $i+1]$j$k[expr 1]; # Pallet node
				
					if { [expr $muSlow] > 0 } { 
						set tolSlider 10.e-6;
						set maxIterSlider 200;
						element flatSliderBearing $nodeI$nodeK $nodeI $nodeK $frnTag $kInit -P $matPyTag -Mz $matPyTag -orient 0 1 0 -1 0 0 -iter $maxIterSlider $tolSlider;
						equalDOF $nodeI $nodeK 3;
						lappend pallet_elem_list $nodeI$nodeK; # Add element to list of pallet elements
					} else {
						# Fixes pallets in place 
						equalDOF $nodeI $nodeK 1 2 3; 
					}			
					
				}
			}
		}
	}

	## Define beam-column connectors
	# Lower level connectors
	for { set i 2 } { $i <= $nL } { incr i } {
		for { set j 1 } { $j < $nc } { incr j } {
			# Left most connector in the bay at level i
			set k 1;
			set nodeI $i$j;
			set nodeJ $i$j$k;
			rotSpring2D $nodeI$nodeJ $nodeI $nodeJ $sprMat;
			lappend connector_elem_list $nodeI$nodeJ; # Add element to list of beam-column connectors
			# Right most connector in the bay at level i
			set k $nb_div;
			set nodeI $i[expr $j+1];
			set nodeJ $i$j$k;
			rotSpring2D $nodeI$nodeJ $nodeI $nodeJ $sprMat;
			lappend connector_elem_list $nodeI$nodeJ; # Add element to list of beam-column connectors
		}
	}
	
	## Top level connectors (for weak connector strong column design, inner top columns must have connectors half as strong)
	for { set j 1 } { $j < $nc } { incr j } {
		if { $j == 1 || $j == $nc } { set Mat $sprMat; } else { set Mat $sprMatTop };
		# Left most connector in the bay at level i
		set i [expr $nL+1]
		set k 1;
		set nodeI $i$j;
		set nodeJ $i$j$k;
		rotSpring2D $nodeI$nodeJ $nodeI $nodeJ $Mat;
		lappend connector_elem_list $nodeI$nodeJ; # Add element to list of beam-column connectors
		# Right most connector in the bay at level i
		set k $nb_div;
		set nodeI $i[expr $j+1];
		set nodeJ $i$j$k;
		rotSpring2D $nodeI$nodeJ $nodeI $nodeJ $Mat;
		lappend connector_elem_list $nodeI$nodeJ; # Add element to list of beam-column connectors
	}
	
	# Define column base nodes and connectors fixities
	if { $bplOption == "fiber" } {
		variable base_node_list {};	
		for { set j 1 } { $j <= $nc } {incr j} {
			set nodeJ 1$j;		
			set c [expr ($j-1)*$L]; # Centre de l'appui	
			use_fiber_bpl $nodeJ $c $bp $a $tp $wp $plateTag $matSup;
		}
	} else {
		for { set j 1 } { $j <= $nc } {incr j} {
			set k 0;
			node $j$k [expr [expr $j-1]*$L] 0;
			set nodeI $j$k;
			set nodeJ 1$j;
			fix $nodeI 1  1  1; # Assign column base fixities
			rotSpring2D $nodeI$nodeJ $nodeI $nodeJ $bplMat;
			lappend bpl_elem_list $nodeI$nodeJ; # Add element to list of base-plate connectors
		}
	}
	
	# Add viscous damping
	add_rayleigh_damping $zeta;
	
	puts "Model created";
	
	# Outputs model information about nodes and elements to the output folder
	print "$output_dir/nodes.out" -node
	print "$output_dir/elements.out" -ele
	
}

proc ::rack_pack::rack_period { {n_modes 1} } {
	
	####################################################################
	# Calculation of 1st fundamental period
	# Prints period in the terminal and returns the period
	#
	# 8 mai 2016, EJ
	####################################################################

	set pi [expr 2.0*asin(1.0)]; # Definition of pi

	set omega [expr pow([eigen $n_modes], 0.5)];
	set Tn [expr double(round([expr 20000 * $pi / $omega]))/10000];
	puts "T1 = $Tn s";
	return $Tn;

}

proc ::rack_pack::rack_gravity_static { Py load_placement } {

	################################################################################################################
	# Adds a point load of Py either on beams or on columns
	# When applied to columns, 1*Py is applied to exterior column nodes and 2*Py is applied to interior columns
	# When applied to beam, 1*Py is applied to each node
	# Then a static analysis is run, loads are held constant and time is reset to zero
	#
	# INPUTS :
	# Py := scalar pallet weight in [N]
	# load_placement := PL will either be placed "on_beams" "on_cols" "beam_left" "beams right"
	#
	# These two examples would give equivalent vertical reactions at base nodes :
	# Example : rack_gravity_static $PL/6 "on_beams"; (assuming the default of six pallet nodes per beam)
	# Example : rack_gravity_static $PL/2 "on_cols";
	#
	# Last mod. : 5 mai 2016, EJ
	################################################################################################################
	
	# Lists of nodes constructed in rackMRF are included in the procedure, theses variables are part of the rack_pack 
	# namespace
	variable right_pallet_node_list;
	variable left_pallet_node_list;
	variable pallet_node_list;
	variable extcol_node_list;
	variable intcol_node_list;
	variable base_node_list;
	
	# Pattern and time series tags are choosen automatically
	set gravityPattern [ get_Tag "Pattern" ];
	set timeSeriesTag [ get_Tag "Time" ];
	timeSeries Linear $timeSeriesTag;
	
	# Depending on the string in load_placement, loads are added to the appropriate nodes
	pattern Plain $gravityPattern $timeSeriesTag {
		if { $load_placement == "on_beams" } { 
			for { set n 0 } { $n < [llength $pallet_node_list] } {incr n} {
				set node [lindex $pallet_node_list $n];
				load $node 0.0 [expr -$Py] 0.0;
				
			}
		}
		if { $load_placement == "on_cols" } { 
			if { [info exists extcol_node_list] == 1} {
				for { set n 0 } { $n < [llength $extcol_node_list] } {incr n} {
					set node [lindex $extcol_node_list $n];
					load $node 0.0 [expr -$Py] 0.0;
				} 
			}
			if { [info exists intcol_node_list] == 1} {
				for { set n 0 } { $n < [llength $intcol_node_list] } {incr n} {
					set node [lindex $intcol_node_list $n];
					load $node 0.0 [expr -$Py*2] 0.0;
				} 
			}
		}
	}

	constraints Plain; # Defines how to handle constraints between DOFs
	numberer RCM; # Specify how OpenSees numbers the DOFs	
	system BandGeneral; # Save the system of equations	

	set tol 1.0e-6; # Convergence tolerance for test
	set iter 6; # The max number of iterations to check before returning failure condition
	test NormDispIncr $tol $iter; # Determine if convergence has been achieved at the end of an iteration step

	# Algo for the resolution of the system of non-linear equations	
	algorithm Newton; # Use Newton's solution algorithm: updates tangent stiffness at every iteration

	# integrator LoadControl $lambda <$numIter $minLambda $maxLambda>
	set Nsteps 10; # Number of steps in which to apply load
	set lambda [expr 1.0/$Nsteps]; # load increment
	integrator LoadControl $lambda;
		
	#Construct the analysis object
	analysis Static;
	
	# Run analysis in Nsteps
	set ok [analyze $Nsteps]; # this will return zero if no convergence problems were encountered
	set time 0.;
	set message "Applied gravity loads... setting time to ";
	if { $ok == 0 } { puts [append message $time " and holding load constant."] };

	loadConst -time $time; # Sets existing loading to be constant and resets time to 0.0
	
}

proc ::rack_pack::rack_TH { dt path accRecord SF {NTINCR 0}} {

	#################################################################################################################
	# Time-history analysis
	#
	# Warning : call graivty load analysis before to include P-delta effects
	# dt := time step in file [s]
	# path := full path to directory of acceleration file
	# accRecord := string, file name of the accelerogram used
	# NTINCR := number (default is the number of lines in the acceleration file), number of steps to be used in analysis
	#
	# Last mod. : 14 mars 2016, EJ
	#################################################################################################################

	# Dynamic Loading
	# ----------------
	variable g;
	set accelSeries "Series -dt $dt -filePath $path$accRecord -factor [expr $g*$SF]";
	set patternTag [get_Tag "Pattern"];
	set dirX 1;
	pattern UniformExcitation $patternTag $dirX -accel $accelSeries;

	# Newmark Algorithm gamma=0.25; beta=0.5
	# ---------------------------------------
	set tol 0.000001;
	set maxNumIter 300;
	test EnergyIncr $tol $maxNumIter;
	constraints Plain;
	numberer Plain;
	system BandGeneral;
	algorithm Newton;
	#integrator Newmark 0.5 0.25;
	analysis Transient;
	
	integrator HHT 0.9 0.3025 0.6; # Use if there are convergence problems with pallet sliding

	# Check if NTINCR has been specified otherwise get number of lines in file
	if {$NTINCR == 0} {
		set infile [open $path$accRecord r];

		while { [gets $infile line] >= 0 } {
		incr NTINCR;
		} 
		close $infile;
	}

	set dursism [expr $NTINCR*$dt];

	# time step size for the analysis
	set dt2 0.001;
	puts_and_save "Time history time step:$dt2"
	puts_and_save  "Starting Time History Analysis with $accRecord at scale factor $SF";
	set startTime [clock clicks -milliseconds]
	set ok [analyze [expr int($dursism/$dt2)] $dt2];
	set endTime [clock clicks -milliseconds]
	set timeElapsed [expr ($endTime - $startTime)/1000]
	if { $ok == 0 } { puts_and_save "TH Analysis Completed in $timeElapsed s" } else {
		puts_and_save "TH Analysis failed after $timeElapsed s"
	};
	
}

proc ::rack_pack::add_rayleigh_damping { zeta } {
	
	variable pi;
	variable col_elem_list; # List of column nodes
	variable beam_elem_list; # List of beam elements
	variable pallet_node_list; # List of nodes to which pallet masses are assigned
	
	## Rayleigh Damping
	set lambdaN [eigen [expr 2]];# eigenvalue analysis for nEigenJ modes
	set w1 [expr pow([lindex $lambdaN 0],0.5)]; # w1 (1st mode circular frequency)
	set w2 [expr pow([lindex $lambdaN 1],0.5)]; # w2 (2nd mode circular frequency)
	set T1 [expr double(round([expr 2000*$pi/$w1]))/1000];
	set T2 [expr double(round([expr 2000*$pi/$w2]))/1000];
	puts_and_save  "T1 = $T1 s";
	puts_and_save  "T2 = $T2 s";
	# Calculate damping parameters (see Chopra p.419)
	set a0 [expr $zeta*2.0*$w1*$w2/($w1 + $w2)]; # mass damping coefficient based on first and second modes
	set a1 [expr $zeta*2.0/($w1 + $w2)]; # stiffness damping coefficient based on first and second modes
	
		# Assign element and node lists to regions and apply damping
	# Calling rayleigh damping : -rayleigh $alphaM $betaK $betaKinit $betaKcomm;
	set cmd "region 1 -ele $col_elem_list -rayleigh $a0 0. $a1 0.";
	eval $cmd; # Because @!&% tcl doesnt evaluate $col_elem_list as a string when I call it normally
	set cmd "region 2 -ele $beam_elem_list -rayleigh $a0 0. $a1 0.";
	eval $cmd;
	set cmd "region 3 -node $pallet_node_list -rayleigh $a0 0. 0. 0."; # Nodes to which pallet masses are assigned
	eval $cmd;
}

proc ::rack_pack::use_linear_springs { sprMat kc } {

	#################################################################################################################
	# Defines material moment-rotation curve to be used for linear zero-length connectors
	#
	# INPUTS:
	# sprMat := material tag [an integer]
	# kc := rotational stiffness [Nm/rad]
	#
	# Last mod. : 15 Juillet 2015, EJ
	#################################################################################################################
	
	puts_and_save  "Using Linear Spring : $kc";
	
	uniaxialMaterial Elastic $sprMat $kc; # Spring is perfectly elastic
	
}

proc ::rack_pack::use_NL_spring { sprMat NLsprName SF } {

	#################################################################################################################
	# Defines material moment rotation curve to be used for non-linear zero-length connectors
	# 
	# INPUTS:
	# sprMat := material tag [an integer]
	# NLsprName := the name of the non-linear connector [text], see the list "Connector_name"
	# 
	# exemple : use_NL_spring 3 "AC175-7"
	#
	# Last mod. : 10 mars 2016, EJ
	#################################################################################################################
	
	set material_picked 0;
	
	if { $NLsprName == "AC175-7" } {
		set ePd1 0.0036;
		set ePf1 [expr $SF*900.];
		set ePd2 0.025; 
		set ePf2 [expr $SF*1120.]; 
		set ePd3 0.1;
		set ePf3 [expr $SF*4700.]; 
		set ePd4 0.4; 
		set ePf4 0.; 
		set rDispP 0.5; 
		set rForceP 0.04; 
		set uForceP 0.07; 
		set eNd1 [expr -$ePd1]; 
		set eNf1 [expr -$ePf1]; 
		set eNd2 [expr -$ePd2];
		set eNf2 [expr -$ePf2]; 
		set eNd3 [expr -$ePd3]; 
		set eNf3 [expr -$ePf3];
		set eNd4 [expr -$ePd4]; 
		set eNf4 [expr -$ePf4]; 
		set rDispN $rDispP;
		set rForceN $rForceP;
		set uForceN $uForceP;
		set gK1 0.; 
		set gK2 0.; 
		set gK3 0.; 
		set gK4 0.; 
		set gKLim -0.1; 
		set gD1 0.; 
		set gD2 0.1;
		set gD3 0.; 
		set gD4 0.1; 
		set gDLim 0.1; 
		set gF1 0.;
		set gF2 0.;
		set gF3 0.;
		set gF4 0.;
		set gFLim 0.01; 
		uniaxialMaterial Pinching4 $sprMat $ePf1 $ePd1 $ePf2 $ePd2 $ePf3 $ePd3 $ePf4 $ePd4 $eNf1 $eNd1 $eNf2 $eNd2 $eNf3 $eNd3 $eNf4 $eNd4 $rDispP $rForceP $uForceP $rDispN $rForceN $uForceN $gK1 $gK2 $gK3 $gK4 $gKLim $gD1 $gD2 $gD3 $gD4 $gDLim $gF1 $gF2 $gF3 $gF4 $gFLim 1. "cycle"; 
		set material_picked 1;
	}	
	if { $NLsprName == "RBS44" } {
		
		# Modification du 18 fevrier
		
		set	ePd1 0.013;
		set	ePf1 [expr $SF*1500.];
		set	ePd2 0.045;
		set	ePf2 [expr $SF*2800.];
		set	ePd3 0.15;
		set	ePf3 [expr $SF*4200.];
		set	ePd4 0.45;
		set	ePf4 0;
		set	rDispP 0.35;
		set	rForceP 0.04;
		set	uForceP 0.04;

		set	eNd1 [expr -1*$ePd1];
		set	eNf1 [expr -1*$ePf1];
		set	eNd2 [expr -1*$ePd2];
		set	eNf2 [expr -1*$ePf2];
		set	eNd3 [expr -1*$ePd3];
		set	eNf3 [expr -1*$ePf3];
		set	eNd4 [expr -1*$ePd4];
		set	eNf4 [expr -1*$ePf4];
		set	rDispN 0.4;
		set	rForceN 0.02;
		set	uForceN 0.02;
				
		set	gK1 0.;
		set	gK2 0.;
		set	gK3	0.;
		set	gK4	0.;
		set	gKLim 0.05;
					
		set	gD1	0.;
		set	gD2	0.;
		set	gD3	0.;
		set	gD4	0.;
		set	gDLim 0.3;
					
		set	gF1	0.2;
		set	gF2	0.2;
		set	gF3	0.2;
		set	gF4	0.2;
		set	gFLim 0.2;

		uniaxialMaterial Pinching4 $sprMat $ePf1 $ePd1 $ePf2 $ePd2 $ePf3 $ePd3 $ePf4 $ePd4 $eNf1 $eNd1 $eNf2 $eNd2 $eNf3 $eNd3 $eNf4 $eNd4 $rDispP $rForceP $uForceP $rDispN $rForceN $uForceN $gK1 $gK2 $gK3 $gK4 $gKLim $gD1 $gD2 $gD3 $gD4 $gDLim $gF1 $gF2 $gF3 $gF4 $gFLim 1. "cycle"; 
		set material_picked 1;
	}
	if { $NLsprName == "fiber" } {
		# Fiber base-plates are largely constructed in the rackMRF procedure
		set material_picked 1;
	}
	if { $NLsprName == "Parallel" } {
		# This is the material that best replicates the behavoir of the base-plate in the F1 test
		# Hysteretic + Parallel Steel 4 Test
		set s1p [expr $SF*0.3*1000]
		set e1p 0.005
		set s2p [expr $SF*0.8*1000]
		set e2p 0.02
		set s3p [expr $SF*1.1*1000]
		set e3p 0.08
		set s1n [expr $SF*-0.3*1000]
		set e1n -0.005
		set s2n [expr $SF*-0.5*1000]
		set e2n -0.02 
		set s3n [expr $SF*-1.0*1000]
		set e3n -0.08
		set pinchX 0.8
		set pinchY 0.2
		set damage1 0.
		set damage2 0.
		set beta 0.01
		
		set hystMat [get_Tag "Material"];
		set stl4_1Mat [get_Tag "Material"];
		set stl4_2Mat [get_Tag "Material"];
		
		uniaxialMaterial Hysteretic $hystMat $s1p $e1p $s2p $e2p $s3p $e3p $s1n $e1n $s2n $e2n $s3n $e3n $pinchX $pinchY $damage1 $damage2 $beta;
		uniaxialMaterial Steel4 $stl4_1Mat [expr $SF*0.4*1000] [expr $SF*5.0*1000] -kin 0.08 15.0 0.925 0.15 -ult [expr 1.0*1000] 20.; # Initial portion of curve
		uniaxialMaterial Steel4 $stl4_2Mat [expr $SF*1.0*1000] [expr $SF*45.0*1000] -kin 0.01 15.0 0.925 0.15 -ult [expr 2.2*1000] 20.; # 2nd portion of curve
		uniaxialMaterial Parallel $sprMat $hystMat $stl4_1Mat $stl4_2Mat;
		
		set material_picked 1;
	}
	if { $NLsprName == "MinMax" } {
		set elMat [get_Tag "Material"];
		uniaxialMaterial Elastic $elMat [expr $SF*10.e3];
		uniaxialMaterial MinMax $sprMat $elMat -min -0.01 -max 0.01;
		set material_picked 1;
	}
	if { $material_picked == 0 } {
		#parse NLsprName to extract custom material string
		set EOL [llength $NLsprName]; #Last element in the list
		set mat_args [lrange $NLsprName 1 $EOL]; # isolate the
		set NLsprName [lindex $NLsprName 0]
		set cmd "uniaxialMaterial $NLsprName $sprMat $mat_args"
		eval $cmd
	}
}

proc ::rack_pack::use_fiber_bpl { originNode c bp a tp wp plateTag matSup } {

	#######################################################################################
	#  	Makes a fibre base-plate assembly as seen below:
	# 	-the node marked "o" is the originNode and is the point of attachement between the column and the base-plate assembly
	# 	-the nodes marked "#" or "bolt#" are fixed nodes
	#	-the elements marked "fibre" are non-linear beam column elements with a fiber section (defined in rackMRF)
	# 	-the elements marked "comp." are compression only zero-length elements
	# 	-the elements marked "el." are very stiff elastic elements which allow a distribution of compression and moment from the column to the base-plate
	#
	#								     |	
	# 	  y                             col.
	#	  |                             |
	# 	  ----x      		x----el.----o----el.----x                         
	#					    = 			            =                        
	#	                                       					   
	#	        	       =	     		          =       					
	#	 bolt#----fibre----x---------fibre--------x----fibre----#bolt       
	#	              	  |                       |
	#				    comp.		   			comp.
	#				     |			             |
	#	        	     #     		             #
	#	INPUTS:
	# 	originNode := the Tag of the node at the start of the adjoining column
	# 	c := the x coordinate of the originNode
	#	bp := the x-distance from originNode to the bolted end of the base-plate
	#	a := the x-distance between the bolt and the column section
	#	tp := the base-plate thickness
	# 	wp : the base-plate width	
	#	plateTag := the section tag of the fibre base-plate section (as defined in rackMRF)
	#	matSup := the material tag of the compression only material (as defined in rackMRF)
	#                          
	#    
	#	OUTPUTS:
	#	The procedure has no outputs HOWEVER the base_node_list and the bpl_elem_list are constructed here.
	#	Since the interesting outputs are the reactions at the base and the uplift of the plate:
	#	- the base_node_list contains ONLY the four fixed nodes of the base-plate assembly
	#	- the bpl_elem_list contains ONLY the two zero-length compression springs 
	#
	#	EJ
	########################################################################################
	
	set hbpl 0.;
	variable E;
	variable base_node_list;
	variable bpl_elem_list;

	# Nodes along the fiber plate element
	set node1 $originNode[expr 0][expr 0][expr 1];
	node $node1 [expr $c-$bp/2] [expr -$hbpl];
	node $originNode[expr 0][expr 0][expr 2] [expr $c-$bp/2+$a] [expr -$hbpl];
	node $originNode[expr 0][expr 0][expr 4] [expr $c+$bp/2-$a] [expr -$hbpl];
	set node5 $originNode[expr 0][expr 0][expr 5];
	node $node5 [expr $c+$bp/2] [expr -$hbpl];
	# End-nodes for tension only springs
	set node6 $originNode[expr 0][expr 0][expr 6];
	node $node6 [expr $c-$bp/2+$a] [expr -$hbpl];
	set node7 $originNode[expr 0][expr 0][expr 7];
	node $node7 [expr $c+$bp/2-$a] [expr -$hbpl];
	lappend base_node_list $node1 $node5 $node6 $node7;
	# Nodes for rigid "cadre" elements
	node $originNode[expr 0][expr 0][expr 8] [expr $c-$bp/2+$a] 0.;
	node $originNode[expr 0][expr 0][expr 9] [expr $c+$bp/2-$a] 0.;

	fix $originNode[expr 0][expr 0][expr 1] 1 1 1;
	fix $originNode[expr 0][expr 0][expr 5] 1 1 1;  
	fix $originNode[expr 0][expr 0][expr 6] 1 1 1;
	fix $originNode[expr 0][expr 0][expr 7] 1 1 1;

	equalDOF $originNode[expr 0][expr 0][expr 2] $originNode[expr 0][expr 0][expr 8] 1 2 3;
	equalDOF $originNode[expr 0][expr 0][expr 4] $originNode[expr 0][expr 0][expr 9] 1 2 3;
	
	# Cadre elements
	set A 0.00102; set Iz 1.61e-6; set Ecadre [expr $E];
	#element elasticBeamColumn $originNode[expr 0][expr 0][expr 4]$originNode[expr 0][expr 0][expr 9] $originNode[expr 0][expr 0][expr 4] $originNode[expr 0][expr 0][expr 9] $A $Ecadre $Iz 1
	#element elasticBeamColumn $originNode[expr 0][expr 0][expr 2]$originNode[expr 0][expr 0][expr 8] $originNode[expr 0][expr 0][expr 2] $originNode[expr 0][expr 0][expr 8] $A $Ecadre $Iz 1
	set Ar [expr $A]; set Izr [expr $Iz]; set Erigid [expr $Ecadre];
	element elasticBeamColumn $originNode[expr 0][expr 0][expr 8]$originNode $originNode[expr 0][expr 0][expr 8] $originNode $Ar $Erigid $Izr 1
	element elasticBeamColumn $originNode[expr 0][expr 0][expr 9]$originNode $originNode[expr 0][expr 0][expr 9] $originNode $Ar $Erigid $Izr 1

	# Compression only support elements
	set eleLeft $originNode[expr 0][expr 0][expr 2]$originNode[expr 0][expr 0][expr 6];
	set eleRight $originNode[expr 0][expr 0][expr 4]$originNode[expr 0][expr 0][expr 7];
	element zeroLength $eleLeft $originNode[expr 0][expr 0][expr 6] $originNode[expr 0][expr 0][expr 2] -mat $matSup -dir 1 -orient 0 1 0 -1 0 0;
	element zeroLength $eleRight $originNode[expr 0][expr 0][expr 7] $originNode[expr 0][expr 0][expr 4] -mat $matSup -dir 1 -orient 0 1 0 -1 0 0;
	lappend bpl_elem_list $eleLeft $eleRight;
	
	# PlateElements
	set npts 4;
	set maxIters 100;
	set tol 0.000001;
	set ele2 $originNode[expr 0][expr 0][expr 1]$originNode[expr 0][expr 0][expr 2];
	set ele4 $originNode[expr 0][expr 0][expr 2]$originNode[expr 0][expr 0][expr 4];
	set ele5 $originNode[expr 0][expr 0][expr 4]$originNode[expr 0][expr 0][expr 5];
	element nonlinearBeamColumn $ele2 $originNode[expr 0][expr 0][expr 1] $originNode[expr 0][expr 0][expr 2] $npts $plateTag 1;# -iter $maxIters $tol;
	element nonlinearBeamColumn $ele4 $originNode[expr 0][expr 0][expr 2] $originNode[expr 0][expr 0][expr 4] $npts $plateTag 1;# -iter $maxIters $tol;
	element nonlinearBeamColumn $ele5 $originNode[expr 0][expr 0][expr 4] $originNode[expr 0][expr 0][expr 5] $npts $plateTag 1;# -iter $maxIters $tol;
		
}

proc ::rack_pack::use_concentrated_plasticity_beamcolumns { s_node e_node elePrfl geoTransf l_pl secTag mD } {

	####################################################################
	# This procedure is called by rackMRF in order to assign concentrated plasticity elements between 2 nodes
	#
	# INPUTS:
	# s_node := the starting local column node, integer (defined in rackMRF)  [unitless]
	# e_node := the ending column node, integer (defined in rackMRF)  [unitless]
	# elePrfl := the name of the profile, string  [unitless]
	# geoTransf := geometric transformation tag, integer [unitless]
	# l_pl := the length of the plastic hinge, double [m]
	# secTag := the tag for the plastic hinge section, integer [unitless]
	# mD := mass per unit length, double [kg/m2]
	#
	# OUTPUT:
	# $s_node$e_node := the tag of the element created
	#
	# Last mod.: 8 aout 2016, EJ
	####################################################################

	variable mSteel;
	
	element forceBeamColumn $s_node$e_node $s_node $e_node $geoTransf "HingeRadau $secTag $l_pl $secTag $l_pl $secTag" -mass $mD; # <-iter $maxIters $tol>

	return $s_node$e_node
	
}

proc ::rack_pack::use_L_beamcolumns { s_node e_node elePrfl boxed geoTransf } {
	
	####################################################################
	# This procedure is called by rackMRF in order to assign elastic columns
	# between 2 nodes
	#
	# INPUTS:
	# s_node := the starting local column node, integer (defined in rackMRF)  [unitless]
	# e_node := the ending column node, integer (defined in rackMRF)  [unitless]
	# elePrfl := the name of the profile, string  [unitless]
	# boxed := used to double the inertia and area of a column profile when it is a boxed column, double [unitless]
	# geoTransf := geometric transformation tag, integer [unitless]
	#
	# OUTPUT:
	# $s_node$e_node := the tag of the element created
	#
	# Last mod.: 16 avril 2016, EJ
	####################################################################
	
	variable mSteel;
	variable E;	
	
	set A [expr $boxed*[secProp $elePrfl "A"]];
	set Ix [expr $boxed*[secProp $elePrfl "Ix"]];
	set mDc [expr $mSteel*$A];  #Mass per unit length
	
	element elasticBeamColumn $s_node$e_node $s_node $e_node $A $E $Ix $geoTransf -mass $mDc;
	
	return $s_node$e_node
}

proc ::rack_pack::use_fiber_beam_ends { s_node e_node beamSecTag mDb BeamTransf } {
	
	####################################################################
	# This procedure is called by rackMRF in order to perform 
	# discretisation of beam end elements into multiple non-linear elements with
	# fiber cross sections
	#
	# INPUTS:
	# s_node := the starting local column node, integer (defined in rackMRF)  [unitless]
	# e_node := the ending column node, integer (defined in rackMRF)  [unitless]
	# beamSecTag := the fiber section specified in rackMRF
	# mDb : = mass per unit length of beam member
	# BeamTransf := geometric transformation tag, integer [unitless]
	#
	# OUTPUT:
	# beam_list := a list of the elements created
	#
	# EJ
	####################################################################
	
	set npoint 3; # must be integer
	set maxIters 1000;
	set tol 0.000001;
	
	# Creation of sub-elements
	set beam_list {}; 
	set nodeI $s_node;
	set nodeJ $e_node;
	element nonlinearBeamColumn $nodeI$nodeJ $nodeI $nodeJ $npoint $beamSecTag $BeamTransf -mass $mDb -iter $maxIters $tol;
	lappend beam_list $nodeI$nodeJ;	
	
	# A list of column elements created is returned in order to add Rayleigh damping to the region in rackMRF
	return $beam_list
}

proc ::rack_pack::use_fiber_columns { s_node e_node NELEM colSecTag mDc ColTransf l_tol } {
	
	####################################################################
	# This procedure is called by rackMRF in order to perform 
	# discretisation of columns into multiple non-linear elements with
	# fiber cross sections
	#
	# INPUTS:
	# s_node := the starting local column node, integer (defined in rackMRF)  [unitless]
	# e_node := the ending column node, integer (defined in rackMRF)  [unitless]
	# NELEM := the number of sub elements, double
	# colPrfl := the name of the profile, string  [unitless]
	# ColTransf := geometric transformation tag, integer [unitless]
	# l_tol := initial local column cambrure, double [m/m]
	#
	# OUTPUT:
	# col_list := a list of the elements created
	#
	# Last mod.: 16 avril 2016, EJ
	####################################################################
	
	# Give a message while creating the first non-linear column
	if { $s_node == "11" } { puts_and_save  "Using non-linear columns with $NELEM sub-elements" };
	
	variable pi;
	variable g;
	variable col_mid_node_list;
	variable col_elem_list;
	
	set npoint 3; # must be integer
	# set maxIters 10;
	# set tol 0.000001;
	
	# Retrieve start and end coordinates of inputed nodes
	set sx [nodeCoord $s_node 1]; set sy [nodeCoord $s_node 2];
	set ex [nodeCoord $e_node 1]; set ey [nodeCoord $e_node 2];	
	set L [expr sqrt([expr pow([expr $sx-$ex],2)+pow([expr $sy-$ey],2)])]; # Distance between column end point and start point
	
	# Angle of the column with respect to the positive x-axis
	set sinT [expr ($sy-$ey)/$L];
	set cosT [expr ($ex-$sx)/$L];
	
	# Creation of nodes
	# A sine curve is created with amplitude l_tol*L and length L and then is rotated and translated to fit the column end-points
	set x(0) 0; set y(0) 0; set x1(0) $sx; set y1(0) $sy;
	for { set i 1 } { $i <= [expr $NELEM-1] } { incr i } {
		if { [expr [expr $s_node/10]%2] != 0 } {
			# Impair levels
			# Sine curve with length L and amplitude l_tol*L
			set x($i) [expr $i*$L/$NELEM];
			set y($i) [expr $l_tol*$L*sin([expr $pi*($x($i)-0)/$L])];
			# Rotation and translation of coordinates
			set x1($i) [expr $x($i)*$cosT+$y($i)*$sinT + $sx]; 
			set y1($i) [expr -$x($i)*$sinT+$y($i)*$cosT + $sy];	
		} else {
			# Pair levels
			# Sine curve with length L and amplitude l_tol*L
			set x($i) [expr $i*$L/$NELEM];
			set y($i) [expr -1*$l_tol*$L*sin([expr $pi*($x($i)-0)/$L])];
			# Rotation and translation of coordinates
			set x1($i) [expr $x($i)*$cosT+$y($i)*$sinT + $sx]; 
			set y1($i) [expr -$x($i)*$sinT+$y($i)*$cosT + $sy];
		}
		node $s_node[expr 0]$i $x1($i) $y1($i); # Add node to model
		if { $i == 4 } { lappend col_mid_node_list $s_node[expr 0]$i }
	}

	# Create first sub element
	set nodeI $s_node;
	set nodeJ $s_node[expr 0][expr 1];
	element nonlinearBeamColumn $nodeI$nodeJ $nodeI $nodeJ $npoint $colSecTag $ColTransf -mass $mDc;# -iter $maxIters $tol;
	lappend col_elem_list $nodeI$nodeJ;
	# Create intermediate sub-elements	
	for { set i 1 } { $i <= [expr $NELEM - 2] } { incr i } {
		set nodeI $s_node[expr 0]$i;
		set nodeJ $s_node[expr 0][expr $i+1];
		element nonlinearBeamColumn $nodeI$nodeJ $nodeI $nodeJ $npoint $colSecTag $ColTransf -mass $mDc;# -iter $maxIters $tol;
		if {$i == 4} { lappend col_elem_list $nodeI$nodeJ; } 	
	}
	
	# Create last sub element
	set nodeI $s_node[expr 0][expr int($NELEM-1)];
	set nodeJ $e_node;
	element nonlinearBeamColumn $nodeI$nodeJ $nodeI $nodeJ $npoint $colSecTag $ColTransf -mass $mDc;# -iter $maxIters $tol;
	lappend col_elem_list $nodeI$nodeJ;	
	
}

proc ::rack_pack::make_fibre_channel_section { matTag d b t w x xo } {
	
	####################################################################
	# This procedure is called by rack_MRF in order to make channel
	# shaped fibre cross-sections
	#
	# INPUTS:
	# (for more detail see tables of CISC Handbook)
	# d := depth, double [m]
	# b := flange width, double [m]
	# t := flange thickness, double [m]
	# w := we thickness, double [m]
	# x := distance from neutral axis to exterior web, double [m]
	# xo := position of loading vis-a-vis neutral axis, double (!!! be careful of sign!) [m]
	#
	# Example : Channel C100x7 section loaded in its shear centre : make_fibre_channel_section 0.102 0.040 0.0075 0.0032 0.0126 -0.0273 
	#
	# OUTPUT : secTag := an identifier for the section, integer [unitless]
	#
	# Last mod.: 16 avril 2016, EJ
	####################################################################

	set secTag [get_Tag "Section"];
	set aggSecTag [get_Tag "Section"];
	
	set nfd 4;	# number of fibers along web depth 
	set nfw 1;	# number of fibers along web thickness
	set nfb 1;	# number of fibers along flange width (you want this many in a bi-directional loading)
	set nft 8;	# number of fibers along flange thickness
	
	# Some coordinates delineating channel web and flange
	set z1 [expr ($d-2*$t)/2]; # distance to nearest flange fibre from the c.g.
	set z2 [expr $d/2]; # distance to the farthest flange fibre from the c.g
	set y1 [expr -$x-$xo]; # distance from the load centre to the web fibre nearest the shear centre
	set y2 [expr -$x-$xo+$w]; # distance from the load centre to the web fibre farthest from the shear centre
	set y3 [expr -$x-$xo+$b]; # distance from the load centre to the flange tip farthest from the shear centre
	
	section fiberSec $secTag {
		# Web
		# patch quad $matTag $nfw $nfd $y1 [expr -$z1] $y2 [expr -$z1] $y2 $z1 $y1 $z1; # Weak axis in the plane
		patch quad $matTag $nfd $nfw [expr -$z1] $y1 $z1 $y1 $z1 $y2 [expr -$z1] $y2; # Strong axis in the plane  
 		# Top Flange
		# patch quad $matTag $nft $nfb $y1 $z1 $y3 $z1 $y3 $z2 $y1 $z2; # Weak axis in the plane
		patch quad $matTag $nfb $nft $z1 $y1 $z2 $y1 $z2 $y3 $z1 $y3; # Strong axis in the plane
 		# Bottom Flange
		# patch quad $matTag $nft $nfb $y1 [expr -$z2] $y3 [expr -$z2] $y3 [expr -$z1] $y1 [expr -$z1]; # Weak axis in the plane
		patch quad $matTag $nfb $nft [expr -$z2] $y1 [expr -$z1] $y1 [expr -$z1] $y3 [expr -$z2] $y3; # Strong axis in the plane
	}
	
	section Aggregator $aggSecTag $matTag My -section $secTag
	
	return $aggSecTag;
}

proc ::rack_pack::make_fibre_boxed_channel_section { matTag d b t w x xo } {
	
	####################################################################
	# This procedure is called by rack_MRF in order to make channel
	# shaped fibre cross-sections
	#
	# INPUTS:
	# (for more detail see tables of CISC Handbook)
	# d := depth, double [m]
	# b := flange width, double [m]
	# t := flange thickness, double [m]
	# w := we thickness, double [m]
	# x := distance from neutral axis to exterior web, double [m]
	# xo := position of loading vis-a-vis neutral axis, double (!!! be careful of sign!) [m]
	#
	# Example : Channel C100x7 section loaded in its shear centre : make_fibre_channel_section 0.102 0.040 0.0075 0.0032 0.0126 -0.0273 
	#
	# OUTPUT : secTag := an identifier for the section, integer [unitless]
	#
	# Last mod.: 16 avril 2016, EJ
	####################################################################

	set secTag [get_Tag "Section"];
	set aggSecTag [get_Tag "Section"];
	
	set nfd 4;	# number of fibers along web depth 
	set nfw 1;	# number of fibers along web thickness
	set nfb 1;	# number of fibers along flange width (you want this many in a bi-directional loading)
	set nft 8;	# number of fibers along flange thickness
	
	#Some coordinates delineating channel web and flange
	set z1 [expr ($d-2*$t)/2]; # distance to nearest flange fibre from the c.g.
	set z2 [expr $d/2]; # distance to the farthest flange fibre from the c.g
	set y1 [expr -$x-$xo]; # distance from the load centre to the web fibre nearest the shear centre
	set y2 [expr -$x-$xo+$w]; # distance from the load centre to the web fibre farthest from the shear centre
	set y3 [expr -$x-$xo+$b]; # distance from the load centre to the flange tip farthest from the shear centre
	
	section fiberSec $secTag {
		# Web1
		patch quad $matTag $nfd $nfw [expr -$z1] $y1 $z1 $y1 $z1 $y2 [expr -$z1] $y2; # Strong axis in the plane  
 		# Top Flange
		patch quad $matTag $nfb $nft $z1 $y1 $z2 $y1 $z2 [expr 2*$y3-$y1] $z1 [expr 2*$y3-$y1]; # Strong axis in the plane
 		# Bottom Flange
		patch quad $matTag $nfb $nft [expr -$z2] $y1 [expr -$z1] $y1 [expr -$z1] [expr 2*$y3-$y1] [expr -$z2] [expr 2*$y3-$y1]; # Strong axis in the plane
		# Web2
		patch quad $matTag $nfd $nfw [expr -$z1] [expr 2*$y3-$y2] $z1 [expr 2*$y3-$y2] $z1 [expr 2*$y3-$y1] [expr -$z1] [expr 2*$y3-$y1]; # Strong axis in the plane
	}
	
	section Aggregator $aggSecTag $matTag My -section $secTag
	
	return $aggSecTag;
}

proc ::rack_pack::secProp { prflName prop } {

	####################################################################
	# Returns a section from data base of properties of a CISC Channel or Miscellaneous Channel
	#
	# INPUTS:
	# (for more detail see tables of CISC Handbook)
	# D := depth, double [m]
	# B := flange width, double [m]
	# T := flange thickness, double [m]
	# W := we thickness, double [m]
	# X := distance from neutral axis to exterior web, double [m]
	# Xo := position of loading vis-a-vis neutral axis, double (!!! be careful of sign!) [m]
	# Ix := Strong axis interia, double [m4]
	# Iy := Weak axis interia, double [m4]
	# A := Area, double [m2]
	#
	# OUTPUT:
	# The value of the requested section property of the specified section, double
	#
	# Example : secProp "C100x7" "D" would return 0.0102
	#
	# Last mod.: 23 mai 2016, EJ
	####################################################################
	
	if { $prflName == "poly-cold-form-beam" } {
		if {$prop == "Ix"} { return [expr 1.8239*pow(25.4/1000,4)] };
		if {$prop == "A"} { return [expr 0.91936*pow(25.4/1000,2)] };
	}
	if { $prflName == "poly-cold-form-col" } {
		if {$prop == "Ix"} { return [expr 1.5606*pow(25.4/1000,4)] };
		if {$prop == "A"} { return [expr 0.9492*pow(25.4/1000,2)] };
	} else {	
	if { file exists StandardAndMiscChannels.db } {
		# Slurp up the data file
		set fp [open "StandardAndMiscChannels.db" r]
		set file_data [read $fp]
		close $fp
		# Process data file
		set data [split $file_data "\n"]
		set counter 0;
		foreach line $data {
			set counter [expr $counter+1];
			if { [string match *$prflName* $line] } {
				set sec_data [split $line "\t"]
				if {$prop == "D"} { return [expr double([lindex $sec_data 6])/1000] };
				if {$prop == "B"} { return [expr double([lindex $sec_data 7])/1000] };
				if {$prop == "T"} { return [expr double([lindex $sec_data 8])/1000] };
				if {$prop == "W"} { return [expr double([lindex $sec_data 9])/1000] };
				if {$prop == "X"} { return [expr double([lindex $sec_data 22])/1000] };
				if {$prop == "Xo"} { return [expr double([lindex $sec_data 25])/1000] };
				if {$prop == "Ix"} { return [expr double([lindex $sec_data 16])/1000000] };
				if {$prop == "Iy"} { return [expr double([lindex $sec_data 19])/1000000] };
				if {$prop == "A"} { return [expr double([lindex $sec_data 15])/1000000] };
			}
		}
	} else {
	# Get inertias from user
		puts "Data-base of sections was not found"
	}
	}
}

proc ::rack_pack::record_rack { thing_to_record what_to_record {record_time 0}  } {
	
	#################################################################################################################################
	# Records nodes or elements according to previously constructed lists of nodes and elements
	#
	# INPUTS:
	# thing_to_record := a string, designating the name of the list of elements or nodes to records ie. "beams" or "connectors" etc.
	# what_to_record := a string, that can be either "force" "deformation" "disp" "reaction" or "section_forces"
	# record_time := double, output is recorded at specified time if greater than 0, default is zero
	#
	################################################################################################################################
	
	variable output_dir
	puts_and_save  "Recorders:"
	
	if { $thing_to_record == "beams" } { variable beam_elem_list; set ele_list $beam_elem_list ; puts_and_save  "Beams Elements :" };
	if { $thing_to_record == "beam_ends"} { variable beam_end_list; set ele_list $beam_end_list ; puts_and_save  "Beam-End Elements :" };
	if { $thing_to_record == "columns" } { variable col_elem_list; set ele_list $col_elem_list ; puts_and_save  "Column Elements :" };
	if { $thing_to_record == "base-plates" } { variable bpl_elem_list; set ele_list $bpl_elem_list ; puts_and_save  "Base-Plate Elements :" };
	if { $thing_to_record == "connectors" } { variable connector_elem_list; set ele_list $connector_elem_list ; puts_and_save  "Connector Elements :" };
	if { $thing_to_record == "pallets" } { variable pallet_elem_list; set ele_list $pallet_elem_list ; puts_and_save  "Pallet Elements :" };
	
	if { $thing_to_record == "palletnodes" } { variable pallet_node_list; set node_list $pallet_node_list ; puts_and_save  "Pallet Nodes :" };
	if { $thing_to_record == "beamnodes"} { variable beam_node_list; set node_list $beam_node_list ; puts_and_save  "Beam Nodes :" };
	if { $thing_to_record == "basenodes" } { variable base_node_list; set node_list $base_node_list ; puts_and_save  "Base Nodes :"}; 
	if { $thing_to_record == "ext_col_nodes" } { variable extcol_node_list; set node_list $extcol_node_list ; puts_and_save  "Exterior Column Nodes @ Joints :"};
	if { $thing_to_record == "int_col_nodes" } { variable intcol_node_list; set node_list $intcol_node_list ; puts_and_save  "Interior Column Nodes @ Joints :"};
	if { $thing_to_record == "col_mid_nodes"} { variable col_mid_node_list; set node_list $col_mid_node_list ; puts_and_save  "Column Nodes @ Mid-Spean :"}
	
	if { $record_time > 0 } { set time_string "-dT $record_time"; } else { set time_string ""; }
	
	if {$what_to_record == "section_forces"} {
	puts_and_save  "Recording all section forces for elements : $ele_list"
		for { set i 0 } { $i <= [expr [llength $ele_list]-1] } { incr i } { 
			set e [lindex $ele_list $i];
			set cmd "recorder Element -file [format $output_dir/ele_force$e.out] $time_string -ele $e force;"; eval $cmd;
		}
	}
	
	if {$what_to_record == "force"} {	
	puts_and_save  "Recording local forces for elements : $ele_list"
		for { set i 0 } { $i <= [expr [llength $ele_list]-1] } { incr i } { 
			set e [lindex $ele_list $i];
			set cmd "recorder Element -file [format $output_dir/ele_force$e.out] $time_string -ele $e localForce;"; eval $cmd;
		}
	}
	
	if {$what_to_record == "deformation"} { 
	puts_and_save  "Recording deformations for elements : $ele_list"
		for { set i 0 } { $i <= [expr [llength $ele_list]-1] } { incr i } { 
			set e [lindex $ele_list $i];
			set cmd "recorder Element -file [format $output_dir/ele_def$e.out] $time_string -ele $e deformation;"; eval $cmd;
		}
	}
	if {$what_to_record == "disp"} {
	puts_and_save  "Recording displacements at nodes : $node_list"
		for { set i 0 } { $i <= [expr [llength $node_list]-1 ] } { incr i } { 
			set n [lindex $node_list $i];
			set cmd "recorder Node -file [format $output_dir/node_disp$n.out] $time_string -node $n -dof 1 2 3 disp;"; eval $cmd;
		}
	}
	
	if {$what_to_record == "reaction"} {
	puts_and_save  "Recording reactions at nodes : $node_list"
		for { set i 0 } { $i <= [expr [llength $node_list]-1 ] } { incr i } {
			set n [lindex $node_list $i];
			set cmd "recorder Node -file [format $output_dir/node_reac$n.out] $time_string -node $n -dof 1 2 3 reaction;"; eval $cmd;
		}
	}
	if {$what_to_record == "accel"} {
	puts_and_save  "Recording reactions at nodes : $node_list"
		for { set i 0 } { $i <= [expr [llength $node_list]-1 ] } { incr i } {
			set n [lindex $node_list $i];
			set cmd "recorder Node -file [format $output_dir/node_acc$n.out] $time_string -node $n -dof 1 2 3 accel;"; eval $cmd;
		}
	}	

}

proc ::rack_pack::check_total_reaction { dof step } {

	################################################################################################################
	# Calculates the total reaction at the base nodes in a specified direction at a specified time step
	#
	# INPUTS :
	# dof := an integer, direction the reaction [unitless] (1=x,2=y,3=z)
	# step := the time step at which the reaction is to be calculated
	#
	# Example :  
	#
	# EJ
	################################################################################################################

	variable base_node_list;
	variable output_dir;
	
	# Base-plate recorders must be closed before reading the files
	remove recorders;
	
	set base_nodes [llength $base_node_list];
	set total_reaction 0.;
	 
	# Read recorded out files node by node
	for { set i 0 } { $i < $base_nodes } { incr i } {
		set nodeI [lindex $base_node_list $i]; # Get item "i" in list of base nodes
		if { [catch {set fp [open $output_dir/node_reac$nodeI.out]} errmsg] } {
			puts_and_save  "Verify that base-node recorders are set: $errmsg";
			return 0.;
		}		
		set file_data [read $fp]; # Slurp up data in file
		close $fp; # Close connection
		file delete $output_dir/node_reac$nodeI.out -force; # Delete the reaction output file
		# Treat slurped up data to get reaction at the specified step
		set data [split $file_data "\n"]; 
		set lastline [lindex $data [expr $step-1]]; # Parse data to get last line which contains the reactions at the end of static analysis		
		set node_reaction [lindex $lastline [expr $dof-1]];
		set total_reaction [expr $total_reaction + $node_reaction];
	}
	
	# Add-up and display total seismic weight
	set message  "Total vertical reaction is : "
	puts_and_save  [append message $total_reaction];
	
	return $total_reaction;

}

proc ::rack_pack::make_output_dir { dirname } {

	#############################################################
	#
	#
	#
	#
	#############################################################
	
	variable output_dir;
	
	set output_dir $dirname;
	set i 1; 
	
	# If an analysis has already been run in the output directory then make the output directory the default with a number appended
	while { [file exists $output_dir$i] == 1 } { set i [expr $i+1] };
	
	set output_dir $output_dir$i;
	file mkdir $output_dir;
	puts "Writing output to $output_dir";

}

proc ::rack_pack::puts_and_save { string_to_save } {
	
	variable output_dir;
	set out_file [open "$output_dir/a_model_and_analysis_info.out" "a+"]
	puts $string_to_save
	puts $out_file $string_to_save
	close $out_file
	
}

proc rotSpring2D {eleID nodeR nodeC matID} {

	# SETS A MULTIPOINT CONSTRAINT ON THE TRANSLATIONAL DEGREES OF FREEDOM,
	# SO DO NOT USE THIS PROCEDURE IF THERE ARE TRANSLATIONAL ZEROLENGTH
	# ELEMENTS ALSO BEING USED BETWEEN THESE TWO NODES
	#
	# Written: MHS
	# Date: Jan 2000
	#
	# Formal arguments
	#	eleID - unique element ID for this zero length rotational spring
	#	nodeR - node ID which will be retained by the multi-point constraint
	#	nodeC - node ID which will be constrained by the multi-point constraint
	#	matID - material ID which represents the moment-rotation relationship
	#	for the spring

	# Create the zero length element
	element zeroLength $eleID $nodeR $nodeC -mat $matID -dir 6
	# Constrain the translational DOF with a multi-point constraint
	# retained constrained DOF_1 DOF_2 ... DOF_n
	equalDOF $nodeR $nodeC 1 2

}
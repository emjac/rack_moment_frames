###
# A collection of useful functions and variables
require(data.table);
options(warn=-1)
proj_dos = "~/Drive/RackProject"; # General project directory
setwd(proj_dos)
my_cols = c("black","blue","red","darkgreen","cornflowerblue","orange","purple", "gray","darkgreen", "darkred","darkorange","pink")

source("read_model_results.r")
source("load_seismic_test_data.r")
source("load_static_test_data.r")
source("SectionDesign.r");
source("spectra.r");
library(signal)

get_secant_stiffness = function(theta, moment){
  
  ## Because of bumps and general fuzziness in the test data sometimes using the maximum moment gets the wrong
  # peak in a cycle (a peak that is not remotely near a theta peak), also getting secant stiffnes too close
  # to zero can result in errroneously large values
  # I calculate the secant stiffness at all points sufficiently far away from zero and also at the max abs moments
  # in each cycle, then I take the least of the calculated values, it is not fool proof but it gives reasonable results
  # Sudden spikes or dips in secant stiffness should be verfied
  
  # Calculate the secant stiffness at every point sufficiently far away from zero
  ksecs = c()
  for (i in c(1:length(theta))) {
    if (abs(theta[i]) > 0.005) { ksecs = c(ksecs,moment[i]/theta[i]) }
  }
  ksec1 = (abs(max(ksecs)) + abs(min(ksecs)))/2;
  
  # Get the values of theta that coorespond to max and min moment (sometimes there is more than one value of theta and in that case the largest/smallest value of theta is taken)
  Mcpeak = (abs(max(moment))+abs(min(moment)))/2;
  thetaCpeakNeg = theta[min(which(moment == min(moment)))];
  thetaCpeakPos = theta[max(which(moment == max(moment)))];
  # Average of min and max
  thetaCpeak = (abs(thetaCpeakPos) + abs(thetaCpeakNeg))/2;
  ksec2 = Mcpeak/thetaCpeak;
  
  thetapeak = (max(theta)+abs(min(theta)))/2
  moment_min_peak = moment[min(which(theta ==min(theta)))];
  moment_max_peak = moment[max(which(theta == max(theta)))];
  ksec3 = (abs(moment_min_peak) + moment_max_peak)/thetapeak
  
  med_ksec = median(c(ksec1,ksec2,ksec3))
  
  # Return the most reasonable of the two values
  if ( abs(ksec2-med_ksec)/med_ksec > 0.2 ) { return(min(ksec1,ksec2,ksec3)) } else { return(ksec2) } 
  
  return( median(c(ksec1,ksec2,ksec3))  )
  
}
get_EDC = function(theta,moment){
  
  EDC = c()
  i = 1;
  
  for (i in c(1:length(theta)-1)) {
    EDC[i] = (theta[i+1]-theta[i])*(moment[i+1]-moment[i])/2 + moment[i]*(theta[i+1]-theta[i]) 
  }
  
  return (sum(EDC))
  
}
get_essai_cycles = function(rot, tol, minpic){
  ## Get cycles from the essai
  cycles = c(1) # Start the first cycle at index 1
  i = 1 # Index counter to loop through drift vector
  passed_a_max = 0 # value to check if the drift has passed by a maximum in the cycle
  passed_a_min = 0
  while (i < length(rot)){
    d = rot[[i]];
    # check if the drift is passing by zero
    if (isTRUE(all.equal(0,d,tol))){ 
      # check if the drift has passed by a maximum in that cycle before arriving at zero
      # this is to prevent recording the pass by zero as a cycle if the drift is essentially
      # stationnary at zero
      if (passed_a_max*passed_a_min == 1){ 
        cycles = c(cycles,i)
        passed_a_max = 0 # reset to zero until next maximum is reached
        passed_a_min = 0
      } 
    } else {
      # Check if rotation is at least
      if ( d < -minpic ) { passed_a_min = 1 }
      if ( d > minpic ) { passed_a_max = 1 }
    }
    i=i+1
  }
  cycles = c(cycles,length(rot)) # Last drift value is the end of the last cycle
  return (cycles)
};
get_EDC_at_targ = function(theta_targ, data){ 
  
  rot = data$AvgRotPeak
  EDC = data$EDC
  
  EDC1st_pass = c(EDC[1])
  rot1st_pass = c(rot[1])
  
  for (i in c(2:(length(rot)-1))) {
    if ( all.equal(rot[i-1],rot[i], tolerance = 0.01) == TRUE ) { 
    }
    else {
      EDC1st_pass = c(EDC1st_pass,EDC[i])
      rot1st_pass =c(rot1st_pass,rot[i])
    }
    
  }
  
  plot(rot1st_pass, EDC1st_pass/1000, col = "black", xlab = expression(paste(theta," [rad]")), ylab = expression(paste(EDC," [kNmrad]")) )
  grid(col = "lightgray", lty = "dotted", equilogs = TRUE)
  k_targ = as.numeric((approx(rot1st_pass, EDC1st_pass, xout=theta_targ)[2]))
  abline(v=theta_targ, col = "blue")
  abline(h=k_targ/1000, col = "cornflowerblue")
  # legend("topleft", "Base-Plate Pseudo Test", pch = 1, col = "blue")
  return (k_targ); # Value of EDC is returned in Nm/rad  
  
}
get_keff_at_targ = function(theta_targ, data) {
  
  rot = data$AvgRotPeak
  ksec = data$ksec
  
  ksec1st_pass = c(ksec[1])
  rot1st_pass = c(rot[1])
  
  for (i in c(2:(length(rot)-1))) {
    # print(abs(rot[i-1]-rot[i]))
    if ( abs(rot[i-1]-rot[i])  > 0.001 ) {
      # if ( all.equal(rot[i-1],rot[i], tolerance = 0.01) == TRUE ) { 
      ksec1st_pass = c(ksec1st_pass,ksec[i])
      rot1st_pass =c(rot1st_pass,rot[i])
    }
    else {
      
    }
  }
  
  plot(rot1st_pass, ksec1st_pass/1000, xlab = expression(paste(theta," [rad]")), ylab = expression(paste(k[sec]," [kNm/rad]")))
  grid(col = "lightgray", lty = "dotted", equilogs = TRUE)
  k_targ = as.numeric((approx(rot1st_pass, ksec1st_pass, xout=theta_targ)[2]))
  abline(v=theta_targ, col = "blue")
  abline(h=k_targ/1000, col = "cornflowerblue")
  return (k_targ); # Value of ksec is returned in Nm/rad
}
run_modal_analysis = function(h1, hx, Lb, beam, Fyb, col, Fyc, boxed_levels, PL, nB, nL, sprOption, bplOption, kc, kctop, kbpl, transOption, g_tol, l_tol) {
  
  if(.Platform$OS.type == "unix") {
    proj_dos = getwd()
    out_dos = paste0(proj_dos,"/rack_nL",nL,"nB",nB,"_1"); # OpenSees output directory
    system("mkdir -p modes");
    system("touch modes/period.out");
    fileConn=file("main.tcl");
    # Edit main.tcl file with model parameters
    writeLines(c("lappend ::auto_path [eval pwd];",
                 "package require rack_pack;",
                 "namespace import rack_pack::*",
                 paste0("rackMRF ", h1," ", hx," ", Lb, ' "',beam,'" ', Fyb, ' "', col,'" ', Fyc," ", boxed_levels, " ", PL," ", nB," ", nL," ",0.," ",1000," ",1000," ", ' "',sprOption, '"', ' "',bplOption, '" ', kc," ",kctop," ",kbpl,";"),
                 "source modal_analysis.tcl;","modal_analysis 1;"), fileConn);
    close(fileConn);
    system("echo source main.tcl | OpenSees"); # Pipe "source main.tcl" command to OpenSees
    df = read.table("modes/period.out");
    Teff = df$V1[1];
    system("rm -rf modes");
    system(paste("rm -rf", basename(out_dos)));
    return(Teff)
  } else {
    shell("mkdir modes"); #Windows
    shell("type nul >modes/period.out"); #Windows
    fileConn=file("main.tcl");
    # Edit main.tcl file with model parameters
    writeLines(c("lappend ::auto_path [eval pwd];",
                 "package require rack_pack;",
                 "namespace import rack_pack::*",
                 paste0("rackMRF ", h1," ", hx," ", Lb, ' "',beam,'" ', Fyb, ' "', col,'" ', Fyc," ", boxed_levels, " ", PL," ", nB," ", nL," ",0.," ",0.," ",0., ' "',sprOption, '"', ' "',bplOption, '" ', kc," ",kctop," ",kbpl,' "',transOption,'" ',g_tol, " ",l_tol, ";"),
                 "source modal_analysis.tcl;","modal_analysis 1;"), fileConn);
    close(fileConn);
    shell("echo source main.tcl | OpenSees"); # Pipe "source main.tcl" command to OpenSees
    shell("del main.tcl /f /q")
    df = read.table("modes/period.out");
    T1 = df$V1[1];
    shell("rmdir modes /s /q");
    return(T1)
  }
}
get_Sd = function (test,Tn) {
  
  spectra = fread(paste0(proj_dos,"/Accelerograms/", test,"_spectra.txt"));
  spectra = spectra[-1,];
  spectra = spectra[-1,];
  
  ymax = as.numeric(spectra$Sd[length(spectra$Sd)])
  
  my_Sd = approx(spectra$Periode, spectra$Sd , xout = Tn)$y;
  plot(0,0, type ="l", xlim = c(0,3), ylim = c(0,ymax), xlab = "Period [s]", ylab = "Sd [m]");
  lines(spectra$Periode, spectra$Sd, lty = 1);
  points(Tn,my_Sd)
  return(my_Sd)
}
get_rot_max = function(theta, moment){
  
  j = 1;
  
  if ( abs(min(theta)) > abs(max(theta)) ) { 
    thetaE = min(theta); 
  } else { 
    thetaE = max(theta); 
  }
  
  while (thetaE != theta[j]) {
    j = j + 1;
  }
  
  return (theta[j])
}

pi = 2*asin(1.0);
prpt = fread("StandardAndMiscChannels.db")
g = 9.80665; # Acceleration of gravity [m/s2]
E = 200e9;
# Load recorded results
#

### Load results
beam_output = function (directory,nB,nL,n_div,output_type){
  
  ####################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when beams is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := of elements the beam in each bay at each level has been divided into
  # output_type := a string, either "def" or "force" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  ####################################################################################################
  
  output_list = vector("list", (n_div-1)*nL*nB); # Initialise the output_list with a certain length
  n_list = 0;
  
  for (i in c(1:nL)) {
    for (j in c(1:nB)) {
      for (k in c(1:(n_div-1))) {
        n_list = n_list+1;
        nodeI = paste0(i+1,j,k);
        nodeJ = paste0(i+1,j,k+1);
        ele = paste0(nodeI,nodeJ)
        if (output_type == "force") {
          df = fread(paste0("ele_force",ele,".out"));
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
        }
        if (output_type == "def"){
          df = fread(paste0("ele_def",ele,".out"));
          colnames(df) = c("dxi","dyi","rzi","dxj","dyj","rzj");
        }
        output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
        names(output_list)[n_list] = paste0(ele); # Name the list item after the element
      }
    }
  }
  
  return(output_list);
  
};
column_output = function (dir,nB,nL,n_div,output_type){
  
  ####################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when columns is passed as a parameter to that procedure
  #
  # Inputs :
  # dir := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := of elements the column in each bay at each level has been divided into
  # steps := the analysis steps that we wish to see outputed
  # output_type := a string, "def", "force", "Fx", "Fy", "Mz" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  ####################################################################################################
  
  output_list = vector("list", n_div*nL*(nB+1));
  n_list = 0;
  
  if (n_div > 1) {
    # Non linear columns ### Not sure if any of the non-linear part works
    # for (i in c(1:nL)) { 
    #   for (j in c(1:(nB+1))) {
    #     for (k in c(1:n_div)) {
    #       n_list = n_list + 1;
    #       # nodeI = paste0(i,j);
    #       # nodeJ = paste0((i+1),j);
    #       ele = paste0(nodeI,nodeJ)
    #       if (output_type == "force") {
    #         df = fread(paste0("ele_force",ele,".out"));
    #         colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
    #       }
    #       if (output_type == "def"){
    #         df = fread(paste0("ele_def",ele,".out"));
    #         colnames(df) = c("dxi","dyi","rzi","dxj","dyj","rzj");
    #       }
    #       output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
    #       names(output_list)[n_list] = paste0(ele); # Name the list item after the element
    #     }
    #   }
    # }
  } else {
    # Linear columns
    output_df = NULL
    ele_names = NULL
    
    for (i in c(1:nL)) { 
      for (j in c(1:(nB+1))) { 
        n_list = n_list + 1;
        nodeI = paste0(i,j);
        nodeJ = paste0((i+1),j);
        ele = paste0(nodeI,nodeJ);
        if (output_type == "Fx") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fyi = df$Fyj = df$Mzi = df$Mzj = NULL
        }
        if (output_type == "Fxi") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fxj = df$Fyi = df$Fyj = df$Mzi = df$Mzj = NULL
        }
        if (output_type == "Fxj") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fxi = df$Fyi = df$Fyj = df$Mzi = df$Mzj = NULL
        }
        if (output_type == "Fy") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fxi = df$Fxj = df$Mzi = df$Mzj = NULL
        }
        if (output_type == "Fyi") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fyj = df$Fxi = df$Fxj = df$Mzi = df$Mzj = NULL;
        }
        if (output_type == "Fyj") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fyi = df$Fxi = df$Fxj = df$Mzi = df$Mzj = NULL
        }
        if (output_type == "Mz") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fxi = df$Fxj = df$Fyi = df$Fyj = NULL
        }
        if (output_type == "Mzi") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fxi = df$Fxj = df$Fyi = df$Fyj = df$Mzj = NULL
        }
        if (output_type == "Mzj") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          df$Fxi = df$Fxj = df$Fyi = df$Fyj = df$Mzi = NULL
        }
        if (output_type == "force") {
          df = fread(paste0(dir,"/ele_force",ele,".out"))[steps];
          colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
          names(output_list)[n_list] = paste0(ele); # Name the list item after the element
        }
        if (output_type == "def"){
          df = fread(paste0(dir,"/ele_def",ele,".out"))[steps];
          colnames(df) = c("dxi","dyi","rzi","dxj","dyj","rzj");
          output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
          names(output_list)[n_list] = paste0(ele); # Name the list item after the element
        }
        output_df = cbind(output_df,c(df));
        ele_names = c(ele_names, ele);
      }
    }
  }
  if (output_type == "def" || output_type == "force") {output_list} else {
    colnames(output_df) = ele_names;
    return(output_df)
  }
};
bpl_output = function (directory,nB,nL,steps,output_type){
  
  ####################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when "baseplates" is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # output_type := a string, either "def" or "force" depending on what output information is requested
  #
  # outputs : 
  # - when fibre base-plates are used the force and deformation outputs are from the 2 zero-length compression springs at the base of each plate
  # - when concentraed plasticity is used the force and deformation outputs are the moment and rotations from the single zero-length spring.
  #
  # EJ
  ####################################################################################################
  
  output_list = NULL;
  n_list = 1;
  output_df = NULL;
  
  for (j in c(1:(nB+1))) {
    nodeI = paste0(j,0) 
    nodeJ = paste0(1,j)
    if (output_type == "force"){
      fiber_base_plate_eles = Sys.glob(paste0(directory,"/ele_force",nodeJ,00,"*",".out")); # Slurp up everything following the naming pattern of fiber base-plate elements
      concentrated_base_plate_eles = Sys.glob(paste0(directory,"/ele_force",nodeI,"*",".out")); # Slurp up everything following the naming pattern of concentrated base-plate elements
      if (length(fiber_base_plate_eles)>0) {
        for ( k in c(1:length(fiber_base_plate_eles))) {
          df = fread(fiber_base_plate_eles[k])[steps];
          colnames(df) = c("dy");
          ele_tag = basename(fiber_base_plate_eles[k]); # Get just the file name
          ele_tag = sub("ele_force","",ele_tag); # Trim off text before node tag
          ele_tag = sub(".out","",ele_tag); # Time off text after node tag
          output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
          names(output_list)[n_list] = paste0(ele_tag); # Name the list item after the element
          n_list = n_list + 1;
        }
      }
      if (length(concentrated_base_plate_eles)>0) {
        for ( k in c(1:length(concentrated_base_plate_eles))) {
          df = fread(concentrated_base_plate_eles[k])[steps];
          colnames(df) = c("dy");
          ele_tag = basename(concentrated_base_plate_eles[k]); # Get just the file name
          ele_tag = sub("ele_force","",ele_tag); # Trim off text before node tag
          ele_tag = sub(".out","",ele_tag); # Time off text after node tag
          output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
          names(output_list)[n_list] = paste0(ele_tag); # Name the list item after the element
          n_list = n_list + 1;
        }
      }
    }
    if (output_type == "def"){
      fiber_base_plate_eles = Sys.glob(paste0(directory,"/ele_def",nodeJ,00,"*",".out")); # Slurp up everything following the naming pattern of fiber base-plate elements
      concentrated_base_plate_eles = Sys.glob(paste0(directory,"/ele_def",nodeI,"*",".out")); # Slurp up everything following the naming pattern of concentrated base-plate elements
      if (length(fiber_base_plate_eles)>0) {
        for ( k in c(1:length(fiber_base_plate_eles))) {
          df = fread(fiber_base_plate_eles[k])[steps];
          colnames(df) = c("dy");
          ele_tag = basename(fiber_base_plate_eles[k]); # Get just the file name
          ele_tag = sub("ele_def","",ele_tag); # Trim off text before node tag
          ele_tag = sub(".out","",ele_tag); # Time off text after node tag
          output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
          names(output_list)[n_list] = paste0(ele_tag); # Name the list item after the element
          n_list = n_list + 1;
        }
      }
      if (length(concentrated_base_plate_eles)>0) {
        for ( k in c(1:length(concentrated_base_plate_eles))) {
          df = fread(concentrated_base_plate_eles[k])[steps];
          colnames(df) = c("dy");
          ele_tag = basename(concentrated_base_plate_eles[k]); # Get just the file name
          ele_tag = sub("ele_def","",ele_tag); # Trim off text before node tag
          ele_tag = sub(".out","",ele_tag); # Time off text after node tag
          output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
          names(output_list)[n_list] = paste0(ele_tag); # Name the list item after the element
          n_list = n_list + 1;
        }
      }
    }
  }
  return(output_list)
}
connector_output = function (directory,nB,nL,n_div,output_type){
  
  ####################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when baseplates is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := the number of beam sub-elements
  # output_type := a string, either "def" or "force" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  ####################################################################################################
  
  output_list = vector("list", nL*(2*nB));
  n_list = 0;
  
  ele_tags = NULL
  output_df = NULL
  
  # Define beam-column connectors
  for ( i in c(2:(nL+1))) {
    for (j in c(1:(nB))) {
      # Leftmost connector in the bay at level i
      n_list = n_list+1;
      k=1;
      nodeI = paste0(i,j);
      nodeJ = paste0(i,j,k);
      ele = paste0(nodeI,nodeJ);
      if (output_type == "Mz") {
        df = fread(paste0(directory,"/ele_force",ele,".out"));
        colnames(df) = c("Mz");
        output_df = cbind(output_df,c(df$Mz));
        ele_tags = c(ele_tags, ele);
      }
      if (output_type == "rz"){
        df = fread(paste0(directory,"/ele_def",ele,".out"));
        colnames(df) = c("rz");
        output_df = cbind(output_df,c(df$rz));
        ele_tags = c(ele_tags, ele);
      }
      output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
      names(output_list)[n_list] = paste0(ele); # Name the list item after the element
      # Rightmost connector in the bay at level i
      n_list = n_list+1;
      k = n_div;
      nodeI = paste0(i,j+1);
      nodeJ = paste0(i,j,k);
      ele = paste0(nodeI,nodeJ);
      if (output_type == "Mz") {
        df = fread(paste0(directory,"/ele_force",ele,".out"));
        colnames(df) = c("Mz");
        output_df = cbind(output_df,c(df$Mz));
        ele_tags = c(ele_tags, ele);
      }
      if (output_type == "rz"){
        df = fread(paste0(directory,"/ele_def",ele,".out"));
        colnames(df) = c("rz");
        output_df = cbind(output_df,c(df$rz));
        ele_tags = c(ele_tags, ele);
      }
      output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
      names(output_list)[n_list] = paste0(ele); # Name the list item after the element
    }
  }
  colnames(output_df) = ele_tags
  return(output_df)
};
pallet_ele_output = function (directory,nB,nL,n_div,output_type){
  ####################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when baseplates is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := of elements the beam in each bay at each level has been divided into
  # output_type := a string, either "def" or "force" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  ####################################################################################################
  
  output_list = vector("list", (n_div-1)*nL*nB); # Initialise the output_list with a certain length
  n_list = 0;
  l = 1; # End node of pallets always end with 1
  
  output_df = NULL
  
  for (i in c(1:nL)) {
    for (j in c(1:nB)) {
      for (k in c(1:(n_div-1))) {
        if (k != 1 && k!=n_div) {
          n_list = n_list+1;
          nodeI = paste0((i+1),j,k); 
          nodeJ = paste0((i+1),j,k,l);
          ele = paste0(nodeI,nodeJ);
          if (output_type == "dx") {
            df = fread(paste0(directory,"/ele_def",ele,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$dx))
          }
          if (output_type == "dy") {
            df = fread(paste0(directory,"/ele_def",ele,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$dy))
          }
          if (output_type == "rz") {
            df = fread(paste0(directory,"/ele_def",ele,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$rz))
          }
          if (output_type == "force") {
            df = fread(paste0(directory,"/ele_force",ele,".out"));
            colnames(df) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
          }
          if (output_type == "deformation"){
            df = fread(paste0(directory,"/ele_def",ele,".out"));
            colnames(df) = c("dx","dy","rz");
          }
          output_list[[n_list]] = df; # Add the dataframe of loaded element output data to the output list
          names(output_list)[n_list] = paste0(ele); # Name the list item after the element 
        }
      }
    }
  }
  
  if (output_type == "deformation" || output_type == "force") {output_list} else {
    return(output_df)
  }
};
beam_node_output = function (directory,nB,nL,n_div,output_type){
  
  #######################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when pallet_nodes is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := of elements the beam in each bay at each level has been divided into
  # output_type := a string, either "reac" or "disp" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  #######################################################################################################
  
  output_list = vector("list", (n_div-2)*nL*nB); # Initialise the output_list with a certain length
  
  output_df = NULL
  node_tags = NULL
  
  for (i in c(1:nL)) {
    for (j in c(1:nB)) {
      for (k in c(2:(n_div-1))) {
        node = paste0((i+1),j,k)
        if (output_type == "dx") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
          output_df = cbind(output_df,c(df$dx))
          node_tags = c(node_tags, node);
        }
        if (output_type == "dy") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
          output_df = cbind(output_df,c(df$dy))
          node_tags = c(node_tags, node);
        }
        if (output_type == "rz") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
          output_df = cbind(output_df,c(df$rz))
          node_tags = c(node_tags, node);
        }
        if (output_type == "ax") {
          df = fread(paste0(directory,"/node_acc",node,".out"));
          colnames(df) = c("ax","ay","az");
          output_df = cbind(output_df,c(df$ax));
          node_tags = c(node_tags, node);
        }
        if (output_type == "ay") {
          df = fread(paste0(directory,"/node_acc",node,".out"));
          colnames(df) = c("ax","ay","az");
          output_df = cbind(output_df,c(df$ay));
          node_tags = c(node_tags, node);
        }
        if (output_type == "az") {
          df = fread(paste0(directory,"/node_acc",node,".out"));
          colnames(df) = c("ax","ay","az");
          output_df = cbind(output_df,c(df$az));
          node_tags = c(node_tags, node);
        }
        if (output_type == "disp") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
        }
        if (output_type == "/reac") {
          df = fread(paste0(directory,"node_reac",node,".out"));
          colnames(df) = c("Fx","Fy","Mz");  
        }
        output_list[[k-1]] = df;
        names(output_list)[k-1] = paste0(node)
      }
    }
  }
  if (output_type == "disp" || output_type == "reac") {output_list} else {
    colnames(output_df) = node_tags
    return(output_df)
  }
  
};
pallet_node_output = function (directory,nB,nL,n_div,output_type){
  
  #######################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when pallet_nodes is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := of elements the beam in each bay at each level has been divided into
  # output_type := a string, either "reac" or "disp" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  #######################################################################################################
  
  output_list = vector("list", (n_div-2)*nL*nB); # Initialise the output_list with a certain length
  
  output_df = NULL
  node_tags = NULL
  
  
  for (i in c(1:nL)) {
    for (j in c(1:nB)) {
      for (k in c(2:(n_div-1))) {
        l = 1; # All pallet node names end with 1
        node = paste0((i+1),j,k,l)
        if (output_type == "dx") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
          output_df = cbind(output_df,c(df$dx))
          node_tags = c(node_tags, node);
        }
        if (output_type == "dy") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
          output_df = cbind(output_df,c(df$dy))
          node_tags = c(node_tags, node);
        }
        if (output_type == "rz") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
          output_df = cbind(output_df,c(df$rz))
          node_tags = c(node_tags, node);
        }
        if (output_type == "ax") {
          df = fread(paste0(directory,"/node_acc",node,".out"));
          colnames(df) = c("ax","ay","az");
          output_df = cbind(output_df,c(df$ax));
          node_tags = c(node_tags, node);
        }
        if (output_type == "ay") {
          df = fread(paste0(directory,"/node_acc",node,".out"));
          colnames(df) = c("ax","ay","az");
          output_df = cbind(output_df,c(df$ay));
          node_tags = c(node_tags, node);
        }
        if (output_type == "az") {
          df = fread(paste0(directory,"/node_acc",node,".out"));
          colnames(df) = c("ax","ay","az");
          output_df = cbind(output_df,c(df$az));
          node_tags = c(node_tags, node);
        }
        if (output_type == "disp") {
          df = fread(paste0(directory,"/node_disp",node,".out"));
          colnames(df) = c("dx","dy","rz");
        }
        if (output_type == "reac") {
          df = fread(paste0(directory,"/node_reac",node,".out"));
          colnames(df) = c("Fx","Fy","Mz");  
        }
        output_list[[k-1]] = df;
        names(output_list)[k-1] = paste0(node)
      }
    }
  }
  if (output_type == "disp" || output_type == "reac") {output_list} else {
    colnames(output_df) = node_tags
    return(output_df)
  }
  
};
base_node_output = function (directory,nB,nL,output_type){
  
  ############
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when basenodes is passed as a parameter
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # output_type := a string, either "reac" or "disp" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  ############
  
  nc = nB + 1; # Number of columns
  output_list = vector("list", nc); # Give the output list the needed length
  i = 1; # Base nodes all start with 1
  
  node_tags = NULL
  output_df = NULL
  
  for (j in c(1:nc)) {
    node = paste0(i,j); # Make a string that contains the nodeTag
    if (output_type == "dx") {
      base_plate_nodes = Sys.glob(paste0(directory,"/node_disp",node,"*",".out")); # Get a vector of everything in the directory that starts with node_reacij
      for ( k in c(1:length(base_plate_nodes))) {
        df = fread(base_plate_nodes[k]);
        colnames(df) = c("dx","dy","rz");
        output_df = cbind(output_df,c(df$dx));
        node_tag = basename(base_plate_nodes[k]); # Get just the file name
        node_tag = sub("node_disp","",node_tag); # Trim off text before node tag
        node_tag = sub(".out","",node_tag); # Time off text after node tag
        node_tags = c(node_tags, node_tag); 
      }
    }
    if (output_type == "dy") {
      base_plate_nodes = Sys.glob(paste0(directory,"/node_disp",node,"*",".out")); # Get a vector of everything in the directory that starts with node_reacij
      for ( k in c(1:length(base_plate_nodes))) {
        df = fread(base_plate_nodes[k]);
        colnames(df) = c("dx","dy","rz");
        output_df = cbind(output_df,c(df$dy));
        node_tag = basename(base_plate_nodes[k]); # Get just the file name
        node_tag = sub("node_disp","",node_tag); # Trim off text before node tag
        node_tag = sub(".out","",node_tag); # Time off text after node tag
        node_tags = c(node_tags, node_tag); 
      }
    }
    if (output_type == "rz") {
      base_plate_nodes = Sys.glob(paste0(directory,"/node_disp",node,"*",".out")); # Get a vector of everything in the directory that starts with node_reacij
      for ( k in c(1:length(base_plate_nodes))) {
        df = fread(base_plate_nodes[k]);
        colnames(df) = c("dx","dy","rz");
        output_df = cbind(output_df,c(df$rz));
        node_tag = basename(base_plate_nodes[k]); # Get just the file name
        node_tag = sub("node_disp","",node_tag); # Trim off text before node tag
        node_tag = sub(".out","",node_tag); # Time off text after node tag
        node_tags = c(node_tags, node_tag); 
      }
    }
    if (output_type == "Fx") {
      # In the case of fiber base-plates there are 4 nodes fixed to the ground whose tags all start with the name of the origin node ij 
      base_plate_nodes = Sys.glob(paste0(directory,"/node_reac",node,"*",".out")); # Get a vector of everything in the directory that starts with node_reacij
      for ( k in c(1:length(base_plate_nodes))) {
        df = fread(base_plate_nodes[k]);
        colnames(df) = c("Fx","Fy","Mz");
        output_df = cbind(output_df,c(df$Fx));
        node_tag = basename(base_plate_nodes[k]); # Get just the file name
        node_tag = sub("node_reac","",node_tag); # Trim off text before node tag
        node_tag = sub(".out","",node_tag); # Time off text after node tag
        node_tags = c(node_tags, node_tag); 
      }
    }
    if (output_type == "Fy") {
      base_plate_nodes = Sys.glob(paste0(directory,"/node_reac",node,"*",".out")); # Get a vector of everything in the directory that starts with node_reacij
      for ( k in c(1:length(base_plate_nodes))) {
        df = fread(base_plate_nodes[k]);
        colnames(df) = c("Fx","Fy","Mz");
        output_df = cbind(output_df,c(df$Fy));
        node_tag = basename(base_plate_nodes[k]); # Get just the file name
        node_tag = sub("node_reac","",node_tag); # Trim off text before node tag
        node_tag = sub(".out","",node_tag); # Time off text after node tag
        node_tags = c(node_tags, node_tag); 
      }
    }
    if (output_type == "Mz") {
      base_plate_nodes = Sys.glob(paste0(directory,"/node_reac",node,"*",".out")); # Get a vector of everything in the directory that starts with node_reacij
      for ( k in c(1:length(base_plate_nodes))) {
        df = fread(base_plate_nodes[k]);
        colnames(df) = c("Fx","Fy","Mz");
        output_df = cbind(output_df,c(df$Mz));
        node_tag = basename(base_plate_nodes[k]); # Get just the file name
        node_tag = sub("node_reac","",node_tag); # Trim off text before node tag
        node_tag = sub(".out","",node_tag); # Time off text after node tag
        node_tags = c(node_tags, node_tag); 
      }
    }
    if (output_type == "disp") {
      df = fread(paste0("/node_disp",node,".out"));
      colnames(df) = c("dx","dy","rz");
    }
    if (output_type == "/reac") {
      df = fread(paste0("node_reac",node,".out"));
      colnames(df) = c("Fx","Fy","Mz");  
    }
    output_list[[j]] = df;
    names(output_list)[j] = paste0(node)
  }
  if (output_type == "disp" || output_type == "reac") {output_list} else {
    colnames(output_df) = node_tags
    return(output_df)
  }
};
ext_col_node_output = function (directory,nB,nL,output_type){
  
  #######################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when ext_col_nodes is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # output_type := a string, either "reac" or "disp" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  #######################################################################################################
  
  nc = nB + 1
  output_list = vector("list", nL*2)
  n_list = 0
  
  node_tags = NULL
  output_df = NULL
  
  for (i in c(1:(nL+1))) {
    for (j in c(1:nc)) {
      if ( i != 1 ) { 
        if ( j == 1 || j == nc ) {
          n_list = n_list+1;
          node =  paste0(i,j);
          if (output_type == "dx") {
            df = fread(paste0(directory,"/node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$dx));
            node_tags = c(node_tags, node);
          }
          if (output_type == "dy") {
            df = fread(paste0(directory,"/node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$dy));
            node_tags = c(node_tags, node);
          }
          if (output_type == "rz") {
            df = fread(paste0(directory,"/node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$rz));
            node_tags = c(node_tags, node);
          }
          if (output_type == "ax") {
            df = fread(paste0(directory,"/node_acc",node,".out"));
            colnames(df) = c("ax","ay","az");
            output_df = cbind(output_df,c(df$ax));
            node_tags = c(node_tags, node);
          }
          if (output_type == "ay") {
            df = fread(paste0(directory,"/node_acc",node,".out"));
            colnames(df) = c("ax","ay","az");
            output_df = cbind(output_df,c(df$ay));
            node_tags = c(node_tags, node);
          }
          if (output_type == "az") {
            df = fread(paste0(directory,"/node_acc",node,".out"));
            colnames(df) = c("ax","ay","az");
            output_df = cbind(output_df,c(df$az));
            node_tags = c(node_tags, node);
          }
          if (output_type == "Fx") {
            df = fread(paste0(directory,"/node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");
            output_df = cbind(output_df,c(df$Fx));
            node_tags = c(node_tags, node);
          }
          if (output_type == "Fy") {
            df = fread(paste0(directory,"/node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");
            output_df = cbind(output_df,c(df$Fy));
            node_tags = c(node_tags, node);
          }
          if (output_type == "Mz") {
            df = fread(paste0(directory,"/node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");
            output_df = cbind(output_df,c(df$Mz));
            node_tags = c(node_tags, node);
          }
          if (output_type == "disp") {
            df = fread(paste0(directory,"/node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
          }
          if (output_type == "reac") {
            df = fread(paste0(directory,"/node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");  
          }
          output_list[[n_list]] = df;
          names(output_list)[n_list] = paste0(node);
        }
      }
    }
  }
  if (output_type == "disp" || output_type == "reac") {output_list} else {
    colnames(output_df) = node_tags
    return(output_df)
  }
};
int_col_node_output = function (directory,nB,nL,output_type){
  
  #######################################################################################################
  # Function that loads OpenSees output from the tcl procedure ::rack_pack::record_rack 
  # when in_col_nodes is passed as a parameter to that procedure
  #
  # Inputs :
  # directory := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # output_type := a string, either "reac" or "disp" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 3 mai 2016, EJ
  #######################################################################################################
  
  nc = nB + 1
  output_list = vector("list", (nB-1)*nL);
  n_list = 0;
  
  node_tags = NULL
  output_df = NULL
  
  for (i in c(1:(nL+1))) {
    for (j in c(1:nc)) {
      if ( i != 1 ) { 
        if ( j != 1 && j != nc ) {
          n_list = n_list+1;
          node =  paste0(i,j);
          if (output_type == "dx") {
            df = fread(paste0(directory,"/node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$dx));
            node_tags = c(node_tags, node);
          }
          if (output_type == "dy") {
            df = fread(paste0(directory,"/node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$dy));
            node_tags = c(node_tags, node);
          }
          if (output_type == "rz") {
            df = fread(paste0(directory,"/node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
            output_df = cbind(output_df,c(df$rz));
            node_tags = c(node_tags, node);
          }
          if (output_type == "ax") {
            df = fread(paste0(directory,"/node_acc",node,".out"));
            colnames(df) = c("ax","ay","az");
            output_df = cbind(output_df,c(df$ax));
            node_tags = c(node_tags, node);
          }
          if (output_type == "ay") {
            df = fread(paste0(directory,"/node_acc",node,".out"));
            colnames(df) = c("ax","ay","az");
            output_df = cbind(output_df,c(df$ay));
            node_tags = c(node_tags, node);
          }
          if (output_type == "az") {
            df = fread(paste0(directory,"/node_acc",node,".out"));
            colnames(df) = c("ax","ay","az");
            output_df = cbind(output_df,c(df$az));
            node_tags = c(node_tags, node);
          }
          if (output_type == "Fx") {
            df = fread(paste0(directory,"/node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");
            output_df = cbind(output_df,c(df$Fx));
            node_tags = c(node_tags, node);
          }
          if (output_type == "Fy") {
            df = fread(paste0(directory,"/node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");
            output_df = cbind(output_df,c(df$Fy));
            node_tags = c(node_tags, node);
          }
          if (output_type == "Mz") {
            df = fread(paste0(directory,"/node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");
            output_df = cbind(output_df,c(df$Mz));
            node_tags = c(node_tags, node);
          }
          if (output_type == "disp") {
            full_path = paste0(directory,"node_disp",node,".out");
            df = fread(paste0(directory,"node_disp",node,".out"));
            colnames(df) = c("dx","dy","rz");
          }
          if (output_type == "reac") {
            df = fread(paste0(directory,"node_reac",node,".out"));
            colnames(df) = c("Fx","Fy","Mz");  
          }
          output_list[[n_list]] = df;
          names(output_list)[n_list] = paste0(node);
        }
      }
    }
  }
  if (output_type == "disp" || output_type == "reac") {output_list} else {
    colnames(output_df) = node_tags
    return(output_df)
  };
};
column_end_output = function (dir,nB,nL,n_div,output_type) {
  
  ####################################################################################################
  # //TODO: Describe what the function does
  #
  # Inputs :
  # dir := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := of elements the column in each bay at each level has been divided into
  # output_type := a string, "Fx", "Fy", "Mz" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 21 nov 2016, EJ
  ####################################################################################################
  
  output_df = NULL
  ele_names = NULL
  
  if ( output_type == "Fx" || output_type == "Fy" || output_type == "Mz") {
    for (i in c(1:(nB+1))) { 
      for (j in c(1:nL)) { 
        eleBottom = paste0(j,i,j,i,0,1);
        eleTop = paste0(j,i,0,n_div-1,j+1,i);
        ele_names = c(ele_names, eleBottom);
        ele_names = c(ele_names, eleTop);
        dfBottom = fread(paste0(dir,"/ele_force",eleBottom,".out"));
        dfTop = fread(paste0(dir,"/ele_force",eleTop,".out"));
        colnames(dfBottom) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
        colnames(dfTop) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
        switch(output_type, 
               Fx={
                 output_df = cbind(output_df, dfBottom$Fxi);
                 output_df = cbind(output_df, dfTop$Fxj);
               },
               Fy={
                 output_df = cbind(output_df, dfBottom$Fyi);
                 output_df = cbind(output_df, dfTop$Fyj);
               },
               Mz={
                 output_df = cbind(output_df, dfBottom$Mzi);
                 output_df = cbind(output_df, dfTop$Mzj);
               }
        )
      }
    }
  }
  
  if ( output_type == "ex" || output_type == "ez") {
    for (i in c(1:(nB+1))) {
      for (j in c(1:nL)) {
        eleBottom = paste0(j,i,j,i,0,1);
        eleTop = paste0(j,i,0,n_div-1,j+1,i);
        ele_names = c(ele_names, eleBottom);
        ele_names = c(ele_names, eleTop);
        dfBottom = fread(paste0(dir,"/ele_def",eleBottom,".out"));
        dfTop = fread(paste0(dir,"/ele_def",eleTop,".out"));
        colnames(dfBottom) = c("exi","ezi","zeroi");
        colnames(dfTop) = c("exj","ezj","zeroj");
        switch(output_type,
               ex={
                 output_df = cbind(output_df, dfBottom$exi);
                 output_df = cbind(output_df, dfTop$exj);
               },
               ez={
                 output_df = cbind(output_df, dfBottom$ezi);
                 output_df = cbind(output_df, dfTop$ezj);
               }
        )
      }
    }
  }
  
  colnames(output_df) = ele_names;
  return(output_df)
  
};
beam_end_output = function (dir,nB,nL,n_div,output_type) {
  
  ####################################################################################################
  # //TODO: Describe what the function does
  #
  # Inputs :
  # dir := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # n_div := of elements the column in each bay at each level has been divided into
  # output_type := a string, "Fx", "Fy", "Mz" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 21 nov 2016, EJ
  ####################################################################################################
  
  output_df = NULL
  ele_names = NULL
  
  for (i in c(2:(nL+1))) { 
    for (j in c(1:nB)) { 
      eleLeft = paste0(i,j,1,i,j,2);
      eleRight = paste0(i,j,n_div,i,j,n_div+1);
      
      fileLeft = paste0(dir,"/ele_force",eleLeft,".out");
      fileRight = paste0(dir,"/ele_force",eleRight,".out");
      
      if (file.exists(fileLeft)) {
        ele_names = c(ele_names, eleLeft);
        dfLeft = fread(fileLeft);
        colnames(dfLeft) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
      }
      
      if (file.exists(fileRight)) {
        ele_names = c(ele_names, eleRight);
        dfRight = fread(fileRight);
        colnames(dfRight) = c("Fxi","Fyi","Mzi","Fxj","Fyj","Mzj");
      }
      
      switch(output_type, 
             Fx={
               if (file.exists(fileLeft)) {
                output_df = cbind(output_df, dfLeft$Fxi);
               }
               if (file.exists(fileRight)) {
                output_df = cbind(output_df, dfRight$Fxj);
               }
             },
             Fy={
               if (file.exists(fileLeft)) {
                output_df = cbind(output_df, dfLeft$Fyi);
               }
               if (file.exists(fileRight)) {
                output_df = cbind(output_df, dfRight$Fyj);
               }
             },
             Mz={
               if (file.exists(fileLeft)) {
                output_df = cbind(output_df, dfLeft$Mzi);
               }
               if (file.exists(fileRight)) {
                output_df = cbind(output_df, dfRight$Mzj);
               }
             }
      )
    }
  }
  colnames(output_df) = ele_names;
  return(output_df)
};
mid_col_node_disp = function (dir,nB,nL,output_type) {
  
  ####################################################################################################
  # //TODO: Describe what the function does
  #
  # Inputs :
  # dir := a string, the directory where the opensees output can be found
  # nB := an integer, the number of bays
  # nL := an integer, the number of levels
  # output_type := a string, "dx", "dy", "rz" depending on what output information is requested
  #
  # outputs : A list of dataframes containing the requested output information
  #
  # Last mod. 21 nov 2016, EJ
  ####################################################################################################
  
  output_df = NULL
  node_names = NULL
  
  for (i in c(1:(nB+1))) { 
    for (j in c(1:nL)) { 
      node = paste0(j,i,"04");
      node_names = c(node_names, node);
      df = fread(paste0(dir,"/node_disp",node,".out"));
      colnames(df) = c("dx", "dy", "rz");
      switch(output_type, 
         dx={
           output_df = cbind(output_df, df$dx);
         },
         dy={
           output_df = cbind(output_df, df$dy);
         },
         rz={
           output_df = cbind(output_df, df$rz);
         }
      )
    }
  }
  colnames(output_df) = node_names;
  return(output_df)
};


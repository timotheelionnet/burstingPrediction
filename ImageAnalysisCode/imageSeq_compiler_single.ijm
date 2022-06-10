run("Image Sequence...", "open=[Image Path Here - Input] file=F02f01 sort");
run("Out [-]");
run("Stack to Hyperstack...", "order=xyzct channels=2 slices=1 frames=120 display=Color");
run("Green"); run("Next Slice [>]"); run("Red");
selectWindow("col1+2");
run("Bleach Correction", "correction=[Exponential Fit]");
selectWindow("DUP_col1+2");
saveAs("Tiff", "Image Path Here - Output");
{ 
      while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
  } 

//********************* FILE NAMES GO HERE **************************************
//root folder where all data is
input = getDirectory("choose your input directory");

//where the analysis images will go 
output = getDirectory("choose your output directory");

//******************************************************************************


//batch mode so nothing is displayed to accelerate execution
setBatchMode(true); 

//create analysis directories if they do not exist

tmplist = getFileList(input);
list = newArray(10000);
ctr = 0;

for (k=0; k<tmplist.length; k++) {
 	
    if (endsWith(tmplist[k], "/") || endsWith(tmplist[k], "\\") ){
       //listFiles(""+dir+list[i]);
       //print(tmplist[k]);
		tmplist2 = getFileList(input+tmplist[k]);	
		for (l=0; l<tmplist2.length; l++) {
			if(indexOf(tmplist2[l], ".tif") != -1){
	        	list[ctr] = input+tmplist[k]+tmplist2[l];
	        	//print(list[ctr]);
	        	ctr++;
        	}
		}  
    }
    else{
    	if(indexOf(tmplist[k], ".tif") != -1){
	        	list[ctr] = input+tmplist[k];
	        	//print(list[ctr]);
	        	ctr++;
        	}
    }     
 }
 list = Array.trim(list,ctr);
 
 print("list of files to treat");
 for (k=0; k<list.length; k++) {
 	print(list[k]);
 }
print("starting file processing...");
//now that file names for the conditions are all combined, lets do the analysis
for (i = 0; i < list.length; i++){
		//print file name in log window
		//print(d2s(indexOf(list[i], ".tif"),1));
		
        if (indexOf(list[i], ".tif") != -1) { // making sure it is an image
	        open(list[i]); 
        	title = getTitle();
        	print("processing file "+ title);
			getDimensions(dummy,dummy,channelcount,dummy,dummy);
			if (channelcount>=0){
				setOption("ScaleConversions", true);
				run("16-bit");
					saveAs(".tif", output+title);
					close();
				}
			}
			else {
				close();
			}
		}

setBatchMode(false); 
        
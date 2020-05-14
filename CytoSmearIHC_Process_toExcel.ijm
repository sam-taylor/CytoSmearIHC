/*
 * Macro to process a folder of scanned .svs files 
 * and output a quantification of the IHC staining in a separate
 * folder. The macro is optimized for whole-slide images scanned
 * at 40X magnification. Different magnifications can be used by
 * making modifications to the Bio-Formats Importer instructions.
 * for files starting at low magnification, for example, change
 * "series_3" to "series_2" or "series_1".
 * 
 * This macro relies on the macro "CytoSmearIHC," which also 
 * needs to be installed in the macros folder of ImageJ.
 * 
 * Output files include a PNG of the total DAB-positive pixels 
 * from the original image, a PNG of the area included in the
 * quantification, a PNG of the quantification histogram, and a
 * .csv file containing the histogram data.
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".svs") suffix

//modify these values as necessary to exclude cell clumps and debris
#@ String (label = "Cell Size Range (pixels)", value = "000-500") cellSize

//modify these values as necessary to completely exclude background in the converted 8-bit images
#@ Integer (label= "Cell Threshold Low", min=0, max=254, value=0) cellThreshLo
#@ Integer (label= "Cell Threshold High (must be > low)", min=0, max=254, value=222) cellThreshHi

// tab variables
var windowTabName="Stat Results Table",nameOfStatTab="["+windowTabName+"]",label="",undoErease="";
// stat variables
var nbPerim=0,TheTotalArea=0,meanPerim=0,meanObject=0;nCells=0;
var v = newArray(256);

setBatchMode(true); //batch mode on
processFolder(input);
setBatchMode(false); //exit batch mode

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	//now save the accumulated results table as an excel file	
}

function processFile(input, output, file) {
	file = file;
	print(file);
	theFile = "open=[" + input + File.separator + file + "]";
	run("Bio-Formats Importer", theFile +  " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_3 series_6");
	saveAs("PNG",  output + File.separator + file + "_Label.png");
	close();
	run("Stack to RGB");

	//these commands create a subsection in the center of the image
  	getDimensions(w, h, channels, slices, frames);
	makeRectangle(w/4, h/4, w/2, h/2);
	run("Crop");
	saveAs("PNG",  output + File.separator + file + "_Input.png");
	cropped = getImageID();

	//here a threshold is applied to remove background pixels from the analysis
	run("Duplicate...", " ");
	run("8-bit");
	setAutoThreshold("Default");
	setThreshold(cellThreshLo, cellThreshHi); 
	run("Convert to Mask");

	//exclude clumped cells
	run("Set Measurements...", "area perimeter display redirect=None decimal=2");
	run("Analyze Particles...", "size=" + cellSize + " show=Masks display clear include record ");
	if (! isOpen("Results")) {exit ("No Results table")}
	makeStatFromResults_AP ();
	//print(meanObject);
	nCells = nResults();
	run("Close");

	//now need to select the area decided as cells
	setAutoThreshold("Default");
	setThreshold(50, 255);//these values should always capture the particles
	//saveAs("PNG",  output + File.separator + file + "_Included.png");
	run("Create Selection");
	
	//the non-background, non-clumped analysis area is applied to the DAB channel
	//of the colour-deconvoluted image
	selectImage(cropped);
	run("Colour Deconvolution", "vectors=[H DAB]");
	close();
	run("Restore Selection");
	saveAs("PNG",  output + File.separator + file + "_TotalDAB.png");

	//DAB intensity of the analysis area is quantified and results are saved
	runMacro("CytoSmearIHC");
	makeStatFromResults_HS ();
	//print("/t" + v[10]);
	getStats();
	
	// function making stats from ImageJ Results Table values. This will work for analyze particles
	//I'll make another function like this for the histogram part

	saveAs("PNG", output + File.separator + file + "_Hist.png");
	saveAs("Results", output + File.separator + file + ".csv");
	print("Processing: " + input + File.separator + file);
	print("Saving to: " + output);
	print(" ");
}

						// --------------- tab functions ---------------//

function makeStatFromResults_AP () {
	
	sumObject=0;sumPerim=0;TheTotalArea=0;meanObject=0;meanPerim=0;
	// extraction from data from the Result Table:
	if (nResults() > 0 && isOpen("Results")) {
		//here instead of the weird image label I'll use the filename as label
		label = file;
		if (columnLabelList ("Area") >= 0) {sumObject=sumColumnTab ("Area"); } else {exit ("No \"Area\" measurements in this Result Table");}
		if (columnLabelList ("Perim.") >= 0) {sumPerim=sumColumnTab ("Perim."); } else {exit ("No \"Perim.\" measurements in this Result Table");}
	} else {exit ("No result table");}
	// Stat calculations:
	TheTotalArea = sumObject;
	if (nResults() != 0) meanObject = TheTotalArea/nResults();
	if (nResults() != 0) meanPerim=sumPerim/nResults();
}

// function returning the number of the column which name is contained in kind. return -1 if doesn't exists
function columnLabelList (kind) {

	columnNumber=-1;
	if (nResults() > 0 && isOpen("Results")) {
		selectWindow("Results");
   		results = getInfo();
   		lines = split(results, "\n");
  		headings = lines[0];
		titlesofcolumns = split(headings, ",\t");
		for (a=0; a<titlesofcolumns.length; a++) {if (titlesofcolumns[a] == kind) columnNumber=a;}
	}
	return columnNumber;
}

// function getting the sum of the value from the column name contained in kind
function sumColumnTab (kind) {
	sum=0;
	if (columnLabelList (kind) >=0) {
		for (a=0; a<nResults(); a++) {
			sum=sum+getResult(kind,a);
		}
	return sum;
}
}
	
function makeStatFromResults_HS () {
	// extraction from data from the Result Table:
	if (nResults() > 0 && isOpen("Results")) {
		if (columnLabelList ("Value") >= 0) {
			recordValues ("Value");
			} 
			else {exit ("No \"Value\" measurements in this Result Table");}
		}
	else {exit ("No result table");}
}

function recordValues (kind) {
	var v = newArray(256);
	if (columnLabelList (kind) >=0) {
		for (a=0; a<nResults(); a++) {
			v[a]=getResult(kind,a);
		}
	}
}	


//takes current results window and extracts the data into the stats table
function getStats() {
	undoErease = ""; windowTabName="Stat Results Table";nameOfStatTab="["+windowTabName+"]";
	if (isOpen(windowTabName)) {addLine=1;} 
	else {addLine=0;}
	//I'll just beef up the below function to add all the variables
	TabWindow (addLine);
}

// function building/managing a table window
function TabWindow (addLine) {
	undoErease = "";
	if (! isOpen(windowTabName)) {	
		run("New... ", "name="+nameOfStatTab+" type=Table");
		valueHeaders = "";
		for (i=0; i<v.length; i++) {
			valueHeaders = valueHeaders + "\t" + i;
		}
		print(nameOfStatTab, "\\Headings:Slice Name\tCount Objects\tTotal Area\tAverage Size\tMean Perim" + valueHeaders);
	}
	if (addLine == 0) {print(nameOfStatTab, "\\Clear");} 
	values = "";
		for (i=0; i<v.length; i++) {
			values = values + "\t" + v[i];
		}
	print(nameOfStatTab,  label+ "\t" + d2s(nCells,2) + "\t" + TheTotalArea + "\t" + d2s(meanObject,2) +  "\t" + d2s(meanPerim,2) + values);
	print(nCells);
	print(d2s(nCells,2));
}

// function removing the last line of the tab
function rmLastLine () {

	if ( isOpen (windowTabName)) {
		selectWindow (windowTabName);
   		tabContent = getInfo();
   		linesInTab = split(tabContent, "\n");
		if (linesInTab[linesInTab.length-1] > 0) {
			print(nameOfStatTab, "\\Clear");
			resteLines="";
			for (i=1; i < (linesInTab.length -1); i++) {
				resteLines=resteLines+linesInTab[i] +"\n";
			}
			if (linesInTab.length > 2) print (nameOfStatTab,resteLines);
			if (linesInTab.length > 1) undoErease=linesInTab[linesInTab.length-1];
		}
	}
}


// function restoring the last ereased line in the table
function undormLastLine () {

	if (undoErease != "") print(nameOfStatTab,undoErease);
	undoErease="";
}

// saving a tab as excel type file
function saveTab (path) {

	if (isOpen(windowTabName)) {
		if (path == "") {
			selectWindow (windowTabName);
			run("Input/Output...", "jpeg=75 gif=-1 file=.xls");
			saveAs("Text");
		}
		if (path != "") {
			selectWindow(windowTabName);
			saveAs("Text", path+windowTabName+".xls");
		}
	}
}

function openTab (path,name) {

	undoErease="";windowTabName=name;
	lines=split(File.openAsString(path), "\n");
	if (lines.length < 2) { exit ("This file doesn't seam to contain data");}
  	headings = lines[0];
	titlesOfColumns = split(headings, ",\t");
	nameOfStatTab="["+windowTabName+"]";
	if (isOpen(windowTabName)) {selectWindow(windowTabName) ;run("Close");}	
	firstLine="";
	for (i=0; i < (titlesOfColumns.length ); i++) {
		firstLine=firstLine+ titlesOfColumns [i];
		if ( i < (titlesOfColumns.length )-1) {firstLine=firstLine+ "\t";}
	}
	toPrint="";
	for (i=1; i < (lines.length ); i++) {
		toPrint=toPrint+lines[i]+"\n";
	}
	run("New... ", "name="+nameOfStatTab+" type=Table");
	print(nameOfStatTab, "\\Headings:"+firstLine+"");
	print(nameOfStatTab,toPrint);
}



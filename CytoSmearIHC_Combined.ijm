/* CytoSmearIHC
 * by Samuel Taylor
 * Weill Cornell Medical College
 * Tri-Institutional MD-PhD Program
 * 
 * Macro to process a folder of scanned .svs files 
 * and output a quantification of the IHC staining in a separate
 * folder. The macro is optimized for whole-slide images scanned
 * at 40X magnification. Different magnifications can be used by
 * making modifications to the Bio-Formats Importer instructions.
 * for files starting at low magnification, for example, change
 * "series_3" to "series_2" or "series_1". "series_6" refers to
 * the scanned slide label--adjust this command as needed.
 * 
 * Output files include a PNG of the total DAB-positive pixels 
 * from the original image, a PNG of the area included in the
 * quantification, a PNG of the quantification histogram, and a
 * tab-delimited text file with statistics for all images 
 * analyzed.
 * 
 * This macro employs code from the following sources:
 * 
 * "CustomTabStatFromResults"
 * by Gilles Carpentier
 * Faculte des Sciences et Technologies,
 * Universite Paris 12 Val de Marne, France.
 * 
 * Varghese F, Bukhari AB, Malhotra R, De A (2014) IHC Profiler: 
 * An Open Source Plugin for the Quantitative Evaluation and 
 * Automated Scoring of Immunohistochemistry Images of Human 
 * Tissue Samples. PLoS ONE 9(5): e96801. 
 * https://doi.org/10.1371/journal.pone.0096801
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
// global variables for hist
var x = 0.025, y = 0.1;

setBatchMode(true); //batch mode on
processFolder(input);
saveTable();
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
	//save the series_6 image as the label
	saveAs("PNG",  output + File.separator + file + "_Label.png");
	close();
	run("Stack to RGB");

	//these commands create a subsection in the center of the image
  	getDimensions(w, h, channels, slices, frames);
	makeRectangle(w/4, h/4, w/2, h/2);
	run("Crop");
	//save the cropped image as the input
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
	nCells = nResults();
	run("Close");

	//now need to select the area decided as cells
	setAutoThreshold("Default");
	setThreshold(50, 255);//these values should always capture the particles
	//this image reflects the area to be analyzed
	saveAs("PNG",  output + File.separator + file + "_Included.png");
	run("Create Selection");
	
	//the non-background, non-clumped analysis area is applied to the DAB channel
	//of the colour-deconvoluted image
	selectImage(cropped);
	run("Colour Deconvolution", "vectors=[H DAB]");
	close();
	run("Restore Selection");
	//this is the total DAB image
	saveAs("PNG",  output + File.separator + file + "_TotalDAB.png");

	//DAB intensity of the analysis area is quantified and results are saved
	getHistogramStats();
	makeStatFromResults_HS ();
	getStats();
	
	// function making stats from ImageJ Results Table values. This will work for analyze particles
	//I'll make another function like this for the histogram part

	saveAs("PNG", output + File.separator + file + "_Hist.png");
	//saveAs("Results", output + File.separator + file + ".csv");
	print(" ");
	closeAllWindows();
}

						// --------------- tab functions ---------------//

 function closeAllWindows () { 
      while (nImages>0) {
          selectImage(nImages);
          close(); 
      } 
  } 

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

// This macro has been developed to calculate the number of pixels of
// different color intensities. Based on these it plots a histogram
// profile and assigns a grade to the image. 

function getHistogramStats(){
	bins = 256;
	maxCount = 0;
	histMin = 0;
	histMax = 0;
	
	if (histMax > 0)
	    getHistogram(values, counts, bins, histMin, histMax);
	else
	    getHistogram(values, counts, bins);
	
	is8bits = bitDepth() == 8 || bitDepth() == 24;
	
	Plot.create("Histogram", "Pixel Value", "Count", values, counts);
	
	if (maxCount > 0)
	    Plot.setLimits(0, 256, 0, maxCount);
		
	n = 0;
	sum = 0;
	min = 9999999;
	max = -9999999;
	Region2 = 0;
	Region3 = 0;
	Region4 = 0;
	Region1 = 0;
	Region0 = 0;
	TotalPixel = 0;
	PercentRegion1 = 0;
	PercentRegion2 = 0;
	PercentRegion4 = 0;
	PercentRegion3 = 0;
	PercentRegion0 = 0;
	Score = 0;
	PixelUnderConsideration = 0;
	
	for (i = 0; i < bins; i++) {
	    count = counts[i];
	    if (count > 0) {
	        n += count;
	        sum += count * i;
	        if (i < min) min = i;
	        if (i > max) max = i;
	    }
	}
	
	print("Pixel Count: " + n);
	
	if (is8bits) {
	    for (i = 0; i < bins; i++) {
	        if (i >= 0 && i < 60)
	            Region4 = Region4 + counts[i];
	        if (i > 60 && i < 120)
	            Region3 = Region3 + counts[i];
	        if (i > 120 && i < 180)
	            Region2 = Region2 + counts[i];
	        if (i > 180 && i < 255)
	            Region1 = Region1 + counts[i];
	        if (i > 255 && i <= 256)
	            Region0 = Region0 + counts[i];
	    }
	}
	
	TotalPixel = TotalPixel + Region1 + Region2 + Region3 + Region4 + Region0;
	PixelUnderConsideration = TotalPixel - Region0;
	
	PercentRegion3 = (Region3 / PixelUnderConsideration) * 100;
	PercentRegion2 = (Region2 / PixelUnderConsideration) * 100;
	PercentRegion1 = (Region1 / PixelUnderConsideration) * 100;
	PercentRegion4 = (Region4 / PixelUnderConsideration) * 100;
	
	print("Percentage contribution of High Positive:  " + PercentRegion4);
	print("Percentage contribution of Positive:  " + PercentRegion3);
	print("Percentage contribution of Low Positive:  " + PercentRegion2);
	print("Percentage contribution of Negative:  " + PercentRegion1);
	
	if ((PercentRegion3 > 66) || (PercentRegion2 > 66) || (PercentRegion1 > 66) || (PercentRegion4 > 66)) {
	    if (PercentRegion4 > 66)
	        print("The score is High Positive  ");
	    if (PercentRegion3 > 66)
	        print("The score is Positive  ");
	    if (PercentRegion2 > 66)
	        print("The score is Low Positive  ");
	    if (PercentRegion1 > 66)
	        print("The score is Negative  ");
	} else {
	    Score = Score + (PercentRegion4 / 100) * 4 + (PercentRegion3 / 100) * 3 + (PercentRegion2 / 100) * 2 + (PercentRegion1 / 100) * 1;
	    if (Score >= 2.95) {
	        print("The score is High Positive  ");
	    }
	    if ((Score >= 1.95) && (Score <= 2.94)) {
	        print("The score is Positive  ");
	    }
	    if ((Score >= 0.95) && (Score <= 1.94)) {
	        print("The score is Low Positive  ");
	    }
	    if ((Score >= 0.0) && (Score <= 0.94)) {
	        print("The score is Negative  ");
	    }
	}
	
	Array.show("Results (row numbers)", counts);
}

function draw (text) {
    Plot.addText(text, x, y);
    y += 0.08;
}

function saveTable () {
	if (! isOpen(windowTabName)) {exit ("No Stat Results Table")}
	selectWindow(windowTabName);
	saveAs("Text", output + File.separator + windowTabName + ".txt");
}



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
}

function processFile(input, output, file) {
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
	run("Analyze Particles...", "size=" + cellSize + " show=Masks exclude include");
	setAutoThreshold("Default");
	setThreshold(50, 255);//these values should always capture the particles
	saveAs("PNG",  output + File.separator + file + "_Included.png");
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
	saveAs("PNG", output + File.separator + file + "_Hist.png");
	saveAs("Results", output + File.separator + file + ".csv");
	print("Processing: " + input + File.separator + file);
	print("Saving to: " + output);
	print(" ");
}

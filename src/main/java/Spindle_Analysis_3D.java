

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


/**
 *
 * @author thomasb and philippem
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.RankFilters;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import ij.util.ThreadUtil;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.FastFilters3D;


public class Spindle_Analysis_3D implements PlugIn {

    private final boolean canceled = false;
    private int nbChannel = 0;
    private Calibration cal = new Calibration();
    private boolean mtocs = true;
    ImageHandler img;

    
// Find channel name in ND file
    public String[] Find_nd_info(String NDdir, String NDfile) throws FileNotFoundException, IOException {
        BufferedReader br = new BufferedReader(new FileReader(NDdir + NDfile));
        int nbline = 0;
        String strLine;
        String[] line;
        nbChannel = 0;
        String[] waveName = {"", "", ""};
        int nbWave = 1;
        while ((strLine = br.readLine()) != null) {
            if (strLine.startsWith("\"NWavelengths\"")) {    // find channel number
                line = strLine.split(", ");
                nbChannel = Integer.valueOf(line[1]);
            }

            if (strLine.startsWith("\"WaveName" + nbWave + "\"")) {
                line = strLine.split("\"");
                waveName[nbWave - 1] = "_w" + nbWave + line[3];
                nbWave++;
            }
        }
        return (waveName);
    }

/*
    Dialogbox Mtocs or Chromosomes
    */ 
    public boolean dialogMtocs() {
        boolean cancel = false;
        String choice;
        GenericDialog gd = new GenericDialog("MTOCS or Chromosomes");
        String[] choices = {"MTOCS", "Chromosomes"};
        gd.addRadioButtonGroup(null, choices, 1, 2, "MTOCS");
        gd.showDialog();
        if(gd.wasCanceled()) cancel = true;
        choice = gd.getNextRadioButton();
        if (choice == "MTOCS") mtocs = true;
        return cancel;
    }
    
   
// find ovocyte boundaries 
    public Roi find_crop(ImagePlus BG) {
        int x, y, w, h;
        ImageProcessor ip = BG.getProcessor();
        RankFilters var = new RankFilters();
        ResultsTable rt = new ResultsTable();
        ParticleAnalyzer measure = new ParticleAnalyzer(ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES + ParticleAnalyzer.INCLUDE_HOLES,
                ParticleAnalyzer.RECT + ParticleAnalyzer.AREA, rt, 200000, 6500000);

        var.rank(ip, 25, RankFilters.VARIANCE);
        ip.setAutoThreshold(AutoThresholder.Method.Triangle, true);
        measure.analyze(BG, ip);
        // crop if object found
        if (rt.getCounter() != 0) {
            x = (int) rt.getValue("BX", 0);
            y = (int) rt.getValue("BY", 0);
            w = (int) rt.getValue("Width", 0);
            h = (int) rt.getValue("Height", 0);
        }
        else {
            x = 0;
            y = 0;
            w = BG.getWidth();
            h = BG.getHeight();
        }
        return (new Roi(x, y, w, h));
    }

// crop GFP or RFP to keep only ovocyte boundaries
    public ImagePlus crop_image(ImagePlus img, Roi roicrop) {

        ImageProcessor ipGfp = img.getProcessor();
        img.setRoi(roicrop);
        ipGfp.crop();
        img.updateAndDraw();
        ImagePlus imgcrop = new Duplicator().run(img);
        return (imgcrop);

    }

// filter image 
    public ImagePlus image_filter(ImagePlus img, String title) {
        // run pure denoise filter
       IJ.run(img, "PureDenoise ...", "parameters='2 3' estimation='Auto Global' ");
       WindowManager.getWindow("Log").setVisible(false);
       img.close();
       img.flush();
       ImagePlus imgDenoised = WindowManager.getCurrentImage();
        GaussianBlur gaussian = new GaussianBlur();
        ImageStack stack = imgDenoised.getStack();
        for (int s = 1; s <= imgDenoised.getNSlices(); s++) {
            imgDenoised.setSlice(s);
            gaussian.blurGaussian(stack.getProcessor(s), 2, 2, 0.02);
        }
        imgDenoised.updateAndDraw();
        imgDenoised.setTitle(title+"_filtered");
        imgDenoised.hide();
        img.flush();
        img.close();
        //IJ.run(imgDenoised,"8-bit","");
        return (imgDenoised);
    }

    /**
     * 
     * @param img
     * @return 
     */
    private void image_threshold(ImagePlus img) {
        img.setZ(img.getNSlices()/2);
        img.updateAndDraw();
        IJ.setAutoThreshold(img, "MaxEntropy dark");
        IJ.run(img, "Convert to Mask", "method=MaxEntropy background=Dark black");
        IJ.run(img, "Options...", "iterations=1 count=1 black do=[Fill Holes] stack");
        img.setCalibration(cal);
    }
    
// Differences of  Gaussian filter RFP image if chromosomes
// 3D adaptative fast filters if mtocs 
    public ImagePlus rfp_filter(ImagePlus img) {
        if (mtocs) {
            ImageInt imgInt = ImageInt.wrap(img);
            IJ.showStatus("filtering RFP ....");
            ImageInt AdaptLocalStack = FastFilters3D.filterIntImage(imgInt, FastFilters3D.ADAPTIVE, 3, 3, 3, ThreadUtil.getNbCpus(), true);
            img = AdaptLocalStack.getImagePlus();
            img.setCalibration(cal);
        }
        else {
            RankFilters median = new RankFilters();
            IJ.run(img, "Difference of Gaussians", "  sigma1=20 sigma2=1 stack");
            for (int s = 1; s <= img.getNSlices(); s++) {
                img.setSlice(s);
                median.rank(img.getProcessor(), 1, RankFilters.MEDIAN);
            }
        }
        img.updateAndDraw();
        img.setTitle("RFP_filtered");
        return (img);
    }
    
    
    // tag spindle number
    void tagsObject(ImageHandler imgObj, Object3D spots, int n) {
        ImagePlus img = imgObj.getImagePlus();
        img.setSlice((int)spots.getCenterZ());
        ImageProcessor ip = img.getProcessor();
        Font tagFont = new Font("SansSerif", Font.PLAIN, 9);
        ip.setFont(tagFont);
        ip.setColor(Color.blue);
        int x = (int)spots.getCenterX();
        int y = (int)spots.getCenterY();
        ip.drawString(Integer.toString(n), x, y);
        img.updateAndDraw();
    }

    @Override
    public void run(String arg) {

        cal.pixelWidth = 0.1613;
        cal.pixelHeight = 0.1613;
        cal.pixelDepth = 1;
        cal.setUnit("microns");

        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            String imageDir = IJ.getDirectory("Choose Directory Containing ND Files...");
            if (imageDir == null) {
                return;
            }
            File inDir = new File(imageDir);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                return;
            }
            if (dialogMtocs() == true) return;
            String[] channelName = null;
            // Write headers results for distance file 
            FileWriter fwDistAnalyze;
            fwDistAnalyze = new FileWriter(imageDir + "Analyze_MTOCS_results.xls", false);
            BufferedWriter outputDistAnalyze = new BufferedWriter(fwDistAnalyze);
            outputDistAnalyze.write("Image\tSpindle Vol\tSpindle Feret\tSphericity\tFlatness\tCompactness\tMtoc Volume\tMtoc Distance to pole\t+"
                    + "Border Mtoc distance to border spindle\tCenter Mtoc distance to border spindle\t"
                    + "Border Mtoc distance to spindle center\tCenter Mtoc distance to spindle center\tEVF Spindle\tEVF Mtocs\tColoc (common voxels)\n");
            outputDistAnalyze.flush();

            // 
            for (int f = 0; f < imageFile.length; f++) {
                if (imageFile[f].endsWith(".nd")) {
                    channelName = Find_nd_info(imageDir, imageFile[f]);
                    String imageName = imageFile[f].substring(0, imageFile[f].indexOf(".nd"));
                    Opener opener = new Opener();
                    ImagePlus brightfield = new ImagePlus();

                    // open brighField image and find Roi to crop images
                    Roi cropRoi;
                    brightfield = opener.openTiff(imageDir, imageName + channelName[0] + ".TIF");
                    cropRoi = find_crop(brightfield);
                    brightfield.close();
                    brightfield.flush();

                    // open GFP channel
                    ImagePlus gfp = new ImagePlus();
                    ImagePlus gfp_crop = new ImagePlus();
                    // check gfp channel
                    if (channelName[1].endsWith("GFP")) {
                        gfp = opener.openTiff(imageDir, imageName + channelName[1] + ".TIF");
                    }
                    else {
                        gfp = opener.openTiff(imageDir, imageName + channelName[2] + ".TIF");
                    }
                    IJ.showStatus("Processing image " + imageName);
                    // crop gfp image to Ovocyte size
                    gfp_crop = crop_image(gfp, cropRoi);
                    gfp.close();
                    gfp.flush();

                    // filter spindle (GFP)
                    gfp_crop = image_filter(gfp_crop, "GFP");
                    image_threshold(gfp_crop);
                    FileSaver gfpSave = new FileSaver(gfp_crop);
                    gfpSave.saveAsTiffStack(imageDir + imageName + channelName[1] + "_mask.tif");

                    //open RFP channel (MTOCS)
                    ImagePlus rfp = new ImagePlus();
                    ImagePlus rfp_crop = new ImagePlus();
                    // check RFP channel
                    if (channelName[1].endsWith("RFP")) {
                        rfp = opener.openTiff(imageDir, imageName + channelName[1] + ".TIF");
                    }
                    else {
                        rfp = opener.openTiff(imageDir, imageName + channelName[2] + ".TIF");
                    }
                    //IJ.run(rfp,"8-bit","");
                    // crop gfp image to Ovocyte size
                    rfp_crop = crop_image(rfp, cropRoi);
                    rfp.close();
                    rfp.flush();

                    // filter spindle (RFP)
                    rfp_crop = image_filter(rfp_crop, "RFP");
                    image_threshold(rfp_crop);
                    FileSaver rfpSave = new FileSaver(rfp_crop);
                    rfpSave.saveAsTiffStack(imageDir + imageName + channelName[2] + "_mask.tif");

                    // ASSUME WE HAVE BINARY IMAGES HERE
                    IJ.showStatus("Computing object population ...");
                    Objects3DPopulation popTmp = getPopFromImage(gfp_crop);
                    Objects3DPopulation spots = getPopFromImage(rfp_crop);
                    // For spindle if more than one object take only the biggest
                    int index = 0;
                    if (popTmp.getNbObjects() > 1) {
                        double volume = 1000;   // minimum size for spindle
                        for (int i = 0; i < popTmp.getNbObjects(); i++) {
                            if (popTmp.getObject(i).getVolumeUnit() > volume) {
                                volume = popTmp.getObject(i).getVolumeUnit();
                                index = i;
                            }
                        }
                    }
                    // Compute infos if some objects exist and size >= 10000       
                    if ((popTmp.getNbObjects() > 0) && (spots.getNbObjects() > 0)) {
                        if (index >= 0) {
                            Object3D spindle = popTmp.getObject(index);
                            IJ.showStatus("computing distances ...");
                            computeInfo(gfp_crop, spindle, spots, outputDistAnalyze, imageDir, imageName);
                        }
                    }
                    gfp_crop.flush();
                    rfp_crop.flush();
                }
            }
            outputDistAnalyze.close();
            IJ.showStatus("Finished");
        } catch (IOException ex) {
            Logger.getLogger(Spindle_Analysis_3D.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(img.getCalibration());
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }

    private void computeInfo(ImagePlus imgSpindle, Object3D spindle, Objects3DPopulation mtocs, BufferedWriter results, String inDir, String image) throws IOException {
        // EVF spindle info 
        ImageInt img = ImageInt.wrap(imgSpindle);
        ImageFloat edtSpindle = EDT.run(img, 128, (float) cal.pixelWidth, (float) cal.pixelDepth, false, 0);
        EDT.normalizeDistanceMap(edtSpindle, img, true);

        // EVF Mtocs info
        Voxel3D Feret1 = spindle.getFeretVoxel1();
        Voxel3D Feret2 = spindle.getFeretVoxel2();
        ImageHandler poles = img.createSameDimensions();
        poles.setPixel(Feret1, 255);
        poles.setPixel(Feret2, 255);
        ImageFloat edtMtocs = EDT.run(poles, 128, (float) cal.pixelWidth, (float) cal.pixelDepth, true, 0);
        EDT.normalizeDistanceMap(edtMtocs, img, true);

        double Feret_length = Feret1.distance(Feret2, cal.pixelWidth,cal.pixelDepth);
        IJ.log("Image: " + image);
        IJ.log("Feret :" + Feret1 + " " + Feret2 + "\n Spindle volume : " + spindle.getVolumeUnit());
        IJ.log("Feret_lenght : " + Feret_length);
        IJ.log("#Mtocs : " + mtocs.getNbObjects());

        ImageHandler imgObjects = img.createSameDimensions();
        imgObjects.set332RGBLut();
        imgObjects.setCalibration(cal);
        spindle.draw(imgObjects, 29);
        int nMtocs = 0;
        double coloc;
        double sphSpindle = spindle.getSphericity(true);
        double flatSpindle = spindle.getMedianElongation();
        double compSpindle = spindle.getCompactness(true);
        for (int i = 0; i < mtocs.getNbObjects(); i++) {
            // nbre de pixel colocalise 
            coloc = spindle.getColoc(mtocs.getObject(i));
            if (spindle.getColoc(mtocs.getObject(i)) > 10) {
                nMtocs++;
                IJ.log("Mtocs volume :" + mtocs.getObject(i).getVolumeUnit());
                
                double distBorder = mtocs.getObject(i).distBorderUnit(spindle);
                IJ.log("Border spot to spindle border " + i + " : " + distBorder);
                
                double distCenterBorder = mtocs.getObject(i).distCenterBorderUnit(spindle);
                IJ.log("Center spot to spindle border " + i + " : " + distCenterBorder);
                
                double distCenter = spindle.distCenterBorderUnit(mtocs.getObject(i));
                IJ.log("Border spot to spindle center " + i + " : " + distCenter);
                
                double distCenterCenter = spindle.distCenterUnit(mtocs.getObject(i));
                IJ.log("Center spot to spindle center " + i + " : " + distCenterCenter);
                
                double dist1 = mtocs.getObject(i).distPixelCenter(Feret1.x, Feret1.y, Feret1.z);
                double dist2 = mtocs.getObject(i).distPixelCenter(Feret2.x, Feret2.y, Feret2.z);
                IJ.log("spot to poles " + i + " : " + dist1 + " " + dist2);
                IJ.log("spot to poles " + i + " : " + Math.min(dist1, dist2));
                
                double EVFSpindle = edtSpindle.getPixel(mtocs.getObject(i).getCenterAsPoint());
                IJ.log("EVF mtocs spindle " + i + " : " + EVFSpindle);
                
                double EVFMtocs = edtMtocs.getPixel(mtocs.getObject(i).getCenterAsPoint());
                IJ.log("EVF mtocs poles " + i + " : " + EVFMtocs);
                results.write(image + "\t" + spindle.getVolumeUnit() + "\t" + Feret_length + "\t" + sphSpindle + "\t" + flatSpindle + "\t"+ compSpindle + "\t"+
                        mtocs.getObject(i).getVolumeUnit() + "\t"
                        + Math.min(dist1, dist2) + "\t" + distBorder + "\t" + distCenterBorder + "\t" + distCenter + "\t" + distCenterCenter + "\t" + EVFMtocs + "\t" +  EVFSpindle + "\t" + coloc + "\n");
                results.flush();
                mtocs.getObject(i).draw(imgObjects, 228);
                // tag object by it number
                tagsObject(imgObjects, mtocs.getObject(i), nMtocs);
            }
        }
        FileSaver objectsFile = new FileSaver(imgObjects.getImagePlus());
        objectsFile.saveAsTiffStack(inDir + image + "_objects.tif");
    }

}

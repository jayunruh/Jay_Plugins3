/*******************************************************************************
 * Copyright (c) 2012 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/
import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import jguis.*;
import jalgs.*;

public class stitch_mosaic_image_jru_v1 implements PlugIn, FrameInterface, gui_interface {
	boolean showrois;
	int tempchan,tempslice,tempframe,slices,channels,frames,trueframes;
	ImageStack tempstack;

	public void run(String arg) {
		ImagePlus imp=WindowManager.getCurrentImage();
		ImageWindow[] iws=jutils.selectPlots(false,1,new String[]{"Position_Plot"});
		if(iws==null) return;
		Object[] retvals=exec(imp,iws[0]);
		ImagePlus retimp=(ImagePlus)retvals[0];
		if(imp.getNChannels()==1) retimp.show();
		else{
			(new CompositeImage(retimp,CompositeImage.COLOR)).show();
		}
	}

	public static float[] getAvgOverlap(float[][] coords,int width,int height){
		float hover=0.0f; int nhpairs=0;
		float vover=0.0f; int nvpairs=0;
		float halfwidth=0.5f*(float)width;
		float halfheight=0.5f*(float)height;
		for(int i=0;i<coords[0].length;i++){
			for(int j=(i+1);j<coords[0].length;j++){
				float xdist=Math.abs(coords[0][j]-coords[0][i]);
				float ydist=Math.abs(coords[1][j]-coords[1][i]);
				if(xdist>halfwidth && xdist<(float)width && ydist<halfheight){
					hover+=xdist; nhpairs++;
				}
				if(xdist<halfwidth && ydist>halfheight && ydist<(float)height){
					vover+=ydist; nvpairs++;
				}
			}
		}
		hover/=(float)nhpairs;
		vover/=(float)nvpairs;
		return new float[]{(float)width-hover,(float)height-vover};
	}

	public Object[] exec(ImagePlus imp,ImageWindow iw){
		float[][] xvals=(float[][])jutils.runPW4VoidMethod(iw,"getXValues");
		float[][] yvals=(float[][])jutils.runPW4VoidMethod(iw,"getYValues");
		int sel=(Integer)jutils.runPW4VoidMethod(iw,"getSelected"); if(sel<0) sel=0;
		float[] overlap=getAvgOverlap(new float[][]{xvals[sel],yvals[sel]},imp.getWidth(),imp.getHeight());
		IJ.log("avg overlap = "+overlap[0]+" , "+overlap[1]);
		return exec(imp,xvals[sel],yvals[sel],overlap[0],overlap[1]);
	}

	public Object[] exec(ImagePlus imp,float[] xvals,float[] yvals,float hoverlap,float voverlap){
		int width=imp.getWidth(); int height=imp.getHeight();
		ImageStack stack=imp.getStack();
		int size=stack.getSize();
		channels=imp.getNChannels();
		slices=imp.getNSlices();
		frames=imp.getNFrames();
		if(frames==1){
			frames=slices; slices=1;
		}
		trueframes=1;
		if(frames>xvals.length){
			trueframes=(int)(frames/xvals.length);
			frames=xvals.length;
		}
		float psize=(float)jutils.get_psize(imp);
		int typeindex=algutils.get_array_type(stack.getPixels(1));
		stitching sclass=new stitching(width,height,xvals,yvals,psize,typeindex,this);
		ImageStack stack2=new ImageStack(sclass.newwidth,sclass.newheight);
		tempstack=stack;
		for(int k=0;k<trueframes;k++){
			for(int i=0;i<slices;i++){
				for(int j=0;j<channels;j++){
					tempchan=j;
					tempslice=i;
					tempframe=k*frames;
					Object stitched=sclass.stitch_frame(this,frames,true,hoverlap,voverlap);
					stack2.addSlice("",stitched);
				}
			}
		}
		if(showrois){
			int[][] rois=sclass.getRois();
			RoiManager rman=RoiManager.getInstance();
			if(rman==null) rman=new RoiManager();
			for(int i=0;i<rois.length;i++){
				rman.addRoi(new Roi(rois[i][0],rois[i][1],rois[i][2],rois[i][3]));
			}
		}
		ImagePlus imp5=new ImagePlus("Stitched Image",stack2);
		imp5.copyScale(imp);
		imp5.setOpenAsHyperStack(true);
		imp5.setDimensions(channels,slices,trueframes);
		//imp5.show();
		return new Object[]{imp5};
	}

	public Object getNextFrame(){
		Object temp=jutils.get3DSlice(tempstack,tempframe,tempslice,tempchan,frames,slices,channels);
		tempframe++;
		return temp;
	}

	public void showMessage(String message){
		IJ.log(message);
	}
	
	public void showProgress(int currpos,int finalpos){
		IJ.showProgress(currpos,finalpos);
	}

}
